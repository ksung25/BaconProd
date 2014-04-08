#include "BaconProd/Ntupler/interface/FillerPhoton.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerPhoton::FillerPhoton(const edm::ParameterSet &iConfig):
  fMinPt       (iConfig.getUntrackedParameter<double>("minPt",10)),
  fPhotonName  (iConfig.getUntrackedParameter<std::string>("edmName","photons")),
  fPFCandName  (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fEleName     (iConfig.getUntrackedParameter<std::string>("edmElectronName","gsfElectrons")),
  fConvName    (iConfig.getUntrackedParameter<std::string>("edmConversionName","allConversions")),
  fEBSCName    (iConfig.getUntrackedParameter<std::string>("edmEBSuperClusterName","correctedHybridSuperClusters")),
  fEESCName    (iConfig.getUntrackedParameter<std::string>("edmEESuperClusterName","correctedMulti5x5SuperClustersWithPreshower")),
  fEBRecHitName(iConfig.getUntrackedParameter<std::string>("edmEBRecHitName","reducedEcalRecHitsEB")),
  fEERecHitName(iConfig.getUntrackedParameter<std::string>("edmEERecHitName","reducedEcalRecHitsEE"))
{}

//--------------------------------------------------------------------------------------------------
FillerPhoton::~FillerPhoton(){}

//--------------------------------------------------------------------------------------------------
void FillerPhoton::fill(TClonesArray *array, 
                        const edm::Event &iEvent, const edm::EventSetup &iSetup,
                        const reco::Vertex &pv,
			const std::vector<const reco::PFCandidate*> &pfNoPU,
		        const std::vector<const reco::PFCandidate*> &pfPU,
		        const std::vector<TriggerRecord> &triggerRecords,
		        const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
  
  // Get photon collection
  edm::Handle<reco::PhotonCollection> hPhotonProduct;
  iEvent.getByLabel(fPhotonName,hPhotonProduct);
  assert(hPhotonProduct.isValid());
  const reco::PhotonCollection *photonCol = hPhotonProduct.product();

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByLabel(fPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();

  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> hEleProduct;
  iEvent.getByLabel(fEleName,hEleProduct);
  assert(hEleProduct.isValid());
  
  // Get conversions collection
  edm::Handle<reco::ConversionCollection> hConvProduct;
  iEvent.getByLabel(fConvName,hConvProduct);
  assert(hConvProduct.isValid());

  // Get SuperCluster collections
  edm::Handle<reco::SuperClusterCollection> hEBSCProduct;
  iEvent.getByLabel(fEBSCName,hEBSCProduct);
  assert(hEBSCProduct.isValid());
  const reco::SuperClusterCollection *ebSCCol = hEBSCProduct.product();
  
  edm::Handle<reco::SuperClusterCollection> hEESCProduct;
  iEvent.getByLabel(fEESCName,hEESCProduct);
  assert(hEESCProduct.isValid());
  const reco::SuperClusterCollection *eeSCCol = hEESCProduct.product();

    
  edm::InputTag ebRecHitTag(fEBRecHitName);
  edm::InputTag eeRecHitTag(fEERecHitName);
  EcalClusterLazyTools lazyTools(iEvent, iSetup, ebRecHitTag, eeRecHitTag);
  
  std::vector<const reco::PFCandidate*> usedPFPhotons;  // keep track of PF photons that are also counted as standard photons
  
  // PF photon cuts for HZZ4l FSR recovery
  const double pfMinPt  = 2;
  const double pfMaxEta = 2.4;
  
  for(reco::PhotonCollection::const_iterator itPho = photonCol->begin(); itPho!=photonCol->end(); ++itPho) {
    
    // Photon cuts
    if(itPho->pt() < fMinPt) continue;
    
    // construct object and place in array
    TClonesArray &rPhotonArr = *array;
    assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
    const int index = rPhotonArr.GetEntries();
    new(rPhotonArr[index]) baconhep::TPhoton();
    baconhep::TPhoton *pPhoton = (baconhep::TPhoton*)rPhotonArr[index];

    const reco::SuperClusterRef sc = itPho->superCluster();
    std::vector<float> vCov = lazyTools.localCovariances(*(sc->seed()));
    
    //
    // Kinematics
    //==============================
    pPhoton->pt  = itPho->pt();
    pPhoton->eta = itPho->eta();
    pPhoton->phi = itPho->phi();
    
    pPhoton->scEt  = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pPhoton->scEta = sc->eta();
    pPhoton->scPhi = sc->phi();
    pPhoton->r9    = lazyTools.e3x3(*(sc->seed())) / (sc->rawEnergy());

    // consider standard photon also to be a PF photon if they share supercluster
    pPhoton->pfPt  = 0;
    pPhoton->pfEta = 0;
    pPhoton->pfPhi = 0;
    bool hasPFMatch=false;
    for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
      if(itPF->particleId() != reco::PFCandidate::gamma) continue;
      if(itPF->pt()        < pfMinPt)  continue;
      if(fabs(itPF->eta()) > pfMaxEta) continue;
      
      const reco::SuperClusterRef pfsc = itPF->superClusterRef();
      if(pfsc.isNull()) continue;
      
      if(pfsc == sc) {
        hasPFMatch = true;
        usedPFPhotons.push_back(&(*itPF));
	pPhoton->pfPt  = itPF->pt();
	pPhoton->pfEta = itPF->eta();
	pPhoton->pfPhi = itPF->phi();
	break;
      }
    }

    //
    // Isolation
    //==============================
    pPhoton->trkIso04  = itPho->trkSumPtHollowConeDR04();
    pPhoton->ecalIso04 = itPho->ecalRecHitSumEtConeDR04();
    pPhoton->hcalIso04 = itPho->hcalTowerSumEtConeDR04();
    
    computeIso(*itPho, pfNoPU,
               pPhoton->chHadIso03,
	       pPhoton->gammaIso03,
	       pPhoton->neuHadIso03);

    pPhoton->isoForFsr03 = -1;
    if(hasPFMatch) {
      const reco::PFCandidate *pfcand = usedPFPhotons.back();
      pPhoton->isoForFsr03 = computeIsoForFSR(pfcand, pfNoPU, pfPU);
      pPhoton->mvaNothingGamma = pfcand->mva_nothing_gamma();
    }
    
    //
    // Identification
    //==============================    
    pPhoton->hovere = itPho->hadronicOverEm();
    pPhoton->sieie  = itPho->sigmaIetaIeta();
    pPhoton->sipip  = isnan(vCov[2]) ? 0 : sqrt(vCov[2]);
    
    pPhoton->fiducialBits=0;
    if(itPho->isEB())        pPhoton->fiducialBits |= kIsEB;
    if(itPho->isEE())        pPhoton->fiducialBits |= kIsEE;
    if(itPho->isEBEEGap())   pPhoton->fiducialBits |= kIsEBEEGap;
    if(itPho->isEBEtaGap())  pPhoton->fiducialBits |= kIsEBEtaGap;
    if(itPho->isEBPhiGap())  pPhoton->fiducialBits |= kIsEBPhiGap;
    if(itPho->isEEDeeGap())  pPhoton->fiducialBits |= kIsEEDeeGap;
    if(itPho->isEERingGap()) pPhoton->fiducialBits |= kIsEERingGap;
    
    pPhoton->typeBits=0;
    if(itPho->isStandardPhoton()) pPhoton->typeBits |= baconhep::kEGamma;    // only standard photons in AOD 'photons' collection?
    if(itPho->isPFlowPhoton())    pPhoton->typeBits |= baconhep::kPFPhoton;  // always 'false' for standard photons?
    if(hasPFMatch)                pPhoton->typeBits |= baconhep::kPFPhoton;  // consider standard photon to be PF if they share supercluster    
   
    pPhoton->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator itSC = ebSCCol->begin(); itSC!=ebSCCol->end(); ++itSC) {
      scIndex++;
      if(itPho->superCluster().get() == &(*itSC)) {
        pPhoton->scID = scIndex;
	break;
      }
    }
    for(reco::SuperClusterCollection::const_iterator itSC = eeSCCol->begin(); itSC!=eeSCCol->end(); ++itSC) {
      scIndex++;
      if(itPho->superCluster().get() == &(*itSC)) {
        pPhoton->scID = scIndex;
	break;
      }
    }
    
    pPhoton->hasPixelSeed = itPho->hasPixelSeed();
    
    pPhoton->isConv = ConversionTools::hasMatchedPromptElectron(itPho->superCluster(), hEleProduct, hConvProduct, pv.position(), 2.0, 1e-6, 0);
    
    pPhoton->hltMatchBits = TriggerTools::matchHLT(pPhoton->eta, pPhoton->phi, triggerRecords, triggerEvent);
  }

  //
  // Include PF photons for HZZ4l FSR recovery
  // Note: ECAL energy associated with PFMuons also considered FSR candidates
  //
  for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
    
    TLorentzVector pfpho(0,0,0,0);    
    if(itPF->particleId() == reco::PFCandidate::mu && itPF->ecalEnergy()>0) {
      pfpho.SetPtEtaPhiM(itPF->pt()*itPF->ecalEnergy()/itPF->p(), itPF->eta(), itPF->phi(), 0);
    } else if(itPF->particleId() == reco::PFCandidate::gamma) {
      pfpho.SetPtEtaPhiM(itPF->pt(), itPF->eta(), itPF->phi(), 0);
    }
    
    if(pfpho.Pt()        <= pfMinPt)  continue;
    if(fabs(pfpho.Eta()) >= pfMaxEta) continue;
    
    bool isUsed=false;
    for(unsigned int iused=0; iused<usedPFPhotons.size(); iused++) {
      if(&(*itPF) == usedPFPhotons[iused]) {
        isUsed=true;
	break;
      }
    }
    if(isUsed) continue;
    
    // construct object and place in array
    TClonesArray &rPhotonArr = *array;
    assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
    const int index = rPhotonArr.GetEntries();
    new(rPhotonArr[index]) baconhep::TPhoton();
    baconhep::TPhoton *pPhoton = (baconhep::TPhoton*)rPhotonArr[index];

    const reco::PhotonRef       pho = itPF->photonRef();
    const reco::SuperClusterRef sc  = itPF->superClusterRef();
    std::vector<float> vCov;
    if(sc.isNonnull()) vCov = lazyTools.localCovariances(*(sc->seed()));
    
    //
    // Kinematics
    //==============================
    pPhoton->pt  = pfpho.Pt();
    pPhoton->eta = pfpho.Eta();
    pPhoton->phi = pfpho.Phi();
    
    pPhoton->scEt  = 0;
    pPhoton->scEta = 0;
    pPhoton->scPhi = 0;
    pPhoton->r9    = 0;    
    if(sc.isNonnull()) {
      pPhoton->scEt  = (sc->energy())*(sc->position().Rho())/(sc->position().R());
      pPhoton->scEta = sc->eta();
      pPhoton->scPhi = sc->phi();
      pPhoton->r9    = lazyTools.e3x3(*(sc->seed())) / (sc->rawEnergy());
    }
    
    // consider standard photon also to be a PF photon if they share supercluster
    pPhoton->pfPt  = pfpho.Pt();
    pPhoton->pfEta = pfpho.Eta();
    pPhoton->pfPhi = pfpho.Phi();
    
    //
    // Isolation
    //==============================
    pPhoton->trkIso04  = -1;
    pPhoton->ecalIso04 = -1;
    pPhoton->hcalIso04 = -1;
    
    pPhoton->chHadIso03  = -1;
    pPhoton->gammaIso03  = -1;
    pPhoton->neuHadIso03 = -1;

    if(pho.isNonnull()) {
      computeIso(*pho, pfNoPU,
                 pPhoton->chHadIso03,
	         pPhoton->gammaIso03,
	         pPhoton->neuHadIso03);
    }
    
    pPhoton->isoForFsr03     = computeIsoForFSR(&(*itPF), pfNoPU, pfPU);
    pPhoton->mvaNothingGamma = itPF->mva_nothing_gamma();
    
    //
    // Identification
    //==============================
    if(pho.isNonnull()) pPhoton->hovere = pho->hadronicOverEm();
    if(pho.isNonnull()) pPhoton->sieie  = pho->sigmaIetaIeta();
    if(sc.isNonnull())  pPhoton->sipip  = isnan(vCov[2]) ? 0 : sqrt(vCov[2]);
    
    pPhoton->fiducialBits=0;
    if(pho.isNonnull()) {
      if(pho->isEB())        pPhoton->fiducialBits |= kIsEB;
      if(pho->isEE())        pPhoton->fiducialBits |= kIsEE;
      if(pho->isEBEEGap())   pPhoton->fiducialBits |= kIsEBEEGap;
      if(pho->isEBEtaGap())  pPhoton->fiducialBits |= kIsEBEtaGap;
      if(pho->isEBPhiGap())  pPhoton->fiducialBits |= kIsEBPhiGap;
      if(pho->isEEDeeGap())  pPhoton->fiducialBits |= kIsEEDeeGap;
      if(pho->isEERingGap()) pPhoton->fiducialBits |= kIsEERingGap;
    }
    
    pPhoton->typeBits = 0;
    if(itPF->particleId() == reco::PFCandidate::mu)    pPhoton->typeBits |= baconhep::kPFMuonPhoton;
    if(itPF->particleId() == reco::PFCandidate::gamma) pPhoton->typeBits |= baconhep::kPFPhoton;
    if(pho.isNonnull() && pho->isStandardPhoton())     pPhoton->typeBits |= baconhep::kEGamma;
       
    pPhoton->scID = -1;
    int scIndex = -1;
    if(sc.isNonnull()) {
      for(reco::SuperClusterCollection::const_iterator itSC = ebSCCol->begin(); itSC!=ebSCCol->end(); ++itSC) {
        scIndex++;
        if(sc.get() == &(*itSC)) {
          pPhoton->scID = scIndex;
  	  break;
        }
      }
      for(reco::SuperClusterCollection::const_iterator itSC = eeSCCol->begin(); itSC!=eeSCCol->end(); ++itSC) {
        scIndex++;
        if(sc.get() == &(*itSC)) {
          pPhoton->scID = scIndex;
	  break;
        }
      }
    }
    
    pPhoton->hasPixelSeed = false;    
    pPhoton->isConv       = false;    
    pPhoton->hltMatchBits = TriggerTools::matchHLT(pPhoton->eta, pPhoton->phi, triggerRecords, triggerEvent);        
  }  
}

//--------------------------------------------------------------------------------------------------
void FillerPhoton::computeIso(const reco::Photon &photon,
                              const std::vector<const reco::PFCandidate*> &pfNoPU,
                              float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const 
{
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double extRadius       = 0.3;
  const double intRadiusChHad  = 0.02;
  const double intRadiusNeuHad = 0.;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<pfNoPU.size(); ipf++) {
    const reco::PFCandidate *pfcand = pfNoPU.at(ipf);
    
    if(pfcand->superClusterRef().isNonnull() && photon.superCluster().isNonnull() &&
       pfcand->superClusterRef() == photon.superCluster())
      continue;
    
    if     (pfcand->particleId() == reco::PFCandidate::h)     { intRadius = intRadiusChHad; }
    else if(pfcand->particleId() == reco::PFCandidate::gamma) { intRadius = photon.isEB() ? 0.015 : 0.00864*fabs(TMath::SinH(photon.superCluster()->eta()))*4.; }
    else if(pfcand->particleId() == reco::PFCandidate::h0)    { intRadius = intRadiusNeuHad; }
    
    double dr = reco::deltaR(pfcand->eta(), pfcand->phi(), photon.eta(), photon.phi());
    if(dr>=extRadius || dr<intRadius) continue;
    
    if     (pfcand->particleId() == reco::PFCandidate::h)     { chHadIso  += pfcand->pt(); }
    else if(pfcand->particleId() == reco::PFCandidate::gamma) { gammaIso  += pfcand->pt(); }
    else if(pfcand->particleId() == reco::PFCandidate::h0)    { neuHadIso += pfcand->pt(); }
  }
  
  out_chHadIso  = chHadIso;
  out_gammaIso  = gammaIso;
  out_neuHadIso = neuHadIso;
}

//--------------------------------------------------------------------------------------------------
float FillerPhoton::computeIsoForFSR(const reco::PFCandidate *photon,
                                     const std::vector<const reco::PFCandidate*> &pfNoPU,
			             const std::vector<const reco::PFCandidate*> &pfPU) const
{ 
  double extRadius = 0.3;
  double intRadius = 0.01;
  
  double chIso  = 0;
  double gamIso = 0;
  double neuIso = 0;
  double puIso  = 0;
  
  // compute isolation contributions from charged hadrons, neutral hadrons, and photons from the PV
  for(unsigned int ipf=0; ipf<pfNoPU.size(); ipf++) {      
    const reco::PFCandidate *pfcand = pfNoPU.at(ipf);
    
    if(pfcand==photon) continue;
    
    // Add p_T to running sum if PFCandidate is close enough
    double dr = reco::deltaR(pfcand->eta(), pfcand->phi(), photon->eta(), photon->phi());
    if(dr >= extRadius || dr < intRadius) continue;
    
    if     (pfcand->particleId() == reco::PFCandidate::h     && pfcand->pt() > 0.2) { chIso  += pfcand->pt(); }
    else if(pfcand->particleId() == reco::PFCandidate::gamma && pfcand->pt() > 0.5) { gamIso += pfcand->pt(); }
    else if(pfcand->particleId() == reco::PFCandidate::h0    && pfcand->pt() > 0.5) { neuIso += pfcand->pt(); }    
  }
  
  // compute isolation contributions from charged particles not from PV
  for(unsigned int ipf=0; ipf<pfPU.size(); ipf++) {
    const reco::PFCandidate *pfcand = pfPU.at(ipf);
    assert(pfcand);

    if(pfcand->trackRef().isNull()) continue;
    if(pfcand->pt() < 0.2)          continue;
    
    // Add p_T to running sum if PFCandidate is within isolation cone
    double dr = reco::deltaR(pfcand->eta(), pfcand->phi(), photon->eta(), photon->phi());
    if(dr > extRadius || dr <= intRadius) continue;
    
    puIso += pfcand->pt();
  }
  
  return (chIso + gamIso + neuIso + puIso);
}
