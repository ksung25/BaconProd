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
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h" 
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h" 
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <map>

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
  fEERecHitName(iConfig.getUntrackedParameter<std::string>("edmEERecHitName","reducedEcalRecHitsEE")),
  fRhoName     (iConfig.getUntrackedParameter<std::string>("edmRhoForEnergyRegression","kt6PFJets")),
  fPVName       (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices"))
{
  fPhotonMVA = new PhotonMVACalculator();
  fPhotonReg = new EGEnergyCorrector();  
}

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
  if(!fPhotonReg->IsInitialized()) { 
      std::string cmssw_base_utils = getenv("CMSSW_BASE"); 
      cmssw_base_utils+="/src/BaconProd/Utils/data/";
      fPhotonReg->Initialize(iSetup,cmssw_base_utils+"gbrv3ph_52x.root");
      fPhotonMVA->initialize(cmssw_base_utils+"2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15.weights.xml",cmssw_base_utils+"/2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15.weights.xml");
  }

  
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

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *vtxCol = hVertexProduct.product();

  //Get SuperCluster Collection
  edm::Handle<reco::SuperClusterCollection> hEESCProduct;
  iEvent.getByLabel(fEESCName,hEESCProduct);
  assert(hEESCProduct.isValid());
  const reco::SuperClusterCollection *eeSCCol = hEESCProduct.product();

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  edm::InputTag rhoTag(fRhoName,"rho","RECO");
  iEvent.getByLabel(rhoTag,hRho);
  assert(hRho.isValid()); 

  //PreShower crap
  edm::Handle<EcalRecHitCollection> hESRecHits;
  iEvent.getByLabel("reducedEcalRecHitsES" , hESRecHits);
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  //const CaloSubdetectorGeometry *geometryES = pGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  EcalPreshowerTopology topology_p(pGeometry);

  //map of preshower rechits for shape calculations                                                                                                                                                   
  std::map<DetId, EcalRecHit> esmap;
  if (hESRecHits.isValid()) {
    EcalRecHitCollection::const_iterator it;
    for (it = hESRecHits->begin(); it != hESRecHits->end(); ++it) {
      // remove bad ES rechits                                                                                                                                                                      
      if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
      //Make the map of DetID, EcalRecHit pairs
      esmap.insert(std::make_pair(it->id(), *it));
    }
  }
    
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

    //Regression corrections
    std::pair<double,double> lCorr = fPhotonReg->CorrectedEnergyWithErrorV3(*itPho,*vtxCol,*hRho,lazyTools,iSetup);
    pPhoton->ptcorr = itPho->pt()*lCorr.first /itPho->energy();
    pPhoton->ptunc  = itPho->pt()*lCorr.second/itPho->energy();
    
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
    //Isolation for Photon MVA
    pPhoton->chHadIso03SelVtx  = -1;
    pPhoton->chHadIso03WstVtx  = -1;
    computeVtxIso(*itPho,*pfCandCol,*vtxCol,
		  pPhoton->chHadIso03SelVtx,pPhoton->chHadIso03WstVtx);
   
    //Preshower
    double lRR = 0; 
    if(fabs(sc->eta()) > 1.5)  { 
      std::vector<float> phoESShape = getESShape(getESHits(sc->x(),sc->y(),sc->z(), esmap, *pGeometry.product(), &topology_p, 0));
      float lXX = phoESShape[0];
      float lYY = phoESShape[1];
      lRR = sqrt(lXX*lXX+lYY*lYY);
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

    //Apply the gg MVA
    pPhoton->mva             = fPhotonMVA->mvaValue((*itPho),lazyTools,*hRho,pPhoton->gammaIso03,pPhoton->chHadIso03SelVtx,pPhoton->chHadIso03WstVtx,lRR);
    
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

    //Regression correction
    //std::pair<double,double> lCorr = fPhotonReg->CorrectedEnergyWithErrorV3(*pPhoton,*vtxCol,*hRho,lazyTools,iSetup);
    //pPhoton->ptcorr = itPho->pt()/itPho->e()*lCorr.first;
    //pPhoton->ptunc  = itPho->pt()/itPho->e()*lCorr.second;

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

    pPhoton->chHadIso03SelVtx  = -1;
    pPhoton->chHadIso03WstVtx  = -1;

    if(pho.isNonnull()) {
      computeVtxIso(*pho,*pfCandCol,*vtxCol,
		    pPhoton->chHadIso03SelVtx,pPhoton->chHadIso03WstVtx);
    }

    pPhoton->isoForFsr03     = computeIsoForFSR(&(*itPF), pfNoPU, pfPU);
    pPhoton->mvaNothingGamma = itPF->mva_nothing_gamma();
    //Preshower
    double lRR = 0; 
    if(sc.isNonnull() && fabs(sc->eta()) > 1.5)  { 
      std::vector<float> phoESShape = getESShape(getESHits(sc->x(),sc->y(),sc->z(), esmap, *pGeometry.product(), &topology_p, 0));
      float lXX = phoESShape[0];
      float lYY = phoESShape[1];
      lRR = sqrt(lXX*lXX+lYY*lYY);
    }
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
    if(pho.isNonnull()) pPhoton->mva          = fPhotonMVA->mvaValue(*pho,lazyTools,*hRho,pPhoton->gammaIso03,pPhoton->chHadIso03SelVtx,pPhoton->chHadIso03WstVtx,lRR);
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
//--------------------------------------------------------------------------------------------------
void FillerPhoton::computeVtxIso(const reco::Photon &photon,
				 const std::vector<reco::PFCandidate>        &pf,
				 const std::vector<reco::Vertex>             &iVertex,
				 float &out_chHadIsoWvtx,float &out_chHadIsoFirstVtx) const 
{
  double extRadius = 0.3;
  double intRadius = 0.02;
  for(unsigned int iVtx=0; iVtx<iVertex.size(); iVtx++) { 
    const reco::Vertex vertex = iVertex.at(iVtx);
    double pIso = 0; 
    for(unsigned int ipf=0; ipf<pf.size(); ipf++) {      
      const reco::PFCandidate pfcand = pf.at(ipf);
      if(pfcand.particleId() != reco::PFCandidate::h) continue;
      // Add p_T to running sum if PFCandidate is close enough
      TVector3 pVec; 
      pVec.SetXYZ(photon.superCluster()->position().x()-vertex.position().x(),
		  photon.superCluster()->position().y()-vertex.position().y(),
		  photon.superCluster()->position().z()-vertex.position().z());
      
      double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), pVec.Eta(), pVec.Phi());
      if(dr >= extRadius || dr < intRadius) continue;
      double dZ = pfcand.trackRef()   ->dz (vertex.position());
      double d0 = pfcand.trackRef()   ->dxy(vertex.position());
      if(fabs(dZ) > 0.2) continue;
      if(fabs(d0) > 0.1) continue;
      pIso  += pfcand.pt(); 
    }
    if(iVtx == 0)               out_chHadIsoFirstVtx = pIso;
    if(out_chHadIsoWvtx < pIso) out_chHadIsoWvtx     = pIso;
  }
}
//--------------------------------------------------------------------------------------------------
//horrible code below copied from globe for preshower cluster shape calculations
std::vector<float> FillerPhoton::getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row) {
  std::vector<float> esHits;

  const GlobalPoint point(X,Y,Z);

  const CaloSubdetectorGeometry *geometry_p ;
  geometry_p = geometry.getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

  std::map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;
    
  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);
    
  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}
//--------------------------------------------------------------------------------------------------
std::vector<float> FillerPhoton::getESShape(std::vector<float> ESHits0)
{
  std::vector<float> esShape;

  const int nBIN = 21;
  float esRH_F[nBIN];
  float esRH_R[nBIN];
  for (int idx=0; idx<nBIN; idx++) {
    esRH_F[idx] = 0.;
    esRH_R[idx] = 0.;
  }

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
    if (ibin==0) {
      esRH_F[(nBIN-1)/2] = ESHits0[ibin];
      esRH_R[(nBIN-1)/2] = ESHits0[ibin+31];
    } else {
      esRH_F[(nBIN-1)/2+ibin] = ESHits0[ibin];
      esRH_F[(nBIN-1)/2-ibin] = ESHits0[ibin+15];
      esRH_R[(nBIN-1)/2+ibin] = ESHits0[ibin+31];
      esRH_R[(nBIN-1)/2-ibin] = ESHits0[ibin+31+15];
    }
  }

  // ---- Effective Energy Deposit Width ---- //                                                                                                                                                     
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;
  for (int id_X=0; id_X<21; id_X++) {
    totalEnergyXX  += esRH_F[id_X];
    EffStatsXX     += esRH_F[id_X]*(id_X-10)*(id_X-10);
    totalEnergyYY  += esRH_R[id_X];
    EffStatsYY     += esRH_R[id_X]*(id_X-10)*(id_X-10);
  }
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;

  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);

  return esShape;
}


