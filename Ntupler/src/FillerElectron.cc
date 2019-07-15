#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <utility>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerElectron::FillerElectron(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt                 (iConfig.getUntrackedParameter<double>("minPt",7)),
  fEleName               (iConfig.getUntrackedParameter<std::string>("edmName","gedGsfElectrons")),
  fBSName                (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot")),
  fPFCandName            (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName             (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fConvName              (iConfig.getUntrackedParameter<std::string>("edmConversionName","allConversions")),
  fSCName                (iConfig.getUntrackedParameter<edm::InputTag>("edmSCName")),
  fPuppiName             (iConfig.getUntrackedParameter<std::string>("edmPuppiName","puppi")),
  fPuppiNoLepName        (iConfig.getUntrackedParameter<std::string>("edmPuppiNoLepName","puppiNoLep")),
  fUsePuppi              (iConfig.getUntrackedParameter<bool>("usePuppi",true)),
  fMVASpring16           (iConfig.getUntrackedParameter<std::string>("edmEleMVASpring16")),
  fMVAFall17V1Iso        (iConfig.getUntrackedParameter<std::string>("edmEleMVAFall17V1Iso")),
  fMVAFall17V1NoIso      (iConfig.getUntrackedParameter<std::string>("edmEleMVAFall17V1NoIso")),
  fMVAFall17V2Iso        (iConfig.getUntrackedParameter<std::string>("edmEleMVAFall17V2Iso")),
  fMVAFall17V2NoIso      (iConfig.getUntrackedParameter<std::string>("edmEleMVAFall17V2NoIso")),
  fMVASpring16HZZ        (iConfig.getUntrackedParameter<std::string>("edmEleMVASpring16HZZ")),
  fMediumMVA             (iConfig.getUntrackedParameter<std::string>("edmEleMediumMVA")),
  fTightMVA              (iConfig.getUntrackedParameter<std::string>("edmEleTightMVA")),
  //fMVA                   (iConfig.getUntrackedParameter<std::string>("edmEleMVA")),
  //fMediumMVAIso          (iConfig.getUntrackedParameter<std::string>("edmEleMediumMVAIso")),
  //fTightMVAIso           (iConfig.getUntrackedParameter<std::string>("edmEleTightMVAIso")),
  //fMVAIso                (iConfig.getUntrackedParameter<std::string>("edmEleMVAIso")),
  //fSecondMVA             (iConfig.getUntrackedParameter<bool>("storeSecondMVA",false)),
  //fStoreHZZMVA           (iConfig.getUntrackedParameter<bool>("storeHZZMVA",false)),
  fUseTO                 (iConfig.getUntrackedParameter<bool>("useTriggerObject",false)),
  fUseAOD                (useAOD)
{
  if(fUseAOD)  fTokEleName        = iC.consumes<reco::GsfElectronCollection>(fEleName);
  if(!fUseAOD) fTokPatEleName     = iC.consumes<pat::ElectronCollection>(fEleName);
  fTokSCName         = iC.consumes<reco::SuperClusterCollection>(fSCName);
  fTokBSName         = iC.consumes<reco::BeamSpot>             (fBSName);
  fTokPFCandName     = iC.consumes<reco::PFCandidateCollection>(fPFCandName);
  if(fUseAOD) fTokPuppiName      = iC.consumes<reco::PFCandidateCollection>(fPuppiName);
  if(fUseAOD) fTokPuppiNoLepName = iC.consumes<reco::PFCandidateCollection>(fPuppiNoLepName);
  if(!fUseAOD) fTokPuppiPATName      = iC.consumes<pat::PackedCandidateCollection>(fPuppiName);
  if(!fUseAOD) fTokPuppiNoLepPATName = iC.consumes<pat::PackedCandidateCollection>(fPuppiNoLepName);
  fTokTrackName      = iC.consumes<reco::TrackCollection>      (fTrackName);
  fTokConvName       = iC.consumes<reco::ConversionCollection>(fConvName);
}

//--------------------------------------------------------------------------------------------------
FillerElectron::~FillerElectron(){}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerElectron::fill(TClonesArray *array,	    
	                  const edm::Event &iEvent, const edm::EventSetup &iSetup,      
	                  const reco::Vertex &pv,
			  const std::vector<TriggerRecord> &triggerRecords,
			  const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
  
  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> hEleProduct;
  iEvent.getByToken(fTokEleName,hEleProduct);
  assert(hEleProduct.isValid());
  const reco::GsfElectronCollection *eleCol = hEleProduct.product();

  // Get beamspot
  edm::Handle<reco::BeamSpot> hBeamSpotProduct;
  iEvent.getByToken(fTokBSName,hBeamSpotProduct);
  assert(hBeamSpotProduct.isValid());
  const reco::BeamSpot *bs = hBeamSpotProduct.product();

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByToken(fTokPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
  
  const reco::PFCandidateCollection *pfPuppi      = 0;
  const reco::PFCandidateCollection *pfPuppiNoLep = 0;
  if(fUsePuppi) { 
    // Get Puppi-candidates collection woof woof
    edm::Handle<reco::PFCandidateCollection> hPuppiProduct;
    iEvent.getByToken(fTokPuppiName,hPuppiProduct);
    assert(hPuppiProduct.isValid());
    pfPuppi = hPuppiProduct.product();
    
    // Get Puppi-no lep candidates collection arf arf
    edm::Handle<reco::PFCandidateCollection> hPuppiNoLepProduct;
    iEvent.getByToken(fTokPuppiNoLepName,hPuppiNoLepProduct);
    assert(hPuppiNoLepProduct.isValid());
    pfPuppiNoLep = hPuppiNoLepProduct.product();
  }
  // Get track collection
  edm::Handle<reco::TrackCollection> hTrackProduct;
  iEvent.getByToken(fTokTrackName,hTrackProduct);
  assert(hTrackProduct.isValid());
  const reco::TrackCollection *trackCol = hTrackProduct.product();
  
  // Get conversions collection
  edm::Handle<reco::ConversionCollection> hConvProduct;
  iEvent.getByToken(fTokConvName,hConvProduct);
  assert(hConvProduct.isValid());

  // Get SuperCluster collection
  edm::Handle<reco::SuperClusterCollection> hSCProduct;
  iEvent.getByToken(fTokSCName,hSCProduct);
  assert(hSCProduct.isValid());
  const reco::SuperClusterCollection *scCol = hSCProduct.product();

  // Track builder for computing 3D impact parameter    
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();


  for(reco::GsfElectronCollection::const_iterator itEle = eleCol->begin(); itEle!=eleCol->end(); ++itEle) {

    const reco::GsfTrackRef gsfTrack = itEle->gsfTrack();
    const reco::SuperClusterRef sc   = itEle->superCluster();

    // ref to access value maps
    edm::RefToBase<reco::GsfElectron> eleBaseRef( edm::Ref<reco::GsfElectronCollection>(hEleProduct, itEle - eleCol->begin()) );

    // electron pT cut
    if(itEle->pt() < fMinPt) continue;

    // construct object and place in array    
    TClonesArray &rElectronArr = *array;
    assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
    const int index = rElectronArr.GetEntries();  
    new(rElectronArr[index]) baconhep::TElectron();
    baconhep::TElectron *pElectron = (baconhep::TElectron*)rElectronArr[index];

    
    //
    // Kinematics
    //==============================    
    pElectron->pt         = itEle->pt();
    pElectron->eta        = itEle->eta();
    pElectron->phi        = itEle->phi();
    pElectron->q          = itEle->charge();
    pElectron->ecalEnergy = itEle->correctedEcalEnergy();
    pElectron->scEt       = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pElectron->scEta      = sc->eta();
    pElectron->scPhi      = sc->phi();

    pElectron->pfPt  = 0;
    pElectron->pfEta = 0;
    pElectron->pfPhi = 0;
    for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
      if( (itEle->closestTrack().isNonnull() && itPF->trackRef().isNonnull() && itEle->closestTrack() == itPF->trackRef()) ||
          (itPF->gsfTrackRef().isNonnull() && itPF->gsfTrackRef() == gsfTrack) ) {

	pElectron->pfPt  = itPF->pt();
	pElectron->pfEta = itPF->eta();
	pElectron->pfPhi = itPF->phi();  
      }
    }


    //
    // Isolation
    //==============================
    pElectron->trkIso        = itEle->dr03TkSumPt();
    pElectron->ecalIso       = itEle->dr03EcalRecHitSumEt();
    pElectron->hcalIso       = itEle->dr03HcalTowerSumEt();
    pElectron->hcalDepth1Iso = itEle->dr03HcalDepth1TowerSumEt();

    pElectron->chHadIso  = itEle->pfIsolationVariables().sumChargedHadronPt;
    pElectron->gammaIso  = itEle->pfIsolationVariables().sumPhotonEt;
    pElectron->neuHadIso = itEle->pfIsolationVariables().sumNeutralHadronEt;
    pElectron->puIso     = itEle->pfIsolationVariables().sumPUPt;

    if(fUsePuppi) { 
      double pEta = pElectron->pfEta;
      double pPhi = pElectron->pfPhi;
      if(pEta == 0) pEta = itEle->eta();
      if(pPhi == 0) pPhi = itEle->phi();
      computeIso(pEta,pPhi, 0.4, (*pfPuppi), 
		 pElectron->puppiChHadIso,
		 pElectron->puppiGammaIso,
		 pElectron->puppiNeuHadIso);
      
      computeIso(pEta,pPhi, 0.4, (*pfPuppiNoLep),
		 pElectron->puppiChHadIsoNoLep,
		 pElectron->puppiGammaIsoNoLep,
		 pElectron->puppiNeuHadIsoNoLep);
    }

    //
    // Impact Parameter
    //==============================
    if(gsfTrack.isNonnull()) { 
      pElectron->d0 = (-1)*(gsfTrack->dxy(pv.position()));  // note: d0 = -dxy
      pElectron->dz = gsfTrack->dz(pv.position());

      //(!) double check recipe
      const reco::TransientTrack &tt = transientTrackBuilder->build(gsfTrack);
      const double gsfsign = (pElectron->d0 >= 0) ? 1. : -1.;
      const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
      pElectron->sip3d = ip3d.first ? gsfsign*ip3d.second.value() / ip3d.second.error() : -999.;
    }


    //
    // Identification
    //==============================
    pElectron->sieie      = itEle->full5x5_sigmaIetaIeta();
    pElectron->e1x5       = itEle->full5x5_e1x5();
    pElectron->e2x5       = itEle->full5x5_e2x5Max();
    pElectron->e5x5       = itEle->full5x5_e5x5();
    pElectron->r9         = itEle->full5x5_r9();
    pElectron->hovere     = itEle->hcalOverEcal();
    pElectron->eoverp     = itEle->eSuperClusterOverP();
    pElectron->fbrem      = itEle->fbrem();
    pElectron->dEtaInSeed = dEtaInSeed(*itEle);
    pElectron->dEtaIn     = itEle->deltaEtaSuperClusterTrackAtVtx();
    pElectron->dPhiIn     = itEle->deltaPhiSuperClusterTrackAtVtx();
    
    //pElectron->mva = -999;
    pElectron->isConv = ConversionTools::hasMatchedConversion(*itEle, hConvProduct, bs->position(), true, 2.0, 1e-6, 0);
    
    if(gsfTrack.isNonnull()) {
      pElectron->nMissingHits = gsfTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    pElectron->typeBits=0;
    if(itEle->ecalDrivenSeed())    pElectron->typeBits |= baconhep::kEcalDriven;
    if(itEle->trackerDrivenSeed()) pElectron->typeBits |= baconhep::kTrackerDriven;
    
    pElectron->fiducialBits=0;
    if(itEle->isEB())        pElectron->fiducialBits |= kIsEB;
    if(itEle->isEE())        pElectron->fiducialBits |= kIsEE;
    if(itEle->isGap())       pElectron->fiducialBits |= kIsGap;
    if(itEle->isEBEEGap())   pElectron->fiducialBits |= kIsEBEEGap;
    if(itEle->isEBGap())     pElectron->fiducialBits |= kIsEBGap;
    if(itEle->isEBEtaGap())  pElectron->fiducialBits |= kIsEBEtaGap;
    if(itEle->isEBPhiGap())  pElectron->fiducialBits |= kIsEBPhiGap;
    if(itEle->isEEGap())     pElectron->fiducialBits |= kIsEEGap;
    if(itEle->isEEDeeGap())  pElectron->fiducialBits |= kIsEEDeeGap;
    if(itEle->isEERingGap()) pElectron->fiducialBits |= kIsEERingGap;

    pElectron->classification = itEle->classification();

    // Obtain a supercluster ID, unique per event. The SC ID is the index in the SC collection.
    pElectron->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator iS = scCol->begin(); iS!=scCol->end(); ++iS) {
      scIndex++;
      if(itEle->superCluster().get() == &(*iS)) {
        pElectron->scID = scIndex;
	break;
      }
    }

    // Obtain a track ID, unique per event. The track ID is the index in the general tracks collection
    pElectron->trkID = -1;
    if(itEle->closestTrack().isNonnull()) {
      int trkIndex = -1;    
      for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
        trkIndex++;
        if(itEle->closestTrack().get() == &(*itTrk)) {
          pElectron->trkID = trkIndex;
	  break;
        }
      }
    }
    
    if(fUseTO) pElectron->hltMatchBits = TriggerTools::matchHLT(pElectron->eta, pElectron->phi, triggerRecords, triggerEvent);
  }
}

// === filler for MINIAOD ===
void FillerElectron::fill(TClonesArray *array,
                          const edm::Event &iEvent, const edm::EventSetup &iSetup,
                          const reco::Vertex &pv,
                          const std::vector<TriggerRecord> &triggerRecords,
                          const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);

  // Get electron collection
  edm::Handle<pat::ElectronCollection> hEleProduct;
  iEvent.getByToken(fTokPatEleName,hEleProduct);
  assert(hEleProduct.isValid());
  const pat::ElectronCollection *eleCol = hEleProduct.product();

  // Get supercluster collection
  edm::Handle<reco::SuperClusterCollection> hSCProduct;
  iEvent.getByToken(fTokSCName,hSCProduct);
  assert(hSCProduct.isValid());
  const reco::SuperClusterCollection *scCol = hSCProduct.product();
  
  const pat::PackedCandidateCollection *pfPuppi      = 0;
  const pat::PackedCandidateCollection *pfPuppiNoLep = 0;
  if(fUsePuppi) { 
    // Get Puppi-candidates collection woof woof
    edm::Handle<pat::PackedCandidateCollection> hPuppiProduct;
    iEvent.getByToken(fTokPuppiPATName,hPuppiProduct);
    assert(hPuppiProduct.isValid());
    pfPuppi = hPuppiProduct.product();
    
    // Get Puppi-no lep candidates collection arf arf
    edm::Handle<pat::PackedCandidateCollection> hPuppiNoLepProduct;
    iEvent.getByToken(fTokPuppiNoLepPATName,hPuppiNoLepProduct);
    assert(hPuppiNoLepProduct.isValid());
    pfPuppiNoLep = hPuppiNoLepProduct.product();
  }
  for(pat::ElectronCollection::const_iterator itEle = eleCol->begin(); itEle!=eleCol->end(); ++itEle) {

    const reco::GsfTrackRef gsfTrack = itEle->gsfTrack();
    const reco::SuperClusterRef sc   = itEle->superCluster();
    edm::RefToBase<reco::GsfElectron> eleBaseRef( edm::Ref<pat::ElectronCollection>(hEleProduct, itEle - eleCol->begin()) );

    // electron pT cut
    if(itEle->pt() < fMinPt) continue;

    // construct object and place in array
    TClonesArray &rElectronArr = *array;
    assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
    const int index = rElectronArr.GetEntries();
    new(rElectronArr[index]) baconhep::TElectron();
    baconhep::TElectron *pElectron = (baconhep::TElectron*)rElectronArr[index];

    //
    // Kinematics
    //==============================
    pElectron->pt         = itEle->pt();
    pElectron->regscale   = itEle->ecalRegressionScale();
    pElectron->regsmear   = itEle->ecalRegressionSmear();
    pElectron->eta        = itEle->eta();
    pElectron->phi        = itEle->phi();
    pElectron->q          = itEle->charge();
    pElectron->ecalEnergy = itEle->ecalEnergy();
    pElectron->scEt       = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pElectron->scEta      = sc->eta();
    pElectron->scPhi      = sc->phi();

    auto corrP4  = itEle->p4() * itEle->userFloat("ecalTrkEnergyPostCorr") / itEle->energy(); // apply scaling/smearing
    pElectron->calibPt = corrP4.Pt();
    pElectron->calibE = corrP4.E();
    pElectron->calibEcalE = itEle->userFloat("ecalEnergyPostCorr");
    pElectron->calibPtErr = itEle->userFloat("ecalTrkEnergyErrPostCorr")*itEle->pt()/itEle->p();
    pElectron->energyScaleStatUp = itEle->userFloat("energyScaleStatUp");
    pElectron->energyScaleStatDown = itEle->userFloat("energyScaleStatDown");
    pElectron->energyScaleSystUp = itEle->userFloat("energyScaleSystUp");
    pElectron->energyScaleSystDown = itEle->userFloat("energyScaleSystDown");
    pElectron->energyScaleGainUp = itEle->userFloat("energyScaleGainUp");
    pElectron->energyScaleGainDown = itEle->userFloat("energyScaleGainDown");
    pElectron->energySigmaRhoUp = itEle->userFloat("energySigmaRhoUp");
    pElectron->energySigmaRhoDown = itEle->userFloat("energySigmaRhoDown");
    pElectron->energySigmaPhiUp = itEle->userFloat("energySigmaPhiUp");
    pElectron->energySigmaPhiDown = itEle->userFloat("energySigmaPhiDown");

    pElectron->pfPt  = 0;
    pElectron->pfEta = 0;
    pElectron->pfPhi = 0;
    double pfDR = 0.05;
    for(unsigned int ipf=0; ipf<itEle->numberOfSourceCandidatePtrs(); ipf++) {
      if(ipf==0 && itEle->pfCandidateRef().isNonnull()) continue;  // PF-candidate collection not in MINIAOD, asking for reference will crash
      const reco::CandidatePtr pfcand = itEle->sourceCandidatePtr(ipf);
      double dR = reco::deltaR(itEle->eta(), itEle->phi(), pfcand->eta(), pfcand->phi());
      if(dR < pfDR) {
        pfDR = dR;
        pElectron->pfPt  = pfcand->pt();
        pElectron->pfEta = pfcand->eta();
        pElectron->pfPhi = pfcand->phi();
      }
    }


    //
    // Isolation
    //==============================
    pElectron->trkIso        = itEle->dr03TkSumPt();
    pElectron->ecalIso       = itEle->dr03EcalRecHitSumEt();
    pElectron->hcalIso       = itEle->dr03HcalTowerSumEt();
    pElectron->hcalDepth1Iso = itEle->dr03HcalDepth1TowerSumEt();

    pElectron->chHadIso  = itEle->pfIsolationVariables().sumChargedHadronPt;
    pElectron->gammaIso  = itEle->pfIsolationVariables().sumPhotonEt;
    pElectron->neuHadIso = itEle->pfIsolationVariables().sumNeutralHadronEt;
    pElectron->puIso     = itEle->pfIsolationVariables().sumPUPt;

    pElectron->ecalPFClusIso = itEle->ecalPFClusterIso();
    pElectron->hcalPFClusIso = itEle->hcalPFClusterIso();

    if(fUsePuppi) { 
      double pEta = pElectron->pfEta;
      double pPhi = pElectron->pfPhi;
      if(pEta == 0) pEta = itEle->eta();
      if(pPhi == 0) pPhi = itEle->phi();
      computeIso(pEta,pPhi, 0.4, (*pfPuppi), 
		 pElectron->puppiChHadIso,
		 pElectron->puppiGammaIso,
		 pElectron->puppiNeuHadIso);
      
      computeIso(pEta,pPhi, 0.4, (*pfPuppiNoLep),
		 pElectron->puppiChHadIsoNoLep,
		 pElectron->puppiGammaIsoNoLep,
		 pElectron->puppiNeuHadIsoNoLep);
    }
    
    //
    // Impact Parameter
    //==============================
    if(gsfTrack.isNonnull()) {
      pElectron->d0    = (-1)*(gsfTrack->dxy(pv.position()));  // note: d0 = -dxy
      pElectron->dz    = gsfTrack->dz(pv.position());
      pElectron->sip3d = (itEle->edB(pat::Electron::PV3D) > 0) ? itEle->dB(pat::Electron::PV3D)/itEle->edB(pat::Electron::PV3D) : -999;
    }

    //
    // Identification
    //==============================
    pElectron->sieie      = itEle->full5x5_sigmaIetaIeta();
    pElectron->e1x5       = itEle->full5x5_e1x5();
    pElectron->e2x5       = itEle->full5x5_e2x5Max();
    pElectron->e5x5       = itEle->full5x5_e5x5();
    pElectron->r9         = itEle->full5x5_r9();
    pElectron->hovere     = itEle->hcalOverEcal();
    pElectron->eoverp     = itEle->eSuperClusterOverP();
    pElectron->fbrem      = itEle->fbrem();
    pElectron->dEtaInSeed = dEtaInSeed(*itEle);
    pElectron->dEtaIn     = itEle->deltaEtaSuperClusterTrackAtVtx();
    pElectron->dPhiIn     = itEle->deltaPhiSuperClusterTrackAtVtx();

    // saving v2 electron MVA IDs
    pElectron->mvaSpring16        = itEle->userFloat(fMVASpring16+"Values");
    pElectron->mvaFall17V1Iso     = itEle->userFloat(fMVAFall17V1Iso+"Values");
    pElectron->mvaFall17V1NoIso   = itEle->userFloat(fMVAFall17V1NoIso+"Values");
    pElectron->mvaFall17V2Iso     = itEle->userFloat(fMVAFall17V2Iso+"Values");
    pElectron->mvaFall17V2NoIso   = itEle->userFloat(fMVAFall17V2NoIso+"Values");
    pElectron->mvaSpring16HZZ     = itEle->userFloat(fMVASpring16HZZ+"Values");
    
    pElectron->mvaSpring16Cat        = itEle->userInt(fMVASpring16+"Categories");
    pElectron->mvaFall17V1IsoCat     = itEle->userInt(fMVAFall17V1Iso+"Categories");
    pElectron->mvaFall17V1NoIsoCat   = itEle->userInt(fMVAFall17V1NoIso+"Categories");
    pElectron->mvaFall17V2IsoCat     = itEle->userInt(fMVAFall17V2Iso+"Categories");
    pElectron->mvaFall17V2NoIsoCat   = itEle->userInt(fMVAFall17V2NoIso+"Categories");
    pElectron->mvaSpring16HZZCat     = itEle->userInt(fMVASpring16HZZ+"Categories");

    //pElectron->mvaV2IsoCat     = itEle->userInt(fMVAV2Iso+"Categories");
    //pElectron->mvaV2NoIsoCat   = itEle->userInt(fMVAV2NoIso+"Categories");

    pElectron->mvaBit     = 0; 
    if (itEle->electronID(fMediumMVA)) pElectron->mvaBit |= baconhep::kEleMVAMedBit;
    if (itEle->electronID(fTightMVA)) pElectron->mvaBit |= baconhep::kEleMVATightBit;
    //pElectron->mva = itEle->userFloat(fMVA+"Values");
    //pElectron->mvaCat = itEle->userInt(fMVA+"Categories");

    /*if (fSecondMVA) {
      pElectron->mvaIsoBit     = 0; 
      if (itEle->electronID(fMediumMVAIso)) pElectron->mvaIsoBit |= baconhep::kEleMVAMedBit;
      if (itEle->electronID(fTightMVAIso)) pElectron->mvaIsoBit |= baconhep::kEleMVATightBit;
      pElectron->mvaIso = itEle->userFloat(fMVAIso+"Values");
      pElectron->mvaIsoCat = itEle->userInt(fMVAIso+"Categories");
    }

    if (fStoreHZZMVA) {
        pElectron->mvaHZZ = itEle->userFloat(fMVAHZZ+"Values");
        pElectron->mvaHZZCat = itEle->userInt(fMVAHZZ+"Categories");
    }*/

    pElectron->isConv     = !itEle->passConversionVeto();

    if(gsfTrack.isNonnull()) {
      pElectron->nMissingHits = gsfTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    pElectron->typeBits=0;
    if(itEle->ecalDrivenSeed())    pElectron->typeBits |= baconhep::kEcalDriven;
    if(itEle->trackerDrivenSeed()) pElectron->typeBits |= baconhep::kTrackerDriven;

    pElectron->fiducialBits=0;
    if(itEle->isEB())        pElectron->fiducialBits |= kIsEB;
    if(itEle->isEE())        pElectron->fiducialBits |= kIsEE;
    if(itEle->isGap())       pElectron->fiducialBits |= kIsGap;
    if(itEle->isEBEEGap())   pElectron->fiducialBits |= kIsEBEEGap;
    if(itEle->isEBGap())     pElectron->fiducialBits |= kIsEBGap;
    if(itEle->isEBEtaGap())  pElectron->fiducialBits |= kIsEBEtaGap;
    if(itEle->isEBPhiGap())  pElectron->fiducialBits |= kIsEBPhiGap;
    if(itEle->isEEGap())     pElectron->fiducialBits |= kIsEEGap;
    if(itEle->isEEDeeGap())  pElectron->fiducialBits |= kIsEEDeeGap;
    if(itEle->isEERingGap()) pElectron->fiducialBits |= kIsEERingGap;

    pElectron->classification = itEle->classification();

    // Obtain a supercluster ID, unique per event. The SC ID is the index in the SC collection.
    pElectron->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator iS = scCol->begin(); iS!=scCol->end(); ++iS) {
      scIndex++;
      if(itEle->superCluster().get() == &(*iS)) {
        pElectron->scID = scIndex;
        break;
      }
    }

    // Obtain a track ID, unique per event. The track ID is the index in the general tracks collection
    pElectron->trkID = -1;  // general tracks not in MINIAOD
    if(fUseTO) pElectron->hltMatchBits = TriggerTools::matchHLT(pElectron->eta, pElectron->phi, triggerRecords, triggerObjects);
  }
}

//--------------------------------------------------------------------------------------------------
double FillerElectron::dEtaInSeed(const reco::GsfElectron& ele) {
  return (ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull()) ? 
         ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : 
         std::numeric_limits<float>::max();
}

double FillerElectron::dEtaInSeed(const pat::Electron& ele) {
  return (ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull()) ?
         ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() :
         std::numeric_limits<float>::max();
}
void FillerElectron::computeIso(double &iEta,double &iPhi, const double extRadius,
				const reco::PFCandidateCollection    &puppi,
				float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{
  // Muon PF isolation with delta-beta PU correction:
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
  
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double ptMin           = 0.1;  
  //const double intRadiusChHad  = 0.0001;
  //const double intRadiusGamma  = 0.01;
  //const double intRadiusNeuHad = 0.01;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<puppi.size(); ipf++) {
    const reco::PFCandidate pfcand = puppi.at(ipf);    
    bool pPass = true;
    double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), iEta, iPhi);
    if(dr < 0.0001) pPass = false; //Use this to avoid float/double bullshit
    if(pPass) { 
      if     (pfcand.particleId() == reco::PFCandidate::h)     { intRadius = 0;}//intRadiusChHad; }
      else if(pfcand.particleId() == reco::PFCandidate::gamma) { intRadius = 0;}//intRadiusGamma;  }
      else if(pfcand.particleId() == reco::PFCandidate::h0)    { intRadius = 0;}//intRadiusNeuHad; }
            
      if(dr>=extRadius || dr<intRadius) continue;
            
      if     (pfcand.particleId() == reco::PFCandidate::h)                             { chHadIso  += pfcand.pt(); }
      else if(pfcand.particleId() == reco::PFCandidate::gamma && pfcand.pt() > ptMin) { gammaIso  += pfcand.pt(); }
      else if(pfcand.particleId() == reco::PFCandidate::h0    && pfcand.pt() > ptMin) { neuHadIso += pfcand.pt(); }
    }
  }
  // compute PU iso
  out_chHadIso  = chHadIso;
  out_gammaIso  = gammaIso;
  out_neuHadIso = neuHadIso;
}
void FillerElectron::computeIso(double &iEta,double &iPhi, const double extRadius,
			    const pat::PackedCandidateCollection    &puppi,
                            float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{
  // Muon PF isolation with delta-beta PU correction:
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
  
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double ptMin           = 0.1;  
  //const double intRadiusChHad  = 0.0001;
  //const double intRadiusGamma  = 0.01;
  //const double intRadiusNeuHad = 0.01;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<puppi.size(); ipf++) {
    const pat::PackedCandidate pfcand = puppi.at(ipf);    
    bool pPass = true;
    double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), iEta, iPhi);
    if(dr < 0.0001) pPass = false; //Use this to avoid float/double bullshit
    if(pPass) { 
      if     (abs(pfcand.pdgId()) == 211)     { intRadius = 0;}//intRadiusChHad; }
      else if(abs(pfcand.pdgId()) == 22)      { intRadius = 0;}//intRadiusGamma;  }
      else if(abs(pfcand.pdgId()) == 130)     { intRadius = 0;}//intRadiusNeuHad; }
            
      if(dr>=extRadius || dr<intRadius) continue;
      if     (abs(pfcand.pdgId()) == 211)                        { chHadIso  += pfcand.pt(); }
      else if(abs(pfcand.pdgId()) == 22  && pfcand.pt() > ptMin) { gammaIso  += pfcand.pt(); }
      else if(abs(pfcand.pdgId()) == 130 && pfcand.pt() > ptMin) { neuHadIso += pfcand.pt(); }      
    }
  }
  // compute PU iso
  out_chHadIso  = chHadIso;
  out_gammaIso  = gammaIso;
  out_neuHadIso = neuHadIso;
}
