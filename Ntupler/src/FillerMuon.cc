#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerMuon::FillerMuon(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt         (iConfig.getUntrackedParameter<double>("minPt",0)),
  fMuonName      (iConfig.getUntrackedParameter<std::string>("edmName","muons")),
  fPFCandName    (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName     (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fSaveTracks    (iConfig.getUntrackedParameter<bool>("doSaveTracks",false)),
  fTrackMinPt    (iConfig.getUntrackedParameter<double>("minTrackPt",20)),
  fPuppiName     (iConfig.getUntrackedParameter<std::string>("edmPuppiName","puppi")),
  fPuppiNoLepName(iConfig.getUntrackedParameter<std::string>("edmPuppiNoLepName","puppiNoLep")),
  fUsePuppi      (iConfig.getUntrackedParameter<bool>("usePuppi",true)),
  fUseTO         (iConfig.getUntrackedParameter<bool>("useTriggerObject",false)),
  fUseAOD        (useAOD)
{
  if(fUseAOD)  fTokMuonName       = iC.consumes<reco::MuonCollection>       (fMuonName);
  if(!fUseAOD) fTokPatMuonName    = iC.consumes<pat::MuonCollection>        (fMuonName);
  fTokPFCandName     = iC.consumes<reco::PFCandidateCollection>(fPFCandName);
  if(fUseAOD)  fTokPuppiName      = iC.consumes<reco::PFCandidateCollection>(fPuppiName);
  if(fUseAOD)  fTokPuppiNoLepName = iC.consumes<reco::PFCandidateCollection>(fPuppiNoLepName);
  if(!fUseAOD) fTokPuppiPATName      = iC.consumes<pat::PackedCandidateCollection>(fPuppiName);
  if(!fUseAOD) fTokPuppiNoLepPATName = iC.consumes<pat::PackedCandidateCollection>(fPuppiNoLepName);
  fTokTrackName      = iC.consumes<reco::TrackCollection>      (fTrackName);
}

//--------------------------------------------------------------------------------------------------
FillerMuon::~FillerMuon(){}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerMuon::fill(TClonesArray *array,
                      const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
		      const std::vector<TriggerRecord> &triggerRecords,
		      const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
  
  // Get muon collection
  edm::Handle<reco::MuonCollection> hMuonProduct;
  iEvent.getByToken(fTokMuonName,hMuonProduct);
  assert(hMuonProduct.isValid());
  const reco::MuonCollection *muonCol = hMuonProduct.product();
  
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
  
  // Track builder for computing 3D impact parameter
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();  

  
  for(reco::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {    
    
    // muon pT cut
    if(itMu->pt() < fMinPt) continue;
    
    // construct object and place in array
    TClonesArray &rArray = *array;
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TMuon();
    baconhep::TMuon *pMuon = (baconhep::TMuon*)rArray[index];

    
    //
    // Kinematics
    //==============================
    pMuon->pt     = itMu->muonBestTrack()->pt();
    pMuon->eta    = itMu->muonBestTrack()->eta();
    pMuon->phi    = itMu->muonBestTrack()->phi();
    pMuon->ptErr  = itMu->muonBestTrack()->ptError();
    pMuon->q      = itMu->muonBestTrack()->charge();
    pMuon->btt    = itMu->muonBestTrackType();
    pMuon->staPt  = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->pt()  : 0;
    pMuon->staEta = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->eta() : 0;
    pMuon->staPhi = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->phi() : 0;   
    
    pMuon->pfPt  = itMu->pfP4().pt();
    pMuon->pfEta = itMu->pfP4().eta();
    pMuon->pfPhi = itMu->pfP4().phi();

    
    //
    // Isolation
    //==============================
    pMuon->trkIso  = itMu->isolationR03().sumPt;
    pMuon->ecalIso = itMu->isolationR03().emEt;
    pMuon->hcalIso = itMu->isolationR03().hadEt;

    pMuon->chHadIso  = itMu->pfIsolationR04().sumChargedHadronPt;
    pMuon->gammaIso  = itMu->pfIsolationR04().sumPhotonEt;
    pMuon->neuHadIso = itMu->pfIsolationR04().sumNeutralHadronEt;
    pMuon->puIso     = itMu->pfIsolationR04().sumPUPt;

    pMuon->chHadIso03  = itMu->pfIsolationR03().sumChargedHadronPt;
    pMuon->gammaIso03  = itMu->pfIsolationR03().sumPhotonEt;
    pMuon->neuHadIso03 = itMu->pfIsolationR03().sumNeutralHadronEt;
    pMuon->puIso03     = itMu->pfIsolationR03().sumPUPt;
    
    if(fUsePuppi) { 
      double pEta = pMuon->pfEta;
      double pPhi = pMuon->pfPhi;
      if(pEta == 0) pEta = itMu->muonBestTrack()->eta();
      if(pPhi == 0) pPhi = itMu->muonBestTrack()->phi();
      computeIso(pEta,pPhi, 0.4, (*pfPuppi), 
		 pMuon->puppiChHadIso,
		 pMuon->puppiGammaIso,
		 pMuon->puppiNeuHadIso);
      
      computeIso(pEta,pPhi, 0.4, (*pfPuppiNoLep),
		 pMuon->puppiChHadIsoNoLep,
		 pMuon->puppiGammaIsoNoLep,
		 pMuon->puppiNeuHadIsoNoLep);
    }
    //
    // Impact Parameter
    //==============================
    pMuon->d0 = (-1)*(itMu->muonBestTrack()->dxy(pv.position()));  // note: d0 = -dxy
    pMuon->dz = itMu->muonBestTrack()->dz(pv.position());
    
    const reco::TransientTrack &tt = transientTrackBuilder->build(itMu->muonBestTrack());
    const double thesign = (pMuon->d0 >= 0) ? 1. : -1.;
    const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
    pMuon->sip3d = ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.;

    
    //
    // Identification
    //==============================
    pMuon->tkNchi2    = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->normalizedChi2()  : -999.;
    pMuon->muNchi2    = itMu->isGlobalMuon()           ? itMu->globalTrack()->normalizedChi2() : -999.;    
    pMuon->trkKink    = itMu->combinedQuality().trkKink;
    pMuon->glbKink    = itMu->combinedQuality().glbKink;
    pMuon->trkHitFrac = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->validFraction() : 0;
    pMuon->chi2LocPos = itMu->combinedQuality().chi2LocalPosition;
    pMuon->segComp    = muon::segmentCompatibility(*itMu);
    pMuon->caloComp   = muon::caloCompatibility(*itMu);
    pMuon->typeBits   = itMu->type();
        
    pMuon->selectorBits=0;
    if(muon::isGoodMuon(*itMu,muon::All))                                    pMuon->selectorBits |= baconhep::kAll;
    if(muon::isGoodMuon(*itMu,muon::AllGlobalMuons))                         pMuon->selectorBits |= baconhep::kAllGlobalMuons;
    if(muon::isGoodMuon(*itMu,muon::AllStandAloneMuons))                     pMuon->selectorBits |= baconhep::kAllStandAloneMuons;
    if(muon::isGoodMuon(*itMu,muon::AllTrackerMuons))                        pMuon->selectorBits |= baconhep::kAllTrackerMuons;
    if(muon::isGoodMuon(*itMu,muon::TrackerMuonArbitrated))                  pMuon->selectorBits |= baconhep::kTrackerMuonArbitrated;
    if(muon::isGoodMuon(*itMu,muon::AllArbitrated))                          pMuon->selectorBits |= baconhep::kAllArbitrated;
    if(muon::isGoodMuon(*itMu,muon::GlobalMuonPromptTight))                  pMuon->selectorBits |= baconhep::kGlobalMuonPromptTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationLoose))                     pMuon->selectorBits |= baconhep::kTMLastStationLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationTight))                     pMuon->selectorBits |= baconhep::kTMLastStationTight;
    if(muon::isGoodMuon(*itMu,muon::TM2DCompatibilityLoose))                 pMuon->selectorBits |= baconhep::kTM2DCompatibilityLoose;
    if(muon::isGoodMuon(*itMu,muon::TM2DCompatibilityTight))                 pMuon->selectorBits |= baconhep::kTM2DCompatibilityTight;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationLoose))                      pMuon->selectorBits |= baconhep::kTMOneStationLoose;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationTight))                      pMuon->selectorBits |= baconhep::kTMOneStationTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedLowPtLoose))       pMuon->selectorBits |= baconhep::kTMLastStationOptimizedLowPtLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedLowPtTight))       pMuon->selectorBits |= baconhep::kTMLastStationOptimizedLowPtTight;
    if(muon::isGoodMuon(*itMu,muon::GMTkChiCompatibility))                   pMuon->selectorBits |= baconhep::kGMTkChiCompatibility;
    if(muon::isGoodMuon(*itMu,muon::GMStaChiCompatibility))                  pMuon->selectorBits |= baconhep::kGMStaChiCompatibility;
    if(muon::isGoodMuon(*itMu,muon::GMTkKinkTight))                          pMuon->selectorBits |= baconhep::kGMTkKinkTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationAngLoose))                  pMuon->selectorBits |= baconhep::kTMLastStationAngLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationAngTight))                  pMuon->selectorBits |= baconhep::kTMLastStationAngTight;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationAngLoose))                   pMuon->selectorBits |= baconhep::kTMOneStationAngLoose;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationAngTight))                   pMuon->selectorBits |= baconhep::kTMOneStationAngTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedBarrelLowPtLoose)) pMuon->selectorBits |= baconhep::kTMLastStationOptimizedBarrelLowPtLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedBarrelLowPtTight)) pMuon->selectorBits |= baconhep::kTMLastStationOptimizedBarrelLowPtTight;
    if(muon::isGoodMuon(*itMu,muon::RPCMuLoose))                             pMuon->selectorBits |= baconhep::kRPCMuLoose;

    pMuon->pogIDBits=0;
    if(muon::isLooseMuon(*itMu))      pMuon->pogIDBits |= baconhep::kPOGLooseMuon;
    if(muon::isMediumMuon(*itMu))     pMuon->pogIDBits |= baconhep::kPOGMediumMuon;
    if(muon::isTightMuon(*itMu, pv))  pMuon->pogIDBits |= baconhep::kPOGTightMuon;
    if(muon::isSoftMuon(*itMu, pv))   pMuon->pogIDBits |= baconhep::kPOGSoftMuon;
    if(muon::isHighPtMuon(*itMu, pv)) pMuon->pogIDBits |= baconhep::kPOGHighPtMuon;

    pMuon->nValidHits = itMu->isGlobalMuon()           ? itMu->globalTrack()->hitPattern().numberOfValidMuonHits()       : 0;        
    pMuon->nTkHits    = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().numberOfValidTrackerHits()     : 0;
    pMuon->nPixHits   = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().numberOfValidPixelHits()       : 0;
    pMuon->nTkLayers  = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0;
    pMuon->nPixLayers = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().pixelLayersWithMeasurement()   : 0;
    pMuon->nMatchStn  = itMu->numberOfMatchedStations();

    // Obtain a track ID, unique per event. The track ID is the index in the general tracks collection    
    pMuon->trkID = -1;
    if(itMu->innerTrack().isNonnull()) {
      int trkIndex = -1;
      for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
        trkIndex++;
	if(itMu->innerTrack().get() == &(*itTrk)) {
	  pMuon->trkID = index;
	  break;
	}
      }
    }
    
    if(fUseTO) pMuon->hltMatchBits = TriggerTools::matchHLT(pMuon->eta, pMuon->phi, triggerRecords, triggerEvent);
  }

  //
  // Save tracks in Bacon muon array (optional)
  // * AOD only
  // * Useful for tag-and-probe for standalone muon efficiency (e.g. isolated track probes)
  //
  if(fSaveTracks) {    
    int trkIndex = -1;
    for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
      trkIndex++;
      
      // check track is not a muon
      bool isMuon = false;
      for(reco::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {
        if(itMu->innerTrack().isNonnull() && itMu->innerTrack().get() == &(*itTrk)) {
	  isMuon = true;
	  break;
	}
      }
      if(isMuon) continue;

      // track pT cut
      if(itTrk->pt() < fTrackMinPt) continue;    
      
      
      TClonesArray &rArray = *array;
      assert(rArray.GetEntries() < rArray.GetSize());
      const int index = rArray.GetEntries();
      new(rArray[index]) baconhep::TMuon();
      baconhep::TMuon *pMuon = (baconhep::TMuon*)rArray[index];

      //
      // Kinematics
      //==============================
      pMuon->pt      = itTrk->pt();
      pMuon->eta     = itTrk->eta();
      pMuon->phi     = itTrk->phi();
      pMuon->ptErr   = itTrk->ptError();
      pMuon->q       = itTrk->charge();
      pMuon->staPt   = 0;
      pMuon->staEta  = 0;
      pMuon->staPhi  = 0;   
    
      pMuon->pfPt  = 0;
      pMuon->pfEta = 0;
      pMuon->pfPhi = 0;
      for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
        if(itPF->trackRef().isNonnull() && &(*itTrk) == itPF->trackRef().get()) {
	  pMuon->pfPt  = itPF->pt();
	  pMuon->pfEta = itPF->eta();
	  pMuon->pfPhi = itPF->phi();
        }
      }

    
      //
      // Isolation
      //==============================
      pMuon->trkIso  = -1;
      pMuon->ecalIso = -1;
      pMuon->hcalIso = -1;
/* (!) need to manually compute isolation
      computeIso(*itTrk, 0.3, pfNoPU, pfPU,
                 pMuon->chHadIso,
                 pMuon->gammaIso,
                 pMuon->neuHadIso,
                 pMuon->puIso);
*/    
            
      //
      // Impact Parameter
      //==============================
      pMuon->d0 = (-1)*(itTrk->dxy(pv.position()));  // note: d0 = -dxy
      pMuon->dz =  itTrk->dz(pv.position());

      const reco::TransientTrack &tt = transientTrackBuilder->build(&(*itTrk));
      const double thesign = (pMuon->d0 >= 0) ? 1. : -1.;
      const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
      pMuon->sip3d = ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.;      


      //
      // Identification
      //==============================
      pMuon->tkNchi2      = itTrk->normalizedChi2();
      pMuon->muNchi2      = -999.;    
      pMuon->trkKink      = 0;
      pMuon->glbKink      = 0;        
      pMuon->typeBits     = 0;
      pMuon->selectorBits = 0;
      pMuon->nValidHits   = 0;        
      pMuon->nTkHits      = itTrk->hitPattern().numberOfValidTrackerHits();
      pMuon->nPixHits     = itTrk->hitPattern().numberOfValidPixelHits();
      pMuon->nTkLayers    = itTrk->hitPattern().trackerLayersWithMeasurement();
      pMuon->nPixLayers   = itTrk->hitPattern().pixelLayersWithMeasurement();
      pMuon->nMatchStn    = 0;
      pMuon->trkID        = trkIndex;
      if(fUseTO) pMuon->hltMatchBits = TriggerTools::matchHLT(pMuon->eta, pMuon->phi, triggerRecords, triggerEvent);
    }    
  } 
}

// === filler for MINIAOD ===
void FillerMuon::fill(TClonesArray *array,
                      const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv,
                      const std::vector<TriggerRecord> &triggerRecords,
                      const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);

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
  // Get muon collection
  edm::Handle<pat::MuonCollection> hMuonProduct;
  iEvent.getByToken(fTokPatMuonName,hMuonProduct);
  assert(hMuonProduct.isValid());
  const pat::MuonCollection *muonCol = hMuonProduct.product();

  for(pat::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {

    // muon pT cut
    if(itMu->pt() < fMinPt) continue;

    // construct object and place in array
    TClonesArray &rArray = *array;
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TMuon();
    baconhep::TMuon *pMuon = (baconhep::TMuon*)rArray[index];


    //
    // Kinematics
    //==============================
    pMuon->pt     = itMu->muonBestTrack()->pt();
    pMuon->eta    = itMu->muonBestTrack()->eta();
    pMuon->phi    = itMu->muonBestTrack()->phi();
    pMuon->ptErr  = itMu->muonBestTrack()->ptError();
    pMuon->q      = itMu->muonBestTrack()->charge();
    pMuon->btt    = itMu->muonBestTrackType();
    pMuon->staPt  = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->pt()  : 0;
    pMuon->staEta = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->eta() : 0;
    pMuon->staPhi = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->phi() : 0;

    pMuon->pfPt  = itMu->pfP4().pt();
    pMuon->pfEta = itMu->pfP4().eta();
    pMuon->pfPhi = itMu->pfP4().phi();

    //
    // Isolation
    //==============================
    pMuon->trkIso  = itMu->isolationR03().sumPt;
    pMuon->ecalIso = itMu->isolationR03().emEt;
    pMuon->hcalIso = itMu->isolationR03().hadEt;

    pMuon->chHadIso  = itMu->pfIsolationR04().sumChargedHadronPt;
    pMuon->gammaIso  = itMu->pfIsolationR04().sumPhotonEt;
    pMuon->neuHadIso = itMu->pfIsolationR04().sumNeutralHadronEt;
    pMuon->puIso     = itMu->pfIsolationR04().sumPUPt;

    pMuon->chHadIso03  = itMu->pfIsolationR03().sumChargedHadronPt;
    pMuon->gammaIso03  = itMu->pfIsolationR03().sumPhotonEt;
    pMuon->neuHadIso03 = itMu->pfIsolationR03().sumNeutralHadronEt;
    pMuon->puIso03     = itMu->pfIsolationR03().sumPUPt;

    if(fUsePuppi) { 
      double pEta = pMuon->pfEta;
      double pPhi = pMuon->pfPhi;
      if(pEta == 0) pEta = itMu->muonBestTrack()->eta();
      if(pPhi == 0) pPhi = itMu->muonBestTrack()->phi();
      computeIso(pEta,pPhi, 0.4, (*pfPuppi), 
		 pMuon->puppiChHadIso,
		 pMuon->puppiGammaIso,
		 pMuon->puppiNeuHadIso);
      
      computeIso(pEta,pPhi, 0.4, (*pfPuppiNoLep),
		 pMuon->puppiChHadIsoNoLep,
		 pMuon->puppiGammaIsoNoLep,
		 pMuon->puppiNeuHadIsoNoLep);
    }

    //
    // Impact Parameter
    //==============================
    pMuon->d0    = (-1)*(itMu->muonBestTrack()->dxy(pv.position()));  // note: d0 = -dxy
    pMuon->dz    = itMu->muonBestTrack()->dz(pv.position());
    pMuon->sip3d = (itMu->edB(pat::Muon::PV3D) > 0) ? itMu->dB(pat::Muon::PV3D)/itMu->edB(pat::Muon::PV3D) : -999;


    //
    // Identification
    //==============================
    pMuon->tkNchi2    = itMu->innerTrack().isNonnull()  ? itMu->innerTrack()->normalizedChi2()  : -999.;
    pMuon->muNchi2    = itMu->globalTrack().isNonnull() ? itMu->globalTrack()->normalizedChi2() : -999.;
    pMuon->trkKink    = itMu->combinedQuality().trkKink;
    pMuon->glbKink    = itMu->combinedQuality().glbKink;
    pMuon->trkHitFrac = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->validFraction() : 0;
    pMuon->chi2LocPos = itMu->combinedQuality().chi2LocalPosition;
    pMuon->segComp    = muon::segmentCompatibility(*itMu);
    pMuon->caloComp   = muon::caloCompatibility(*itMu);
    pMuon->typeBits   = itMu->type();

    pMuon->selectorBits=0;
    if(itMu->muonID("All"))                                    pMuon->selectorBits |= baconhep::kAll;
    if(itMu->muonID("AllGlobalMuons"))                         pMuon->selectorBits |= baconhep::kAllGlobalMuons;
    if(itMu->muonID("AllStandAloneMuons"))                     pMuon->selectorBits |= baconhep::kAllStandAloneMuons;
    if(itMu->muonID("AllTrackerMuons"))                        pMuon->selectorBits |= baconhep::kAllTrackerMuons;
    if(itMu->muonID("TrackerMuonArbitrated"))                  pMuon->selectorBits |= baconhep::kTrackerMuonArbitrated;
    if(itMu->muonID("AllArbitrated"))                          pMuon->selectorBits |= baconhep::kAllArbitrated;
    if(itMu->muonID("GlobalMuonPromptTight"))                  pMuon->selectorBits |= baconhep::kGlobalMuonPromptTight;
    if(itMu->muonID("TMLastStationLoose"))                     pMuon->selectorBits |= baconhep::kTMLastStationLoose;
    if(itMu->muonID("TMLastStationTight"))                     pMuon->selectorBits |= baconhep::kTMLastStationTight;
    if(itMu->muonID("TM2DCompatibilityLoose"))                 pMuon->selectorBits |= baconhep::kTM2DCompatibilityLoose;
    if(itMu->muonID("TM2DCompatibilityTight"))                 pMuon->selectorBits |= baconhep::kTM2DCompatibilityTight;
    if(itMu->muonID("TMOneStationLoose"))                      pMuon->selectorBits |= baconhep::kTMOneStationLoose;
    if(itMu->muonID("TMOneStationTight"))                      pMuon->selectorBits |= baconhep::kTMOneStationTight;
    if(itMu->muonID("TMLastStationOptimizedLowPtLoose"))       pMuon->selectorBits |= baconhep::kTMLastStationOptimizedLowPtLoose;
    if(itMu->muonID("TMLastStationOptimizedLowPtTight"))       pMuon->selectorBits |= baconhep::kTMLastStationOptimizedLowPtTight;
    if(itMu->muonID("GMTkChiCompatibility"))                   pMuon->selectorBits |= baconhep::kGMTkChiCompatibility;
    if(itMu->muonID("GMStaChiCompatibility"))                  pMuon->selectorBits |= baconhep::kGMStaChiCompatibility;
    if(itMu->muonID("GMTkKinkTight"))                          pMuon->selectorBits |= baconhep::kGMTkKinkTight;
    if(itMu->muonID("TMLastStationAngLoose"))                  pMuon->selectorBits |= baconhep::kTMLastStationAngLoose;
    if(itMu->muonID("TMLastStationAngTight"))                  pMuon->selectorBits |= baconhep::kTMLastStationAngTight;
    if(itMu->muonID("TMOneStationAngLoose"))                   pMuon->selectorBits |= baconhep::kTMOneStationAngLoose;
    if(itMu->muonID("TMOneStationAngTight"))                   pMuon->selectorBits |= baconhep::kTMOneStationAngTight;
    if(itMu->muonID("TMLastStationOptimizedBarrelLowPtLoose")) pMuon->selectorBits |= baconhep::kTMLastStationOptimizedBarrelLowPtLoose;
    if(itMu->muonID("TMLastStationOptimizedBarrelLowPtTight")) pMuon->selectorBits |= baconhep::kTMLastStationOptimizedBarrelLowPtTight;
    if(itMu->muonID("RPCMuLoose"))                             pMuon->selectorBits |= baconhep::kRPCMuLoose;

    pMuon->pogIDBits=0;
    if(itMu->isLooseMuon())    pMuon->pogIDBits |= baconhep::kPOGLooseMuon;
    if(itMu->isMediumMuon())   pMuon->pogIDBits |= baconhep::kPOGMediumMuon;
    if(itMu->isTightMuon(pv))  pMuon->pogIDBits |= baconhep::kPOGTightMuon;
    if(itMu->isSoftMuon(pv))   pMuon->pogIDBits |= baconhep::kPOGSoftMuon;
    if(itMu->isHighPtMuon(pv)) pMuon->pogIDBits |= baconhep::kPOGHighPtMuon;

    pMuon->nValidHits = itMu->globalTrack().isNonnull() ? itMu->globalTrack()->hitPattern().numberOfValidMuonHits()       : 0;
    pMuon->nTkHits    = itMu->innerTrack().isNonnull()  ? itMu->innerTrack()->hitPattern().numberOfValidTrackerHits()     : 0;
    pMuon->nPixHits   = itMu->innerTrack().isNonnull()  ? itMu->innerTrack()->hitPattern().numberOfValidPixelHits()       : 0;
    pMuon->nTkLayers  = itMu->innerTrack().isNonnull()  ? itMu->innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0;
    pMuon->nPixLayers = itMu->innerTrack().isNonnull()  ? itMu->innerTrack()->hitPattern().pixelLayersWithMeasurement()   : 0;
    pMuon->nMatchStn  = itMu->numberOfMatchedStations();

    // Obtain a track ID, unique per event. The track ID is the index in the general tracks collection
    pMuon->trkID = -1;  // general tracks not in MINIAOD

    if(fUseTO) pMuon->hltMatchBits = TriggerTools::matchHLT(pMuon->eta, pMuon->phi, triggerRecords, triggerObjects);
  }
}
void FillerMuon::computeIso(double &iEta,double &iPhi, const double extRadius,
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
void FillerMuon::computeIso(double &iEta,double &iPhi, const double extRadius,
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
