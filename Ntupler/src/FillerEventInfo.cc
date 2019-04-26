#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/Handle.h"
#include <TLorentzVector.h>
#include <string>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerEventInfo::FillerEventInfo(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fPFCandName    (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fPUInfoName    (iConfig.getUntrackedParameter<std::string>("edmPileupInfoName","addPileupInfo")),
  fPVName        (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fBSName        (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot")),
  fCaloMETName   (iConfig.getUntrackedParameter<std::string>("edmCaloMETName","caloMet")),
  fMETName       (iConfig.getUntrackedParameter<std::string>("edmMETName","slimmedMETs")),
  fPFMETName     (iConfig.getUntrackedParameter<edm::InputTag>("edmPFMETName")),
  fPFMETCName    (iConfig.getUntrackedParameter<std::string>("edmPFMETCorrName","pfType1CorrectedMet")),
  fMVAMETName    (iConfig.getUntrackedParameter<std::string>("edmMVAMETName","pfMEtMVA")),
  fPUPPETName    (iConfig.getUntrackedParameter<edm::InputTag>("edmPuppETName")),//fPUPPETName    (iConfig.getUntrackedParameter<std::string>("edmPuppETName")),
  fPUPPETCName   (iConfig.getUntrackedParameter<std::string>("edmPuppETCorrName","pfType1CorrectedMetPuppi")),
  fALPACAMETName (iConfig.getUntrackedParameter<std::string>("edmAlpacaMETName"    ,"pfMetAlpacaMC")),
  fPALPACAMETName(iConfig.getUntrackedParameter<std::string>("edmPupAlpacaMETName","pfMetPuppiAlpacaMC")),
  fRhoIsoName    (iConfig.getUntrackedParameter<std::string>("edmRhoForIsoName","fixedGridRhoFastjetAll")),
  fRhoJetName    (iConfig.getUntrackedParameter<std::string>("edmRhoForJetEnergy","fixedGridRhoFastjetAll")),
  fUseFilters    (iConfig.getUntrackedParameter<bool>("doFillMETFilters",true)),
  fUseAOD        (useAOD)
{
  fTokPUInfoName     = iC.consumes< std::vector<PileupSummaryInfo> > (fPUInfoName);
  fTokBSName         = iC.consumes<reco::BeamSpot>                   (fBSName);
  fTokCaloMETName    = iC.consumes<reco::CaloMETCollection>          (fCaloMETName);
  fTokCaloMETPATName = iC.consumes<pat::METCollection>               (fMETName);
  fTokPFMETName      = iC.consumes<reco::PFMETCollection>            (fPFMETName);
  fTokPFMETPATName   = iC.consumes<pat::METCollection>               (fPFMETName);
  fTokPFMETCName     = iC.consumes<reco::PFMETCollection>            (fPFMETCName);
  fTokMVAMETName     = iC.consumes<reco::PFMETCollection>            (fMVAMETName);
  fTokPUPPETName     = iC.consumes<reco::PFMETCollection>            (fPUPPETName);
  fTokPUPPETCName    = iC.consumes<reco::PFMETCollection>            (fPUPPETCName);
  fTokPUPPETPATName  = iC.consumes<pat::METCollection>               (fPUPPETName);
  fTokALPACAMETName  = iC.consumes<reco::PFMETCollection>            (fALPACAMETName);
  fTokPALPACAMETName = iC.consumes<reco::PFMETCollection>            (fPALPACAMETName);
  fTokPVName         = iC.consumes<reco::VertexCollection>           (fPVName);
  if(fUseAOD)  fTokPFCandName     = iC.consumes<reco::PFCandidateCollection>    (fPFCandName);
  if(!fUseAOD) fTokPackCandName   = iC.consumes<pat::PackedCandidateCollection> (fPFCandName);
  fTokRhoIso         = iC.consumes<double>                            (fRhoIsoName);
  fTokRhoJet         = iC.consumes<double>                            (fRhoJetName);
  edm::InputTag lBeamHaloSummary("BeamHaloSummary");
  edm::InputTag lHBHENoiseFilter("HBHENoiseFilterResultProducer","HBHENoiseFilterResult");
  edm::InputTag lLaserEvtFilter ("hcalLaserEventFilter");
  edm::InputTag lEEBadScFilter  ("eeBadScFilter");
  edm::InputTag lEcalDeadCell   ("EcalDeadCellTriggerPrimitiveFilter");
  edm::InputTag lTrkFailure     ("trackingFailureFilter");
  edm::InputTag lManyStrip53X   ("manystripclus53X");
  edm::InputTag lTooMany53X     ("toomanystripclus53X");
  edm::InputTag lLogError       ("logErrorTooManyClusters");
  edm::InputTag lBadChCand      ("BadChargedCandidateFilter");
  edm::InputTag lBadPFMuon      ("BadPFMuonFilter");
  edm::InputTag lMetFilters     ("TriggerResults","");
  edm::InputTag lPrefiring      ("prefiringweight:nonPrefiringProb");
  edm::InputTag lPrefiringUp    ("prefiringweight:nonPrefiringProbUp");
  edm::InputTag lPrefiringDown  ("prefiringweight:nonPrefiringProbDown");
  fTokBeamHaloSummary                    = iC.consumes<reco::BeamHaloSummary> (lBeamHaloSummary);
  fTokHBHENoiseFilterResultProducer      = iC.consumes<bool>                  (lHBHENoiseFilter);
  fTokHcalLaserEventFilter               = iC.consumes<bool>                  (lLaserEvtFilter);
  fTokEEBadScFilter                      = iC.consumes<bool>                  (lEEBadScFilter);
  fTokEcalDeadCellTriggerPrimitiveFilter = iC.consumes<bool>                  (lEcalDeadCell);
  fToktrackingFailureFilter              = iC.consumes<bool>                  (lTrkFailure);
  fTokManystripClus53X                   = iC.consumes<bool>                  (lManyStrip53X);
  fTokTooManyStripClus53X                = iC.consumes<bool>                  (lTooMany53X);
  fToklogErrorTooManyClusters            = iC.consumes<bool>                  (lLogError);
  fTokBadChCand                          = iC.consumes<bool>                  (lBadChCand);
  fTokBadPFMuon                          = iC.consumes<bool>                  (lBadPFMuon);
  fTokMetFiltersTag                      = iC.consumes<edm::TriggerResults>   (lMetFilters);
  fTokPrefWeight                         = iC.consumes<double>                (lPrefiring);
  fTokPrefWeightUp                       = iC.consumes<double>                (lPrefiringUp);
  fTokPrefWeightDown                     = iC.consumes<double>                (lPrefiringDown);
}
//--------------------------------------------------------------------------------------------------
FillerEventInfo::~FillerEventInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerEventInfo::fill(TEventInfo *evtInfo,
                           const edm::Event &iEvent, const reco::Vertex &pv,
                           const bool hasGoodPV,
			   const TriggerBits triggerBits)
{
  assert(evtInfo);
  
  evtInfo->runNum  = iEvent.id().run();
  evtInfo->lumiSec = iEvent.luminosityBlock();
  evtInfo->evtNum  = iEvent.id().event();
  //
  // Pile-up info
  //==============================
  if(!iEvent.isRealData()) {
    edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct;
    iEvent.getByToken(fTokPUInfoName,hPileupInfoProduct);
    assert(hPileupInfoProduct.isValid());
    const std::vector<PileupSummaryInfo> *inPUInfos = hPileupInfoProduct.product();
    for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo!=inPUInfos->end(); ++itPUInfo) {
      if(itPUInfo->getBunchCrossing()==0) {
        evtInfo->nPU      = itPUInfo->getPU_NumInteractions();
        evtInfo->nPUmean  = itPUInfo->getTrueNumInteractions();
      } else if(itPUInfo->getBunchCrossing()==-1) { 
        evtInfo->nPUm     = itPUInfo->getPU_NumInteractions();
        evtInfo->nPUmeanm = itPUInfo->getTrueNumInteractions();
      } else if(itPUInfo->getBunchCrossing()==1) {
        evtInfo->nPUp     = itPUInfo->getPU_NumInteractions();
        evtInfo->nPUmeanp = itPUInfo->getTrueNumInteractions();
      }
    }
  }

  
  //
  // primary vertex info
  //==============================
  evtInfo->pvx = pv.x();
  evtInfo->pvy = pv.y();
  evtInfo->pvz = pv.z();
  evtInfo->hasGoodPV = hasGoodPV;
 
  
  //
  // beam spot info
  //==============================
  edm::Handle<reco::BeamSpot> hBeamSpotProduct;
  iEvent.getByToken(fTokBSName,hBeamSpotProduct);
  assert(hBeamSpotProduct.isValid());
  const reco::BeamSpot *bs = hBeamSpotProduct.product();
  evtInfo->bsx = bs->x0();
  evtInfo->bsy = bs->y0();
  evtInfo->bsz = bs->z0();

  // Egamma prefiring weights
  edm::Handle<double> prefweight;
  iEvent.getByToken(fTokPrefWeight, prefweight);
  evtInfo->prefweight = (*prefweight);

  edm::Handle<double> prefweightUp;
  iEvent.getByToken(fTokPrefWeightUp, prefweightUp);
  evtInfo->prefweightUp = (*prefweightUp);
  
  edm::Handle<double> prefweightDown;
  iEvent.getByToken(fTokPrefWeightDown, prefweightDown);
  evtInfo->prefweightDown = (*prefweightDown);
  	

  //
  // MET filter
  //==============================
  evtInfo->metFilterFailBits=0;
  if(fUseFilters) { 
    if(fUseAOD) {  // === AOD ===
      // beam halo filter using CSCs
      edm::Handle<reco::BeamHaloSummary> hBeamHaloSummary;
      iEvent.getByToken(fTokBeamHaloSummary,hBeamHaloSummary);
      assert(hBeamHaloSummary.isValid());
      const reco::BeamHaloSummary *beamHaloSummary = hBeamHaloSummary.product();
      if(beamHaloSummary->CSCTightHaloId()) {  // if true, then event has identified beam halo
	evtInfo->metFilterFailBits |= kCSCTightHaloFilter;
      }
      
      // HB,HE anomalous noise filter
      edm::Handle<bool> hHBHENoiseFilterResult;
      iEvent.getByToken(fTokHBHENoiseFilterResultProducer,hHBHENoiseFilterResult);
      assert(hHBHENoiseFilterResult.isValid());
      if(!(*hHBHENoiseFilterResult)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHBHENoiseFilter;
      }
      
      // HCAL laser filter
      edm::Handle<bool> hHCALLaserEventFilter;
      iEvent.getByToken(fTokHcalLaserEventFilter,hHCALLaserEventFilter);
      assert(hHCALLaserEventFilter.isValid());
      if(!(*hHCALLaserEventFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHCALLaserEventFilter;
      }
      
      // bad EE SuperCrystal filter
      edm::Handle<bool> hEEBadScFilter;
      iEvent.getByToken(fTokEEBadScFilter,hEEBadScFilter);
      assert(hEEBadScFilter.isValid());
      if(!(*hEEBadScFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kEEBadScFilter;
      }
      
      // ECAL dead cell filter using trigger primitives
      edm::Handle<bool> hECALDeadCellTriggerPrimitiveFilter;
      iEvent.getByToken(fTokEcalDeadCellTriggerPrimitiveFilter,hECALDeadCellTriggerPrimitiveFilter);
      assert(hECALDeadCellTriggerPrimitiveFilter.isValid());
      if(!(*hECALDeadCellTriggerPrimitiveFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kECALDeadCellTriggerPrimitiveFilter;
      }

      /*
      // ECAL bad laser correction filter
      edm::Handle<bool> hECALLaserCorrFilter;
      //iEvent.getByLabel("ecalLaserCorrFilter",hECALLaserCorrFilter);
      assert(hECALLaserCorrFilter.isValid());
      if(!(*hECALLaserCorrFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kECALLaserCorrFilter;
      }
      */
      // tracking failure filter
      edm::Handle<bool> hTrackingFailureFilter;
      iEvent.getByToken(fToktrackingFailureFilter,hTrackingFailureFilter);
      assert(hTrackingFailureFilter.isValid());
      if(!(*hTrackingFailureFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrackingFailureFilter;
      }
      
      // tracking POG filters
      edm::Handle<bool> hTrkPOGFilter_manystripclus53X;
      iEvent.getByToken(fTokManystripClus53X,hTrkPOGFilter_manystripclus53X);
      assert(hTrkPOGFilter_manystripclus53X.isValid());
      if(*hTrkPOGFilter_manystripclus53X) {  // if result is "true", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrkPOGFilter_manystripclus53X;
      }
      
      edm::Handle<bool> hTrkPOGFilter_toomanystripclus53X;
      iEvent.getByToken(fTokTooManyStripClus53X,hTrkPOGFilter_toomanystripclus53X);
      assert(hTrkPOGFilter_toomanystripclus53X.isValid());
      if(*hTrkPOGFilter_toomanystripclus53X) {  // if result is "true", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrkPOGFilter_toomanystripclus53X;
      }
      
      edm::Handle<bool> hTrkPOGFilter_logErrorTooManyClusters;
      iEvent.getByToken(fToklogErrorTooManyClusters,hTrkPOGFilter_logErrorTooManyClusters);
      assert(hTrkPOGFilter_logErrorTooManyClusters.isValid());
      if(*hTrkPOGFilter_logErrorTooManyClusters) {  // if result is "true", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrkPOGFilter_logErrorTooManyClusters;
      }


      
    } else {  // === MINIAOD ===
      /*
      // beam halo filter using CSCs
      edm::Handle<reco::BeamHaloSummary> hBeamHaloSummary;
      //iEvent.getByLabel("BeamHaloSummary",hBeamHaloSummary);
      assert(hBeamHaloSummary.isValid());
      const reco::BeamHaloSummary *beamHaloSummary = hBeamHaloSummary.product();
      if(beamHaloSummary->CSCTightHaloId()) {  // if true, then event has identified beam halo
	evtInfo->metFilterFailBits |= kCSCTightHaloFilter;
      }
      
      // HB,HE anomalous noise filter
      edm::Handle<bool> hHBHENoiseFilterResult;
      //iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult",hHBHENoiseFilterResult);
      assert(hHBHENoiseFilterResult.isValid());
      if(!(*hHBHENoiseFilterResult)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHBHENoiseFilter;
      }
      
      // HCAL laser filter
      edm::Handle<bool> hHCALLaserEventFilter;
      //iEvent.getByLabel("hcalLaserEventFilter",hHCALLaserEventFilter);
      assert(hHCALLaserEventFilter.isValid());
      if(!(*hHCALLaserEventFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHCALLaserEventFilter;
      }
      
      // tracking failure filter
      edm::Handle<bool> hTrackingFailureFilter;
      //iEvent.getByLabel("trackingFailureFilter",hTrackingFailureFilter);
      assert(hTrackingFailureFilter.isValid());
      if(!(*hTrackingFailureFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrackingFailureFilter;
      }
      */

      //New MET Filters
      edm::Handle<bool> hFilterBadChCand;
      iEvent.getByToken(fTokBadChCand, hFilterBadChCand);
      assert(hFilterBadChCand.isValid());
      if(!(*hFilterBadChCand)) { 
	evtInfo->metFilterFailBits |= kBadChCandFilter;
      }
      
      edm::Handle<bool> hFilterBadPFMuon;
      iEvent.getByToken(fTokBadPFMuon, hFilterBadPFMuon);
      assert(hFilterBadPFMuon.isValid());
      if(!(*hFilterBadPFMuon)) { 
	evtInfo->metFilterFailBits |= kBadPFMuonFilter;
      }

      //edm::InputTag metFiltersTag("TriggerResults","","PAT");
      edm::InputTag metFiltersTag("TriggerResults","","HLT");
      edm::Handle<edm::TriggerResults> hMETFilters;
      iEvent.getByToken(fTokMetFiltersTag,hMETFilters);
      assert(hMETFilters.isValid());
      const edm::TriggerNames &metFilterNames = iEvent.triggerNames(*hMETFilters);
      
      unsigned int index;

      // beam halo filter using CSCs
      index = metFilterNames.triggerIndex("Flag_CSCTightHalo2015Filter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kCSCTightHaloFilter;
	}
      }

      index = metFilterNames.triggerIndex("Flag_BadPFMuonFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kBadPFMuonFilter;
	}
      }

      index = metFilterNames.triggerIndex("Flag_BadChargedCandidateFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kBadChCandFilter;
	}
      }
      
      // HB,HE anomalous noise filter
      index = metFilterNames.triggerIndex("Flag_HBHENoiseFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kHBHENoiseFilter;
	}
      }

      // HB,HE anomalous noise filter
      index = metFilterNames.triggerIndex("Flag_HBHENoiseIsoFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kHBHENoiseIsoFilter;
	}
      }
      
      // HCAL laser filter
      index = metFilterNames.triggerIndex("Flag_hcalLaserEventFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kHCALLaserEventFilter;
	}
      }
      // bad EE SuperCrystal filter
      index = metFilterNames.triggerIndex("Flag_eeBadScFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kEEBadScFilter;
	}
      }
      
      // ECAL dead cell filter using trigger primitives
      index = metFilterNames.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kECALDeadCellTriggerPrimitiveFilter;
	}
      }
      // Good vertex MET Filter
      index = metFilterNames.triggerIndex("Flag_goodVertices");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kGoodVerticesFilter;
	}
      }
      // Tight beam Halo filter
      index = metFilterNames.triggerIndex("Flag_globalTightHalo2016Filter");
      if(index < hMETFilters->size()) {  // check for valid index                                                                                                                                       
        if(!hMETFilters->accept(index)) {
          evtInfo->metFilterFailBits |= kGlobalTightHalo2016Filter;
        }
      }
      // Global Super Tight beam Halo filter
      index = metFilterNames.triggerIndex("Flag_globalSuperTightHalo2016Filter");
      if(index < hMETFilters->size()) {  // check for valid index                                                                                                                                       
        if(!hMETFilters->accept(index)) {
          evtInfo->metFilterFailBits |= kGlobalSuperTightHalo2016Filter;
        }
      }
      // Good vertex MET Filter
      index = metFilterNames.triggerIndex("Flag_chargedHadronTrackResolutionFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kChargedHadronTrackResolutionFilter;
	}
      }
      index = metFilterNames.triggerIndex("Flag_muonBadTrackFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kMuonBadTrackFilter;
	}
      }
      index = metFilterNames.triggerIndex("Flag_ecalBadCalibFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kEcalBadCalibFilter;
	}
      }
      /*
      // ECAL bad laser correction filter
      index = metFilterNames.triggerIndex("Flag_ecalLaserCorrFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kECALLaserCorrFilter;
	}
      }
      */
      // tracking failure filter
      index = metFilterNames.triggerIndex("Flag_trackingFailureFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kTrackingFailureFilter;
	}
      }
      
      // tracking POG filters
      index = metFilterNames.triggerIndex("Flag_trkPOG_manystripclus53X");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {  //(!) is this correct?
	  evtInfo->metFilterFailBits |= kTrkPOGFilter_manystripclus53X;
	}
      }
      
      index = metFilterNames.triggerIndex("Flag_trkPOG_toomanystripclus53X");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {  //(!) is this correct?
	  evtInfo->metFilterFailBits |= kTrkPOGFilter_toomanystripclus53X;
	}
      }
      
      index = metFilterNames.triggerIndex("Flag_trkPOG_logErrorTooManyClusters");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {  //(!) is this correct?
	  evtInfo->metFilterFailBits |= kTrkPOGFilter_logErrorTooManyClusters;
	}
      }
    }
  }
  //
  // MET info
  //==============================
  if(fUseAOD) {  // === AOD ===
    // CaloMET
    edm::Handle<reco::CaloMETCollection> hCaloMETProduct;
    iEvent.getByToken(fTokCaloMETName,hCaloMETProduct);
    assert(hCaloMETProduct.isValid());
    const reco::CaloMET &inCaloMET = hCaloMETProduct.product()->front();
    evtInfo->caloMET      = inCaloMET.pt();
    evtInfo->caloMETphi   = inCaloMET.phi();
    //evtInfo->caloMETCov00 = inCaloMET.getSignificanceMatrix()(0,0);
    //evtInfo->caloMETCov01 = inCaloMET.getSignificanceMatrix()(0,1);
    //evtInfo->caloMETCov11 = inCaloMET.getSignificanceMatrix()(1,1);
    
    // Raw PF MET
    edm::Handle<reco::PFMETCollection> hPFMETProduct;
    iEvent.getByToken(fTokPFMETName,hPFMETProduct);
    assert(hPFMETProduct.isValid());
    const reco::PFMET &inPFMET = hPFMETProduct.product()->front();
    evtInfo->pfMET      = inPFMET.pt();
    evtInfo->pfMETphi   = inPFMET.phi();
    evtInfo->pfMETCov00 = inPFMET.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCov01 = inPFMET.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCov11 = inPFMET.getSignificanceMatrix()(1,1);

    // Corrected PF MET
    edm::Handle<reco::PFMETCollection> hPFMETCProduct;
    iEvent.getByToken(fTokPFMETCName,hPFMETCProduct);
    assert(hPFMETCProduct.isValid());
    const reco::PFMET &inPFMETC = hPFMETCProduct.product()->front();
    evtInfo->pfMETC      = inPFMETC.pt();
    evtInfo->pfMETCphi   = inPFMETC.phi();
    evtInfo->pfMETCCov00 = inPFMETC.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCCov01 = inPFMETC.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCCov11 = inPFMETC.getSignificanceMatrix()(1,1);

    // MVA MET
    /*
    if(fMVAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hMVAMETProduct;
      iEvent.getByToken(fTokMVAMETName,hMVAMETProduct);
      assert(hMVAMETProduct.isValid());
      const reco::PFMET &inMVAMET = hMVAMETProduct.product()->front();
      evtInfo->mvaMET      = inMVAMET.pt();
      evtInfo->mvaMETphi   = inMVAMET.phi();
      evtInfo->mvaMETCov00 = inMVAMET.getSignificanceMatrix()(0,0);
      evtInfo->mvaMETCov01 = inMVAMET.getSignificanceMatrix()(0,1);
      evtInfo->mvaMETCov11 = inMVAMET.getSignificanceMatrix()(1,1);
    }
    // MVA MET with unity response
    edm::Handle<reco::PFMETCollection> hMVAMETUProduct;
    iEvent.getByToken(fTokMVAMETUName,hMVAMETUProduct);
    assert(hMVAMETUProduct.isValid());
    const reco::PFMET &inMVAMETU = hMVAMETUProduct.product()->front();
    evtInfo->mvaMETU      = inMVAMETU.pt();
    evtInfo->mvaMETUphi   = inMVAMETU.phi();
    //    evtInfo->mvaMETUCov00 = inMVAMETU.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMETUCov01 = inMVAMETU.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMETUCov11 = inMVAMETU.getSignificanceMatrix()(1,1);

    // MVA MET without jet smearing (relevant only for MC)
    edm::Handle<reco::PFMETCollection> hMVAMET0Product;
    iEvent.getByToken(fTokMVAMET0Name,hMVAMET0Product);
    assert(hMVAMET0Product.isValid());
    const reco::PFMET &inMVAMET0 = hMVAMET0Product.product()->front();
    evtInfo->mvaMET0      = inMVAMET0.pt();
    evtInfo->mvaMET0phi   = inMVAMET0.phi();
    //    evtInfo->mvaMET0Cov00 = inMVAMET0.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMET0Cov01 = inMVAMET0.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMET0Cov11 = inMVAMET0.getSignificanceMatrix()(1,1);
    */
    // PUPPI MET
    edm::Handle<reco::PFMETCollection> hPuppET;
    iEvent.getByToken(fTokPUPPETName,hPuppET);
    assert(hPuppET.isValid());
    const reco::PFMET &inPuppET = hPuppET.product()->front();
    evtInfo->puppET      = inPuppET.pt();
    evtInfo->puppETphi   = inPuppET.phi();
    evtInfo->puppETCov00 = inPuppET.getSignificanceMatrix()(0,0);
    evtInfo->puppETCov01 = inPuppET.getSignificanceMatrix()(0,1);
    evtInfo->puppETCov11 = inPuppET.getSignificanceMatrix()(1,1);

    // Type 1 PUPPI MET
    edm::Handle<reco::PFMETCollection> hPuppETC;
    iEvent.getByToken(fTokPUPPETCName,hPuppETC);
    assert(hPuppETC.isValid());
    const reco::PFMET &inPuppETC = hPuppETC.product()->front();
    evtInfo->puppETC      = inPuppETC.pt();
    evtInfo->puppETCphi   = inPuppETC.phi();
    evtInfo->puppETCCov00 = inPuppETC.getSignificanceMatrix()(0,0);
    evtInfo->puppETCCov01 = inPuppETC.getSignificanceMatrix()(0,1);
    evtInfo->puppETCCov11 = inPuppETC.getSignificanceMatrix()(1,1);

    // Alpaca MET
    if(fALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hAlpacaMET;
      iEvent.getByToken(fTokALPACAMETName,hAlpacaMET);
      assert(hAlpacaMET.isValid());
      const reco::PFMET &inAlpacaMET = hAlpacaMET.product()->front();
      evtInfo->alpacaMET      = inAlpacaMET.pt();
      evtInfo->alpacaMETphi   = inAlpacaMET.phi();
    }
    if(fPALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPAlpacaMET;
      iEvent.getByToken(fTokPALPACAMETName,hPAlpacaMET);
      assert(hPAlpacaMET.isValid());
      const reco::PFMET &inPAlpacaMET = hPAlpacaMET.product()->front();
      evtInfo->pcpMET      = inPAlpacaMET.pt();
      evtInfo->pcpMETphi   = inPAlpacaMET.phi();
    }

    // Track MET
    /*
    edm::Handle<reco::VertexCollection> hVertexProduct;
    iEvent.getByToken(fTokPVName,hVertexProduct);
    assert(hVertexProduct.isValid());
    const reco::VertexCollection *pvCol = hVertexProduct.product();

    edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
    iEvent.getByToken(fTokPFCandName,hPFCandProduct);
    assert(hPFCandProduct.isValid());
    const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
    computeTrackMET(0, pvCol, pfCandCol, evtInfo->trkMET, evtInfo->trkMETphi);
    */
  } else {  // === MINIAOD ===

    edm::Handle<pat::METCollection> hMETProduct;
    iEvent.getByToken(fTokPFMETPATName,hMETProduct);
    assert(hMETProduct.isValid());
    const pat::MET &inMET = hMETProduct->front();
    // Raw PF MET
    evtInfo->pfMET      = inMET.uncorPt();
    evtInfo->pfMETphi   = inMET.uncorPhi();
    evtInfo->pfMETCov00 = inMET.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCov01 = inMET.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCov11 = inMET.getSignificanceMatrix()(1,1);

    // Corrected PF MET
    evtInfo->pfMETC      = inMET.pt();
    evtInfo->pfMETCphi   = inMET.phi();
    evtInfo->pfMETCCov00 = inMET.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCCov01 = inMET.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCCov11 = inMET.getSignificanceMatrix()(1,1);
    evtInfo->pfMETCjerup = inMET.shiftedP4(pat::MET::JetResUp).pt();
    evtInfo->pfMETCjerdn = inMET.shiftedP4(pat::MET::JetResDown).pt();
    evtInfo->pfMETCjenup = inMET.shiftedP4(pat::MET::JetEnUp).pt();
    evtInfo->pfMETCjendn = inMET.shiftedP4(pat::MET::JetEnDown).pt();
    evtInfo->pfMETCuncup = inMET.shiftedP4(pat::MET::UnclusteredEnUp).pt();
    evtInfo->pfMETCuncdn = inMET.shiftedP4(pat::MET::UnclusteredEnDown).pt();
    evtInfo->pfMETCjrsup = inMET.shiftedP4(pat::MET::PhotonEnUp).pt();
    evtInfo->pfMETCjrsdn = inMET.shiftedP4(pat::MET::PhotonEnDown).pt();
    evtInfo->pfMETCphijerup = inMET.shiftedP4(pat::MET::JetResUp).phi();
    evtInfo->pfMETCphijerdn = inMET.shiftedP4(pat::MET::JetResDown).phi();
    evtInfo->pfMETCphijenup = inMET.shiftedP4(pat::MET::JetEnUp).phi();
    evtInfo->pfMETCphijendn = inMET.shiftedP4(pat::MET::JetEnDown).phi();
    evtInfo->pfMETCphiuncup = inMET.shiftedP4(pat::MET::UnclusteredEnUp).phi();
    evtInfo->pfMETCphiuncdn = inMET.shiftedP4(pat::MET::UnclusteredEnDown).phi();
    evtInfo->pfMETCphijrsup = inMET.shiftedP4(pat::MET::PhotonEnUp).phi();
    evtInfo->pfMETCphijrsdn = inMET.shiftedP4(pat::MET::PhotonEnDown).phi();

    //Calo MET
    edm::Handle<pat::METCollection> hCaloMETProduct;
    iEvent.getByToken(fTokPFMETPATName,hCaloMETProduct);
    assert(hCaloMETProduct.isValid());
    const pat::MET &inCaloMET = hCaloMETProduct->front();
    evtInfo->caloMET      = inCaloMET.caloMETPt();
    evtInfo->caloMETphi   = inCaloMET.caloMETPhi();
    //evtInfo->caloMETCov00 = inCaloMET.getSignificanceMatrix()(0,0);
    //evtInfo->caloMETCov01 = inCaloMET.getSignificanceMatrix()(0,1);
    //evtInfo->caloMETCov11 = inCaloMET.getSignificanceMatrix()(1,1);
    // MVA MET
    /*
    if(fMVAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hMVAMETProduct;
      iEvent.getByToken(fTokMVAMETName,hMVAMETProduct);
      assert(hMVAMETProduct.isValid());
      const reco::PFMET &inMVAMET = hMVAMETProduct.product()->front();
      evtInfo->mvaMET      = inMVAMET.pt();
      evtInfo->mvaMETphi   = inMVAMET.phi();
    }
    */
    //    evtInfo->mvaMETCov00 = inMVAMET.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMETCov01 = inMVAMET.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMETCov11 = inMVAMET.getSignificanceMatrix()(1,1);

    // MVA MET with unity response
    //edm::Handle<reco::PFMETCollection> hMVAMETUProduct;
    //iEvent.getByToken(fTokMVAMETUName,hMVAMETUProduct);
    //assert(hMVAMETUProduct.isValid());
    //const reco::PFMET &inMVAMETU = hMVAMETUProduct.product()->front();
    //evtInfo->mvaMETU      = inMVAMETU.pt();
    //evtInfo->mvaMETUphi   = inMVAMETU.phi();
    //    evtInfo->mvaMETUCov00 = inMVAMETU.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMETUCov01 = inMVAMETU.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMETUCov11 = inMVAMETU.getSignificanceMatrix()(1,1);

    // MVA MET without jet smearing (relevant only for MC)
    //edm::Handle<reco::PFMETCollection> hMVAMET0Product;
    //iEvent.getByToken(fTokMVAMET0Name,hMVAMET0Product);
    //assert(hMVAMET0Product.isValid());
    //const reco::PFMET &inMVAMET0 = hMVAMET0Product.product()->front();
    //evtInfo->mvaMET0      = inMVAMET0.pt();
    //evtInfo->mvaMET0phi   = inMVAMET0.phi();
    //    evtInfo->mvaMET0Cov00 = inMVAMET0.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMET0Cov01 = inMVAMET0.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMET0Cov11 = inMVAMET0.getSignificanceMatrix()(1,1);

    // PUPPI MET
    edm::Handle<pat::METCollection> hPuppET;
    iEvent.getByToken(fTokPUPPETPATName,hPuppET);
    assert(hPuppET.isValid());
    const pat::MET &inPuppET = hPuppET.product()->front();
    evtInfo->puppET      = inPuppET.uncorPt();
    evtInfo->puppETphi   = inPuppET.uncorPhi();
    //Type1 PUPPI MET
    evtInfo->puppETC      = inPuppET.pt();
    evtInfo->puppETCphi   = inPuppET.phi();
    evtInfo->puppETCov00  = inPuppET.getSignificanceMatrix()(0,0);
    evtInfo->puppETCov01  = inPuppET.getSignificanceMatrix()(0,1);
    evtInfo->puppETCov11  = inPuppET.getSignificanceMatrix()(1,1);
    /*
    evtInfo->puppETCjerup = inPuppET.shiftedPt(pat::MET::METUncertainty::JetResUp);
    evtInfo->puppETCjerdn = inPuppET.shiftedPt(pat::MET::METUncertainty::JetResDown);
    evtInfo->puppETCjenup = inPuppET.shiftedPt(pat::MET::METUncertainty::JetEnUp);
    evtInfo->puppETCjendn = inPuppET.shiftedPt(pat::MET::METUncertainty::JetEnDown);
    evtInfo->puppETCuncup = inPuppET.shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp);
    evtInfo->puppETCuncdn = inPuppET.shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown);
    evtInfo->puppETCjrsup = inPuppET.shiftedPt(pat::MET::METUncertainty::JetResUpSmear);
    evtInfo->puppETCjrsdn = inPuppET.shiftedPt(pat::MET::METUncertainty::JetResDownSmear);
    evtInfo->puppETCphijerup = inPuppET.shiftedPhi(pat::MET::METUncertainty::JetResUp);
    evtInfo->puppETCphijerdn = inPuppET.shiftedPhi(pat::MET::METUncertainty::JetResDown);
    evtInfo->puppETCphijenup = inPuppET.shiftedPhi(pat::MET::METUncertainty::JetEnUp);
    evtInfo->puppETCphijendn = inPuppET.shiftedPhi(pat::MET::METUncertainty::JetEnDown);
    evtInfo->puppETCphiuncup = inPuppET.shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp);
    evtInfo->puppETCphiuncdn = inPuppET.shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown);
    evtInfo->puppETCphijrsup = inPuppET.shiftedPhi(pat::MET::METUncertainty::JetResUpSmear);
    evtInfo->puppETCphijrsdn = inPuppET.shiftedPhi(pat::MET::METUncertainty::JetResDownSmear);
    */
    // Alpaca MET
    if(fALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hAlpacaMET;
      iEvent.getByToken(fTokALPACAMETName,hAlpacaMET);
      assert(hAlpacaMET.isValid());
      const reco::PFMET &inAlpacaMET = hAlpacaMET.product()->front();
      evtInfo->alpacaMET      = inAlpacaMET.pt();
      evtInfo->alpacaMETphi   = inAlpacaMET.phi();
    }
    if(fPALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPAlpacaMET;
      iEvent.getByToken(fTokPALPACAMETName,hPAlpacaMET);
      assert(hPAlpacaMET.isValid());
      const reco::PFMET &inPAlpacaMET = hPAlpacaMET.product()->front();
      evtInfo->pcpMET      = inPAlpacaMET.pt();
      evtInfo->pcpMETphi   = inPAlpacaMET.phi();
    }
    // Track MET
    edm::Handle<pat::PackedCandidateCollection> hPFCandProduct;
    iEvent.getByToken(fTokPackCandName,hPFCandProduct);
    assert(hPFCandProduct.isValid());
    const pat::PackedCandidateCollection *pfCandCol = hPFCandProduct.product();
    computeTrackMET(pfCandCol, evtInfo->trkMET, evtInfo->trkMETphi);
  }
  
  
  //
  // event energy density
  //==============================
  
  // Rho for isolation correction
  edm::Handle<double> hRhoIso;
  iEvent.getByToken(fTokRhoIso,hRhoIso);
  assert(hRhoIso.isValid());
  evtInfo->rhoIso = *hRhoIso;
  
  // Rho for jet energy correction
  edm::Handle<double> hRhoJet;
  iEvent.getByToken(fTokRhoJet,hRhoJet);
  assert(hRhoJet.isValid());
  evtInfo->rhoJet = *hRhoJet;


  //
  // fired triggers
  //==============================
  evtInfo->triggerBits = triggerBits;
}

void FillerEventInfo::computeTrackMET(const unsigned int ipv,
                                      const reco::VertexCollection *pvCol,
                                      const reco::PFCandidateCollection *pfCandCol,
                                      float &out_met, float &out_metphi)
{
  out_met    = 0;
  out_metphi = 0;

  const reco::Vertex pv = pvCol->at(ipv);

  double metx=0, mety=0;
  for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
    if(itPF->bestTrack()!=0) {
      if(itPF->trackRef().isNonnull() && pv.trackWeight(itPF->trackRef())>0) {
        metx -= itPF->px();
        mety -= itPF->py();

      } else if(itPF->particleId()==reco::PFCandidate::e || itPF->particleId()==reco::PFCandidate::mu) {
        metx -= itPF->px();
        mety -= itPF->py();

      } else {
        unsigned int iclosest = pvCol->size();  // note: initialized to an index not within bounds!                                                                                                         
        double minDZ = 999.;
        for(unsigned int iv=0; iv<pvCol->size(); iv++) {
          const reco::Vertex vtx = pvCol->at(iv);

          // if track is used in some other vertex fit, don't use it for track MET                                                                                                                          
          if(itPF->trackRef().isNonnull() && vtx.trackWeight(itPF->trackRef())>0) {
            iclosest = pvCol->size();  // note: set to an index not within bounds!                                                                                                                          
            break;
          }

          double dz = fabs(itPF->vertex().z() - vtx.z());
          if(dz < minDZ) {
            minDZ    = dz;
            iclosest = iv;
          }
        }

        if(iclosest==ipv) {
          metx -= itPF->px();
          mety -= itPF->py();
        }
      }
    }
  }
  TLorentzVector met;
  met.SetPxPyPzE(metx,mety,0,0);
  out_met    = met.Pt();
  out_metphi = met.Phi();
}


void FillerEventInfo::computeTrackMET(const pat::PackedCandidateCollection *pfCandCol,
                                      float &out_met, float &out_metphi)
{
  out_met    = 0;
  out_metphi = 0;

  double metx=0, mety=0;
  for(pat::PackedCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
    // track is:
    // 1) used in the fit of the PV (status=3)
    // 2) not used in fit of any PV but closest in z to the PV (status=2)
    if(itPF->bestTrack()!=0 && itPF->fromPV()>1) {  // (!) MINIAOD: with respect to PV[0]
      metx  -= itPF->px();
      mety  -= itPF->py();
    }
  }

  TLorentzVector met;
  met.SetPxPyPzE(metx,mety,0,0);
  out_met    = met.Pt();
  out_metphi = met.Phi();
}
