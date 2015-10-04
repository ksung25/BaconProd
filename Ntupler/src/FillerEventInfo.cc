#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include <TLorentzVector.h>
#include <string>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerEventInfo::FillerEventInfo(const edm::ParameterSet &iConfig, const bool useAOD):
  fPFCandName    (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fPUInfoName    (iConfig.getUntrackedParameter<std::string>("edmPileupInfoName","addPileupInfo")),
  fPVName        (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fBSName        (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot")),
  fCaloMETName   (iConfig.getUntrackedParameter<std::string>("edmCaloMETName","caloMet")),
  fMETName       (iConfig.getUntrackedParameter<std::string>("edmPFMETName","slimmedMETs")),
  fPFMETName     (iConfig.getUntrackedParameter<std::string>("edmPFMETName","pfMet")),
  fPFMETCName    (iConfig.getUntrackedParameter<std::string>("edmPFMETCorrName","pfType1CorrectedMet")),
  fMVAMETName    (iConfig.getUntrackedParameter<std::string>("edmMVAMETName","pfMEtMVA")),
  fPUPPETName    (iConfig.getUntrackedParameter<std::string>("edmPuppETName","pfMetPuppi")),
  fPUPPETCName   (iConfig.getUntrackedParameter<std::string>("edmPuppETCorrName","pfType1CorrectedMetPuppi")),
  fPFMET30Name   (iConfig.getUntrackedParameter<std::string>("edmPFMET30Name","pfMet30")),
  fPFMETC30Name  (iConfig.getUntrackedParameter<std::string>("edmPFMET30CorrName","pfType1CorrectedMet30")),
  fMVAMET30Name  (iConfig.getUntrackedParameter<std::string>("edmMVAMET30Name","pfMEtMVA30")),
  fPUPPET30Name  (iConfig.getUntrackedParameter<std::string>("edmPuppET30Name","pfMetPuppi30")),
  fPUPPETC30Name (iConfig.getUntrackedParameter<std::string>("edmPuppET30CorrName","pfType1CorrectedMetPuppi30")),
  fALPACAMETName (iConfig.getUntrackedParameter<std::string>("edmAlpacaMETName"    ,"pfMetAlpacaMC")),
  fPALPACAMETName(iConfig.getUntrackedParameter<std::string>("edmPupAlpacaMETName","pfMetPuppiAlpacaMC")),
  fRhoIsoName    (iConfig.getUntrackedParameter<std::string>("edmRhoForIsoName","fixedGridRhoFastjetAll")),
  fRhoJetName    (iConfig.getUntrackedParameter<std::string>("edmRhoForJetEnergy","fixedGridRhoFastjetAll")),
  fUseFilters    (iConfig.getUntrackedParameter<bool>("doFillMETFilters",true)),
  fUseAOD        (useAOD)
{}

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
    iEvent.getByLabel(fPUInfoName,hPileupInfoProduct);
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
  iEvent.getByLabel(fBSName,hBeamSpotProduct);
  assert(hBeamSpotProduct.isValid());
  const reco::BeamSpot *bs = hBeamSpotProduct.product();
  evtInfo->bsx = bs->x0();
  evtInfo->bsy = bs->y0();
  evtInfo->bsz = bs->z0();
  	

  //
  // MET filter
  //==============================
  evtInfo->metFilterFailBits=0;
  if(fUseFilters) { 
    if(fUseAOD) {  // === AOD ===
      // beam halo filter using CSCs
      edm::Handle<reco::BeamHaloSummary> hBeamHaloSummary;
      iEvent.getByLabel("BeamHaloSummary",hBeamHaloSummary);
      assert(hBeamHaloSummary.isValid());
      const reco::BeamHaloSummary *beamHaloSummary = hBeamHaloSummary.product();
      if(beamHaloSummary->CSCTightHaloId()) {  // if true, then event has identified beam halo
	evtInfo->metFilterFailBits |= kCSCTightHaloFilter;
      }
      
      // HB,HE anomalous noise filter
      edm::Handle<bool> hHBHENoiseFilterResult;
      iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult",hHBHENoiseFilterResult);
      assert(hHBHENoiseFilterResult.isValid());
      if(!(*hHBHENoiseFilterResult)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHBHENoiseFilter;
      }
      
      // HCAL laser filter
      edm::Handle<bool> hHCALLaserEventFilter;
      iEvent.getByLabel("hcalLaserEventFilter",hHCALLaserEventFilter);
      assert(hHCALLaserEventFilter.isValid());
      if(!(*hHCALLaserEventFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHCALLaserEventFilter;
      }
      
      // bad EE SuperCrystal filter
      edm::Handle<bool> hEEBadScFilter;
      iEvent.getByLabel("eeBadScFilter",hEEBadScFilter);
      assert(hEEBadScFilter.isValid());
      if(!(*hEEBadScFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kEEBadScFilter;
      }
      
      // ECAL dead cell filter using trigger primitives
      edm::Handle<bool> hECALDeadCellTriggerPrimitiveFilter;
      iEvent.getByLabel("EcalDeadCellTriggerPrimitiveFilter",hECALDeadCellTriggerPrimitiveFilter);
      assert(hECALDeadCellTriggerPrimitiveFilter.isValid());
      if(!(*hECALDeadCellTriggerPrimitiveFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kECALDeadCellTriggerPrimitiveFilter;
      }
      
      // ECAL bad laser correction filter
      edm::Handle<bool> hECALLaserCorrFilter;
      iEvent.getByLabel("ecalLaserCorrFilter",hECALLaserCorrFilter);
      assert(hECALLaserCorrFilter.isValid());
      if(!(*hECALLaserCorrFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kECALLaserCorrFilter;
      }
      
      // tracking failure filter
      edm::Handle<bool> hTrackingFailureFilter;
      iEvent.getByLabel("trackingFailureFilter",hTrackingFailureFilter);
      assert(hTrackingFailureFilter.isValid());
      if(!(*hTrackingFailureFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrackingFailureFilter;
      }
      
      // tracking POG filters
      edm::Handle<bool> hTrkPOGFilter_manystripclus53X;
      iEvent.getByLabel("manystripclus53X",hTrkPOGFilter_manystripclus53X);
      assert(hTrkPOGFilter_manystripclus53X.isValid());
      if(*hTrkPOGFilter_manystripclus53X) {  // if result is "true", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrkPOGFilter_manystripclus53X;
      }
      
      edm::Handle<bool> hTrkPOGFilter_toomanystripclus53X;
      iEvent.getByLabel("toomanystripclus53X",hTrkPOGFilter_toomanystripclus53X);
      assert(hTrkPOGFilter_toomanystripclus53X.isValid());
      if(*hTrkPOGFilter_toomanystripclus53X) {  // if result is "true", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrkPOGFilter_toomanystripclus53X;
      }
      
      edm::Handle<bool> hTrkPOGFilter_logErrorTooManyClusters;
      iEvent.getByLabel("logErrorTooManyClusters",hTrkPOGFilter_logErrorTooManyClusters);
      assert(hTrkPOGFilter_logErrorTooManyClusters.isValid());
      if(*hTrkPOGFilter_logErrorTooManyClusters) {  // if result is "true", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrkPOGFilter_logErrorTooManyClusters;
      }
      
    } else {  // === MINIAOD ===

      //edm::InputTag metFiltersTag("TriggerResults","","PAT");
      edm::InputTag metFiltersTag("TriggerResults","","RECO");
      edm::Handle<edm::TriggerResults> hMETFilters;
      iEvent.getByLabel(metFiltersTag,hMETFilters);
      assert(hMETFilters.isValid());
      const edm::TriggerNames &metFilterNames = iEvent.triggerNames(*hMETFilters);
      
      unsigned int index;
      
      // beam halo filter using CSCs
      index = metFilterNames.triggerIndex("Flag_CSCTightHaloFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kCSCTightHaloFilter;
	}
      }
      
      // HB,HE anomalous noise filter
      index = metFilterNames.triggerIndex("Flag_HBHENoiseFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kHBHENoiseFilter;
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
      
      // ECAL bad laser correction filter
      index = metFilterNames.triggerIndex("Flag_ecalLaserCorrFilter");
      if(index < hMETFilters->size()) {  // check for valid index
	if(!hMETFilters->accept(index)) {
	  evtInfo->metFilterFailBits |= kECALLaserCorrFilter;
	}
      }
      
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
    ////!!!!! Temporary
    // Raw PF MET
    edm::Handle<reco::CaloMETCollection> hCaloMETProduct;
    iEvent.getByLabel(fCaloMETName,hCaloMETProduct);
    assert(hCaloMETProduct.isValid());
    const reco::CaloMET &inCaloMET = hCaloMETProduct.product()->front();
    evtInfo->caloMET      = inCaloMET.pt();
    evtInfo->caloMETphi   = inCaloMET.phi();
    evtInfo->caloMETCov00 = inCaloMET.getSignificanceMatrix()(0,0);
    evtInfo->caloMETCov01 = inCaloMET.getSignificanceMatrix()(0,1);
    evtInfo->caloMETCov11 = inCaloMET.getSignificanceMatrix()(1,1);

    // Raw PF MET
    edm::Handle<reco::PFMETCollection> hPFMETProduct;
    iEvent.getByLabel(fPFMETName,hPFMETProduct);
    assert(hPFMETProduct.isValid());
    const reco::PFMET &inPFMET = hPFMETProduct.product()->front();
    evtInfo->pfMET      = inPFMET.pt();
    evtInfo->pfMETphi   = inPFMET.phi();
    evtInfo->pfMETCov00 = inPFMET.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCov01 = inPFMET.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCov11 = inPFMET.getSignificanceMatrix()(1,1);

    // Corrected PF MET
    edm::Handle<reco::PFMETCollection> hPFMETCProduct;
    iEvent.getByLabel(fPFMETCName,hPFMETCProduct);
    assert(hPFMETCProduct.isValid());
    const reco::PFMET &inPFMETC = hPFMETCProduct.product()->front();
    evtInfo->pfMETC      = inPFMETC.pt();
    evtInfo->pfMETCphi   = inPFMETC.phi();
    evtInfo->pfMETCCov00 = inPFMETC.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCCov01 = inPFMETC.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCCov11 = inPFMETC.getSignificanceMatrix()(1,1);

    // MVA MET
    edm::Handle<reco::PFMETCollection> hMVAMETProduct;
    iEvent.getByLabel(fMVAMETName,hMVAMETProduct);
    assert(hMVAMETProduct.isValid());
    const reco::PFMET &inMVAMET = hMVAMETProduct.product()->front();
    evtInfo->mvaMET      = inMVAMET.pt();
    evtInfo->mvaMETphi   = inMVAMET.phi();
    evtInfo->mvaMETCov00 = inMVAMET.getSignificanceMatrix()(0,0);
    evtInfo->mvaMETCov01 = inMVAMET.getSignificanceMatrix()(0,1);
    evtInfo->mvaMETCov11 = inMVAMET.getSignificanceMatrix()(1,1);
/*
    // MVA MET with unity response
    edm::Handle<reco::PFMETCollection> hMVAMETUProduct;
    iEvent.getByLabel(fMVAMETUName,hMVAMETUProduct);
    assert(hMVAMETUProduct.isValid());
    const reco::PFMET &inMVAMETU = hMVAMETUProduct.product()->front();
    evtInfo->mvaMETU      = inMVAMETU.pt();
    evtInfo->mvaMETUphi   = inMVAMETU.phi();
//    evtInfo->mvaMETUCov00 = inMVAMETU.getSignificanceMatrix()(0,0);
//    evtInfo->mvaMETUCov01 = inMVAMETU.getSignificanceMatrix()(0,1);
//    evtInfo->mvaMETUCov11 = inMVAMETU.getSignificanceMatrix()(1,1);

    // MVA MET without jet smearing (relevant only for MC)
    edm::Handle<reco::PFMETCollection> hMVAMET0Product;
    iEvent.getByLabel(fMVAMET0Name,hMVAMET0Product);
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
    iEvent.getByLabel(fPUPPETName,hPuppET);
    assert(hPuppET.isValid());
    const reco::PFMET &inPuppET = hPuppET.product()->front();
    evtInfo->puppET      = inPuppET.pt();
    evtInfo->puppETphi   = inPuppET.phi();
    evtInfo->puppETCov00 = inPuppET.getSignificanceMatrix()(0,0);
    evtInfo->puppETCov01 = inPuppET.getSignificanceMatrix()(0,1);
    evtInfo->puppETCov11 = inPuppET.getSignificanceMatrix()(1,1);

    // Type 1 PUPPI MET
    edm::Handle<reco::PFMETCollection> hPuppETC;
    iEvent.getByLabel(fPUPPETCName,hPuppETC);
    assert(hPuppETC.isValid());
    const reco::PFMET &inPuppETC = hPuppETC.product()->front();
    evtInfo->puppETC      = inPuppETC.pt();
    evtInfo->puppETCphi   = inPuppETC.phi();
    evtInfo->puppETCCov00 = inPuppETC.getSignificanceMatrix()(0,0);
    evtInfo->puppETCCov01 = inPuppETC.getSignificanceMatrix()(0,1);
    evtInfo->puppETCCov11 = inPuppETC.getSignificanceMatrix()(1,1);

    // Raw PF MET
    if(fPFMET30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPFMET30Product;
      iEvent.getByLabel(fPFMET30Name,hPFMET30Product);
      assert(hPFMET30Product.isValid());
      const reco::PFMET &inPFMET30 = hPFMET30Product.product()->front();
      evtInfo->pfMET30      = inPFMET30.pt();
      evtInfo->pfMET30phi   = inPFMET30.phi();
    }
    // Corrected PF MET
    if(fPFMETC30Name.size() > 0) {
      edm::Handle<reco::PFMETCollection> hPFMETC30Product;
      iEvent.getByLabel(fPFMETC30Name,hPFMETC30Product);
      assert(hPFMETC30Product.isValid());
      const reco::PFMET &inPFMETC30 = hPFMETC30Product.product()->front();
      evtInfo->pfMETC30      = inPFMETC30.pt();
      evtInfo->pfMETC30phi   = inPFMETC30.phi();
    }
    // MVA MET
    if(fMVAMET30Name.size() > 0) {
      edm::Handle<reco::PFMETCollection> hMVAMET30Product;
      iEvent.getByLabel(fMVAMET30Name,hMVAMET30Product);
      assert(hMVAMET30Product.isValid());
      const reco::PFMET &inMVAMET30 = hMVAMET30Product.product()->front();
      evtInfo->mvaMET30      = inMVAMET30.pt();
      evtInfo->mvaMET30phi   = inMVAMET30.phi();
    }
    // PUPPI MET
    if(fPUPPET30Name.size() > 0) {
      edm::Handle<reco::PFMETCollection> hPuppET30;
      iEvent.getByLabel(fPUPPET30Name,hPuppET30);
      assert(hPuppET30.isValid());
      const reco::PFMET &inPuppET30 = hPuppET30.product()->front();
      evtInfo->puppET30      = inPuppET30.pt();
      evtInfo->puppET30phi   = inPuppET30.phi();
    }
    // Type 1 PUPPI MET
    if(fPUPPETC30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPuppETC30;
      iEvent.getByLabel(fPUPPETC30Name,hPuppETC30);
      assert(hPuppETC30.isValid());
      const reco::PFMET &inPuppETC30 = hPuppETC30.product()->front();
      evtInfo->puppETC30      = inPuppETC30.pt();
      evtInfo->puppETC30phi   = inPuppETC30.phi();
    }
    // Alpaca MET
    if(fALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hAlpacaMET;
      iEvent.getByLabel(fALPACAMETName,hAlpacaMET);
      assert(hAlpacaMET.isValid());
      const reco::PFMET &inAlpacaMET = hAlpacaMET.product()->front();
      evtInfo->alpacaMET      = inAlpacaMET.pt();
      evtInfo->alpacaMETphi   = inAlpacaMET.phi();
    }
    if(fPALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPAlpacaMET;
      iEvent.getByLabel(fPALPACAMETName,hPAlpacaMET);
      assert(hPAlpacaMET.isValid());
      const reco::PFMET &inPAlpacaMET = hPAlpacaMET.product()->front();
      evtInfo->pcpMET      = inPAlpacaMET.pt();
      evtInfo->pcpMETphi   = inPAlpacaMET.phi();
    }

    // Track MET
    edm::Handle<reco::VertexCollection> hVertexProduct;
    iEvent.getByLabel(fPVName,hVertexProduct);
    assert(hVertexProduct.isValid());
    const reco::VertexCollection *pvCol = hVertexProduct.product();

    edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
    iEvent.getByLabel(fPFCandName,hPFCandProduct);
    assert(hPFCandProduct.isValid());
    const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
    computeTrackMET(0, pvCol, pfCandCol, evtInfo->trkMET, evtInfo->trkMETphi);

  } else {  // === MINIAOD ===

    edm::Handle<pat::METCollection> hMETProduct;
    iEvent.getByLabel(fMETName,hMETProduct);
    assert(hMETProduct.isValid());
    const pat::MET &inMET = hMETProduct->front();
    // Raw PF MET
    evtInfo->pfMET      = inMET.pt();//shiftedPt (pat::MET::NoShift, pat::MET::Raw);
    evtInfo->pfMETphi   = inMET.phi();//shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
    //evtInfo->pfMETCov00 = inPFMET.getSignificanceMatrix()(0,0);
    //evtInfo->pfMETCov01 = inPFMET.getSignificanceMatrix()(0,1);
    //evtInfo->pfMETCov11 = inPFMET.getSignificanceMatrix()(1,1);
    
    // Corrected PF MET
    edm::Handle<reco::PFMETCollection> hPFMETCProduct;
    iEvent.getByLabel(fPFMETCName,hPFMETCProduct);
    assert(hPFMETCProduct.isValid());
    const reco::PFMET &inPFMETC = hPFMETCProduct.product()->front();
    evtInfo->pfMETC      = inPFMETC.pt();
    evtInfo->pfMETCphi   = inPFMETC.phi();
    evtInfo->pfMETCCov00 = inPFMETC.getSignificanceMatrix()(0,0);
    evtInfo->pfMETCCov01 = inPFMETC.getSignificanceMatrix()(0,1);
    evtInfo->pfMETCCov11 = inPFMETC.getSignificanceMatrix()(1,1);

    // MVA MET
    edm::Handle<reco::PFMETCollection> hMVAMETProduct;
    iEvent.getByLabel(fMVAMETName,hMVAMETProduct);
    assert(hMVAMETProduct.isValid());
    const reco::PFMET &inMVAMET = hMVAMETProduct.product()->front();
    evtInfo->mvaMET      = inMVAMET.pt();
    evtInfo->mvaMETphi   = inMVAMET.phi();
    //    evtInfo->mvaMETCov00 = inMVAMET.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMETCov01 = inMVAMET.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMETCov11 = inMVAMET.getSignificanceMatrix()(1,1);

    // MVA MET with unity response
    //edm::Handle<reco::PFMETCollection> hMVAMETUProduct;
    //iEvent.getByLabel(fMVAMETUName,hMVAMETUProduct);
    //assert(hMVAMETUProduct.isValid());
    //const reco::PFMET &inMVAMETU = hMVAMETUProduct.product()->front();
    //evtInfo->mvaMETU      = inMVAMETU.pt();
    //evtInfo->mvaMETUphi   = inMVAMETU.phi();
    //    evtInfo->mvaMETUCov00 = inMVAMETU.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMETUCov01 = inMVAMETU.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMETUCov11 = inMVAMETU.getSignificanceMatrix()(1,1);

    // MVA MET without jet smearing (relevant only for MC)
    //edm::Handle<reco::PFMETCollection> hMVAMET0Product;
    //iEvent.getByLabel(fMVAMET0Name,hMVAMET0Product);
    //assert(hMVAMET0Product.isValid());
    //const reco::PFMET &inMVAMET0 = hMVAMET0Product.product()->front();
    //evtInfo->mvaMET0      = inMVAMET0.pt();
    //evtInfo->mvaMET0phi   = inMVAMET0.phi();
    //    evtInfo->mvaMET0Cov00 = inMVAMET0.getSignificanceMatrix()(0,0);
    //    evtInfo->mvaMET0Cov01 = inMVAMET0.getSignificanceMatrix()(0,1);
    //    evtInfo->mvaMET0Cov11 = inMVAMET0.getSignificanceMatrix()(1,1);

    // PUPPI MET
    edm::Handle<reco::PFMETCollection> hPuppET;
    iEvent.getByLabel(fPUPPETName,hPuppET);
    assert(hPuppET.isValid());
    const reco::PFMET &inPuppET = hPuppET.product()->front();
    evtInfo->puppET      = inPuppET.pt();
    evtInfo->puppETphi   = inPuppET.phi();
    //evtInfo->puppETCov00 = inPuppET.getSignificanceMatrix()(0,0);
    //evtInfo->puppETCov01 = inPuppET.getSignificanceMatrix()(0,1);
    //evtInfo->puppETCov11 = inPuppET.getSignificanceMatrix()(1,1);

    //Type1 PUPPI MET
    edm::Handle<reco::PFMETCollection> hPuppETC;
    iEvent.getByLabel(fPUPPETCName,hPuppETC);
    assert(hPuppETC.isValid());
    const reco::PFMET &inPuppETC = hPuppETC.product()->front();
    evtInfo->puppETC      = inPuppETC.pt();
    evtInfo->puppETCphi   = inPuppETC.phi();
    //evtInfo->puppETCCov00 = inPuppETC.getSignificanceMatrix()(0,0);
    //evtInfo->puppETCCov01 = inPuppETC.getSignificanceMatrix()(0,1);
    //evtInfo->puppETCCov11 = inPuppETC.getSignificanceMatrix()(1,1);

    // Corrected PF MET
    if(fPFMET30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPFMET30Product;
      iEvent.getByLabel(fPFMET30Name,hPFMET30Product);
      assert(hPFMET30Product.isValid());
      const reco::PFMET &inPFMET30 = hPFMET30Product.product()->front();
      evtInfo->pfMET30      = inPFMET30.pt();
      evtInfo->pfMET30phi   = inPFMET30.phi();
    }
    // Corrected PF MET
    if(fPFMETC30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPFMETC30Product;
      iEvent.getByLabel(fPFMETC30Name,hPFMETC30Product);
      assert(hPFMETC30Product.isValid());
      const reco::PFMET &inPFMETC30 = hPFMETC30Product.product()->front();
      evtInfo->pfMETC30      = inPFMETC30.pt();
      evtInfo->pfMETC30phi   = inPFMETC30.phi();
    }
    // MVA MET
    if(fMVAMET30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hMVAMET30Product;
      iEvent.getByLabel(fMVAMET30Name,hMVAMET30Product);
      assert(hMVAMET30Product.isValid());
      const reco::PFMET &inMVAMET30 = hMVAMET30Product.product()->front();
      evtInfo->mvaMET30      = inMVAMET30.pt();
      evtInfo->mvaMET30phi   = inMVAMET30.phi();
    }
    // PUPPI MET
    if(fPUPPET30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPuppET30;
      iEvent.getByLabel(fPUPPET30Name,hPuppET30);
      assert(hPuppET30.isValid());
      const reco::PFMET &inPuppET30 = hPuppET30.product()->front();
      evtInfo->puppET30      = inPuppET30.pt();
      evtInfo->puppET30phi   = inPuppET30.phi();
    }
    // Type 1 PUPPI MET
    if(fPUPPETC30Name.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPuppETC30;
      iEvent.getByLabel(fPUPPETC30Name,hPuppETC30);
      assert(hPuppETC30.isValid());
      const reco::PFMET &inPuppETC30 = hPuppETC30.product()->front();
      evtInfo->puppETC30      = inPuppETC30.pt();
      evtInfo->puppETC30phi   = inPuppETC30.phi();
    }
    // Alpaca MET
    if(fALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hAlpacaMET;
      iEvent.getByLabel(fALPACAMETName,hAlpacaMET);
      assert(hAlpacaMET.isValid());
      const reco::PFMET &inAlpacaMET = hAlpacaMET.product()->front();
      evtInfo->alpacaMET      = inAlpacaMET.pt();
      evtInfo->alpacaMETphi   = inAlpacaMET.phi();
    }
    if(fPALPACAMETName.size() > 0) { 
      edm::Handle<reco::PFMETCollection> hPAlpacaMET;
      iEvent.getByLabel(fPALPACAMETName,hPAlpacaMET);
      assert(hPAlpacaMET.isValid());
      const reco::PFMET &inPAlpacaMET = hPAlpacaMET.product()->front();
      evtInfo->pcpMET      = inPAlpacaMET.pt();
      evtInfo->pcpMETphi   = inPAlpacaMET.phi();
    }
    // Track MET
    edm::Handle<pat::PackedCandidateCollection> hPFCandProduct;
    iEvent.getByLabel(fPFCandName,hPFCandProduct);
    assert(hPFCandProduct.isValid());
    const pat::PackedCandidateCollection *pfCandCol = hPFCandProduct.product();
    computeTrackMET(pfCandCol, evtInfo->trkMET, evtInfo->trkMETphi);
  }
  
  
  //
  // event energy density
  //==============================
  
  // Rho for isolation correction
  edm::Handle<double> hRhoIso;
  edm::InputTag rhoIsoTag(fRhoIsoName,"");
  iEvent.getByLabel(rhoIsoTag,hRhoIso);
  assert(hRhoIso.isValid());
  evtInfo->rhoIso = *hRhoIso;
  
  // Rho for jet energy correction
  edm::Handle<double> hRhoJet;
  edm::InputTag rhoJetTag(fRhoJetName,"");
  iEvent.getByLabel(rhoJetTag,hRhoJet);
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
