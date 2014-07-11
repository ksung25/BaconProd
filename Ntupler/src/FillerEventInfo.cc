#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TSusyGen.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include <TLorentzVector.h>
#include <string>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerEventInfo::FillerEventInfo(const edm::ParameterSet &iConfig):
  fPFCandName (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fPUInfoName (iConfig.getUntrackedParameter<std::string>("edmPileupInfoName","addPileupInfo")),
  fBSName     (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot")),
  fPFMETName  (iConfig.getUntrackedParameter<std::string>("edmPFMETName","pfMet")),
  fPFMETCName (iConfig.getUntrackedParameter<std::string>("edmPFMETCorrName","pfType0p1CorrectedMet")),
  fMVAMETName (iConfig.getUntrackedParameter<std::string>("edmMVAMETName","pfMEtMVA")),
  fMVAMETUName(iConfig.getUntrackedParameter<std::string>("edmMVAMETUnityName","pfMEtMVAUnity")),
  fMVAMET0Name(iConfig.getUntrackedParameter<std::string>("edmMVAMETNoSmearName","pfMEtMVANoSmear")),
  fRhoIsoName (iConfig.getUntrackedParameter<std::string>("edmRhoForIsoName","kt6PFJets")),
  fRhoJetName (iConfig.getUntrackedParameter<std::string>("edmRhoForJetEnergy","kt6PFJets")),
  fFillMET    (iConfig.getUntrackedParameter<bool>("doFillMET",true)),
  fFillMETFilters(iConfig.getUntrackedParameter<bool>("doFillMETFilters",true)),
  fAddSusyGen (iConfig.getUntrackedParameter<bool>("addSusyGen",false))
{}

//--------------------------------------------------------------------------------------------------
FillerEventInfo::~FillerEventInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerEventInfo::fill(TEventInfo *evtInfo,
                           const edm::Event &iEvent, const reco::Vertex &pv, const bool hasGoodPV,
			   const TriggerBits triggerBits,TSusyGen *iSusyGen)
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
  	
  
  evtInfo->metFilterFailBits=0;

  if(fAddSusyGen) { 
    edm::Handle<LHEEventProduct> comments;
    iEvent.getByLabel("source", comments);
    
    std::string pSusy;
    int lId = 0;
    for(LHEEventProduct::comments_const_iterator pComment = comments->comments_begin(); pComment != comments->comments_end(); pComment++) { 
      if(lId == 1) pSusy = *pComment;
      lId++;
    }
    std::string delimeter  = "_";
    std::string delimeter1 = " ";
    int lSpace = pSusy.find(delimeter1,2)+1;
    std::string pSubSusy   = pSusy.substr(lSpace,pSusy.find(delimeter1,lSpace+2)-lSpace);
    lSpace     = pSubSusy.find(delimeter)+1;
    int m1 = atoi((pSubSusy.substr(lSpace,(pSubSusy.rfind(delimeter)-lSpace))).c_str());
    int m2 = atoi((pSubSusy.substr(pSubSusy.rfind(delimeter)+1)).c_str());
    iSusyGen->id = pSusy;
    iSusyGen->m1 = m1;
    iSusyGen->m2 = m2;
  }
  
  if(fFillMETFilters) { 
    //
    // MET filter tags
    //==============================
    
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
  }
   
  if(fFillMET) { 
    //
    // MET info
    //==============================
    
    // PF MET
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
    
    // MVA MET with unity response
    edm::Handle<reco::PFMETCollection> hMVAMETUProduct;
    iEvent.getByLabel(fMVAMETUName,hMVAMETUProduct);
    assert(hMVAMETUProduct.isValid());
    const reco::PFMET &inMVAMETU = hMVAMETUProduct.product()->front();
    evtInfo->mvaMETU      = inMVAMETU.pt();
    evtInfo->mvaMETUphi   = inMVAMETU.phi();
    evtInfo->mvaMETUCov00 = inMVAMETU.getSignificanceMatrix()(0,0);
    evtInfo->mvaMETUCov01 = inMVAMETU.getSignificanceMatrix()(0,1);
    evtInfo->mvaMETUCov11 = inMVAMETU.getSignificanceMatrix()(1,1);
    
    // MVA MET without jet smearing (relevant only for MC)
    edm::Handle<reco::PFMETCollection> hMVAMET0Product;
    iEvent.getByLabel(fMVAMET0Name,hMVAMET0Product);
    assert(hMVAMET0Product.isValid());
    const reco::PFMET &inMVAMET0 = hMVAMET0Product.product()->front();
    evtInfo->mvaMET0      = inMVAMET0.pt();
    evtInfo->mvaMET0phi   = inMVAMET0.phi();
    evtInfo->mvaMET0Cov00 = inMVAMET0.getSignificanceMatrix()(0,0);
    evtInfo->mvaMET0Cov01 = inMVAMET0.getSignificanceMatrix()(0,1);
    evtInfo->mvaMET0Cov11 = inMVAMET0.getSignificanceMatrix()(1,1);
    
    // Track MET
    edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
    iEvent.getByLabel(fPFCandName,hPFCandProduct);
    assert(hPFCandProduct.isValid());
    const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
    computeTrackMET(pv, pfCandCol, evtInfo->trkMET, evtInfo->trkMETphi);    
  }
  
  
  //
  // event energy density
  //==============================
  
  // Rho for isolation correction
  edm::Handle<double> hRhoIso;
  edm::InputTag rhoIsoTag(fRhoIsoName,"rho");
  iEvent.getByLabel(rhoIsoTag,hRhoIso);
  assert(hRhoIso.isValid());
  evtInfo->rhoIso = *hRhoIso;
  
  // Rho for jet energy correction
  edm::Handle<double> hRhoJet;
  edm::InputTag rhoJetTag(fRhoJetName,"rho");
  iEvent.getByLabel(rhoJetTag,hRhoJet);
  assert(hRhoJet.isValid());
  evtInfo->rhoJet = *hRhoJet;


  //
  // fired triggers
  //==============================
  evtInfo->triggerBits = triggerBits;
}


//--------------------------------------------------------------------------------------------------
void FillerEventInfo::computeTrackMET(const reco::Vertex &pv, const reco::PFCandidateCollection *pfCandCol,
                                      float &out_met, float &out_metphi)
{  
  out_met    = 0;
  out_metphi = 0;


  
  double metx=0, mety=0;
  for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
    if(itPF->trackRef().isNonnull() && pv.trackWeight(itPF->trackRef())>0) {
      metx  -= itPF->px();
      mety  -= itPF->py();
    }
  }
  
  TLorentzVector met;
  met.SetPxPyPzE(metx,mety,0,0);
  out_met    = met.Pt();
  out_metphi = met.Phi();
}
