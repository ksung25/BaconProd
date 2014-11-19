#include "ExpertMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenWeight.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"

#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "BaconProd/Ntupler/interface/FillerPF.hh"

// tools to parse HLT name patterns
#include <boost/foreach.hpp>
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// data format classes
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>


//--------------------------------------------------------------------------------------------------
ExpertMod::ExpertMod(const edm::ParameterSet &iConfig):
  fSkipOnHLTFail  (iConfig.getUntrackedParameter<bool>("skipOnHLTFail",false)),
  fHLTTag         ("TriggerResults","","HLT"),
  fHLTObjTag      ("hltTriggerSummaryAOD","","HLT"),
  fHLTFile        (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fFillerEvtInfo  (0),
  fFillerGenInfo  (0),
  fFillerPV       (0),
  fFillerPF       (0),
  fTrigger        (0),
  fIsActiveEvtInfo(false),
  fIsActiveGenInfo(false),
  fIsActivePV     (false),
  fIsActivePF     (false),
  fOutputName     (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile     (0),
  fTotalEvents    (0),
  fEventTree      (0),
  fEvtInfo        (0),
  fGenEvtInfo     (0),
  fGenWeight      (0),
  fGenParArr      (0),
  fPFParArr       (0),
  fPVArr          (0)
{
  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenWeight::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  baconhep::TPFPart::Class()->IgnoreTObjectStreamer();
  
  //
  // Set up bacon objects and configure fillers
  // 
  if(iConfig.existsAs<edm::ParameterSet>("Info",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Info"));
    fIsActiveEvtInfo = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveEvtInfo) {
      fEvtInfo       = new baconhep::TEventInfo();         assert(fEvtInfo);
      fFillerEvtInfo = new baconhep::FillerEventInfo(cfg); assert(fFillerEvtInfo);
    }
  }
  
  if(iConfig.existsAs<edm::ParameterSet>("GenInfo",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("GenInfo"));
    fIsActiveGenInfo = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveGenInfo) {
      fGenEvtInfo    = new baconhep::TGenEventInfo();                   assert(fGenEvtInfo);
      fGenWeight     = new baconhep::TGenWeight();                      assert(fGenWeight);
      fGenParArr     = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
      fFillerGenInfo = new baconhep::FillerGenInfo(cfg);                assert(fFillerGenInfo);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("PFCand",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PFCand"));
    fIsActivePF = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActivePF) {
      fPFParArr = new TClonesArray("baconhep::TPFPart",5000); assert(fPFParArr);
      fFillerPF = new baconhep::FillerPF(cfg);                assert(fFillerPF);
    }
  }
  
  if(iConfig.existsAs<edm::ParameterSet>("PV",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PV"));
    fIsActivePV = cfg.getUntrackedParameter<bool>("isActive");

    // create array and filler even if vertices won't be saved to output (i.e. fIsActivePV == false),
    // because FillerVertex::fill(...) is used to find the event primary vertex
    // (not elegant, but I suppose a dedicated PV finding function can be implemented somewhere...)
    fPVArr    = new TClonesArray("baconhep::TVertex",5000); assert(fPVArr);
    fFillerPV = new baconhep::FillerVertex(cfg);            assert(fFillerPV);
  }  
}

//--------------------------------------------------------------------------------------------------
ExpertMod::~ExpertMod()
{
  delete fFillerEvtInfo;
  delete fFillerGenInfo;
  delete fFillerPV;
  delete fFillerPF;

  delete fEvtInfo;
  delete fGenEvtInfo;
  delete fGenWeight;
  delete fGenParArr;
  delete fPFParArr;
  delete fPVArr;
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::beginJob()
{
  //
  // Create output file, trees, and histograms
  //
  fOutputFile  = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fEventTree   = new TTree("Events","Events");
  
  if(fIsActiveEvtInfo) { fEventTree->Branch("Info",fEvtInfo); }
  if(fIsActiveGenInfo) {
    fEventTree->Branch("GenEvtInfo" ,fGenEvtInfo);
    fEventTree->Branch("GenWeight"  ,fGenWeight);
    fEventTree->Branch("GenParticle",&fGenParArr);
  }
  if(fIsActivePF) { fEventTree->Branch("PFPart",&fPFParArr); }
  if(fIsActivePV) { fEventTree->Branch("PV",    &fPVArr); }

  //
  // Triggers
  //
  setTriggers();
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::endJob() 
{
  //
  // Save to ROOT file
  //
  //fEventTree->Print();
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::setTriggers()
{
  std::string cmssw_base_src = getenv("CMSSW_BASE"); cmssw_base_src+="/src/";
  fTrigger = new baconhep::TTrigger(cmssw_base_src + fHLTFile);
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fTotalEvents->Fill(1);
  
  edm::Handle<edm::TriggerResults> hTrgRes;
  iEvent.getByLabel(fHLTTag,hTrgRes);
  assert(hTrgRes.isValid());  
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);
  Bool_t config_changed = false;
  if(fTriggerNamesID != triggerNames.parameterSetID()) {
    fTriggerNamesID = triggerNames.parameterSetID();
    config_changed  = true;
  }
  if(config_changed) {
    initHLT(*hTrgRes, triggerNames);
  }
  TriggerBits triggerBits;
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    if(fTrigger->fRecords[irec].hltPathIndex == (unsigned int)-1) continue;
    if(hTrgRes->accept(fTrigger->fRecords[irec].hltPathIndex)) {
      triggerBits [fTrigger->fRecords[irec].baconTrigBit] = 1;
    }
  }
  if(fSkipOnHLTFail && triggerBits == 0) return;  

  if(fIsActiveGenInfo) {
    fGenParArr->Clear();
    float lTmp = 1.;
    fFillerGenInfo->fill(fGenEvtInfo, fGenWeight, fGenParArr, iEvent,lTmp);
  }
  fPVArr->Clear();
  int nvertices = 0;
  const reco::Vertex *pv = fFillerPV->fill(fPVArr, nvertices, iEvent);
  assert(pv);
  
//  separatePileUp(iEvent, *pv);
  
  if(fIsActiveEvtInfo) {
    fFillerEvtInfo->fill(fEvtInfo, iEvent, *pv, (nvertices>0), triggerBits);
  }
  
  if(fIsActivePF) {
    fPFParArr->Clear();  
    fFillerPF->fill(fPFParArr,fPVArr,iEvent);
  }
  
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames)
{
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    fTrigger->fRecords[irec].hltPathName  = "";
    fTrigger->fRecords[irec].hltPathIndex = (unsigned int)-1;
    const std::string pattern = fTrigger->fRecords[irec].hltPattern;
    if(edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
      std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
      if(matches.empty()) {
        std::cout << "requested pattern [" << pattern << "] does not match any HLT paths" << std::endl;
      } else {
        BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) {
          fTrigger->fRecords[irec].hltPathName = *match;
       }
      }
    } else {  // take full HLT path name given
      fTrigger->fRecords[irec].hltPathName = pattern;
    }
    // Retrieve index in trigger menu corresponding to HLT path
    unsigned int index = triggerNames.triggerIndex(fTrigger->fRecords[irec].hltPathName);
    if(index < result.size()) {  // check for valid index
      fTrigger->fRecords[irec].hltPathIndex = index;
    }
  }
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void ExpertMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}
void ExpertMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void ExpertMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void ExpertMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(ExpertMod);
