#include "NtuplerMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TCaloJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TRHPart.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Ntupler/interface/FillerPhoton.hh"
#include "BaconProd/Ntupler/interface/FillerTau.hh"
#include "BaconProd/Ntupler/interface/FillerCaloJet.hh"
#include "BaconProd/Ntupler/interface/FillerJet.hh"
#include "BaconProd/Ntupler/interface/FillerGenJets.hh"
#include "BaconProd/Ntupler/interface/FillerPF.hh"
#include "BaconProd/Ntupler/interface/FillerRH.hh"

// tools to parse HLT name patterns
#include <boost/foreach.hpp>
#include "FWCore/Utilities/interface/RegexMatch.h"

// data format classes
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>


//--------------------------------------------------------------------------------------------------
NtuplerMod::NtuplerMod(const edm::ParameterSet &iConfig):
  fSkipOnHLTFail     (iConfig.getUntrackedParameter<bool>("skipOnHLTFail",false)),
  fUseAOD            (iConfig.getUntrackedParameter<bool>("useAOD",true)),
  fHLTTag            ("TriggerResults","","HLT"),
  fHLTObjTag         (iConfig.getUntrackedParameter<std::string>("TriggerObject","hltTriggerSummaryAOD")),
  fHLTFile           (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fPVName            (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fGenRunInfoName    (iConfig.getUntrackedParameter<std::string>("edmGenRunInfoName","generator")),
  fFillLHEWgt        (false),
  fUseAODJet           (true),
  fUseAODFatJet        (true),
  fUseAODFatterJet     (true),
  fUseAODPuppiJet      (true),
  fUseAODFatPuppiJet   (true),
  fUseAODFatterPuppiJet(true),
  fComputeFullJetInfo(false),
  fComputeFullFatJetInfo(false),
  fComputeFullFatterJetInfo(false),
  fComputeFullPuppiJetInfo(false),
  fComputeFullFatPuppiJetInfo(false),
  fComputeFullFatterPuppiJetInfo(false),
  fFillerEvtInfo     (0),
  fFillerGenInfo     (0),
  fFillerGenJet      (0),
  fFillerGenFatJet   (0),
  fFillerPV          (0),
  fFillerEle         (0),
  fFillerMuon        (0),
  fFillerPhoton      (0),
  fFillerTau         (0),
  //fFillerCaloJet     (0),
  fFillerJet         (0),
  fFillerFatJet      (0),
  fFillerFatterJet   (0),
  fFillerPuppiJet      (0),
  fFillerFatPuppiJet   (0),
  fFillerFatterPuppiJet(0),
  fFillerPF          (0),
  fFillerRH          (0),
  fTrigger           (0),
  fIsActiveEvtInfo   (false),
  fIsActiveGenInfo   (false),
  fIsActiveGenJet    (false),
  fIsActiveGenFatJet (false),
  fIsActivePV        (false),
  fIsActiveEle       (false),
  fIsActiveMuon      (false),
  fIsActivePhoton    (false),
  fIsActiveTau       (false),
  fIsActiveJet       (false),
  fIsActiveFatJet    (false),
  fIsActiveFatterJet (false),
  fIsActivePuppiJet       (false),
  fIsActiveFatPuppiJet    (false),
  fIsActiveFatterPuppiJet (false),
  fIsActivePF        (false),
  fIsActiveRH        (false),
  fUseTrigger        (false),
  fOutputName        (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile        (0),
  fTotalEvents       (0),
  fXS                (0),
  fEventTree         (0),
  fEvtInfo           (0),
  fGenEvtInfo        (0),
  fLHEWgtArr         (0),
  fGenParArr         (0),
  fGenJetArr         (0),
  fGenFatJetArr      (0),
  fEleArr            (0),
  fMuonArr           (0),
  fTauArr            (0),
  fJetArr            (0),
  fFatJetArr         (0),
  fFatterJetArr      (0),
  fPuppiJetArr       (0),
  fFatPuppiJetArr    (0),
  fFatterPuppiJetArr (0),
  fPhotonArr         (0),
  fPVArr             (0),
  fAddJetArr         (0),
  fAddFatJetArr      (0),
  fAddFatterJetArr   (0),
  fAddPuppiJetArr         (0),
  fAddFatPuppiJetArr      (0),
  fAddFatterPuppiJetArr   (0),
  fPFParArr          (0),
  fRHParArr          (0)
{
  fUseTrigger          = iConfig.getUntrackedParameter<bool>("useTrigger",true);
  fTokGenRunInfo       = consumes<GenRunInfoProduct>(edm::InputTag(fGenRunInfoName)); 
  fTokTrgRes           = consumes<edm::TriggerResults>(edm::InputTag(fHLTTag)); 
  fTokTrgEvt           = consumes<trigger::TriggerEvent>(edm::InputTag(fHLTTag)); 
  fTokTrgObj           = consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag(fHLTObjTag)); 
  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TLHEWeight::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  baconhep::TGenJet::Class()->IgnoreTObjectStreamer();
  baconhep::TMuon::Class()->IgnoreTObjectStreamer();
  baconhep::TElectron::Class()->IgnoreTObjectStreamer();
  baconhep::TTau::Class()->IgnoreTObjectStreamer();
  baconhep::TJet::Class()->IgnoreTObjectStreamer();
  baconhep::TPhoton::Class()->IgnoreTObjectStreamer();
  baconhep::TVertex::Class()->IgnoreTObjectStreamer();
  baconhep::TAddJet::Class()->IgnoreTObjectStreamer();
  baconhep::TPFPart::Class()->IgnoreTObjectStreamer();
  baconhep::TRHPart::Class()->IgnoreTObjectStreamer();

  // trigger object information
  if(iConfig.existsAs<std::string>("TriggerFile",false) && fUseTrigger) {
    fUseTrigger = true;
    if(fUseAOD) {
      fHLTObjTag = edm::InputTag("hltTriggerSummaryAOD","","HLT");
    } else {
      //fHLTObjTag = edm::InputTag("selectedPatTrigger","","PAT");
      fHLTObjTag = edm::InputTag("selectedPatTrigger");
    }
  }

  //
  // Set up bacon objects and configure fillers
  // 
  if(iConfig.existsAs<edm::ParameterSet>("Info",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Info"));
    fIsActiveEvtInfo = cfg.getUntrackedParameter<bool>("isActive");

    if(fIsActiveEvtInfo) {
      fEvtInfo       = new baconhep::TEventInfo();                  assert(fEvtInfo);
      fFillerEvtInfo = new baconhep::FillerEventInfo(cfg, fUseAOD,consumesCollector()); assert(fFillerEvtInfo);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("GenInfo",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("GenInfo"));
    fIsActiveGenInfo = cfg.getUntrackedParameter<bool>("isActive");

    fFillLHEWgt = cfg.getUntrackedParameter<bool>("fillLHEWeights");
    if(fIsActiveGenInfo) {
      fGenEvtInfo = new baconhep::TGenEventInfo();              assert(fGenEvtInfo);
      fGenParArr  = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
      if(fFillLHEWgt) {
        fLHEWgtArr = new TClonesArray("baconhep::TLHEWeight",5000); assert(fLHEWgtArr);
      }
      fFillerGenInfo = new baconhep::FillerGenInfo(cfg,consumesCollector()); assert(fFillerGenInfo);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("PV",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PV"));
    fIsActivePV = cfg.getUntrackedParameter<bool>("isActive");
    
    // create array and filler even if vertices won't be saved to output (i.e. fIsActivePV == false),
    // because FillerVertex::fill(...) is used to find the event primary vertex
    // (not elegant, but I suppose a dedicated PV finding function can be implemented somewhere...)
    fPVArr    = new TClonesArray("baconhep::TVertex");    assert(fPVArr);
    fFillerPV = new baconhep::FillerVertex(cfg, fUseAOD,consumesCollector()); assert(fFillerPV);
  }
    
  if(iConfig.existsAs<edm::ParameterSet>("Electron",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Electron"));
    fIsActiveEle = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveEle) {
      fEleArr    = new TClonesArray("baconhep::TElectron");    assert(fEleArr);
      fFillerEle = new baconhep::FillerElectron(cfg, fUseAOD,consumesCollector()); assert(fFillerEle);
    }
  }  

  if(iConfig.existsAs<edm::ParameterSet>("Muon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Muon"));
    fIsActiveMuon = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveMuon) {
      fMuonArr    = new TClonesArray("baconhep::TMuon");    assert(fMuonArr);
      fFillerMuon = new baconhep::FillerMuon(cfg, fUseAOD,consumesCollector()); assert(fFillerMuon);
    }
  }  

  if(iConfig.existsAs<edm::ParameterSet>("Photon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Photon"));
    fIsActivePhoton = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActivePhoton) {
      fPhotonArr    = new TClonesArray("baconhep::TPhoton");    assert(fPhotonArr);
      fFillerPhoton = new baconhep::FillerPhoton(cfg, fUseAOD,consumesCollector()); assert(fFillerPhoton);
    }
  } 

  if(iConfig.existsAs<edm::ParameterSet>("Tau",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Tau"));
    fIsActiveTau = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveTau) {
      fTauArr    = new TClonesArray("baconhep::TTau");    assert(fTauArr);
      fFillerTau = new baconhep::FillerTau(cfg, fUseAOD,consumesCollector()); assert(fFillerTau);
    }
  }
  //!!!! Temporary Calo Jets
  //if(iConfig.existsAs<edm::ParameterSet>("AK4Calo",false)) {
  //  edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("AK4Calo"));
  //  fIsActiveCaloJet = cfg.getUntrackedParameter<bool>("isActive");
  //  if(fIsActiveCaloJet) {
  //    fFillerCaloJet = new baconhep::FillerCaloJet(cfg); assert(fFillerCaloJet);
  //    fCaloJetArr = new TClonesArray("baconhep::TCaloJet");       assert(fCaloJetArr);
  //  }
  //}
  if(iConfig.existsAs<edm::ParameterSet>("AK4CHS",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("AK4CHS"));
    fIsActiveJet = cfg.getUntrackedParameter<bool>("isActive");
    fUseAODJet   = cfg.getUntrackedParameter<bool>("useAOD");
    fComputeFullJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActiveJet) {
      fFillerJet = new baconhep::FillerJet(cfg, fUseAODJet,consumesCollector()); assert(fFillerJet);
      fJetArr = new TClonesArray("baconhep::TJet");       assert(fJetArr);
      if(fComputeFullJetInfo) {
        fAddJetArr = new TClonesArray("baconhep::TAddJet"); assert(fAddJetArr);
      }
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("AK8CHS",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("AK8CHS"));
    fIsActiveFatJet = cfg.getUntrackedParameter<bool>("isActive");
    fUseAODFatJet   = cfg.getUntrackedParameter<bool>("useAOD");
    fComputeFullFatJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActiveFatJet) {
      fFillerFatJet = new baconhep::FillerJet(cfg, fUseAODFatJet,consumesCollector()); assert(fFillerFatJet);
      fFatJetArr = new TClonesArray("baconhep::TJet");             assert(fFatJetArr);
      if(fComputeFullFatJetInfo) {
        fAddFatJetArr = new TClonesArray("baconhep::TAddJet"); assert(fAddFatJetArr);
      }
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("CA15CHS",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("CA15CHS"));
    fIsActiveFatterJet = cfg.getUntrackedParameter<bool>("isActive");
    fUseAODFatterJet   = cfg.getUntrackedParameter<bool>("useAOD");
    fComputeFullFatterJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActiveFatterJet) {
      fFillerFatterJet = new baconhep::FillerJet(cfg, fUseAODFatterJet,consumesCollector()); assert(fFillerFatterJet);
      fFatterJetArr = new TClonesArray("baconhep::TJet");                assert(fFatterJetArr);
      if(fComputeFullFatterJetInfo) {
        fAddFatterJetArr = new TClonesArray("baconhep::TAddJet"); assert(fAddFatterJetArr);
      }
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("AK4Puppi",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("AK4Puppi"));
    fIsActivePuppiJet = cfg.getUntrackedParameter<bool>("isActive");
    fUseAODPuppiJet   = cfg.getUntrackedParameter<bool>("useAOD");    
    fComputeFullPuppiJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActivePuppiJet) {
      fFillerPuppiJet = new baconhep::FillerJet(cfg, fUseAODPuppiJet,consumesCollector()); assert(fFillerPuppiJet);
      fPuppiJetArr = new TClonesArray("baconhep::TJet");               assert(fPuppiJetArr);
      if(fComputeFullPuppiJetInfo) {
        fAddPuppiJetArr = new TClonesArray("baconhep::TAddJet"); assert(fAddPuppiJetArr);
      }
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("CA8Puppi",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("CA8Puppi"));
    fIsActiveFatPuppiJet = cfg.getUntrackedParameter<bool>("isActive");
    fUseAODFatPuppiJet   = cfg.getUntrackedParameter<bool>("useAOD");
    fComputeFullFatPuppiJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActiveFatPuppiJet) {
      fFillerFatPuppiJet = new baconhep::FillerJet(cfg, fUseAODFatPuppiJet,consumesCollector()); assert(fFillerFatPuppiJet);
      fFatPuppiJetArr = new TClonesArray("baconhep::TJet");                 assert(fFatPuppiJetArr);
      if(fComputeFullFatPuppiJetInfo) {
        fAddFatPuppiJetArr = new TClonesArray("baconhep::TAddJet"); assert(fAddFatPuppiJetArr);
      }
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("CA15Puppi",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("CA15Puppi"));
    fIsActiveFatterPuppiJet = cfg.getUntrackedParameter<bool>("isActive");
    fUseAODFatterPuppiJet   = cfg.getUntrackedParameter<bool>("useAOD");
    fComputeFullFatterPuppiJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActiveFatterPuppiJet) {
      fFillerFatterPuppiJet = new baconhep::FillerJet(cfg, fUseAODFatterPuppiJet,consumesCollector()); assert(fFillerFatterPuppiJet);
      fFatterPuppiJetArr    = new TClonesArray("baconhep::TJet");                  assert(fFatterPuppiJetArr);
      if(fComputeFullFatterPuppiJetInfo) {
        fAddFatterPuppiJetArr = new TClonesArray("baconhep::TAddJet"); assert(fAddFatterPuppiJetArr);
      }
    }
  }
  if(iConfig.existsAs<edm::ParameterSet>("PFCand",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PFCand"));
    fIsActivePF = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActivePF) {
      fPFParArr = new TClonesArray("baconhep::TPFPart",5000); assert(fPFParArr);
      fFillerPF = new baconhep::FillerPF(cfg,consumesCollector());                assert(fFillerPF);
    }
  } 

  if(iConfig.existsAs<edm::ParameterSet>("RecHit",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("RecHit"));
    fIsActiveRH = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveRH) {
      fRHParArr = new TClonesArray("baconhep::TRHPart",50000); assert(fRHParArr);
      fFillerRH = new baconhep::FillerRH(cfg,consumesCollector());                 assert(fFillerRH);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("GenJet",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("GenJet"));
    fIsActiveGenJet     = cfg.getUntrackedParameter<bool>("isActive");
    fIsActiveGenFatJet  = cfg.getUntrackedParameter<bool>("isActiveFatJet");
    if(fIsActiveGenJet) {
      fGenJetArr     = new TClonesArray("baconhep::TGenJet");                  assert(fGenJetArr);
      fFillerGenJet  = new baconhep::FillerGenJets(cfg,consumesCollector());  assert(fFillerGenJet);
    }
    if(fIsActiveGenFatJet) {
      fGenFatJetArr  = new TClonesArray("baconhep::TGenJet");           assert(fGenFatJetArr);
    }
  }
}

//--------------------------------------------------------------------------------------------------
NtuplerMod::~NtuplerMod()
{
  delete fFillerEvtInfo;
  delete fFillerGenInfo;
  delete fFillerPV;
  delete fFillerEle;
  delete fFillerMuon;
  delete fFillerPhoton;
  delete fFillerTau;
  //delete fFillerCaloJet;
  delete fFillerJet;
  delete fFillerFatJet;
  delete fFillerFatterJet;
  delete fFillerPuppiJet;
  delete fFillerFatPuppiJet;
  delete fFillerFatterPuppiJet;

  delete fFillerPF;
  delete fFillerRH;

  delete fTrigger;
  delete fEvtInfo;
  delete fGenEvtInfo;
  delete fLHEWgtArr;
  delete fGenParArr;
  delete fGenJetArr;
  delete fEleArr;
  delete fMuonArr;
  delete fTauArr;
  //delete fCaloJetArr;
  delete fJetArr;
  delete fFatJetArr;
  delete fFatterJetArr;
  delete fPuppiJetArr;
  delete fFatPuppiJetArr;
  delete fFatterPuppiJetArr;
  delete fPhotonArr;
  delete fPVArr;
  delete fPFParArr;
  delete fRHParArr;
}
//--------------------------------------------------------------------------------------------------
void NtuplerMod::respondToOpenInputFile(edm::FileBlock const&) 
{  
  if(fFillerJet            != 0) fFillerJet           ->initPUJetId();
  if(fFillerFatJet         != 0) fFillerFatJet        ->initPUJetId();
  if(fFillerFatterJet      != 0) fFillerFatterJet     ->initPUJetId();
  if(fFillerPuppiJet       != 0) fFillerPuppiJet      ->initPUJetId();
  if(fFillerFatPuppiJet    != 0) fFillerFatPuppiJet   ->initPUJetId();
  if(fFillerFatterPuppiJet != 0) fFillerFatterPuppiJet->initPUJetId();
}
//--------------------------------------------------------------------------------------------------
void NtuplerMod::beginJob()
{  
  //
  // Create output file, trees, and histograms
  //
  fOutputFile  = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fEventTree   = new TTree("Events","Events");
  if(fIsActiveEvtInfo) { fEventTree->Branch("Info",fEvtInfo); }
  if(fIsActiveGenInfo) {
    fEventTree->Branch("GenEvtInfo",fGenEvtInfo);
    fEventTree->Branch("GenParticle",&fGenParArr);
    if(fFillLHEWgt) { fEventTree->Branch("LHEWeight",&fLHEWgtArr); }
  }
  if(fIsActiveGenJet)    { fEventTree->Branch("GenJet"     ,&fGenJetArr);}
  if(fIsActiveGenFatJet) { fEventTree->Branch("GenFatJet"  ,&fGenFatJetArr);}
  if(fIsActiveEle)    { fEventTree->Branch("Electron", &fEleArr); }
  if(fIsActiveMuon)   { fEventTree->Branch("Muon",     &fMuonArr); }
  if(fIsActiveTau)    { fEventTree->Branch("Tau",      &fTauArr); }
  if(fIsActivePhoton) { fEventTree->Branch("Photon",   &fPhotonArr); }
  if(fIsActivePV)     { fEventTree->Branch("PV",       &fPVArr); }
  //if(fIsActiveCaloJet){ fEventTree->Branch("CaloJet",  &fCaloJetArr); }
  if(fIsActiveJet) {
    fEventTree->Branch("AK4CHS", &fJetArr);
    if(fComputeFullJetInfo) { fEventTree->Branch("AddAK4CHS", &fAddJetArr); }
  }
  if(fIsActiveFatJet) {
    fEventTree->Branch("AK8CHS", &fFatJetArr);
    if(fComputeFullFatJetInfo) { fEventTree->Branch("AddAK8CHS", &fAddFatJetArr); }
  }
  if(fIsActiveFatterJet) {
    fEventTree->Branch("CA15CHS", &fFatterJetArr);
    if(fComputeFullFatterJetInfo) { fEventTree->Branch("AddCA15CHS", &fAddFatterJetArr); }
  }
  if(fIsActivePuppiJet) {
    fEventTree->Branch("AK4Puppi", &fPuppiJetArr);
    if(fComputeFullPuppiJetInfo) { fEventTree->Branch("AddAK4Puppi", &fAddPuppiJetArr); }
  }
  if(fIsActiveFatPuppiJet) {
    fEventTree->Branch("CA8Puppi", &fFatPuppiJetArr);
    if(fComputeFullFatPuppiJetInfo) { fEventTree->Branch("AddCA8Puppi", &fAddFatPuppiJetArr); }
  }
  if(fIsActiveFatterPuppiJet) {
    fEventTree->Branch("CA15Puppi", &fFatterPuppiJetArr);
    if(fComputeFullFatterPuppiJetInfo) { fEventTree->Branch("AddCA15Puppi", &fAddFatterPuppiJetArr); }
  }
  if(fIsActivePF) { fEventTree->Branch("PFPart", &fPFParArr); }
  if(fIsActiveRH) { fEventTree->Branch("RHPart", &fRHParArr); }
  // Triggers
  if(fUseTrigger) setTriggers();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::endJob() 
{
  //
  // Save to ROOT file
  //
//  fEventTree->Print();
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::setTriggers()
{
  std::string cmssw_base_src = getenv("CMSSW_BASE"); cmssw_base_src+="/src/";
  fTrigger = new baconhep::TTrigger(cmssw_base_src + fHLTFile);
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  fTotalEvents->Fill(1);
  TriggerBits triggerBits;
  if(fUseTrigger) { 
    edm::Handle<edm::TriggerResults> hTrgRes;
    iEvent.getByToken(fTokTrgRes,hTrgRes);
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
    for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
      if(fTrigger->fRecords[irec].hltPathIndex == (unsigned int)-1) continue;
      if(hTrgRes->accept(fTrigger->fRecords[irec].hltPathIndex)) {
	triggerBits [fTrigger->fRecords[irec].baconTrigBit] = 1;
      }
    }
    if(fSkipOnHLTFail && triggerBits == 0) return;  
  }
  if(fIsActiveGenInfo) {
    fGenParArr->Clear();
    if(fFillLHEWgt) {
      fLHEWgtArr->Clear();
      fFillerGenInfo->fill(fGenEvtInfo, fGenParArr, fLHEWgtArr, iEvent, fXS);
    } else {
      fFillerGenInfo->fill(fGenEvtInfo, fGenParArr,          0, iEvent, fXS);
    }
  }
  if(fIsActiveGenJet) {
    fGenJetArr->Clear();
    if(fGenFatJetArr != 0) fGenFatJetArr->Clear();
    fFillerGenJet->fill(fGenJetArr, fGenFatJetArr, iEvent);
  }
  const reco::Vertex *pv = 0;
  int nvertices = 0;
  if(fIsActiveEvtInfo) {
    fPVArr->Clear();
    pv = fFillerPV->fill(fPVArr, nvertices, iEvent);
    assert(pv);
  }
  if(fIsActiveEvtInfo) {
    fFillerEvtInfo->fill(fEvtInfo, iEvent, *pv, (nvertices>0), triggerBits);//,fSusyGen);
  }
  edm::Handle<trigger::TriggerEvent> hTrgEvt;
  edm::Handle<pat::TriggerObjectStandAloneCollection> hTrgObjs;
  const trigger::TriggerEvent*                  hTrgEvtDummy      = 0; 
  const pat::TriggerObjectStandAloneCollection* hTrgObjsDummy     = 0; 
  if(fUseTrigger) { 
    iEvent.getByToken(fTokTrgEvt,hTrgEvt);
    iEvent.getByToken(fTokTrgObj,hTrgObjs);
    if(fUseAOD) {hTrgEvtDummy  = &(*hTrgEvt); }
    else        {hTrgObjsDummy = &(*hTrgObjs); }
  }
  if(fIsActiveEle) {
    fEleArr->Clear();
    if(fUseAOD) { fFillerEle->fill(fEleArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);  }
    else        { fFillerEle->fill(fEleArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
  }
  if(fIsActiveMuon) {
    fMuonArr->Clear();  
    if(fUseAOD) { fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);  }
    else        { fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
  }
  if(fIsActivePhoton) {
    fPhotonArr->Clear();
    if(fUseAOD) { fFillerPhoton->fill(fPhotonArr, iEvent, iSetup, fTrigger->fRecords, *hTrgEvt);  }
    else        { fFillerPhoton->fill(fPhotonArr, iEvent, iSetup, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
  }
  if(fIsActiveTau) {
    fTauArr->Clear();
    if(fUseAOD) { fFillerTau->fill(fTauArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);  }
    else        { fFillerTau->fill(fTauArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
  }
  //  if(fIsActiveCaloJet) {
  //  fCaloJetArr->Clear();
  //  if(fUseAOD) { fFillerCaloJet->fill(fCaloJetArr,iEvent, iSetup,fTrigger->fRecords, *hTrgEvt);  }
  //}
  if(fIsActiveJet) {
    fJetArr->Clear();
    if(fComputeFullJetInfo) {
      fAddJetArr->Clear();      
      if(fUseAODJet) { fFillerJet->fill(fJetArr, fAddJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy ,hTrgObjsDummy);  }
      else           { fFillerJet->fill(fJetArr, fAddJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
    } else {
      if(fUseAODJet) { fFillerJet->fill(fJetArr,          0, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy ,hTrgObjsDummy);  }
      else           { fFillerJet->fill(fJetArr,          0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
    }
  }
  if(fIsActiveFatJet) {
    fFatJetArr->Clear();
    if(fComputeFullFatJetInfo) {
      fAddFatJetArr->Clear();      
      if(fUseAODFatJet) { fFillerFatJet->fill(fFatJetArr, fAddFatJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy ,hTrgObjsDummy);  }
      else              { fFillerFatJet->fill(fFatJetArr, fAddFatJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
    } else {
      if(fUseAODFatJet) { fFillerFatJet->fill(fFatJetArr,             0, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy ,hTrgObjsDummy);  }
      else              { fFillerFatJet->fill(fFatJetArr,             0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
    }
  }
  if(fIsActiveFatterJet) {
    fFatterJetArr->Clear();
    if(fComputeFullFatterJetInfo) {
      fAddFatterJetArr->Clear();      
      if(fUseAODFatterJet) { fFillerFatterJet->fill(fFatterJetArr, fAddFatterJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy ,hTrgObjsDummy);  }
      else                 { fFillerFatterJet->fill(fFatterJetArr, fAddFatterJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); } 
    } else {
      if(fUseAODFatterJet) { fFillerFatterJet->fill(fFatterJetArr,                0, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy ,hTrgObjsDummy); }
      else                 { fFillerFatterJet->fill(fFatterJetArr,                0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs);  } 
    }
  }
  if(fIsActivePuppiJet) {
    fPuppiJetArr->Clear();
    if(fComputeFullPuppiJetInfo) {
      fAddPuppiJetArr->Clear();      
      if(fUseAODPuppiJet) { fFillerPuppiJet->fill(fPuppiJetArr, fAddPuppiJetArr, iEvent, iSetup, *pv, fTrigger->fRecords,  hTrgEvtDummy, hTrgObjsDummy);  }
      else                { fFillerPuppiJet->fill(fPuppiJetArr, fAddPuppiJetArr, iEvent, iSetup, *pv, fTrigger->fRecords,  *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
    } else {
      if(fUseAODPuppiJet) { fFillerPuppiJet->fill(fPuppiJetArr,          0, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy, hTrgObjsDummy);  }
      else                { fFillerPuppiJet->fill(fPuppiJetArr,          0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }  // (!) consolidate fillers for AOD and MINIAOD
    }
  }
  if(fIsActiveFatPuppiJet) {
    fFatPuppiJetArr->Clear();
    if(fComputeFullFatPuppiJetInfo) {
      fAddFatPuppiJetArr->Clear();      
      if(fUseAODFatPuppiJet) { fFillerFatPuppiJet->fill(fFatPuppiJetArr, fAddFatPuppiJetArr, iEvent, iSetup, *pv, fTrigger->fRecords,hTrgEvtDummy,hTrgObjsDummy);  }
      else                   { fFillerFatPuppiJet->fill(fFatPuppiJetArr, fAddFatPuppiJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); } 
    } else {
      if(fUseAODFatPuppiJet) { fFillerFatPuppiJet->fill(fFatPuppiJetArr,             0, iEvent, iSetup, *pv, fTrigger->fRecords, hTrgEvtDummy,hTrgObjsDummy);  }
      else                   { fFillerFatPuppiJet->fill(fFatPuppiJetArr,             0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }
    }
  }
  if(fIsActiveFatterPuppiJet) {
    fFatterPuppiJetArr->Clear();
    if(fComputeFullFatterPuppiJetInfo) {
      fAddFatterPuppiJetArr->Clear();      
      if(fUseAODFatterPuppiJet) { fFillerFatterPuppiJet->fill(fFatterPuppiJetArr, fAddFatterPuppiJetArr, iEvent, iSetup, *pv, fTrigger->fRecords,hTrgEvtDummy,hTrgObjsDummy);  }
      else                      { fFillerFatterPuppiJet->fill(fFatterPuppiJetArr, fAddFatterPuppiJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); } 
    } else {
      if(fUseAODFatterPuppiJet) { fFillerFatterPuppiJet->fill(fFatterPuppiJetArr,                     0, iEvent, iSetup, *pv, fTrigger->fRecords,hTrgEvtDummy,hTrgObjsDummy);  }
      else                      { fFillerFatterPuppiJet->fill(fFatterPuppiJetArr,                     0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); } 
    }
  }
  if(fIsActivePF) { 
    fPFParArr->Clear();
    if(fUseAOD) { fFillerPF->fill(fPFParArr,fPVArr,iEvent); }
    else        { fFillerPF->fillMiniAOD(fPFParArr,fPVArr,iEvent); }
  }
  if(fIsActiveRH) { 
    fRHParArr->Clear();
    fFillerRH->fill(fRHParArr,iEvent,iSetup);
  }
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames)
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
void NtuplerMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  /*
  if(fIsActiveGenInfo) { 
    // Get generator event information
    edm::Handle<GenRunInfoProduct> hGenRunInfoProduct;
    iRun.getByToken(fTokGenRunInfo,hGenRunInfoProduct);
    assert(hGenRunInfoProduct.isValid());
    fXS = float(hGenRunInfoProduct->crossSection());
  }
  */
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void NtuplerMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void NtuplerMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerMod);
