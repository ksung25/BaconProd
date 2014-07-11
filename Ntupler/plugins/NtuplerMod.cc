#include "NtuplerMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TSusyGen.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/DataFormats/interface/TTopJet.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TRHPart.hh"

#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Ntupler/interface/FillerPhoton.hh"
#include "BaconProd/Ntupler/interface/FillerTau.hh"
#include "BaconProd/Ntupler/interface/FillerJet.hh"
#include "BaconProd/Ntupler/interface/FillerPF.hh"
#include "BaconProd/Ntupler/interface/FillerRH.hh"

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
NtuplerMod::NtuplerMod(const edm::ParameterSet &iConfig):
  fSkipOnHLTFail     (iConfig.getUntrackedParameter<bool>("skipOnHLTFail",false)),
  fHLTTag            ("TriggerResults","","HLT"),
  fHLTObjTag         ("hltTriggerSummaryAOD","","HLT"),
  fHLTFile           (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fPVName            (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fPFCandName        (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fComputeFullJetInfo(false),
  fFillerEvtInfo     (0),
  fFillerGenInfo     (0),
  fFillerPV          (0),
  fFillerEle         (0),
  fFillerMuon        (0),
  fFillerPhoton      (0),
  fFillerTau         (0),
  fFillerJet         (0),
  fFillerPF          (0),
  fFillerRH          (0),
  fTrigger           (0),
  fIsActiveEvtInfo   (false),
  fIsActiveGenInfo   (false),
  fIsActivePV        (false),
  fIsActiveEle       (false),
  fIsActiveMuon      (false),
  fIsActivePhoton    (false),
  fIsActiveTau       (false),
  fIsActiveJet       (false),
  fIsActivePF        (false),
  fIsActiveRH        (false),
  fOutputName        (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile        (0),
  fTotalEvents       (0),
  fEventTree         (0),
  fEvtInfo           (0),
  fGenEvtInfo        (0),
  fGenParArr         (0),
  fEleArr            (0),
  fMuonArr           (0),
  fTauArr            (0),
  fJetArr            (0),
  fPhotonArr         (0),
  fPVArr             (0),
  fAddJetArr         (0),
  fTopJetArr         (0),
  fPFParArr          (0),
  fRHParArr          (0)
{
  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TSusyGen::Class()->IgnoreTObjectStreamer();
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  baconhep::TMuon::Class()->IgnoreTObjectStreamer();
  baconhep::TElectron::Class()->IgnoreTObjectStreamer();
  baconhep::TTau::Class()->IgnoreTObjectStreamer();
  baconhep::TJet::Class()->IgnoreTObjectStreamer();
  baconhep::TPhoton::Class()->IgnoreTObjectStreamer();
  baconhep::TVertex::Class()->IgnoreTObjectStreamer();
  baconhep::TAddJet::Class()->IgnoreTObjectStreamer();
  baconhep::TPFPart::Class()->IgnoreTObjectStreamer();
  baconhep::TRHPart::Class()->IgnoreTObjectStreamer();
  
  //
  // Set up bacon objects and configure fillers
  // 
  if(iConfig.existsAs<edm::ParameterSet>("Info",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Info"));
    fIsActiveEvtInfo = cfg.getUntrackedParameter<bool>("isActive");
    fAddSusyGen      = cfg.getUntrackedParameter<bool>("addSusyGen",false);

    if(fIsActiveEvtInfo) {
      fEvtInfo       = new baconhep::TEventInfo();                        assert(fEvtInfo);
      fFillerEvtInfo = new baconhep::FillerEventInfo(cfg);                assert(fFillerEvtInfo);
      if(fAddSusyGen)  fSusyGen       = new baconhep::TSusyGen();         assert(fSusyGen);
    }
  }
  
  if(iConfig.existsAs<edm::ParameterSet>("GenInfo",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("GenInfo"));
    fIsActiveGenInfo = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveGenInfo) {
      fGenEvtInfo    = new baconhep::TGenEventInfo();                   assert(fGenEvtInfo);
      fGenParArr     = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
      fFillerGenInfo = new baconhep::FillerGenInfo(cfg);                assert(fFillerGenInfo);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("PV",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PV"));
    fIsActivePV = cfg.getUntrackedParameter<bool>("isActive");
    
    // create array and filler even if vertices won't be saved to output (i.e. fIsActivePV == false),
    // because FillerVertex::fill(...) is used to find the event primary vertex
    // (not elegant, but I suppose a dedicated PV finding function can be implemented somewhere...)
    fPVArr    = new TClonesArray("baconhep::TVertex"); assert(fPVArr);
    fFillerPV = new baconhep::FillerVertex(cfg);       assert(fFillerPV);
  }
    
  if(iConfig.existsAs<edm::ParameterSet>("Electron",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Electron"));
    fIsActiveEle = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveEle) {
      fEleArr    = new TClonesArray("baconhep::TElectron"); assert(fEleArr);
      fFillerEle = new baconhep::FillerElectron(cfg);       assert(fFillerEle);
    }
  }  

  if(iConfig.existsAs<edm::ParameterSet>("Muon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Muon"));
    fIsActiveMuon = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveMuon) {
      fMuonArr    = new TClonesArray("baconhep::TMuon"); assert(fMuonArr);
      fFillerMuon = new baconhep::FillerMuon(cfg);       assert(fFillerMuon);
    }
  }  

  if(iConfig.existsAs<edm::ParameterSet>("Photon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Photon"));
    fIsActivePhoton = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActivePhoton) {
      fPhotonArr    = new TClonesArray("baconhep::TPhoton"); assert(fPhotonArr);
      fFillerPhoton = new baconhep::FillerPhoton(cfg);       assert(fFillerPhoton);
    }
  } 

  if(iConfig.existsAs<edm::ParameterSet>("Tau",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Tau"));
    fIsActiveTau = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveTau) {
      fTauArr    = new TClonesArray("baconhep::TTau"); assert(fTauArr);
      fFillerTau = new baconhep::FillerTau(cfg);       assert(fFillerTau);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("Jet",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Jet"));
    fIsActiveJet = cfg.getUntrackedParameter<bool>("isActive");
    
    std::vector<double> empty_vdouble;
    fConeSizes = cfg.getUntrackedParameter< std::vector<double> >("coneSizes",empty_vdouble);
    const unsigned int kNCones = fConeSizes.size();
    if(kNCones==0) { fIsActiveJet = false; }

    std::vector<std::string> empty_vstring;
    fJetPostFix = cfg.getUntrackedParameter< std::vector<std::string> >("postFix",empty_vstring);
  
    fComputeFullJetInfo = cfg.getUntrackedParameter<bool>("doComputeFullJetInfo");
    if(fIsActiveJet) {
      fFillerJet = new baconhep::FillerJet*[kNCones*(fJetPostFix.size())];
      fJetArr    = new TClonesArray        *[kNCones*(fJetPostFix.size())];
      if(fComputeFullJetInfo) {
	fAddJetArr = new TClonesArray*[kNCones*(fJetPostFix.size())]; 
	fTopJetArr = new TClonesArray*[kNCones*(fJetPostFix.size())];
      }
      
      for(unsigned int i0=0; i0<kNCones; i0++) {
	for(unsigned int i1 = 0; i1 < fJetPostFix.size(); i1++) { 
	  int iId = i0*fJetPostFix.size()+i1;
	  fJetArr[iId] = new TClonesArray("baconhep::TJet");
	  assert(fJetArr[iId]);
	  
	  if(fComputeFullJetInfo) {
	    fAddJetArr[iId] = new TClonesArray("baconhep::TAddJet");
	    assert(fAddJetArr[iId]);
	    fTopJetArr[iId] = new TClonesArray("baconhep::TTopJet");
	    assert(fTopJetArr[iId]);
	  }
	  
	  std::stringstream pSS; pSS << "AK" << int(fConeSizes[i0]*10);  // assumes cone sizes are multiples of 0.1 
	  fFillerJet[iId] = new baconhep::FillerJet(cfg, fConeSizes[i0], pSS.str(),fJetPostFix[i1]);
	  assert(fFillerJet);
	}
      }
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
  if(iConfig.existsAs<edm::ParameterSet>("RecHit",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("RecHit"));
    fIsActiveRH = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveRH) {
      fRHParArr = new TClonesArray("baconhep::TRHPart",50000); assert(fRHParArr);
      fFillerRH = new baconhep::FillerRH(cfg);                 assert(fFillerRH);
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
  delete fFillerPF;
  delete fFillerRH;
  
  delete fTrigger;
  delete fSusyGen;
  delete fEvtInfo;
  delete fGenEvtInfo;
  delete fGenParArr;
  delete fEleArr;
  delete fMuonArr;
  delete fTauArr;
  delete fPhotonArr;
  delete fPVArr;
  delete fPFParArr;
  delete fRHParArr;
  
  if(fIsActiveJet) {
    for(unsigned int i0=0; i0<fConeSizes.size(); i0++) { delete fFillerJet[i0]; }
    delete [] fFillerJet;
    
    for(unsigned int i0=0; i0<fConeSizes.size(); i0++) { delete fJetArr[i0]; }
    delete [] fJetArr;
    if(fComputeFullJetInfo) {
      for(unsigned int i0=0; i0<fConeSizes.size(); i0++) { delete fAddJetArr[i0]; }
      delete [] fAddJetArr;
      for(unsigned int i0=0; i0<fConeSizes.size(); i0++) { delete fTopJetArr[i0]; }
      delete [] fTopJetArr;
    }
  }
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
  
  if(fIsActiveEvtInfo) { 
    fEventTree->Branch("Info",fEvtInfo); 
    if(fAddSusyGen)     fEventTree->Branch("SusyGen",fSusyGen); 
  }
  if(fIsActiveGenInfo) {
    fEventTree->Branch("GenEvtInfo",fGenEvtInfo);
    fEventTree->Branch("GenParticle",&fGenParArr);
  }
  if(fIsActiveEle)    { fEventTree->Branch("Electron", &fEleArr); }
  if(fIsActiveMuon)   { fEventTree->Branch("Muon",     &fMuonArr); }
  if(fIsActiveTau)    { fEventTree->Branch("Tau",      &fTauArr); }
  if(fIsActivePhoton) { fEventTree->Branch("Photon",   &fPhotonArr); }
  if(fIsActivePV)     { fEventTree->Branch("PV",       &fPVArr); }
  if(fIsActiveJet) {
    for(unsigned int i0=0; i0 < fConeSizes.size(); i0++) { 
      for(unsigned int i1 = 0; i1 < fJetPostFix.size(); i1++) { 
	std::stringstream pSS; pSS << "Jet0" << int(fConeSizes[i0]*10) << fJetPostFix[i1];  // assumes cone sizes are multiples of 0.1 
	int iId = i0*fJetPostFix.size() + i1;
	fEventTree->Branch(pSS.str().c_str(), &fJetArr[iId]);
	if(fComputeFullJetInfo) fEventTree->Branch(("Add"+pSS.str()).c_str(), &fAddJetArr[iId]);
	if(fComputeFullJetInfo) fEventTree->Branch(("Top"+pSS.str()).c_str(), &fTopJetArr[iId]);
      }
    }
  }
  if(fIsActivePF) fEventTree->Branch("PFPart", &fPFParArr);
  if(fIsActiveRH) fEventTree->Branch("RHPart", &fRHParArr);
  
  //
  // Triggers
  //
  setTriggers();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::endJob() 
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
void NtuplerMod::setTriggers()
{
  std::string cmssw_base_src = getenv("CMSSW_BASE"); cmssw_base_src+="/src/";
  fTrigger = new baconhep::TTrigger(cmssw_base_src + fHLTFile);
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    fFillerGenInfo->fill(fGenEvtInfo, fGenParArr, iEvent);
  }
    
  fPVArr->Clear();
  int nvertices = 0;
  const reco::Vertex *pv = fFillerPV->fill(fPVArr, nvertices, iEvent);
  assert(pv);
  
  separatePileUp(iEvent, *pv);
  
  if(fIsActiveEvtInfo) {
    fFillerEvtInfo->fill(fEvtInfo, iEvent, *pv, (nvertices>0), triggerBits,fSusyGen);
  }
  
  edm::Handle<trigger::TriggerEvent> hTrgEvt;
  iEvent.getByLabel(fHLTObjTag,hTrgEvt);
  
  if(fIsActiveEle) {
    fEleArr->Clear();
    fFillerEle->fill(fEleArr, iEvent, iSetup, *pv, nvertices, fPFNoPU, fTrigger->fRecords, *hTrgEvt);
  }

  if(fIsActiveMuon) {
    fMuonArr->Clear();  
    fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fPFNoPU, fPFPU, fTrigger->fRecords, *hTrgEvt);
  }

  if(fIsActivePhoton) {
    fPhotonArr->Clear();  
    fFillerPhoton->fill(fPhotonArr, iEvent, iSetup, *pv, fPFNoPU, fPFPU, fTrigger->fRecords, *hTrgEvt);
  }

  if(fIsActiveTau) {
    fTauArr->Clear();
    fFillerTau->fill(fTauArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);
  }
  
  if(fIsActiveJet) {
    for(unsigned int i0=0; i0<fConeSizes.size()*(fJetPostFix.size()); i0++) { 
      fJetArr[i0]->Clear();
      if(fComputeFullJetInfo) {
        fAddJetArr[i0]->Clear();      
        fTopJetArr[i0]->Clear();      
        fFillerJet[i0]->fill(fJetArr[i0], fAddJetArr[i0], fTopJetArr[i0], iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);
      } else {
        fFillerJet[i0]->fill(fJetArr[i0], 0, 0, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);
      }
    }
  }
  
  if(fIsActivePF) { 
    fPFParArr->Clear();
    fFillerPF->fill(fPFParArr,fPVArr,iEvent);
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
void NtuplerMod::separatePileUp(const edm::Event &iEvent, const reco::Vertex &pv)
{
  // recipe from Matthew Chan

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByLabel(fPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();  
  
  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();
  
  fPFNoPU.clear();
  fPFPU.clear();
  
  for(reco::PFCandidateCollection::const_iterator iP = pfCandCol->begin(); iP!=pfCandCol->end(); ++iP) {
    if(iP->particleId() == reco::PFCandidate::h) {  // charged hadrons
      if(iP->trackRef().isNonnull() && pv.trackWeight(iP->trackRef())>0) {
        // charged hadrons with track used to compute PV
	fPFNoPU.push_back(&(*iP)); 
      
      } else {
        // Find closest vertex to charged hadron's vertex source
	bool vertexFound = false;
	const reco::Vertex *closestVtx = 0;
	double dzmin = 10000;
	
	for(reco::VertexCollection::const_iterator iV = pvCol->begin(); iV!=pvCol->end(); ++iV) {
	  if(iP->trackRef().isNonnull() && iV->trackWeight(iP->trackRef())>0) {
	    vertexFound = true;
	    closestVtx  = &(*iV);
	    break;
	  }
	  
	  double dz = fabs(iP->vertex().z() - iV->z());
	  if(dz < dzmin) {
	    closestVtx = &(*iV);
	    dzmin      = dz;
	  }
	}
	
	if(vertexFound || closestVtx != &pv) {
	  fPFPU.push_back(&(*iP));
	} else {
	  fPFNoPU.push_back(&(*iP));  // Note: when no associated vertex found, assume to come from PV
	}
      }
      
    } else {  // all non-charged-hadron PFCandidates are considered to be from PV
      fPFNoPU.push_back(&(*iP));
    }
  }
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void NtuplerMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void NtuplerMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerMod);
