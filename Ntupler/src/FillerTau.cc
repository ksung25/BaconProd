#include "BaconProd/Ntupler/interface/FillerTau.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerTau::FillerTau(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt  (iConfig.getUntrackedParameter<double>("minPt",10)),
  fTauName(iConfig.getUntrackedParameter<std::string>("edmName","hpsPFTauProducer")),
  fPuppiName     (iConfig.getUntrackedParameter<std::string>("edmPuppiName","puppi")),
  fPuppiNoLepName(iConfig.getUntrackedParameter<std::string>("edmPuppiNoLepName","puppiNoLep")),
  fUsePuppi      (iConfig.getUntrackedParameter<bool>("usePuppi",true)),
  fUseTO         (iConfig.getUntrackedParameter<bool>("useTriggerObject",false)),
  fUseAOD (useAOD)
{
  if(fUseAOD) {
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseElectronRejection",kByLooseElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumElectronRejection",kByMediumElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightElectronRejection",kByTightElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA6VLooseElectronRejection",kByMVA6VLooseElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA6LooseElectronRejection",kByMVA6LooseElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA6MediumElectronRejection",kByMVA6MediumElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA6TightElectronRejection",kByMVA6TightElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA6VTightElectronRejection",kByMVA6VTightElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseMuonRejection",kByLooseMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumMuonRejection",kByMediumMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightMuonRejection",kByTightMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseMuonRejection3",kByLooseMuonRejection3));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightMuonRejection3",kByTightMuonRejection3));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVALooseMuonRejection",kByMVALooseMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVAMediumMuonRejection",kByMVAMediumMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVATightMuonRejection",kByMVATightMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByDecayModeFinding",kByDecayModeFinding));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolation",kByVLooseIsolation));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr",kByVLooseCombinedIsolationDBSumPtCorr));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr",kByLooseCombinedIsolationDBSumPtCorr));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr",kByMediumCombinedIsolationDBSumPtCorr));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr",kByTightCombinedIsolationDBSumPtCorr));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits",kByLooseCombinedIsolationDBSumPtCorr3Hits));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits",kByMediumCombinedIsolationDBSumPtCorr3Hits));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits",kByTightCombinedIsolationDBSumPtCorr3Hits));
    /*
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT",kByVLooseIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwoLT",kByLooseIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwoLT",kByMediumIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwoLT",kByTightIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT",kByVTightIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwoLT",kByVVTightIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwLT",kByVLooseIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT",kByLooseIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT",kByMediumIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT",kByTightIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwLT",kByVTightIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwLT",kByVVTightIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwoLT",kByVLooseIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationMVA3newDMwoLT",kByLooseIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationMVA3newDMwoLT",kByMediumIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationMVA3newDMwoLT",kByTightIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVTightIsolationMVA3newDMwoLT",kByVTightIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwoLT",kByVVTightIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwLT",kByVLooseIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationMVA3newDMwLT",kByLooseIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationMVA3newDMwLT",kByMediumIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationMVA3newDMwLT",kByTightIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVTightIsolationMVA3newDMwLT",kByVTightIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwLT",kByVVTightIsolationMVA3newDMwLT));
    */
  } else {
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronLoose",kByLooseElectronRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronMedium",kByMediumElectronRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronTight",kByTightElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronVLooseMVA6",kByMVA6VLooseElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronLooseMVA6",kByMVA6LooseElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronMediumMVA6",kByMVA6MediumElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronTightMVA6",kByMVA6TightElectronRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstElectronVTightMVA6",kByMVA6VTightElectronRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonLoose",kByLooseMuonRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonMedium",kByMediumMuonRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonTight",kByTightMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonLoose3",kByLooseMuonRejection3));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonTight3",kByTightMuonRejection3));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonLooseMVA",kByMVALooseMuonRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonMediumMVA",kByMVAMediumMuonRejection));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("againstMuonTightMVA",kByMVATightMuonRejection));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("decayModeFinding",kByDecayModeFinding));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVLooseIsolation",kByVLooseIsolation));                                       (!) not in MINIAOD
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVLooseCombinedIsolationDBSumPtCorr",kByVLooseCombinedIsolationDBSumPtCorr)); (!) not in MINIAOD
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byLooseCombinedIsolationDBSumPtCorr",kByLooseCombinedIsolationDBSumPtCorr));   (!) not in MINIAOD
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byMediumCombinedIsolationDBSumPtCorr",kByMediumCombinedIsolationDBSumPtCorr)); (!) not in MINIAOD
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byTightCombinedIsolationDBSumPtCorr",kByTightCombinedIsolationDBSumPtCorr));   (!) not in MINIAOD
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byLooseCombinedIsolationDeltaBetaCorr3Hits",kByLooseCombinedIsolationDBSumPtCorr3Hits));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byMediumCombinedIsolationDeltaBetaCorr3Hits",kByMediumCombinedIsolationDBSumPtCorr3Hits));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byTightCombinedIsolationDeltaBetaCorr3Hits",kByTightCombinedIsolationDBSumPtCorr3Hits));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVLooseIsolationMVA3oldDMwoLT",kByVLooseIsolationMVA3oldDMwoLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byLooseIsolationMVA3oldDMwoLT",kByLooseIsolationMVA3oldDMwoLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byMediumIsolationMVA3oldDMwoLT",kByMediumIsolationMVA3oldDMwoLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byTightIsolationMVA3oldDMwoLT",kByTightIsolationMVA3oldDMwoLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVTightIsolationMVA3oldDMwoLT",kByVTightIsolationMVA3oldDMwoLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVVTightIsolationMVA3oldDMwoLT",kByVVTightIsolationMVA3oldDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVLooseIsolationMVArun2v1DBoldDMwLT",kByVLooseIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byLooseIsolationMVArun2v1DBoldDMwLT",kByLooseIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byMediumIsolationMVArun2v1DBoldDMwLT",kByMediumIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byTightIsolationMVArun2v1DBoldDMwLT",kByTightIsolationMVA3oldDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVTightIsolationMVArun2v1DBoldDMwLT",kByVTightIsolationMVA3oldDMwLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVVTightIsolationMVA3oldDMwLT",kByVVTightIsolationMVA3oldDMwLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVLooseIsolationMVA3newDMwoLT",kByVLooseIsolationMVA3newDMwoLT));
    //    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byLooseIsolationMVA3newDMwoLT",kByLooseIsolationMVA3newDMwoLT));
    //    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byMediumIsolationMVA3newDMwoLT",kByMediumIsolationMVA3newDMwoLT));
    //    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byTightIsolationMVA3newDMwoLT",kByTightIsolationMVA3newDMwoLT));
    //    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVTightIsolationMVA3newDMwoLT",kByVTightIsolationMVA3newDMwoLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVVTightIsolationMVA3newDMwoLT",kByVVTightIsolationMVA3newDMwoLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVLooseIsolationMVArun2v1DBnewDMwLT",kByVLooseIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byLooseIsolationMVArun2v1DBnewDMwLT",kByLooseIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byMediumIsolationMVArun2v1DBnewDMwLT",kByMediumIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byTightIsolationMVArun2v1DBnewDMwLT",kByTightIsolationMVA3newDMwLT));
    fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVTightIsolationMVArun2v1DBnewDMwLT",kByVTightIsolationMVA3newDMwLT));
    //fMyTauDiscHandles.push_back(new MyTauDiscHandle("byVVTightIsolationMVA3newDMwLT",kByVVTightIsolationMVA3newDMwLT));
  }
  if(fUseAOD)  fTokTauName      = iC.consumes<reco::PFTauCollection>(fTauName);
  if(!fUseAOD) fTokPatTauName   = iC.consumes<pat::TauCollection>   (fTauName);
  for(unsigned int idisc = 0; idisc < fTokTauHandles.size(); idisc++) { 
    edm::EDGetTokenT<reco::PFTauDiscriminator>  lTokTauHandles = iC.consumes<reco::PFTauDiscriminator>(fMyTauDiscHandles[idisc]->name);
    fTokTauHandles.push_back(lTokTauHandles);
  }
  // Get raw value and category of "MVA6" electron rejection discriminator
  std::string lTokMVA6EleRejRaw              = "hpsPFTauDiscriminationByMVA6rawElectronRejection";
  std::string lTokMVA6EleRejCat              = "hpsPFTauDiscriminationByMVA6rawElectronRejection:category";
  std::string lTokMVAMuonRejRaw              = "hpsPFTauDiscriminationByMVArawMuonRejection";
  std::string lTokCombIsoDBSumPtCorr3HitsRaw = "hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits";
  std::string lTokIsoMVA3newwRaw             = "hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw";
  std::string lTokIsoMVA3oldwRaw             = "hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw";


  fTokMVA6EleRejRaw              = iC.consumes<reco::PFTauDiscriminator>(lTokMVA6EleRejRaw);
  fTokMVA6EleRejCat              = iC.consumes<reco::PFTauDiscriminator>(lTokMVA6EleRejCat);
  fTokMVAMuonRejRaw              = iC.consumes<reco::PFTauDiscriminator>(lTokMVAMuonRejRaw);
  fTokCombIsoDBSumPtCorr3HitsRaw = iC.consumes<reco::PFTauDiscriminator>(lTokCombIsoDBSumPtCorr3HitsRaw);
  fTokIsoMVA3newwRaw             = iC.consumes<reco::PFTauDiscriminator>(lTokIsoMVA3newwRaw);
  fTokIsoMVA3oldwRaw             = iC.consumes<reco::PFTauDiscriminator>(lTokIsoMVA3oldwRaw);
  if(fUseAOD)  fTokPuppiName      = iC.consumes<reco::PFCandidateCollection>(fPuppiName);
  if(fUseAOD)  fTokPuppiNoLepName = iC.consumes<reco::PFCandidateCollection>(fPuppiNoLepName);
  if(!fUseAOD) fTokPuppiPATName      = iC.consumes<pat::PackedCandidateCollection>(fPuppiName);
  if(!fUseAOD) fTokPuppiNoLepPATName = iC.consumes<pat::PackedCandidateCollection>(fPuppiNoLepName);
  // event energy density for ring isolation
  //edm::Handle<double> hRho;

  //std::string cmssw_base_src = getenv("CMSSW_BASE");
  //cmssw_base_src += "/src/";
  
  //fRingIso.initialize (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("ringIsoFile", ""),true);
  //fRingIso2.initialize(cmssw_base_src + iConfig.getUntrackedParameter<std::string>("ringIso2File",""),true);
}

//--------------------------------------------------------------------------------------------------
FillerTau::~FillerTau(){}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerTau::fill(TClonesArray *array,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
		     const std::vector<TriggerRecord> &triggerRecords,
		     const trigger::TriggerEvent &triggerEvent) 
{
  assert(array);

  // Get tau collection
  edm::Handle<reco::PFTauCollection> hTauProduct;
  iEvent.getByToken(fTokTauName,hTauProduct);
  assert(hTauProduct.isValid());
  const reco::PFTauCollection *tauCol = hTauProduct.product();
  
  // Get HPS tau discriminators
  for(unsigned int idisc=0; idisc<fMyTauDiscHandles.size(); idisc++) {
    iEvent.getByToken(fTokTauHandles[idisc], fMyTauDiscHandles[idisc]->handle);
  }

   // Get raw value and category of "MVA6" electron rejection discriminator
  edm::Handle<reco::PFTauDiscriminator> hMVA6EleRejRaw;
  iEvent.getByToken(fTokMVA6EleRejRaw,hMVA6EleRejRaw);
  edm::Handle<reco::PFTauDiscriminator> hMVA6EleRejCat;
  iEvent.getByToken(fTokMVA6EleRejCat,hMVA6EleRejCat);
  
  // Get raw value of new mva anti-muon raw value
  edm::Handle<reco::PFTauDiscriminator> hMVAMuonRejRaw;
  iEvent.getByToken(fTokMVAMuonRejRaw,hMVAMuonRejRaw);

  // Get raw value for isolation discriminators
  edm::Handle<reco::PFTauDiscriminator> hCombIsoDBSumPtCorr3HitsRaw;
  iEvent.getByToken(fTokCombIsoDBSumPtCorr3HitsRaw,hCombIsoDBSumPtCorr3HitsRaw);
  //edm::Handle<reco::PFTauDiscriminator> hIsoMVA3oldwoRaw;
  //iEvent.getByLabel("hpsPFTauDiscriminationByIsolationMVA3oldDMwoLTraw",hIsoMVA3oldwoRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3oldwRaw;
  iEvent.getByToken(fTokIsoMVA3oldwRaw,hIsoMVA3oldwRaw);
  //edm::Handle<reco::PFTauDiscriminator> hIsoMVA3newwoRaw;
  //iEvent.getByToken(fTokIsoMVA3newwoRaw,hIsoMVA3newwoRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3newwRaw;
  iEvent.getByToken(fTokIsoMVA3newwRaw,hIsoMVA3newwRaw);

  // event energy density for ring isolation
  //edm::Handle<double> hRho;
  //edm::InputTag rhoTag(fRhoName,"","RECO");
  //iEvent.getByLabel(rhoTag,hRho);

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

  for(reco::PFTauCollection::const_iterator itTau = tauCol->begin(); itTau!=tauCol->end(); ++itTau) {
    
    reco::PFTauRef tauRef(hTauProduct, itTau - tauCol->begin());
    
    // tau pT cut
    if(itTau->pt() < fMinPt) continue;
    
    // construct object and place in array
    TClonesArray &rArray = *array;
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TTau();
    baconhep::TTau *pTau = (baconhep::TTau*)rArray[index];

  
    //
    // Kinematics
    //==============================    
    pTau->pt  = itTau->pt();
    pTau->eta = itTau->eta();
    pTau->phi = itTau->phi();
    pTau->m   = itTau->mass();
    pTau->e   = itTau->energy();
    pTau->q   = itTau->charge();

    
    //
    // Impact Parameter
    //==============================
    const reco::PFCandidatePtr& leadChHad = itTau->leadPFChargedHadrCand();
    pTau->dzLeadChHad = (leadChHad.isNonnull() && leadChHad->trackRef().isNonnull()) ? leadChHad->trackRef()->dz(pv.position()) : -999.;    
    pTau->d0LeadChHad = (leadChHad.isNonnull() && leadChHad->trackRef().isNonnull()) ? -leadChHad->trackRef()->dxy(pv.position()) : -999.;


   
    //
    // Isolation
    //==============================
    //pTau->ringIso     = fRingIso.isInitialized()  ? fRingIso.mvaValue (*itTau, *hRho) : -1;
    //pTau->ringIso2    = fRingIso2.isInitialized() ? fRingIso2.mvaValue(*itTau, *hRho) : -1;
    pTau->rawIso3Hits = hCombIsoDBSumPtCorr3HitsRaw.isValid() ? (*hCombIsoDBSumPtCorr3HitsRaw)[tauRef] : 0;
    //pTau->rawIsoMVA3oldDMwoLT = hIsoMVA3oldwoRaw.isValid() ? (*hIsoMVA3oldwoRaw)[tauRef] : 0;
    pTau->rawIsoMVA3oldDMwLT  = hIsoMVA3oldwRaw.isValid()  ? (*hIsoMVA3oldwRaw)[tauRef]  : 0;
    //pTau->rawIsoMVA3newDMwoLT = hIsoMVA3newwoRaw.isValid() ? (*hIsoMVA3newwoRaw)[tauRef] : 0;
    pTau->rawIsoMVA3newDMwLT  = hIsoMVA3newwRaw.isValid()  ? (*hIsoMVA3newwRaw)[tauRef]  : 0;
    if(fUsePuppi) { 
      const std::vector<reco::PFCandidatePtr> pCands = itTau->signalPFCands();
      double pEta = itTau->eta();
      double pPhi = itTau->phi();
      computeIso(pEta,pPhi,pCands, 0.4, (*pfPuppi), 
		 pTau->puppiChHadIso,
		 pTau->puppiGammaIso,
		 pTau->puppiNeuHadIso);
      
      computeIso(pEta,pPhi,pCands, 0.4, (*pfPuppiNoLep),
		 pTau->puppiChHadIsoNoLep,
		 pTau->puppiGammaIsoNoLep,
		 pTau->puppiNeuHadIsoNoLep);
    }

    //
    // Identification
    //==============================               
    pTau->nSignalChHad = itTau->signalPFChargedHadrCands().size();
    pTau->nSignalGamma = itTau->signalPFGammaCands().size();
    pTau->decaymode    = itTau->decayMode();
    pTau->hpsDisc=0;
    for(unsigned int idisc=0; idisc<fMyTauDiscHandles.size(); idisc++) {
      if(fMyTauDiscHandles[idisc]->value(tauRef)>0) {
        pTau->hpsDisc |= fMyTauDiscHandles[idisc]->flag;
      }
    }
    
    pTau->antiEleMVA6      = hMVA6EleRejRaw.isValid() ? (*hMVA6EleRejRaw)[tauRef] : 0;
    pTau->antiEleMVA6Cat   = hMVA6EleRejCat.isValid() ? (*hMVA6EleRejCat)[tauRef] : 0;
    pTau->rawMuonRejection = hMVAMuonRejRaw.isValid() ? (*hMVAMuonRejRaw)[tauRef] : 0;
    
    if(fUseTO) pTau->hltMatchBits = TriggerTools::matchHLT(pTau->eta, pTau->phi, triggerRecords, triggerEvent);
  } 
}

// === filler for MINIAOD ===
void FillerTau::fill(TClonesArray *array,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv,
                     const std::vector<TriggerRecord> &triggerRecords,
                     const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);

  // Get tau collection
  edm::Handle<pat::TauCollection> hTauProduct;
  iEvent.getByToken(fTokPatTauName,hTauProduct);
  assert(hTauProduct.isValid());
  const pat::TauCollection *tauCol = hTauProduct.product();

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

  for(pat::TauCollection::const_iterator itTau = tauCol->begin(); itTau!=tauCol->end(); ++itTau) {

    // tau pT cut
    if(itTau->pt() < fMinPt) continue;

    // construct object and place in array
    TClonesArray &rArray = *array;
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TTau();
    baconhep::TTau *pTau = (baconhep::TTau*)rArray[index];


    //
    // Kinematics
    //==============================
    pTau->pt  = itTau->pt();
    pTau->eta = itTau->eta();
    pTau->phi = itTau->phi();
    pTau->m   = itTau->mass();
    pTau->e   = itTau->energy();
    pTau->q   = itTau->charge();

    //
    // Impact Parameter (!) how to get in MINIAOD?
    //==============================
    //const reco::CandidatePtr& leadChHad = itTau->leadChargedHadrCand();
    //pTau->dzLeadChHad = (leadChHad.isNonnull() && leadChHad->trackRef().isNonnull()) ? leadChHad->trackRef()->dz(pv.position()) : -999.;
    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(itTau->leadChargedHadrCand().get());
    pTau->dzLeadChHad = (packedLeadTauCand)?packedLeadTauCand->dz():-999;
    pTau->d0LeadChHad = (packedLeadTauCand)?-packedLeadTauCand->dxy():-999;
    
    //
    // Isolation
    //==============================
    //pTau->ringIso     = fRingIso.isInitialized()  ? fRingIso.mvaValue (*itTau, *hRho) : -1;
    //pTau->ringIso2    = fRingIso2.isInitialized() ? fRingIso2.mvaValue(*itTau, *hRho) : -1;
    pTau->rawIso3Hits         = itTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    //pTau->rawIsoMVA3oldDMwoLT = itTau->tauID("byIsolationMVA3oldDMwoLTraw");
    //pTau->rawIsoMVA3oldDMwLT  = itTau->tauID("byIsolationMVA3oldDMwLTraw");
    //pTau->rawIsoMVA3newDMwoLT = itTau->tauID("byIsolationMVA3newDMwoLTraw");
    pTau->rawIsoMVA3newDMwLT  = itTau->tauID("byIsolationMVArun2v1DBnewDMwLTraw");
    pTau->rawIsoMVA3oldDMwLT  = itTau->tauID("byIsolationMVArun2v1DBoldDMwLTraw");

    
    if(fUsePuppi) { 
      const std::vector<reco::PFCandidatePtr> pCands = itTau->signalPFCands();
      double pEta = itTau->eta();
      double pPhi = itTau->phi();
      computeIso(pEta,pPhi,pCands, 0.4, (*pfPuppi), 
		 pTau->puppiChHadIso,
		 pTau->puppiGammaIso,
		 pTau->puppiNeuHadIso);
      
      computeIso(pEta,pPhi,pCands, 0.4, (*pfPuppiNoLep),
		 pTau->puppiChHadIsoNoLep,
		 pTau->puppiGammaIsoNoLep,
		 pTau->puppiNeuHadIsoNoLep);
    }
    
    //
    // Identification
    //==============================
    pTau->nSignalChHad = itTau->signalPFChargedHadrCands().size();
    pTau->nSignalGamma = itTau->signalPFGammaCands().size();
    pTau->decaymode    = itTau->decayMode();
    pTau->hpsDisc=0;
    for(unsigned int idisc=0; idisc<fMyTauDiscHandles.size(); idisc++) {
      if(itTau->tauID(fMyTauDiscHandles[idisc]->name)>0.0) {
        pTau->hpsDisc |= fMyTauDiscHandles[idisc]->flag;
      }
    }
    pTau->antiEleMVA6      = itTau->tauID("againstElectronMVA6Raw");
    pTau->antiEleMVA6Cat   = itTau->tauID("againstElectronMVA6category");
    //pTau->rawMuonRejection = itTau->tauID("againstMuonMVAraw");
    if(fUseTO) pTau->hltMatchBits = TriggerTools::matchHLT(pTau->eta, pTau->phi, triggerRecords, triggerObjects);
  }
}
void FillerTau::computeIso(double &iEta,double &iPhi,const std::vector<reco::PFCandidatePtr>& iCands, const double extRadius,
			   const reco::PFCandidateCollection    &puppi,
			   float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{

  
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double ptMin           = 0.1;  
  //const double intRadiusChHad  = 0.0001;
  //const double intRadiusGamma  = 0.01;
  //const double intRadiusNeuHad = 0.01;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<puppi.size(); ipf++) {
    const reco::PFCandidate pfcand = puppi.at(ipf);    
    bool pPass = true;
    for(unsigned int iTau=0; iTau<iCands.size(); iTau++) { 
      const reco::PFCandidatePtr tauCand = iCands[iTau];
      double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), tauCand->eta(), tauCand->phi());
      if(dr < 0.0001) pPass = false;  
    }
    if(pPass) { 
      if     (pfcand.particleId() == reco::PFCandidate::h)     { intRadius = 0;}//intRadiusChHad; }
      else if(pfcand.particleId() == reco::PFCandidate::gamma) { intRadius = 0;}//intRadiusGamma;  }
      else if(pfcand.particleId() == reco::PFCandidate::h0)    { intRadius = 0;}//intRadiusNeuHad; }
            
      double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), iEta, iPhi);
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
void FillerTau::computeIso(double &iEta,double &iPhi,const std::vector<reco::PFCandidatePtr>& iCands, const double extRadius,
			   const pat::PackedCandidateCollection    &puppi,
			   float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{

  
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double ptMin           = 0.1;  
  //const double intRadiusChHad  = 0.0001;
  //const double intRadiusGamma  = 0.01;
  //const double intRadiusNeuHad = 0.01;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<puppi.size(); ipf++) {
    const pat::PackedCandidate pfcand = puppi.at(ipf);    
    bool pPass = true;
    for(unsigned int iTau=0; iTau<iCands.size(); iTau++) { 
      const reco::PFCandidatePtr tauCand = iCands[iTau];
      double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), tauCand->eta(), tauCand->phi());
      if(dr < 0.0001) pPass = false;  
    }
    if(pPass) { 
      if     (abs(pfcand.pdgId()) == 211)     { intRadius = 0;}//intRadiusChHad; }
      else if(abs(pfcand.pdgId()) == 22)      { intRadius = 0;}//intRadiusGamma;  }
      else if(abs(pfcand.pdgId()) == 130)     { intRadius = 0;}//intRadiusNeuHad; }
            
      double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), iEta, iPhi);
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
