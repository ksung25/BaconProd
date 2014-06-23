#include "BaconProd/Ntupler/interface/FillerTau.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerTau::FillerTau(const edm::ParameterSet &iConfig):
  fMinPt  (iConfig.getUntrackedParameter<double>("minPt",15)),
  fTauName(iConfig.getUntrackedParameter<std::string>("edmName","hpsPFTauProducer"))
{
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByDecayModeFinding",kByDecayModeFinding));

  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolation",kByVLooseIsolation));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolation", kByLooseIsolation));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolation",kByMediumIsolation));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolation", kByTightIsolation));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr",kByVLooseIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr", kByLooseIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr",kByMediumIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationDBSumPtCorr", kByTightIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits",kByLooseCombinedIsolationDBSumPtCorr3Hits));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits",kByMediumCombinedIsolationDBSumPtCorr3Hits));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", kByTightCombinedIsolationDBSumPtCorr3Hits));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseElectronRejection", kByLooseElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumElectronRejection",kByMediumElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightElectronRejection", kByTightElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA3LooseElectronRejection", kByMVA3LooseElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA3MediumElectronRejection",kByMVA3MediumElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA3TightElectronRejection", kByMVA3TightElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA3VTightElectronRejection",kByMVA3VTightElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseMuonRejection", kByLooseMuonRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumMuonRejection",kByMediumMuonRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightMuonRejection", kByTightMuonRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseMuonRejection2", kByLooseMuonRejection2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumMuonRejection2",kByMediumMuonRejection2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightMuonRejection2", kByTightMuonRejection2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseMuonRejection3",kByLooseMuonRejection3));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightMuonRejection3",kByTightMuonRejection3));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA5VLooseElectronRejection",kByMVA5VLooseElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA5LooseElectronRejection",kByMVA5LooseElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA5MediumElectronRejection",kByMVA5MediumElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA5TightElectronRejection",kByMVA5TightElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVA5VTightElectronRejection",kByMVA5VTightElectronRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVALooseMuonRejection",kByMVALooseMuonRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVAMediumMuonRejection",kByMVAMediumMuonRejection));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMVATightMuonRejection",kByMVATightMuonRejection));
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
   
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  
  //fRingIso.initialize (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("ringIsoFile", ""),true);
  //fRingIso2.initialize(cmssw_base_src + iConfig.getUntrackedParameter<std::string>("ringIso2File",""),true);
}

//--------------------------------------------------------------------------------------------------
FillerTau::~FillerTau(){}

//--------------------------------------------------------------------------------------------------
void FillerTau::fill(TClonesArray *array,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
		     const std::vector<TriggerRecord> &triggerRecords,
		     const trigger::TriggerEvent &triggerEvent) 
{
  assert(array);

  // Get tau collection
  edm::Handle<reco::PFTauCollection> hTauProduct;
  iEvent.getByLabel(fTauName,hTauProduct);
  assert(hTauProduct.isValid());
  const reco::PFTauCollection *tauCol = hTauProduct.product();
  
  // Get HPS tau discriminators
  for(unsigned int idisc=0; idisc<fMyTauDiscHandles.size(); idisc++) {
    iEvent.getByLabel(fMyTauDiscHandles[idisc]->name, fMyTauDiscHandles[idisc]->handle);
  }

  // Get raw value and category of "MVA3" electron rejection discriminator
  edm::Handle<reco::PFTauDiscriminator> hMVA3EleRejRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVA3rawElectronRejection",hMVA3EleRejRaw);
  edm::Handle<reco::PFTauDiscriminator> hMVA3EleRejCat;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVA3rawElectronRejection:category",hMVA3EleRejCat);

   // Get raw value and category of "MVA5" electron rejection discriminator
  edm::Handle<reco::PFTauDiscriminator> hMVA5EleRejRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVA5rawElectronRejection",hMVA5EleRejRaw);
  edm::Handle<reco::PFTauDiscriminator> hMVA5EleRejCat;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVA5rawElectronRejection:category",hMVA5EleRejCat);
  
  // Get raw value of new mva anti-muon raw value
  edm::Handle<reco::PFTauDiscriminator> hMVAMuonRejRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVArawMuonRejection",hMVAMuonRejRaw);

  // Get raw value for isolation discriminators
  edm::Handle<reco::PFTauDiscriminator> hCombIsoDBSumPtCorr3HitsRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits",hCombIsoDBSumPtCorr3HitsRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3oldwoRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByIsolationMVA3oldDMwoLTraw",hIsoMVA3oldwoRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3oldwRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByIsolationMVA3oldDMwLTraw",hIsoMVA3oldwRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3newwoRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw",hIsoMVA3newwoRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3newwRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByIsolationMVA3newDMwLTraw",hIsoMVA3newwRaw);

  // event energy density for ring isolation
  //edm::Handle<double> hRho;
  //edm::InputTag rhoTag(fRhoName,"","RECO");
  //iEvent.getByLabel(rhoTag,hRho);

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
   
    //
    // Isolation
    //==============================
    //pTau->ringIso     = fRingIso.isInitialized()  ? fRingIso.mvaValue (*itTau, *hRho) : -1;
    //pTau->ringIso2    = fRingIso2.isInitialized() ? fRingIso2.mvaValue(*itTau, *hRho) : -1;
    pTau->rawIso3Hits = hCombIsoDBSumPtCorr3HitsRaw.isValid() ? (*hCombIsoDBSumPtCorr3HitsRaw)[tauRef] : 0;
    pTau->rawIsoMVA3oldDMwoLT = hIsoMVA3oldwoRaw.isValid() ? (*hIsoMVA3oldwoRaw)[tauRef] : 0;
    pTau->rawIsoMVA3oldDMwLT  = hIsoMVA3oldwRaw.isValid() ? (*hIsoMVA3oldwRaw)[tauRef] : 0;
    pTau->rawIsoMVA3newDMwoLT = hIsoMVA3newwoRaw.isValid() ? (*hIsoMVA3newwoRaw)[tauRef] : 0;
    pTau->rawIsoMVA3newDMwLT  = hIsoMVA3newwRaw.isValid() ? (*hIsoMVA3newwRaw)[tauRef] : 0;

    //
    // Identification
    //==============================               
    pTau->nSignalChHad = itTau->signalPFChargedHadrCands().size();
    pTau->nSignalGamma = itTau->signalPFGammaCands().size();

    pTau->hpsDisc=0;
    for(unsigned int idisc=0; idisc<fMyTauDiscHandles.size(); idisc++) {
      if(fMyTauDiscHandles[idisc]->value(tauRef)>0) {
        pTau->hpsDisc |= fMyTauDiscHandles[idisc]->flag;
      }
    }
    
    pTau->antiEleMVA3    = hMVA3EleRejRaw.isValid() ? (*hMVA3EleRejRaw)[tauRef] : 0;
    pTau->antiEleMVA3Cat = hMVA3EleRejCat.isValid() ? (*hMVA3EleRejCat)[tauRef] : 0;
    pTau->antiEleMVA5    = hMVA5EleRejRaw.isValid() ? (*hMVA5EleRejRaw)[tauRef] : 0;
    pTau->antiEleMVA5Cat = hMVA5EleRejCat.isValid() ? (*hMVA5EleRejCat)[tauRef] : 0;
    pTau->rawMuonRejection = hMVAMuonRejRaw.isValid() ? (*hMVAMuonRejRaw)[tauRef] : 0;
    
    pTau->hltMatchBits = 0.0; //TriggerTools::matchHLT(pTau->eta, pTau->phi, triggerRecords, triggerEvent);
  } 
}
