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
  fMinPt  (iConfig.getUntrackedParameter<double>("minPt",20)),
  fTauName(iConfig.getUntrackedParameter<std::string>("edmName","hpsPFTauProducer")),
  fRhoName(iConfig.getUntrackedParameter<std::string>("edmRhoForRingIso","kt6PFJets"))
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
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr",kByVLooseCombinedIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr", kByLooseCombinedIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr",kByMediumCombinedIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr", kByTightCombinedIsolationDBSumPtCorr));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationMVA", kByLooseIsolationMVA));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationMVA",kByMediumIsolationMVA));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationMVA", kByTightIsolationMVA));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByLooseIsolationMVA2", kByLooseIsolationMVA2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByMediumIsolationMVA2",kByMediumIsolationMVA2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationByTightIsolationMVA2", kByTightIsolationMVA2));

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
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationkByLooseMuonRejection2", kByLooseMuonRejection2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationkByMediumMuonRejection2",kByMediumMuonRejection2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationkByTightMuonRejection2", kByTightMuonRejection2));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationkByLooseMuonRejection3",kByLooseMuonRejection3));
  fMyTauDiscHandles.push_back(new MyTauDiscHandle("hpsPFTauDiscriminationkByTightMuonRejection3",kByTightMuonRejection3));
   
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  
  fRingIso.initialize (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("ringIsoFile", ""),true);
  fRingIso2.initialize(cmssw_base_src + iConfig.getUntrackedParameter<std::string>("ringIso2File",""),true);
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
  
  // Get raw value for isolation discriminator
  edm::Handle<reco::PFTauDiscriminator> hCombIsoDBSumPtCorr3HitsRaw;
  iEvent.getByLabel("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits",hCombIsoDBSumPtCorr3HitsRaw);
  edm::Handle<reco::PFTauDiscriminator> hIsoMVA3Raw;
  iEvent.getByLabel("hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw",hIsoMVA3Raw);
  
  // event energy density for ring isolation
  edm::Handle<double> hRho;
  edm::InputTag rhoTag(fRhoName,"rho","RECO");
  iEvent.getByLabel(rhoTag,hRho);


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
    const reco::PFCandidateRef& leadChHad = itTau->leadPFChargedHadrCand();
    pTau->dzLeadChHad = (leadChHad.isNonnull() && leadChHad->trackRef().isNonnull()) ? leadChHad->trackRef()->dz(pv.position()) : -999.;    
   
    //
    // Isolation
    //==============================
    pTau->ringIso     = fRingIso.isInitialized()  ? fRingIso.mvaValue (*itTau, *hRho) : -1;
    pTau->ringIso2    = fRingIso2.isInitialized() ? fRingIso2.mvaValue(*itTau, *hRho) : -1;
    pTau->rawIso3Hits = hCombIsoDBSumPtCorr3HitsRaw.isValid() ? (*hCombIsoDBSumPtCorr3HitsRaw)[tauRef] : 0;
    pTau->rawIsoMVA3  = hIsoMVA3Raw.isValid() ? (*hIsoMVA3Raw)[tauRef] : 0;

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
    
    pTau->hltMatchBits = TriggerTools::matchHLT(pTau->eta, pTau->phi, triggerRecords, triggerEvent);
  } 
}
