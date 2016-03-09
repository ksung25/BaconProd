#ifndef BACONPROD_NTUPLER_FILLERCALOJET_HH
#define BACONPROD_NTUPLER_FILLERCALOJET_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "TRandom2.h"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
class FactorizedJetCorrector;
class JetCorrectionUncertainty;
namespace trigger {
  class TriggerEvent;
}

namespace baconhep
{
  class FillerCaloJet
  {
    public:
       FillerCaloJet(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerCaloJet();

      // === filler for AOD ===
      void fill(TClonesArray                     *array,                        // output array to be filled
                const edm::Event                 &iEvent,                       // event info
		const edm::EventSetup            &iSetup,                       // event setup info
		const std::vector<TriggerRecord> &triggerRecords,               // list of trigger names and objects
		const trigger::TriggerEvent      &triggerEvent);                // event trigger objects

    protected:
      void initJetCorr(const std::vector<std::string> &jecFiles, 
                       const std::vector<std::string> &jecUncFiles);
      
      const reco::BasicJet* match(const reco::CaloJet *jet, const reco::BasicJetCollection *jets);
      const reco::GenJet*   match(const reco::CaloJet *jet, const reco::GenJetCollection *jets);
      
      // Jet cuts
      double fMinPt;
 
      // Do matching to GenJets?
      bool fUseGen;
      
      // EDM object collection names
      std::string fPVName;
      std::string fRhoName;
      std::string fJetName;
      std::string fGenJetName;
      std::string fJetFlavorName;
      double      fConeSize;

      // JEC corrector
      FactorizedJetCorrector   *fJetCorr;
      JetCorrectionUncertainty *fJetUnc;
      edm::EDGetTokenT<reco::CaloJetCollection> fTokJetName;
      edm::EDGetTokenT<reco::GenJetCollection>  fTokGenJetName;
      edm::EDGetTokenT<double>                  fTokRhoTag;
  };
}
#endif
