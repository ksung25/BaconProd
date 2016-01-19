#ifndef BACONPROD_NTUPLER_FILLERJET_HH
#define BACONPROD_NTUPLER_FILLERJET_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/JetPUIDMVACalculator.hh"
#include "BaconProd/Utils/interface/ShowerDeco.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
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
  class FillerJet
  {
    public:
       FillerJet(const edm::ParameterSet &iConfig, const bool useAOD);
      ~FillerJet();

      // === filler for AOD ===
      void fill(TClonesArray                     *array,                        // output array to be filled
		TClonesArray                     *iExtraArray,                  // Extra Array to be filled
                const edm::Event                 &iEvent,                       // event info
		const edm::EventSetup            &iSetup,                       // event setup info
	        const reco::Vertex		 &pv,	                        // event primary vertex
		const std::vector<TriggerRecord> &triggerRecords,               // list of trigger names and objects
		const trigger::TriggerEvent      *triggerEvent,                 // event trigger objects
		const pat::TriggerObjectStandAloneCollection *patTriggerObjects);  // hack for AOD type filler with miniAOD

      // === filler for MINIAOD ===
      void fill(TClonesArray                                 *array,            // output array to be filled
                TClonesArray                                 *iExtraArray,      // Extra Array to be filled
                const edm::Event                             &iEvent,           // event info
                const edm::EventSetup                        &iSetup,           // event setup info
                const reco::Vertex                           &pv,               // event primary vertex
                const std::vector<TriggerRecord>             &triggerRecords,   // list of trigger names and objects
                const pat::TriggerObjectStandAloneCollection &triggerObjects);  // event trigger objects
     void  initPUJetId();

    protected:
      void initJetCorr(const std::vector<std::string> &jecFiles, 
                       const std::vector<std::string> &jecUncFiles);
    
      void  addJet(TAddJet *pAddJet, const edm::Event &iEvent, const reco::PFJet &itJet, const reco::JetBaseRef &jetBaseRef);
      void  addJet(baconhep::TAddJet *pAddJet, const edm::Event &iEvent, const pat::Jet &itJet);

      const reco::BasicJet* match(const reco::PFJet *jet, const reco::BasicJetCollection *jets);
      const reco::BasicJet* match(const pat::Jet *iJet,   const reco::BasicJetCollection *jets);
      const reco::GenJet*   match(const reco::PFJet *jet, const reco::GenJetCollection *jets);
      const reco::GenJet*   match(const pat::Jet *iJet,   const reco::GenJetCollection *jets);
      
      // Jet cuts
      double fMinPt;
 
      // Do matching to GenJets?
      bool fUseGen;
      bool fApplyJEC;
      
      // EDM object collection names
      std::string fPVName;
      std::string fRhoName;
      std::string fJetName;
      std::string fGenJetName;
      std::string fJetFlavorName;
      std::string fPrunedJetName;
      std::string fTrimmedJetName;
      std::string fSoftDropJetName;
      std::string fSubJetName;
      std::string fCSVbtagName;
      std::string fCSVbtagSubJetName;
      std::string fCSVDoubleBtagName;
      std::string fJettinessName;
      std::string fQGLikelihood;
      std::string fQGLikelihoodSubJets;
      std::string fTopTaggerName;
      std::string fLowPtWeightFile;
      std::string fHighPtWeightFile;
      std::string fShowerDecoConf;
      double      fConeSize;
      bool        fComputeFullJetInfo;
      
      // Jet ID MVA
      JetPUIDMVACalculator fJetPUIDMVACalc;
      ShowerDeco*          fShowerDeco;

      // Random number generator for Q-jet volatility
      TRandom2* fRand;

      // JEC corrector
      FactorizedJetCorrector   *fJetCorr;
      JetCorrectionUncertainty *fJetUnc;

      bool fUseAOD;
  };
}
#endif
