#ifndef BACONPROD_NTUPLER_FILLERELECTRON_HH
#define BACONPROD_NTUPLER_FILLERELECTRON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/ElectronMomentumCorrector.hh"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
class EcalClusterLazyTools;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerElectron
  {
    public:
      FillerElectron(const edm::ParameterSet &iConfig);
      ~FillerElectron();
      
      void fill(TClonesArray                                *array,           // output array to be filled
                const edm::Event                            &iEvent,          // event info
		const edm::EventSetup                       &iSetup,          // event setup info
		const reco::Vertex                          &pv,              // event primary vertex
		const int                                    nvtx,            // number of primary vertices
		const std::vector<const reco::PFCandidate*> &pfNoPU,          // PFNoPU candidates
		const std::vector<TriggerRecord>            &triggerRecords,  // list of trigger names and objects
		const trigger::TriggerEvent                 &triggerEvent);   // event trigger objects
  
    protected:
      void computeIso(const reco::GsfElectron &ele, const double extRadius,
                      const std::vector<const reco::PFCandidate*> &pfNoPU,
                      float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;
      
      double evalEleIDMVA(const reco::GsfElectron &ele, EcalClusterLazyTools &lazyTools);
      
      
      // Electron cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fEleName;
      std::string fPFCandName;
      std::string fTrackName;
      std::string fConvName;
      std::string fRhoName;
      std::string fEBSCName;
      std::string fEESCName;
      std::string fEBRecHitName;
      std::string fEERecHitName;

      // Electron ID MVA
      EGammaMvaEleEstimator fEleIDMVA;
    
      // Electron momentum corrector
      bool fDoEleCorr;
      baconhep::ElectronMomentumCorrector fEleCorr;
  };
}
#endif
