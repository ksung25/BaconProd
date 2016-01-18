#ifndef BACONPROD_NTUPLER_FILLERELECTRON_HH
#define BACONPROD_NTUPLER_FILLERELECTRON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
namespace trigger {
  class TriggerEvent;
}
namespace pat {
  class Electron;
}


namespace baconhep
{
  class FillerElectron
  {
    public:
      FillerElectron(const edm::ParameterSet &iConfig, const bool useAOD);
      ~FillerElectron();

      // === filler for AOD ===
      void fill(TClonesArray                     *array,                        // output array to be filled
                const edm::Event                 &iEvent,                       // event info
		const edm::EventSetup            &iSetup,                       // event setup info
		const reco::Vertex               &pv,                           // event primary vertex
		const std::vector<TriggerRecord> &triggerRecords,               // list of trigger names and objects
		const trigger::TriggerEvent      &triggerEvent);                // event trigger objects

      // === filler for MINIAOD ===
      void fill(TClonesArray                                 *array,            // output array to be filled
                const edm::Event                             &iEvent,           // event info
                const edm::EventSetup                        &iSetup,           // event setup info
                const reco::Vertex                           &pv,               // event primary vertex
                const std::vector<TriggerRecord>             &triggerRecords,   // list of trigger names and objects
                const pat::TriggerObjectStandAloneCollection &triggerObjects);  // event trigger objects

    protected:
    double dEtaInSeed(const reco::GsfElectron& ele);
    double dEtaInSeed(const pat::Electron& ele);      
    void computeIso(double &iEta,double &iPhi, const double extRadius,
		    const reco::PFCandidateCollection    &puppi,
		    float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;    
  
      // Electron cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fEleName;
      std::string fBSName;
      std::string fPFCandName;
      std::string fTrackName;
      std::string fConvName;
      std::string fSCName;
      // Puppi
      std::string fPuppiName; 
      std::string fPuppiNoLepName; 
      bool fUsePuppi;

      // PF cluster isolation info (not in AOD)
      edm::InputTag fEcalPFClusterIsoMapTag;
      edm::InputTag fHcalPFClusterIsoMapTag;

      bool fUseAOD;
  };
}
#endif
