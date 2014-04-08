#ifndef BACONPROD_NTUPLER_FILLERMUON_HH
#define BACONPROD_NTUPLER_FILLERMUON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/MuonMomentumCorrector.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerMuon
  {
    public:
      FillerMuon(const edm::ParameterSet &iConfig);
      ~FillerMuon();
      
      void fill(TClonesArray				    *array,	      // output array to be filled
                const edm::Event			    &iEvent,	      // event info
	        const edm::EventSetup			    &iSetup,	      // event setup info
	        const reco::Vertex			    &pv,	      // event primary vertex
	        const std::vector<const reco::PFCandidate*> &pfNoPU,	      // PFNoPU candidates
	        const std::vector<const reco::PFCandidate*> &pfPU,	      // PFPU candidates
	        const std::vector<TriggerRecord>	    &triggerRecords,  // list of trigger names and objects
	        const trigger::TriggerEvent		    &triggerEvent);   // event trigger objects
     
      void computeIso(const reco::Track &track, const double extRadius,
                      const std::vector<const reco::PFCandidate*> &pfNoPU,
		      const std::vector<const reco::PFCandidate*> &pfPU,
                      float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso, float &out_puIso) const;
      
      
      // Muon cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fMuonName;
      std::string fPFCandName;
      std::string fTrackName;

      // Muon momentum corrector
      bool fDoMuCorr;
      MuonMomentumCorrector fMuCorr;

      // general tracks cuts
      bool   fSaveTracks;
      double fTrackMinPt;
  };
}
#endif
