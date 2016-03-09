#ifndef BACONPROD_NTUPLER_FILLERMUON_HH
#define BACONPROD_NTUPLER_FILLERMUON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerMuon
  {
    public:
       FillerMuon(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC);
      ~FillerMuon();

      // === filler for AOD ===
      void fill(TClonesArray			 *array,	                // output array to be filled
                const edm::Event		 &iEvent,	                // event info
	        const edm::EventSetup		 &iSetup,	                // event setup info
	        const reco::Vertex		 &pv,	                        // event primary vertex
	        const std::vector<TriggerRecord> &triggerRecords,               // list of trigger names and objects
	        const trigger::TriggerEvent	 &triggerEvent);                // event trigger objects

      // === filler for MINIAOD ===
      void fill(TClonesArray                                 *array,            // output array to be filled
                const edm::Event                             &iEvent,           // event info
                const edm::EventSetup                        &iSetup,           // event setup info
                const reco::Vertex                           &pv,               // event primary vertex
                const std::vector<TriggerRecord>             &triggerRecords,   // list of trigger names and objects
                const pat::TriggerObjectStandAloneCollection &triggerObjects);  // event trigger objects

    protected:

    void computeIso(double &iEta,double &iPhi, const double extRadius,
		    const reco::PFCandidateCollection    &puppi,
		    float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;

      // Muon cuts
      double fMinPt;

      // EDM object collection names
      std::string fMuonName;
      std::string fPFCandName;
      std::string fTrackName;

      // general tracks cuts
      bool   fSaveTracks;
      double fTrackMinPt;
      //Puppi
      std::string fPuppiName; 
      std::string fPuppiNoLepName; 
      bool fUsePuppi;
      bool fUseAOD;
      edm::EDGetTokenT<reco::MuonCollection>         fTokMuonName;
      edm::EDGetTokenT<pat::MuonCollection>          fTokPatMuonName;
      edm::EDGetTokenT<reco::PFCandidateCollection>  fTokPFCandName;
      edm::EDGetTokenT<reco::PFCandidateCollection>  fTokPuppiName;
      edm::EDGetTokenT<reco::PFCandidateCollection>  fTokPuppiNoLepName;
      edm::EDGetTokenT<reco::TrackCollection>        fTokTrackName;
  };
}
#endif
