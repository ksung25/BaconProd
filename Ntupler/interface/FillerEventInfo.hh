#ifndef BACONPROD_NTUPLER_FILLEREVENTINFO_HH
#define BACONPROD_NTUPLER_FILLEREVENTINFO_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class TEventInfo;  // foward declaration
  class TSusyGen;    // ditto
  class FillerEventInfo
  {
    public:
      FillerEventInfo(const edm::ParameterSet &iConfig);
      ~FillerEventInfo();
      
      void fill(TEventInfo         *evtInfo,       // output object to be filled
                const edm::Event   &iEvent,        // EDM event info
		const reco::Vertex &pv,            // event primary vertex
		const bool          hasGoodPV,     // flag for if PV passing cuts is found
		const TriggerBits   triggerBits,   // bits for corresponding fired triggers
		TSusyGen           *susyGen=0);      // output for SUSY objects
	       
    protected:
      void computeTrackMET(const reco::Vertex &pv, 
                           const reco::PFCandidateCollection *pfCandCol,
                           float &out_met, float &out_metphi);
    
    
      // EDM object collection names
      std::string fPFCandName;
      std::string fPUInfoName;
      std::string fBSName;
      std::string fPFMETName;
      std::string fPFMETCName;
      std::string fMVAMETName;
      std::string fMVAMETUName;
      std::string fMVAMET0Name;
      std::string fRhoIsoName;
      std::string fRhoJetName;
      bool        fFillMET;
      bool        fFillMETFilters;
      bool        fAddSusyGen;
  };
}
#endif
