#ifndef BACONPROD_NTUPLER_FILLEREVENTINFO_HH
#define BACONPROD_NTUPLER_FILLEREVENTINFO_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class TEventInfo;  // foward declaration

  class FillerEventInfo
  {
    public:
      FillerEventInfo(const edm::ParameterSet &iConfig, const bool useAOD);
      ~FillerEventInfo();
      
      void fill(TEventInfo         *evtInfo,       // output object to be filled
                const edm::Event   &iEvent,        // EDM event info
		const reco::Vertex &pv,            // event primary vertex
		const bool          hasGoodPV,     // flag for if PV passing cuts is found
		const TriggerBits   triggerBits);  // bits for corresponding fired triggers
	       
    protected:
      void computeTrackMET(const unsigned int ipv,
			   const reco::VertexCollection *pvCol,
			   const reco::PFCandidateCollection *pfCandCol,
			   float &out_met, float &out_metphi);
    
      void computeTrackMET(const pat::PackedCandidateCollection *pfCandCol,
                           float &out_met, float &out_metphi);
    
      // EDM object collection names
      std::string fPFCandName;
      std::string fPUInfoName;
      std::string fPVName;
      std::string fBSName;
      std::string fCaloMETName;
      std::string fMETName;
      std::string fPFMETName;
      std::string fPFMETCName;
      std::string fMVAMETName;
      std::string fPUPPETName;
      std::string fPUPPETCName;
      std::string fPFMET30Name;
      std::string fPFMETC30Name;
      std::string fMVAMET30Name;
      std::string fPUPPET30Name;
      std::string fPUPPETC30Name;
      std::string fALPACAMETName;
      std::string fPALPACAMETName;
      std::string fRhoIsoName;
      std::string fRhoJetName;
      bool fUseFilters;
      bool fUseAOD;
  };
}
#endif
