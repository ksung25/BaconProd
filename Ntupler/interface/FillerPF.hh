#ifndef BACONPROD_NTUPLER_FILLERPF_HH
#define BACONPROD_NTUPLER_FILLERPF_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"

class TClonesArray;


namespace baconhep
{
  class FillerPF
  {
    public:
      FillerPF(const edm::ParameterSet &iConfig);
      ~FillerPF();
      
       void fill(TClonesArray       *array,    // output array to be filled
		 TClonesArray       *iVtxCol,
		 const edm::Event   &iEvent);  // event info

       void fillMiniAOD(TClonesArray       *array,    // output array to be filled
			TClonesArray       *iVtxCol,
			const edm::Event   &iEvent);  // event info
   
    protected:
      //Useful tools
      float depthDeltaR(const reco::PFCandidate *iPF,const reco::PFRecHitCollection &iPFCol,double iDR=0.08);
      float timeDeltaR (const reco::PFCandidate *iPF,const reco::PFRecHitCollection &iPFCol,double iDR=0.08);
      float depth(const reco::PFCandidate *iPF);
      float time (const reco::PFCandidate *iPF);
      
      // EDM object collection names
      std::string fPFName;
      std::string fPVName;
      bool        fAddDepthTime;
  };
}
#endif
