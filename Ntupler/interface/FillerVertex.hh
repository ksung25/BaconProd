#ifndef BACONPROD_NTUPLER_FILLERVERTEX_HH
#define BACONPROD_NTUPLER_FILLERVERTEX_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;


namespace baconhep
{
  class FillerVertex
  {
    public:
      FillerVertex(const edm::ParameterSet &iConfig);
      ~FillerVertex();
      
      const reco::Vertex* fill(TClonesArray       *array,    // output array to be filled
		               int                &nvtx,     // output number of primary vertices
		               const edm::Event   &iEvent);  // event info

      
      // EDM object collection names
      std::string fPVName;

      // Primary vertex selection
      unsigned int fMinNTracksFit;
      double       fMinNdof;
      double       fMaxAbsZ;
      double       fMaxRho;
  };
}
#endif
