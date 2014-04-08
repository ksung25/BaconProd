#ifndef BACONPROD_NTUPLER_FILLERGENINFO_HH
#define BACONPROD_NTUPLER_FILLERGENINFO_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
class TClonesArray;


namespace baconhep
{
  class TGenEventInfo;  // foward declaration

  class FillerGenInfo
  {
    public:
      FillerGenInfo(const edm::ParameterSet &iConfig);
      ~FillerGenInfo();
      
      void fill(TGenEventInfo    *genEvtInfo,     // output object to be filled
                TClonesArray     *particlesArr,   // output array of particles to be filled
		const edm::Event &iEvent);        // EDM event info
  
    protected:  
      
      // EDM object collection names
      std::string fGenEvtInfoName;
      std::string fGenParName;
      bool        fFillAll;
  };
}
#endif
