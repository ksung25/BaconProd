//--------------------------------------------------------------------------------------------------
//
// ElectronEnergySmearingScaling
//
// Helper Class for applying electron energy scale and resolution corrections
//
//--------------------------------------------------------------------------------------------------

#ifndef BACONPROD_UTILS_ELECTRONEPCOMBINATION_HH
#define BACONPROD_UTILS_ELECTRONEPCOMBINATION_HH

#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include <utility>

namespace reco { class GsfElectron; }

namespace baconhep {

  class ElectronEpCombination {
    public:
      ElectronEpCombination():fForest(0){}
      ~ElectronEpCombination() { if(fForest) delete fForest; }
      
      void initialize(const char *regressionFilename);
      
      bool isInitialized() const {return fIsInitialized;}
      
      std::pair<double,double> evaluate(const reco::GsfElectron *ele,           // electron object
                                        const double             energy,        // electron energy
		                        const double             energyError,   // electron energy uncertainty
					const bool               printDebug=false);
    
    protected:
      bool fIsInitialized;
      GBRForest *fForest;
  };
}
#endif
