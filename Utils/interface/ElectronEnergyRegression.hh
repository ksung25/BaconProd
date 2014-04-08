//--------------------------------------------------------------------------------------------------
//
// ElectronEnergyRegression
//
// Helper Class for applying electron energy regression calculation
//
//--------------------------------------------------------------------------------------------------

#ifndef BACONPROD_UTILS_ELECTRONENERGYREGRESSION_HH
#define BACONPROD_UTILS_ELECTRONENERGYREGRESSION_HH

#include <utility>

class GBRForest;
class EcalClusterLazyTools;
namespace reco { class GsfElectron; }
namespace edm  { class EventSetup; }

namespace baconhep {      
  
  class ElectronEnergyRegression {
    public:
    
      enum RegressionType {
        kNoTrkVar,      // without tracker variables
        kNoTrkVarV1,    // without tracker variables V1
        kWithTrkVar,    // with tracker variables
        kWithTrkVarV1,  // with tracker variables V1
        kWithTrkVarV2,  // with tracker variables V2
        kWithSubCluVar  // with sub-cluster variables and without tracker variables
      };      
    
      ElectronEnergyRegression();
      ~ElectronEnergyRegression();

      void initialize(const char *regressionFilename, RegressionType type);

      bool isInitialized() const {return fIsInitialized;}
      
      std::pair<double,double> evaluate(const reco::GsfElectron *ele,                // electron object
                                        const double             rho,                // event energy density
					const int                nvertices,          // number of primary vertices
	                                const edm::EventSetup   &iSetup,             // event setup handle
					EcalClusterLazyTools    &lazyTools,          // class to compute ECAL cluster quantities
					const bool               printDebug=false);

    private:
      bool fIsInitialized;
      RegressionType fVersionType;
      
      GBRForest *forestCorrection_eb;   // pointer to energy correction forest for barrel
      GBRForest *forestCorrection_ee;   // pointer to energy correction forest for endcap
      GBRForest *forestUncertainty_eb;  // pointer to energy uncertainty forest for barrel
      GBRForest *forestUncertainty_ee;  // pointer to energy uncertainty forest for endcap
  };
}
#endif
