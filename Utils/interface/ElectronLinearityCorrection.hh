//--------------------------------------------------------------------------------------------------
//
// ElectronLinearityCorrection
//
// Apply electron energy scale and resolution corrections
//
// Based on ElectronEnergyCalibrator class in the package:
//   V00-00-08 /CMSSW/EgammaAnalysis/ElectronTools 
//
//--------------------------------------------------------------------------------------------------

#ifndef BACONPROD_UTILS_ELECTRONLINEARITYCORRECTION_HH
#define BACONPROD_UTILS_ELECTRONLINEARITYCORRECTION_HH

#include "BaconProd/Utils/interface/ElectronEnergySmearingScaling.hh"
#include <vector>

#define NLINCAT 6

namespace reco { class GsfElectron; }

namespace baconhep {
  
  class LinearityRecord {
    public:      
      LinearityRecord():ptMin(0),ptMax(999999) {  
        for(int i=0; i<NLINCAT; i++)
	  corr[i]=1.;
      }
      ~LinearityRecord(){}
      
      unsigned int ptMin;  // pT range minimum
      unsigned int ptMax;  // pT range maximum
      
      // linearity corrections by category
      //  0: barrel, classification<2, |eta|<1
      //  1: barrel, classification<2, |eta|>1
      //  2: endcap, classification<2
      //  3: barrel, classification>2, |eta|<1
      //  4: barrel, classification>2, |eta|>1
      //  5: endcap, classification>2
      double corr[NLINCAT]; 
  };  
  
  
  class ElectronLinearityCorrection {
    public:
      ElectronLinearityCorrection();
      ~ElectronLinearityCorrection();
      
      void initialize(const ElectronEnergySmearingScaling::DatasetType dataset,  // dataset type
		      const char *infilename);                                   // linearity correction data
      
      bool isInitialized() const {return fIsInitialized;}
      
      double corrScale(const reco::GsfElectron *ele,	            // pointer to electron object
                       const double             momentum,           // electron momentum
	               const bool               printDebug=false);
  
    protected:
      bool fIsInitialized;
      ElectronEnergySmearingScaling::DatasetType fDatasetType;  // dataset type
      std::vector<LinearityRecord> fLinearityRecords;           // data for linearity corrections
  };
}
#endif
