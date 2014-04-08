//--------------------------------------------------------------------------------------------------
//
// ElectronEnergySmearingScaling
//
// Apply electron energy scale and resolution corrections
//
// Based on ElectronEnergyCalibrator class in the package:
//   V00-00-08 /CMSSW/EgammaAnalysis/ElectronTools 
//
//--------------------------------------------------------------------------------------------------

#ifndef BACONPROD_UTILS_ELECTRONENERGYSMEARINGSCALING_HH
#define BACONPROD_UTILS_ELECTRONENERGYSMEARINGSCALING_HH

#include <vector>
#include <utility>

// number of categories
#define NCAT 8

class TRandom3;
class EcalClusterLazyTools;
namespace reco { class GsfElectron; }

namespace baconhep {
  
  class CorrRecord {
    public:      
      CorrRecord():runNumMin(0),runNumMax(999999) {  
        for(int i=0; i<NCAT; i++)
	  corr[i]=1.;
      }
      ~CorrRecord(){}
      
      unsigned int runNumMin;  // run range minimum
      unsigned int runNumMax;  // run range maximum
      
      // corrections by category
      //  0: barrel, |eta|<1, R9< 0.94
      //  1: barrel, |eta|<1, R9>=0.94
      //  2: barrel, |eta|>1, R9< 0.94
      //  3: barrel, |eta|>1, R9>=0.94
      //  4: endcap, |eta|<2, R9< 0.94
      //  5: endcap, |eta|<2, R9>=0.94
      //  6: endcap, |eta|>2, R9< 0.94
      //  7: endcap, |eta|>2, R9>=0.94
      double corr[NCAT]; 
  };  
  
  class ElectronEnergySmearingScaling {
    public:
      ElectronEnergySmearingScaling();
      ~ElectronEnergySmearingScaling();
      
      enum DatasetType {
        kSummer11,                // 2011 MC
	kFall11,                  // 2011 MC
	kSummer12,                // 2012 MC	
	kSummer12_DR53X_HCP2012,  // 2012 MC
	kSummer12_LegacyPaper,    // 2012 MC
	
	kReReco,                  // 2011 Data
	kJan16ReReco,	          // 2011 Data        
        kICHEP2012,               // 2012 Data
	kMoriond2013,	          // 2012 Data
	k22Jan2013ReReco          // 2012 Data
      };

      void initialize(const DatasetType  dataset,                // dataset type
		      const int          corrType,               // correction type
                      const char        *scalesFilename,         // scale correction
                      const char        *smearsType1Filename,    // resolution correction type 1
                      const char        *smearsType2Filename,    // resolution correction type 2
                      const char        *smearsType3Filename,    // resolution correction type 3
		      const bool         doRand    = true,       // flag to toggle randomization
		      const int          seed      = 0xDEADBEEF, // seed for randomization
		      const double       lumiRatio = 0);         // fraction of total luminosity from 2012D
      
      bool isInitialized() const {return fIsInitialized;}
      
      std::pair<double,double> evaluate(
        const reco::GsfElectron *ele,        // pointer to electron object
        const double	         energy,     // electron energy
	const double	         error,      // eletron energy uncertainty
	const unsigned int       runNum,     // run number
	EcalClusterLazyTools    &lazyTools,  // class to compute ECAL cluster quantities
	const bool               printDebug=false);
  
    protected:
      void loadCorrections(const char *infilename, std::vector<CorrRecord> &records);
      
      bool fIsInitialized;
      
      bool      fDoRand;  	                    // flag to do randomization based on smearing factor
      TRandom3 *fRand;                              // random number generator
      
      int         fCorrType;                        // correction type
      DatasetType fDatasetType;                     // dataset type
      double      fLumiRatio;                       // ratio of luminosity for 2012A,B,C to 2012D
      
      std::vector<CorrRecord> fScalesRecords;       // data for scale corrections
      std::vector<CorrRecord> fSmearsType1Records;  // data for type 1 resolution corrections
      std::vector<CorrRecord> fSmearsType2Records;  // data for type 2 resolution corrections
      std::vector<CorrRecord> fSmearsType3Records;  // data for type 3 resolution corrections
  };
}
#endif
