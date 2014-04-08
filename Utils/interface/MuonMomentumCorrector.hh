//--------------------------------------------------------------------------------------------------
//
// MuonMomentumCorrector
//
// correct the four-momentum of the muon due to misalignment of the tracker
//
//--------------------------------------------------------------------------------------------------

#ifndef BACONPROD_UTILS_MUONMOMENTUMCORRECTOR_HH
#define BACONPROD_UTILS_MUONMOMENTUMCORRECTOR_HH
#include <TLorentzVector.h>

class MuScleFitCorrector;

namespace baconhep {

  class MuonMomentumCorrector
  {
    public:
      
      enum MuCorrType {
	kMuScleFall11_START42,              // MuScleFit corrections for 2011 MC
        kMuScleData2011_42X,                // MuScleFit corrections for 2011 Data
	kMuScleSummer12_DR53X_smearReReco,  // MuScleFit corrections for Summer12 MC
	kMuScleData53X_ReReco               // MuScleFit corrections for 2012 Data in 53X ReReco
      }; 

      MuonMomentumCorrector();
      ~MuonMomentumCorrector();
      
      void initialize(const MuCorrType type, const char *corrDataDir="", const bool doRand=true);
      
      bool isInitialized() const {return fIsInitialized;}
      
      TLorentzVector evaluate(const TLorentzVector &mu,                 // muon kinematics
                              const int             charge,             // muon charge
			      const unsigned int    runNum,             // run number
			      const bool            printDebug=false);
          
    protected:
      bool                fIsInitialized;      
      bool                fDoRand;              // flag to perform randomization
      MuCorrType          fType;                // correction type
      MuScleFitCorrector *fMuScleFitCorr;       // MuScleFit muon corrector
      MuScleFitCorrector *fMuScleFitCorr2012D;  // MuScleFit muon corrector for 2012D data
  };
}
#endif
