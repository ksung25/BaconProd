#ifndef BACONPROD_UTILS_JETPUIDMVACALCULATOR_HH
#define BACONPROD_UTILS_JETPUIDMVACALCULATOR_HH

#include <string>

// forward class declarations
namespace TMVA {
  class Reader;
}

namespace baconhep {

  class JetPUIDMVACalculator
  {
    public:
      enum JetIDType {
        kBaseline,
	k42,
	k52,
	kCut,
	k53,
	k53MET,
	k53METFull
      };
      
      JetPUIDMVACalculator();
      ~JetPUIDMVACalculator();
      
      void initialize(const JetPUIDMVACalculator::JetIDType jetIdType,
                      const std::string lowPtMethodTag, const std::string lowPtWeightFile,
                      const std::string highPtMethodTag, const std::string highPtWeightFile);
      
      bool isInitialized() const {return fIsInitialized;}
      
      float mvaValue(const float nvtx,
                     const float jetPt, const float jetEta, const float jetPhi,
		     const float d0, const float dZ,
		     const float beta, const float betaStar,
		     const float nCharged, const float nNeutrals,
		     const float dRMean, const float dR2Mean, const float ptD,
		     const float frac01, const float frac02, const float frac03, const float frac04, const float frac05,
		     const bool printDebug=false);
     
    
    private:
      void initReader(TMVA::Reader *reader, const std::string filename);
      
      bool fIsInitialized;
      
      TMVA::Reader *fLowPtReader;
      TMVA::Reader *fHighPtReader;
      std::string fLowPtMethodTag;
      std::string fHighPtMethodTag;
      
      // input variables to compute MVA values
      float _nvtx,
            _jetPt, _jetEta, _jetPhi,
	    _d0, _dZ,
            _beta, _betaStar,
            _nCharged, _nNeutrals,
            _dRMean, _dR2Mean, _ptD,
            _frac01, _frac02, _frac03, _frac04, _frac05;    
  };
}
#endif
