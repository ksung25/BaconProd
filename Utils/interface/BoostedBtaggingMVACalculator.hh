#ifndef BACONPROD_UTILS_BOOSTEDBTAGGINGMVACALCULATOR_HH
#define BACONPROD_UTILS_BOOSTEDBTAGGINGMVACALCULATOR_HH

#include <string>

// forward class declarations
namespace TMVA {
  class Reader;
}

namespace baconhep {

  class BoostedBtaggingMVACalculator
  {
    public:
      
      BoostedBtaggingMVACalculator();
      ~BoostedBtaggingMVACalculator();
      
      void initialize(
                      const std::string MethodTag, const std::string WeightFile);
      
      bool isInitialized() const {return fIsInitialized;}
      
      float mvaValue(
	 	     		     const float massPruned, const float flavour, const int nbHadrons, const float ptPruned, const float etaPruned,
                                     const float SubJet_csv,const float z_ratio, const float trackSipdSig_3, const float trackSipdSig_2, const float trackSipdSig_1,
                                     const float trackSipdSig_0, const float trackSipdSig_1_0, const float trackSipdSig_0_0, const float trackSipdSig_1_1,
                                     const float trackSipdSig_0_1, const float trackSip2dSigAboveCharm_0, const float trackSip2dSigAboveBottom_0,
                                     const float trackSip2dSigAboveBottom_1, const float tau0_trackEtaRel_0, const float tau0_trackEtaRel_1, const float tau0_trackEtaRel_2,
                                     const float tau1_trackEtaRel_0, const float tau1_trackEtaRel_1, const float tau1_trackEtaRel_2, const float tau_vertexMass_0,
                                     const float tau_vertexEnergyRatio_0, const float tau_vertexDeltaR_0, const float tau_flightDistance2dSig_0, const float tau_vertexMass_1,
                                     const float tau_vertexEnergyRatio_1, const float tau_flightDistance2dSig_1, const int jetNTracks, const int nSV,
		     		     const bool printDebug=false);
     
    
    private:
      void initReader(TMVA::Reader *reader, const std::string filename);
      
      bool fIsInitialized;
      
      TMVA::Reader *fReader;
      std::string fMethodTag;
      // input variables to compute MVA value
      //
      float _SubJet_csv, _z_ratio , _trackSipdSig_3 ,_trackSipdSig_2,_trackSipdSig_1,_trackSipdSig_0,_trackSipdSig_1_0,_trackSipdSig_0_0,_trackSipdSig_1_1,_trackSipdSig_0_1,_trackSip2dSigAboveCharm_0,_trackSip2dSigAboveBottom_0,_trackSip2dSigAboveBottom_1,_tau0_trackEtaRel_0,_tau0_trackEtaRel_1,_tau0_trackEtaRel_2,_tau1_trackEtaRel_0,_tau1_trackEtaRel_1,_tau1_trackEtaRel_2,_tau_vertexMass_0,_tau_vertexEnergyRatio_0,_tau_vertexDeltaR_0,_tau_flightDistance2dSig_0,_tau_vertexMass_1,_tau_vertexEnergyRatio_1,_tau_flightDistance2dSig_1,_massPruned,_flavour,_ptPruned,_etaPruned;
      int _jetNTracks,_nSV, _nbHadrons;
  };
}
#endif
