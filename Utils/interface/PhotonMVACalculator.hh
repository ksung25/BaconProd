#ifndef BACONPROD_UTILS_PHOTONMVACALCULATOR_H
#define BACONPROD_UTILS_PHOTONMVACALCULATOR_H

#include <vector>
#include <string>

// forward class declarations
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
namespace TMVA {
  class Reader;
}

namespace baconhep {   
  class PhotonMVACalculator 
  {
    public:      
      PhotonMVACalculator();
      ~PhotonMVACalculator(); 
    
      void initialize(const std::string weightFileB,const std::string weightFileE);
      void initialize(TMVA::Reader* iReader, const std::string iWeightFileB,bool iEndcap);
      
      bool isInitialized() const {return fIsInitialized;}
      
      float mvaValue(const reco::Photon &iPhoton,EcalClusterLazyTools &iLazyTools, const double &rho, const float &iGammaIso,const float &iCHadIso,const float &iBadIso,const float &iRR);
      
    
    protected:  
      bool fIsInitialized;
  
      TMVA::Reader *fReaderB;
      TMVA::Reader *fReaderE;
  };
}
#endif
