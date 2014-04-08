#ifndef BACONPROD_UTILS_TAUISOMVACALCULATOR_H
#define BACONPROD_UTILS_TAUISOMVACALCULATOR_H

#include <vector>
#include <string>

// forward class declarations
#include "DataFormats/TauReco/interface/PFTauFwd.h"
namespace TMVA {
  class Reader;
}
class GBRForest;

namespace baconhep {

  class TauIsoRings
  {
    public:
      TauIsoRings(){}
      ~TauIsoRings(){}
      
      std::vector<int> niso;
      std::vector< std::vector<float> > rings;
      std::vector< std::vector<float> > shapes;
      
      std::vector<float> getVector() {
	std::vector<float> all;
	all.reserve(33);
      
	for(unsigned int i=0; i<niso.size(); i++)   { all.push_back(niso[i]); }       
	for(unsigned int i=0; i<rings.size(); i++)  { all.insert(all.end(), rings[i].begin(), rings[i].end()); }      
	for(unsigned int i=0; i<shapes.size(); i++) { all.insert(all.end(), shapes[i].begin(), shapes[i].end()); }
	
	return all;
      }
  };
  
  class TauIsoMVACalculator 
  {
    public:      
      TauIsoMVACalculator();
      ~TauIsoMVACalculator(); 
      
      void initialize(const std::string weightFile, const bool useGBR=false);
      
      bool isInitialized() const {return fIsInitialized;}
      
      float mvaValue(const reco::PFTau &tau, const double rho);
      
      baconhep::TauIsoRings computeIsoRings(const reco::PFTau &tau);
    
    protected:  
      bool fIsInitialized;
      bool fUseGBR;
  
      TMVA::Reader *fReader;
      GBRForest    *fGBRReader;
  };
}
#endif
