#include "BaconProd/Utils/interface/ElectronLinearityCorrection.hh"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <fstream>
#include <string>
#include <sstream>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
ElectronLinearityCorrection::ElectronLinearityCorrection():fIsInitialized(false)
{}

//--------------------------------------------------------------------------------------------------
ElectronLinearityCorrection::~ElectronLinearityCorrection()
{}

//--------------------------------------------------------------------------------------------------
void ElectronLinearityCorrection::initialize(const ElectronEnergySmearingScaling::DatasetType dataset,  // dataset type
                                             const char *infilename)                                    // linearity correction data
{      
  fDatasetType = dataset;
  
  //
  // Parse .csv file of run dependent corrections
  // Expected .csv file format:
  //    ptMin, ptMax, corr[0], corr[1], corr[2], corr[3], corr[4], corr[5]
  //
  std::ifstream ifs;
  std::string line;  
  LinearityRecord rec;
  char c0,c1,c2,c3,c4,c5,c6;  // eat commas in CSV file
  
  // load lineary corrections
  ifs.open(infilename);
  assert(ifs.is_open());
  while(std::getline(ifs,line)) {
    std::istringstream ss(line);
    ss >> rec.ptMin   >> c0 >> rec.ptMax   >> c1 
       >> rec.corr[0] >> c2 >> rec.corr[1] >> c3 >> rec.corr[2] >> c4
       >> rec.corr[3] >> c5 >> rec.corr[4] >> c6 >> rec.corr[5];
    fLinearityRecords.push_back(rec);
  }
  ifs.close();
  
  fIsInitialized=true;
}

//--------------------------------------------------------------------------------------------------
double ElectronLinearityCorrection::corrScale(const reco::GsfElectron *ele,         // pointer to electron object
                                              const double             momentum,    // electron momentum
                                              const bool               printDebug)  
{
  bool   isEB           = ele->isEB();
  double scEta		= ele->superCluster()->eta();
  double theta  	= 2*atan(exp(-scEta));
  double pt		= momentum * fabs(sin(theta));
  int	 classification = ele->classification();
  bool   isMC  = (fDatasetType==ElectronEnergySmearingScaling::kSummer11) ||
                 (fDatasetType==ElectronEnergySmearingScaling::kFall11) ||
		 (fDatasetType==ElectronEnergySmearingScaling::kSummer12) ||
		 (fDatasetType==ElectronEnergySmearingScaling::kSummer12_DR53X_HCP2012) ||
		 (fDatasetType==ElectronEnergySmearingScaling::kSummer12_LegacyPaper);

  double corr = 0;
  
  if(!isMC) {
    for(unsigned int irec=0; irec < fLinearityRecords.size(); irec++) {
      if( (pt >= fLinearityRecords[irec].ptMin) && (pt <= fLinearityRecords[irec].ptMax) ) {
        if(isEB) {
          if(fabs(scEta) < 1) {
            if(classification<2) { corr = fLinearityRecords[irec].corr[0]; } 
	    else                 { corr = fLinearityRecords[irec].corr[3]; }        
	  } else {
            if(classification<2) { corr = fLinearityRecords[irec].corr[1]; }  
	    else                 { corr = fLinearityRecords[irec].corr[4]; }      
          }
      
        } else { // !isEB
          if(classification<2) { corr = fLinearityRecords[irec].corr[2]; } 
	  else                 { corr = fLinearityRecords[irec].corr[5]; }
        }
      }
    }
  }
  
  double scale = 1./(1.+corr);	  
  if(printDebug) 
  {
      std::cout << "[ElectronLinearityCorrection]" << std::endl;
      std::cout << "correction scale = " << scale << std::endl; 
  }
  
  return scale;
}
