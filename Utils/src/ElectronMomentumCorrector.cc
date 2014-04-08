#include "BaconProd/Utils/interface/ElectronMomentumCorrector.hh"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <iostream>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
void ElectronMomentumCorrector::initialize(
  const char   *regressionFilename,			     // regression weights file name
  const ElectronEnergyRegression::RegressionType type,       // regression type
  const ElectronEnergySmearingScaling::DatasetType dataset,  // dataset type
  const int	corrType,  			             // correction type
  const char   *scalesFilename, 			     // scale correction
  const char   *smearsType1Filename,			     // resolution correction type 1
  const char   *smearsType2Filename,			     // resolution correction type 2
  const char   *smearsType3Filename,			     // resolution correction type 3
  const char   *linearityFilename,                           // linearity correction data
  const bool	doRand,			                     // flag to toggle randomization
  const int	seed, 		                             // seed for randomization
  const double  lumiRatio)  			             // fraction of total luminosity from 2012D
{
  fRegression.initialize(regressionFilename, type);
  fSmearScale.initialize(dataset, corrType, scalesFilename, smearsType1Filename, smearsType2Filename, smearsType3Filename, doRand, seed, lumiRatio);
  fEpCombine.initialize(regressionFilename);
  fLinearity.initialize(dataset, linearityFilename);
  
  fIsInitialized = fRegression.isInitialized() && fSmearScale.isInitialized() && fEpCombine.isInitialized() && fLinearity.isInitialized();
}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> ElectronMomentumCorrector::ElectronMomentumCorrector::evaluate(
  const reco::GsfElectron *ele,
  const double             rho,
  const int	           nvertices,
  const unsigned int       runNum,
  const edm::EventSetup   &iSetup,
  EcalClusterLazyTools    &lazyTools,
  const bool               printDebug)
{
  if(printDebug) {
    std::cout << "[ElectronMomentumCorrector]" << std::endl;
    std::cout << " Electron pt = " << ele->pt() << " eta = " << ele->eta() << " phi = " << ele->phi() << std::endl;
  }
  
  std::pair<double,double>
  result = fRegression.evaluate(ele,
                                rho,
				nvertices,
				iSetup,
				lazyTools,
				printDebug);  
  
  result = fSmearScale.evaluate(ele,
                                result.first,
			        result.second,
			        runNum,
				lazyTools,
				printDebug);
  
  result = fEpCombine.evaluate(ele,
                               result.first,
			       result.second,
			       printDebug);
  
  result.first *= fLinearity.corrScale(ele, result.first, printDebug);
  
  return result;
}
