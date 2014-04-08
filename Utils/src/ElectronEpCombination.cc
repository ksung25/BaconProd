#include "BaconProd/Utils/interface/ElectronEpCombination.hh"
#include "Cintex/Cintex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <TFile.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
void ElectronEpCombination::initialize(const char *regressionFilename)
{
  ROOT::Cintex::Cintex::Enable();  // Needed to read non-TObject classes from a ROOT file!
  
  TFile file(regressionFilename);
  fForest = (GBRForest*)file.Get("CombinationWeight");
  assert(fForest);
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> ElectronEpCombination::evaluate(
  const reco::GsfElectron *ele,
  const double             energy,
  const double             energyError,
  const bool               printDebug)
{
  assert(ele);
  assert(fIsInitialized);
  
  double momentum      = ele->trackMomentumAtVtx().R();
  double momentumError = ele->trackMomentumError();
  
  // compute relative errors and ratio of errors
  double energyRelError   = energyError/energy;
  double momentumRelError = momentumError/momentum;
  double errorRatio       = energyRelError/momentumRelError;
  
  // calculate E/p and corresponding error
  double eOverP = energy/momentum;
  double eOverPerror = sqrt( (energyError/momentum)*(energyError/momentum) +
                             (energy*momentumError/momentum/momentum)*(energy*momentumError/momentum/momentum) );
  
  // fill input variables
  float regressionInputs[11];
  regressionInputs[0]  = energy;
  regressionInputs[1]  = energyRelError;
  regressionInputs[2]  = momentum;
  regressionInputs[3]  = momentumRelError;
  regressionInputs[4]  = errorRatio;
  regressionInputs[5]  = eOverP;
  regressionInputs[6]  = eOverPerror;
  regressionInputs[7]  = (ele->ecalDrivenSeed()) ? 1. : 0.;
  regressionInputs[8]  = (ele->trackerDrivenSeed()) ? 1. : 0.;
  regressionInputs[9]  = (float)(ele->classification());
  regressionInputs[10] = (ele->isEB()) ? 1. : 0.;
  
  // retrieve combination weight
  double weight = 0.;
  if( eOverP>0.025 && 
      fabs(momentum-energy)<15.*sqrt(momentumError*momentumError + energyError*energyError) )  // protect against crazy track measurement
  {  
    weight = fForest->GetResponse(regressionInputs);
    if(weight>1.) weight = 1.;
    else if(weight<0.) weight = 0.;
  }
  
  double combMom = -1;
  double combMomError = -1;
  if(!ele->ecalDrivenSeed()) {
    combMom = ele->p();
    
    // track momentum error not properly determined for pure tracker electrons (error = 999)
    // Recompute using parametrization from RecoEgamma/EgammaElectronAlgos/src/ElectronEnergyCorrector.cc::simpleParameterizationUncertainty()
    double eleMom = ele->p()<15. ? 15. : ele->p();
    if(ele->isEB()) {
      float parEB[3] = { 5.24e-02, 2.01e-01, 1.00e-02 };
      combMomError = eleMom * sqrt( pow(parEB[0]/sqrt(eleMom),2) + pow(parEB[1]/eleMom,2) + pow(parEB[2],2) ); 
    } else {
      float parEE[3] = { 1.46e-01, 9.21e-01, 1.94e-03} ;
      combMomError = eleMom * sqrt( pow(parEE[0]/sqrt(eleMom),2) + pow(parEE[1]/eleMom,2) + pow(parEE[2],2) );
    }
  } 
  
  if(momentumError!=999 || weight==0) {
    combMom      = weight*momentum + (1.-weight)*energy;
    combMomError = sqrt( weight*weight*momentumError*momentumError + (1.-weight)*(1.-weight)*energyError*energyError );    
  }
  
  if(printDebug) {    
    std::cout << "[ElectronEpCombination]" << std::endl;
    for(uint v=0; v < 11; ++v) std::cout << v << "=" << regressionInputs[v] << ", ";
    std::cout << std::endl;

    std::cout << "input energy   = " << energy   << " +/- " << energyError   << std::endl;
    std::cout << "input momentum = " << momentum << " +/- " << momentumError << std::endl;
    std::cout << "weight = " << weight << std::endl;
    std::cout << "combined momentum = " << combMom << std::endl;
    std::cout << "combined momentum uncertainty = " << combMomError << std::endl;
  }
    
  return std::pair<double,double>(combMom,combMomError);
}
