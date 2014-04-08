#include "BaconProd/Utils/interface/ElectronEnergySmearingScaling.hh"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include <TRandom3.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
ElectronEnergySmearingScaling::ElectronEnergySmearingScaling():
fIsInitialized(false),fDoRand(true),fRand(0),fCorrType(1),fLumiRatio(0)
{}

//--------------------------------------------------------------------------------------------------
ElectronEnergySmearingScaling::~ElectronEnergySmearingScaling()
{
  if(fRand) delete fRand;
}

//--------------------------------------------------------------------------------------------------
void ElectronEnergySmearingScaling::initialize(const DatasetType  dataset,
					       const int          corrType,
                                               const char        *scalesFilename,
                                               const char        *smearsType1Filename,
                                               const char        *smearsType2Filename,
                                               const char        *smearsType3Filename,
					       const bool         doRand,
					       const int          seed,
					       const double       lumiRatio)
{
  fDatasetType = dataset;
  fCorrType    = corrType;
    
  loadCorrections(scalesFilename,      fScalesRecords);       // load scale corrections
  loadCorrections(smearsType1Filename, fSmearsType1Records);  // load type 1 resolution corrections
  loadCorrections(smearsType2Filename, fSmearsType2Records);  // load type 2 resolution corrections
  loadCorrections(smearsType3Filename, fSmearsType3Records);  // load type 3 resolution corrections
  
  fDoRand = doRand;
  if(fDoRand)
    fRand = new TRandom3(seed);
  
  fLumiRatio = lumiRatio;
  
  fIsInitialized = true;
}


//--------------------------------------------------------------------------------------------------
std::pair<double,double> ElectronEnergySmearingScaling::evaluate(
  const reco::GsfElectron *ele,	        // pointer to electron object
  const double             energy,	// electron energy
  const double             error,	// eletron energy uncertainty
  const unsigned int       runNum,	// run number
  EcalClusterLazyTools    &lazyTools,   // class to compute ECAL cluster quantities
  const bool               printDebug)  
{
  double scEta = ele->superCluster()->eta();
  double r9    = lazyTools.e3x3(*(ele->superCluster()->seed())) / ele->superCluster()->rawEnergy();
  bool   isEB  = ele->isEB();
  bool   isMC  = (fDatasetType==kSummer11) ||
                 (fDatasetType==kFall11) ||
		 (fDatasetType==kSummer12) ||
		 (fDatasetType==kSummer12_DR53X_HCP2012) ||
		 (fDatasetType==kSummer12_LegacyPaper);
  
  int icat = -1;   
  if(isEB) {
    if(fabs(scEta)<1) { icat = (r9<0.94) ? 0 : 1; }
    else              { icat = (r9<0.94) ? 2 : 3; }
  } else {
    if(fabs(scEta)<2) { icat = (r9<0.94) ? 4 : 5; }
    else              { icat = (r9<0.94) ? 6 : 7; }
  }
  assert(icat>-1);

  //
  // scale correction
  //
  int irun = -1;
  for(unsigned int irec=0; irec<fScalesRecords.size(); irec++) {
    if(fScalesRecords[irec].runNumMin <= runNum && runNum <= fScalesRecords[irec].runNumMax)
      irun = irec;
  }
  if(!isMC && irun<0) { // if data event is a run with unspecified corrections, return without corrections
    if(printDebug) {
      std::cout << "[ElectronEnergySmearingScaling]" << std::endl;
      std::cout << " * category = " << icat << std::endl;
      std::cout << " *   scale = 1" << std::endl;
      std::cout << "old energy = " << energy << " +/- " << error << std::endl;
      std::cout << "new energy = " << energy << " +/- " << error << std::endl;
    }

    return std::pair<double,double>(energy, error);
  }

  //
  // resolution correction
  //
  // (NOTE 1: Indices for fSmearsType*Records[] are hard-coded; make sure indices correspond to entries in input .csv files)
  // (NOTE 2: While smearing values are assigned for data and MC events, smearing is done for MC only;
  //          it's superfluous and may be a bit confusing, but staying close to 
  //          /CMSSW/EgammaAnalysis/ElectronTools/src/ElectronEnergyCalibrator.cc for now...)
  //
  double dsigMC=0;
  if(fCorrType==1) {
    if(fDatasetType==kFall11 || fDatasetType==kJan16ReReco) {
      dsigMC = fSmearsType1Records[0].corr[icat];
      
    } else if(fDatasetType==kSummer12_DR53X_HCP2012 || fDatasetType==kMoriond2013) {
      if(isMC) {
        if(fLumiRatio==0) {
          dsigMC = fSmearsType1Records[1].corr[icat];
        } else if(fLumiRatio==1) {
          dsigMC = fSmearsType1Records[2].corr[icat];
        } else {
          double rn = gRandom->Uniform();
          // NOTE: In /CMSSW/EgammaAnalysis/ElectronTools/src/ElectronEnergyCalibrator.cc the boundary
          //       cases fLumiRatio=0 and fLumiRatio=1 contradict the meaning of fLumiRatio for 
          //       intermediate values. The implementation below makes the meaning is consistent.
	  //       In CMSSW, the 2012D values are applied when rn > fLumiRatio ...
	  if(rn<fLumiRatio) dsigMC = fSmearsType1Records[1].corr[icat]; 
          else              dsigMC = fSmearsType1Records[2].corr[icat];
        }
      } else {
        dsigMC = (runNum<=203002) ? fSmearsType1Records[0].corr[icat] : fSmearsType1Records[1].corr[icat];
      }
    } else {
      assert(0);
    }
  
  } else if(fCorrType==2) {
    if(fDatasetType==kFall11 || fDatasetType==kJan16ReReco)
      dsigMC = fSmearsType2Records[0].corr[icat];    
    else if(fDatasetType==kSummer12_LegacyPaper || fDatasetType==k22Jan2013ReReco)
      dsigMC = fSmearsType2Records[1].corr[icat];    
    else
      assert(0);
    
  } else if(fCorrType==3) {
    if     (fDatasetType==kSummer11               || fDatasetType==kReReco)      { dsigMC = fSmearsType3Records[0].corr[icat]; }
    else if(fDatasetType==kFall11		  || fDatasetType==kJan16ReReco) { dsigMC = fSmearsType3Records[1].corr[icat]; }
    else if(fDatasetType==kSummer12  	          || fDatasetType==kICHEP2012)   { dsigMC = fSmearsType3Records[2].corr[icat]; }
    else if(fDatasetType==kSummer12_DR53X_HCP2012 || fDatasetType==kMoriond2013) { dsigMC = fSmearsType3Records[3].corr[icat]; }
    else
      assert(0);
  
  } else {
    assert(0);
  }
  
  double newEnergy = energy;  
  if(isMC) {
    double corrMC = fDoRand ? fRand->Gaus(1.,dsigMC) : (1.+dsigMC);
    newEnergy *= corrMC;
  } else {
    double scale = fScalesRecords[irun].corr[icat];
    newEnergy *= scale;
  }
  
  double newError = sqrt(error*error + dsigMC*dsigMC*energy*energy);

  if(printDebug) {    
    std::cout << "[ElectronEnergySmearingScaling]" << std::endl;
    std::cout << " * category = " << icat << std::endl;
    if(isMC) std::cout << " *   dsigMC = " << dsigMC << std::endl;
    else     std::cout << " *   scale = " << newEnergy/energy << std::endl;

    std::cout << "old energy = " << energy    << " +/- " << error    << std::endl;
    std::cout << "new energy = " << newEnergy << " +/- " << newError << std::endl;
  }
    
  return std::pair<double,double>(newEnergy, newError);
}

//--------------------------------------------------------------------------------------------------
void ElectronEnergySmearingScaling::loadCorrections(const char *infilename, std::vector<CorrRecord> &records)
{  
  //
  // Parse .csv file of run dependent corrections
  // Expected .csv file format:
  //    runNumMin,runNumMax,corr[0],corr[1],corr[2],corr[3],corr[4],corr[5],corr[6],corr[7]
  //
  std::ifstream ifs;
  std::string line;  
  CorrRecord rec;
  char c0,c1,c2,c3,c4,c5,c6,c7,c8;  // eat commas in CSV file

  ifs.open(infilename);
  if(!ifs.is_open()) { std::cout << "[ElectronEnergySmearingScaling] " << infilename << " not found!" << std::endl; assert(0); }
  while(std::getline(ifs,line)) {
    std::istringstream ss(line);
    ss >> rec.runNumMin >> c0 >> rec.runNumMax >> c1 
       >> rec.corr[0]   >> c2 >> rec.corr[1]   >> c3 >> rec.corr[2] >> c4 >> rec.corr[3] >> c5
       >> rec.corr[4]   >> c6 >> rec.corr[5]   >> c7 >> rec.corr[6] >> c8 >> rec.corr[7];
    records.push_back(rec);
  }
  ifs.close();
}
