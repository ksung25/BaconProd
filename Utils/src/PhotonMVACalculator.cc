#include "BaconProd/Utils/interface/PhotonMVACalculator.hh"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "TMVA/Reader.h"
#include <TFile.h>

using namespace baconhep;



//--------------------------------------------------------------------------------------------------
PhotonMVACalculator::PhotonMVACalculator():
  fIsInitialized(false),
  fReaderB       (0),
  fReaderE       (0)
{}

//--------------------------------------------------------------------------------------------------
PhotonMVACalculator::~PhotonMVACalculator()
{
  delete fReaderB;
  delete fReaderE;
}

//--------------------------------------------------------------------------------------------------
void PhotonMVACalculator::initialize(const std::string weightFileB,const std::string weightFileE) 
{
  fIsInitialized = true;
  fReaderB = new TMVA::Reader("!Color:!Silent");
  initialize(fReaderB,weightFileB,false);
  fReaderE = new TMVA::Reader("!Color:!Silent");
  initialize(fReaderE,weightFileE,true);
}
//--------------------------------------------------------------------------------------------------
void PhotonMVACalculator::initialize(TMVA::Reader* iReader, const std::string iWeightFile,bool iEndcap) {
  iReader->AddVariable("ph.scrawe",                   (float *)0);
  iReader->AddVariable("ph.r9",                       (float *)0);
  iReader->AddVariable("ph.sigietaieta",              (float *)0);
  iReader->AddVariable("ph.scetawidth",               (float *)0);
  iReader->AddVariable("ph.scphiwidth",               (float *)0);
  iReader->AddVariable("ph.idmva_CoviEtaiPhi",        (float *)0);
  iReader->AddVariable("ph.idmva_s4ratio",            (float *)0);
  iReader->AddVariable("ph.idmva_GammaIso",           (float *)0);
  iReader->AddVariable("ph.idmva_ChargedIso_selvtx"  ,(float *)0);
  iReader->AddVariable("ph.idmva_ChargedIso_worstvtx",(float *)0);
  iReader->AddVariable("ph.sceta",                    (float *)0);
  iReader->AddVariable("rho",                         (float *)0);
  if(iEndcap) iReader->AddVariable("ph.idmva_PsEffWidthSigmaRR",(float *)0);
  iReader->BookMVA("BDTG", iWeightFile);
}
//--------------------------------------------------------------------------------------------------
float PhotonMVACalculator::mvaValue(const reco::Photon &iPhoton,EcalClusterLazyTools &iLazyTools, const double &rho, const float &iGammaIso,const float &iCHadIso,const float &iBadIso,const float &iRR) {
  std::vector<float> mvaInput;
  const reco::SuperClusterRef sc = iPhoton.superCluster();
  std::vector<float> vCov = iLazyTools.localCovariances(*(sc->seed()));
  float lSCE   = sc->energy();
  float lR9    = iLazyTools.e3x3(*(sc->seed())) / (sc->rawEnergy());
  float lSIE   = iPhoton.sigmaIetaIeta();
  float lSCEW  = sc->etaWidth();
  float lSCPW  = sc->phiWidth();
  float lSEP   = vCov[1];
  float lS4R   = iLazyTools.e2x2(*(sc->seed())) / iLazyTools.e5x5(*(sc->seed()));
  float lPFIso = iGammaIso;
  float lQIsoG = iCHadIso;
  float lQIsoB = iBadIso;
  float lEta   = sc->eta();
  float lRho   = rho;
  mvaInput.push_back(lSCE);
  mvaInput.push_back(lR9);
  mvaInput.push_back(lSIE);
  mvaInput.push_back(lSCEW);
  mvaInput.push_back(lSCPW);
  mvaInput.push_back(lSEP);
  mvaInput.push_back(lS4R);
  mvaInput.push_back(lPFIso);
  mvaInput.push_back(lQIsoG);
  mvaInput.push_back(lQIsoB);
  mvaInput.push_back(lEta);
  mvaInput.push_back(lRho);
  if(fabs(lEta) > 1.5)  mvaInput.push_back(iRR);
  if(fabs(lEta) > 1.5)  return fReaderE->EvaluateMVA(mvaInput,"BDTG");
  return fReaderB->EvaluateMVA(mvaInput,"BDTG");
}
