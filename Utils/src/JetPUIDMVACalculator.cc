#include "BaconProd/Utils/interface/JetPUIDMVACalculator.hh"
#include "TMVA/Reader.h"
#include <iostream>
#include <cmath>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
JetPUIDMVACalculator::JetPUIDMVACalculator():
  fIsInitialized(false),
  fLowPtReader(0),
  fHighPtReader(0),
  fLowPtMethodTag(""),
  fHighPtMethodTag("")
{}

//--------------------------------------------------------------------------------------------------
JetPUIDMVACalculator::~JetPUIDMVACalculator()
{
  delete fLowPtReader;
  delete fHighPtReader;
}

//--------------------------------------------------------------------------------------------------
void JetPUIDMVACalculator::initialize(const JetPUIDMVACalculator::JetIDType jetIDType,
                                      const std::string lowPtMethodTag, const std::string lowPtWeightFile,
                                      const std::string highPtMethodTag, const std::string highPtWeightFile)
{
  fLowPtMethodTag  = lowPtMethodTag;
  fHighPtMethodTag = highPtMethodTag;
  
  if(lowPtWeightFile.length()>0) {
    fLowPtReader = new TMVA::Reader();
    fLowPtReader->AddVariable("nvtx", &_nvtx); 
    if(jetIDType != k53) fLowPtReader->AddVariable("jetPt",  &_jetPt);  
    if(jetIDType != k53) fLowPtReader->AddVariable("jetEta", &_jetEta);
    if(jetIDType != k53) fLowPtReader->AddVariable("jetPhi", &_jetPhi);	     
    fLowPtReader->AddVariable("dZ", &_dZ);
    if(jetIDType != k53 && jetIDType != k53MET) fLowPtReader->AddVariable("d0", &_d0);
    fLowPtReader->AddVariable("beta",      &_beta);
    fLowPtReader->AddVariable("betaStar",  &_betaStar);
    fLowPtReader->AddVariable("nCharged",  &_nCharged);
    fLowPtReader->AddVariable("nNeutrals", &_nNeutrals);
    if(jetIDType != k53 && jetIDType != k53MET) fLowPtReader->AddVariable("dRMean",  &_dRMean);
    if(jetIDType == k53 || jetIDType == k53MET) fLowPtReader->AddVariable("dR2Mean", &_dR2Mean);
    if(jetIDType == k53 || jetIDType == k53MET) fLowPtReader->AddVariable("ptD",     &_ptD);
    fLowPtReader->AddVariable("frac01", &_frac01);
    fLowPtReader->AddVariable("frac02", &_frac02);
    fLowPtReader->AddVariable("frac03", &_frac03);
    fLowPtReader->AddVariable("frac04", &_frac04);
    fLowPtReader->AddVariable("frac05", &_frac05);
    
    if(jetIDType == k53) fLowPtReader->AddSpectator("jetPt",  &_jetPt);  
    if(jetIDType == k53) fLowPtReader->AddSpectator("jetEta", &_jetEta);
    if(jetIDType == k53) fLowPtReader->AddSpectator("jetPhi", &_jetPhi);  
    
    fLowPtReader->BookMVA(fLowPtMethodTag, lowPtWeightFile);
  }
  
  if(highPtWeightFile.length()>0) {
    fHighPtReader = new TMVA::Reader();
    if(jetIDType == kBaseline) {
      fHighPtReader->AddVariable("nvtx",     &_nvtx);
      fHighPtReader->AddVariable("jetPt",    &_jetPt);
      fHighPtReader->AddVariable("jetEta",   &_jetEta);
      fHighPtReader->AddVariable("jetPhi",   &_jetPhi);
      fHighPtReader->AddVariable("dZ",       &_dZ);
      fHighPtReader->AddVariable("d0",       &_d0);
      fHighPtReader->AddVariable("beta",     &_beta);
      fHighPtReader->AddVariable("betaStar", &_betaStar);
      fHighPtReader->AddVariable("nCharged", &_nCharged);
      fHighPtReader->AddVariable("nNeutrals",&_nNeutrals);
      fHighPtReader->AddVariable("dRMean",   &_dRMean);
      fHighPtReader->AddVariable("frac01",   &_frac01);
      fHighPtReader->AddVariable("frac02",   &_frac02);
      fHighPtReader->AddVariable("frac03",   &_frac03);
      fHighPtReader->AddVariable("frac04",   &_frac04);
      fHighPtReader->AddVariable("frac05",   &_frac05);
      
    } else if(jetIDType == k42) {
      fHighPtReader->AddVariable("nvtx",     &_nvtx);
      fHighPtReader->AddVariable("dZ",       &_dZ);
      fHighPtReader->AddVariable("beta",     &_beta);
      fHighPtReader->AddVariable("betaStar", &_betaStar);
      fHighPtReader->AddVariable("nCharged", &_nCharged);
      fHighPtReader->AddVariable("nNeutrals",&_nNeutrals);
      fHighPtReader->AddVariable("frac01",   &_frac01);
      fHighPtReader->AddVariable("frac02",   &_frac02);
      fHighPtReader->AddVariable("frac03",   &_frac03);
      fHighPtReader->AddVariable("frac04",   &_frac04);
      fHighPtReader->AddVariable("frac05",   &_frac05);

      fHighPtReader->AddSpectator("jetPt", &_jetPt);
      fHighPtReader->AddSpectator("jetEta",&_jetEta);     
 
    } else if(jetIDType == k52) {
      fHighPtReader->AddVariable("nvtx",     &_nvtx);
      fHighPtReader->AddVariable("dZ",       &_dZ);
      fHighPtReader->AddVariable("beta",     &_beta);
      fHighPtReader->AddVariable("betaStar", &_betaStar);
      fHighPtReader->AddVariable("nCharged", &_nCharged);
      fHighPtReader->AddVariable("nNeutrals",&_nNeutrals);
      fHighPtReader->AddVariable("dR2Mean",  &_dR2Mean);
      fHighPtReader->AddVariable("frac01",   &_frac01);
      fHighPtReader->AddVariable("frac02",   &_frac02);
      fHighPtReader->AddVariable("frac03",   &_frac03);
      fHighPtReader->AddVariable("frac04",   &_frac04);
      fHighPtReader->AddVariable("frac05",   &_frac05);

      fHighPtReader->AddSpectator("jetPt", &_jetPt);
      fHighPtReader->AddSpectator("jetEta",&_jetEta);
      
    } else if(jetIDType == k53) {
      fHighPtReader->AddVariable("nvtx",     &_nvtx);
      fHighPtReader->AddVariable("dZ",       &_dZ);
      fHighPtReader->AddVariable("beta",     &_beta);
      fHighPtReader->AddVariable("betaStar", &_betaStar);
      fHighPtReader->AddVariable("nCharged", &_nCharged);
      fHighPtReader->AddVariable("nNeutrals",&_nNeutrals);
      fHighPtReader->AddVariable("dR2Mean",  &_dR2Mean);
      fHighPtReader->AddVariable("ptD",      &_ptD);
      fHighPtReader->AddVariable("frac01",   &_frac01);
      fHighPtReader->AddVariable("frac02",   &_frac02);
      fHighPtReader->AddVariable("frac03",   &_frac03);
      fHighPtReader->AddVariable("frac04",   &_frac04);
      fHighPtReader->AddVariable("frac05",   &_frac05);
      
      fHighPtReader->AddSpectator("jetPt", &_jetPt);
      fHighPtReader->AddSpectator("jetEta",&_jetEta);
      fHighPtReader->AddSpectator("jetPhi",&_jetPhi);

    } else if(jetIDType == k53MET) {
      fHighPtReader->AddVariable("nvtx",     &_nvtx);
      fHighPtReader->AddVariable("jetPt",    &_jetPt);
      fHighPtReader->AddVariable("jetEta",   &_jetEta);
      fHighPtReader->AddVariable("jetPhi",   &_jetPhi);
      fHighPtReader->AddVariable("dZ",       &_dZ);
      fHighPtReader->AddVariable("beta",     &_beta);
      fHighPtReader->AddVariable("betaStar", &_betaStar);
      fHighPtReader->AddVariable("nCharged", &_nCharged);
      fHighPtReader->AddVariable("nNeutrals",&_nNeutrals);
      fHighPtReader->AddVariable("dR2Mean",  &_dR2Mean);
      fHighPtReader->AddVariable("ptD",      &_ptD);
      fHighPtReader->AddVariable("frac01",   &_frac01);
      fHighPtReader->AddVariable("frac02",   &_frac02);
      fHighPtReader->AddVariable("frac03",   &_frac03);
      fHighPtReader->AddVariable("frac04",   &_frac04);
      fHighPtReader->AddVariable("frac05",   &_frac05);
    }
    
    fHighPtReader->BookMVA(fHighPtMethodTag, highPtWeightFile);
  } 
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
float JetPUIDMVACalculator::mvaValue(const float nvtx, const float jetPt, const float jetEta, const float jetPhi,
		                     const float d0, const float dZ, const float beta, const float betaStar,
		                     const float nCharged, const float nNeutrals, const float dRMean, const float dR2Mean, const float ptD,
		                     const float frac01, const float frac02, const float frac03, const float frac04, const float frac05,
				     const bool printDebug)
{
  _nvtx      = nvtx;
  _jetPt     = jetPt;
  _jetEta    = jetEta;
  _jetPhi    = jetPhi;
  _d0        = fabs(d0);
  _dZ        = fabs(dZ);
  _beta      = beta;
  _betaStar  = betaStar;
  _nCharged  = nCharged;
  _nNeutrals = nNeutrals;
  _dRMean    = dRMean;
  _dR2Mean   = dR2Mean;
  _ptD       = ptD;
  _frac01    = frac01;
  _frac02    = frac02;
  _frac03    = frac03;
  _frac04    = frac04;
  _frac05    = frac05;
  
  double val = -2;
  if(jetPt < 10) { val = (fLowPtReader  ? fLowPtReader->EvaluateMVA(fLowPtMethodTag)   : -2); } 
  else           { val = (fHighPtReader ? fHighPtReader->EvaluateMVA(fHighPtMethodTag) : -2); }
  
  if(printDebug) {
    std::cout << "[JetPUIDMVACalculator]" << std::endl;
    std::cout << "Inputs: nvtx= " << _nvtx;
    std::cout << "  jetPt= " << _jetPt << "  jetEta= " << _jetEta << "  jetPhi= " << _jetPhi;
    std::cout << "  |d0|= " << _d0 << "  |dZ|= " << _dZ;
    std::cout << "  beta= " << _beta << "  betaStar= " << _betaStar;
    std::cout << "  nCharged= " << _nCharged << "  nNeutrals= " << nNeutrals;
    std::cout << "  dRMean= " << _dRMean << "  dR2Mean= " << _dR2Mean << "  ptD= " << _ptD;
    std::cout << "  frac01= " << _frac01 << "  frac02= " << _frac02 << "  frac03= " << _frac03 << "  frac04= " << _frac04 << "  frac05= " << _frac05;
    std::cout << std::endl;
    std::cout << " > MVA value = " << val << std::endl;
  }
  
  return val;
}
