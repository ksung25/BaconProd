#include "BaconProd/Utils/interface/MuonMomentumCorrector.hh"
#include "MuScleFit/Calibration/interface/MuScleFitCorrector.h"
#include <iostream>
#include <cassert>

using namespace baconhep;

MuonMomentumCorrector::MuonMomentumCorrector():
fIsInitialized     (false),
fDoRand            (true),
fType              (kMuScleSummer12_DR53X_smearReReco),
fMuScleFitCorr     (0),
fMuScleFitCorr2012D(0)
{}

MuonMomentumCorrector::~MuonMomentumCorrector() {
  if(fMuScleFitCorr)      delete fMuScleFitCorr;
  if(fMuScleFitCorr2012D) delete fMuScleFitCorr2012D;
}

void MuonMomentumCorrector::initialize(const MuCorrType type, const char *corrDataDir, const bool doRand)
{  
  fDoRand = doRand;
  fType   = type;
  
  char fname[100];
  if(type == kMuScleFall11_START42) {
    sprintf(fname,"%s/MuScleFit_2011_MC_42X.txt",corrDataDir);
    fMuScleFitCorr = new MuScleFitCorrector(fname);
  
  } else if(type == kMuScleData2011_42X) {
    sprintf(fname,"%s/MuScleFit_2011_DATA_42X.txt",corrDataDir);
    fMuScleFitCorr = new MuScleFitCorrector(fname);
  
  } else if(type == kMuScleSummer12_DR53X_smearReReco) {
    sprintf(fname,"%s/MuScleFit_2012_MC_53X_smearReReco.txt",corrDataDir);
    fMuScleFitCorr = new MuScleFitCorrector(fname);
  
  } else if(type == kMuScleData53X_ReReco) {
    sprintf(fname,"%s/MuScleFit_2012ABC_DATA_ReReco_53X.txt",corrDataDir);
    fMuScleFitCorr = new MuScleFitCorrector(fname);
    sprintf(fname,"%s/MuScleFit_2012D_DATA_ReReco_53X.txt",corrDataDir);
    fMuScleFitCorr2012D = new MuScleFitCorrector(fname);
  
  } else {
    assert(0);
  }
  
  fIsInitialized = true;
}

TLorentzVector MuonMomentumCorrector::evaluate(const TLorentzVector &mu, const int charge, 
                                               const unsigned int runNum, const bool printDebug)
{
  assert(fIsInitialized);
  
  TLorentzVector p4(mu);
  
  if(fType == kMuScleFall11_START42) {
    assert(fMuScleFitCorr);
    fMuScleFitCorr->applyPtCorrection(p4, charge);
    fMuScleFitCorr->applyPtSmearing(p4, charge, !fDoRand);
  
  } else if(fType == kMuScleData2011_42X) {
    assert(fMuScleFitCorr);
    fMuScleFitCorr->applyPtCorrection(p4, charge);
  
  } else if(fType == kMuScleSummer12_DR53X_smearReReco) {
    assert(fMuScleFitCorr);
    fMuScleFitCorr->applyPtCorrection(p4, charge);
    fMuScleFitCorr->applyPtSmearing(p4, charge, !fDoRand);
  
  } else if(fType == kMuScleData53X_ReReco) {
    if(runNum < 203773) {
      assert(fMuScleFitCorr);
      fMuScleFitCorr->applyPtCorrection(p4, charge);
    } else {
      assert(fMuScleFitCorr2012D);
      fMuScleFitCorr2012D->applyPtCorrection(p4, charge);
    }
  }
    
  if(printDebug) {
    std::cout << "[MuonMomentumCorrector]" << std::endl;
    std::cout << "  Type: ";
    if     (fType == kMuScleFall11_START42)             { std::cout << "MuScleFit Fall11 MC" << std::endl; }
    else if(fType == kMuScleData2011_42X)               { std::cout << "MuScleFit 2011 Data" << std::endl; }
    else if(fType == kMuScleSummer12_DR53X_smearReReco) { std::cout << "MuScleFit Summer12 MC" << std::endl; }
    else if(fType == kMuScleData53X_ReReco)             { std::cout << "MuScleFit 2012 Data 53X ReReco" << std::endl; }
    std::cout << "  Run: " << runNum << std::endl;
    std::cout << "  >> Before: pT = " << mu.Pt() << ", eta = " << mu.Eta() << ", phi = " << mu.Phi() << std::endl;
    std::cout << "  <<  After: pT = " << p4.Pt() << ", eta = " << p4.Eta() << ", phi = " << p4.Phi() << std::endl;
  }
  
  return p4;
}
