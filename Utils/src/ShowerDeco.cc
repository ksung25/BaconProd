#include "BaconProd/Utils/interface/ShowerDeco.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <iostream>

ShowerDeco::ShowerDeco(std::string iConfig,double iMicrojetConeSize,int iNMaxMicroJets) {
  // set up shower deconstruction stuff
  fMicrojetConeSize = iMicrojetConeSize;
  fNMaxMicrojets    = iNMaxMicroJets;
  fParam       = new AnalysisParameters(iConfig);
  fSignal      = new Deconstruction::TopGluonModel  (*fParam);
  fBackground  = new Deconstruction::BackgroundModel(*fParam);
  fISR         = new Deconstruction::ISRModel       (*fParam);
  fDeconstruct = new Deconstruction::Deconstruct    (*fParam, *fSignal, *fBackground, *fISR);
}
ShowerDeco::~ShowerDeco(){}

double ShowerDeco::chi(double iPt, std::vector<fastjet::PseudoJet> &iParts) { 
  // shower deconstruction
  double microconesize = fMicrojetConeSize;
  if  (fMicrojetConeSize<0){
    if      (iPt < 500) microconesize=0.3;
    else if (iPt < 700) microconesize=0.2;
    else                microconesize=0.1;
  }
  fastjet::JetDefinition reclustering(fastjet::JetAlgorithm::kt_algorithm, microconesize);
  fastjet::ClusterSequence * cs_micro = new fastjet::ClusterSequence(iParts, reclustering);
  std::vector<fastjet::PseudoJet> microjets = fastjet::sorted_by_pt(cs_micro->inclusive_jets(10.));
  if (int(microjets.size())>fNMaxMicrojets) microjets.erase(microjets.begin()+fNMaxMicrojets,microjets.end());
  double lChi = -1;
  try {
    double Psignal     = 0.0;
    double Pbackground = 0.0;
    lChi        = fDeconstruct->deconstruct(microjets, Psignal, Pbackground);
  } catch(Deconstruction::Exception &e) {
    std::cout << "Exception while running SD: " << e.what() << std::endl;
  }
  if (cs_micro->inclusive_jets(0.).size()>0) cs_micro->delete_self_when_unused();
  delete cs_micro;
  return lChi;
}

