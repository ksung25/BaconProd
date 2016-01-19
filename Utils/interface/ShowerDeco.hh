#ifndef BaconProd_Utils_ShowerDeco_hh
#define BaconProd_Utils_ShowerDeco_hh
#include <vector>
#include "fastjet/PseudoJet.hh"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Message.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/TopGluonModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/BackgroundModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/ISRModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Deconstruct.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/ParseUtils.h"


class ShowerDeco { 
 public:
  ShowerDeco(std::string iConfig="",double iMicrojetConeSize=-1,int iNMaxMicroJets=10);
  ~ShowerDeco();
  double chi(double iPt, std::vector<fastjet::PseudoJet> &iParts);
  
 private:
  double fMicrojetConeSize;
  int    fNMaxMicrojets;
  Deconstruction::Deconstruct*     fDeconstruct;
  AnalysisParameters*              fParam;
  Deconstruction::TopGluonModel*   fSignal;
  Deconstruction::BackgroundModel* fBackground;
  Deconstruction::ISRModel*        fISR;
};

#endif
