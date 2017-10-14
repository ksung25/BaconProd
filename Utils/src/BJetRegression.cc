#include "BaconProd/Utils/interface/BJetRegression.hh"
#include "BaconProd/Utils/interface/JetTools.hh"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <cmath>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
BJetRegression::BJetRegression():
  fIsInitialized(false),
  fReader(0),
  fMethodTag("")
{}

//--------------------------------------------------------------------------------------------------
BJetRegression::~BJetRegression() {
  delete fReader;
  fIsInitialized = false;
}

//--------------------------------------------------------------------------------------------------
void BJetRegression::initialize(const std::string iMethodTag, const std::string iPtWeightFile) {
  fMethodTag  = iMethodTag;
  if(iPtWeightFile.length()>0) {
    if(fReader !=0) delete fReader;
    fReader = new TMVA::Reader();
    fReader->AddVariable("Jet_pt", &_jetPt); 
    fReader->AddVariable("nPVs",    &_nvtx); 
    fReader->AddVariable("Jet_eta", &_jetEta); 
    fReader->AddVariable("Jet_mt",  &_jetMass);
    fReader->AddVariable("Jet_leadTrackPt",&_leadTrackPt);
    fReader->AddVariable("Jet_leptonPtRel",&_leptonPtRel);
    fReader->AddVariable("Jet_leptonPt",&_leptonPt);
    fReader->AddVariable("Jet_leptonDeltaR",&_leptonDR);
    fReader->AddVariable("Jet_neHEF",&_nHEF);
    fReader->AddVariable("Jet_neEmEF",&_nEmEF);
    fReader->AddVariable("Jet_vtxPt",&_vtxPt);
    fReader->AddVariable("Jet_vtxMass",&_vtxMass);
    fReader->AddVariable("Jet_vtx3dL",&_vtx3dL);
    fReader->AddVariable("Jet_vtxNtrk",&_vtxNtrk);
    fReader->AddVariable("Jet_vtx3deL",&_vtx3deL);
    fReader->BookMVA(fMethodTag, iPtWeightFile);
  }
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
float BJetRegression::mvaValue(const float nvtx,      const float jetPt,    const float jetEta, const float jetMass,
			       const float leadTrack, const float lepPtRel, const float lepPt,  const float lepDR,
			       const float nHEF,      const float nEmEF,    const float vtxMass,
			       const float vtxPt,     const float vtx3dL,   const float vtxNtrk, const float vtx3deL,
			       const bool printDebug) {
  _nvtx        = nvtx;
  _jetPt       = jetPt;
  _jetEta      = jetEta;
  _jetMass     = jetMass;
  _leadTrackPt = leadTrack;
  _leptonPtRel = lepPtRel;
  _leptonPt    = lepPt;
  _leptonDR    = lepDR;
  _nHEF        = nHEF;
  _nEmEF       = nEmEF;
  _vtxMass     = vtxMass;
  _vtxPt       = vtxPt;
  _vtx3dL      = vtx3dL;
  _vtxNtrk     = vtxNtrk;
  _vtx3deL     = vtx3deL;
  
  double val = -2;
  if(fReader != 0) {
    val = fReader->EvaluateRegression(fMethodTag)[0];
  }
  if(printDebug) {
    std::cout << "[BJetRegression]" << std::endl;
    std::cout << "Inputs: nvtx= " << _nvtx;
    std::cout << "  jetPt= " << _jetPt << "  jetEta= " << _jetEta << "  jetMass= " << _jetMass;
    std::cout << "  leadTrackPt= " << _leadTrackPt << "  leptonPtRel= " << _leptonPtRel;
    std::cout << "  leptonPt= "    << _leptonPt    << "  leptonDR= "    << _leptonDR;
    std::cout << "  nHEF= "        << _nHEF        << "  nEMEF= "       << _nEmEF;
    std::cout << "  vtxPt= "       << _vtxPt       << "  nvtx3dL= "     << _vtx3dL;
    std::cout << "  vtxNtrk= "     << _vtxNtrk     << "  nvtx3deL= "    << _vtx3deL;
    std::cout << std::endl;
    std::cout << " > MVA value = " << val << std::endl;
  }
  return val;
}
float BJetRegression::mvaValue(int iNPV,double iCorr,const pat::Jet &iJet) { 
  float vtxPt      = sqrt(iJet.userFloat("vtxPx")*iJet.userFloat("vtxPx") + iJet.userFloat("vtxPy")*iJet.userFloat("vtxPy"));
  float vtxMass    = iJet.userFloat("vtxMass");
  float vtx3dL     = std::max(float(0.),iJet.userFloat("vtx3DVal"));
  float vtxNtracks = iJet.userFloat("vtxNtracks");
  float vtx3DSig   = iJet.userFloat("vtx3DSig");
  if(vtx3DSig > 0) vtx3DSig = vtx3dL/vtx3DSig;
  return mvaValue(iNPV,iJet.pt(),iJet.eta(),iCorr*iJet.mass(),
		  JetTools::leadPt(iJet),JetTools::leptons(iJet,1),JetTools::leptons(iJet,0),JetTools::leptons(iJet,2),
		  iJet.neutralEmEnergyFraction(),iJet.neutralHadronEnergyFraction(),    vtxMass,
		  vtxPt,     vtx3dL,   vtxNtracks, vtx3DSig);
} 
float BJetRegression::mvaValue(int iNPV,double iCorr,const reco::PFJet &iJet,float iVtxPt,float iVtxMass,float iVtx3DVal,float iVtxNtracks,float iVtx3DeL) { 
  float vtxPt      = iVtxPt;
  float vtxMass    = iVtxMass;
  float vtx3dL     = iVtx3DVal;
  float vtxNtracks = iVtxNtracks;
  float vtx3DeL    = iVtx3DeL;
  return mvaValue(iNPV,iJet.pt(),iJet.eta(),iCorr*iJet.mass(),
		  JetTools::leadPt(iJet),JetTools::leptons(iJet,1),JetTools::leptons(iJet,0),JetTools::leptons(iJet,2),
		  iJet.neutralEmEnergyFraction(),iJet.neutralHadronEnergyFraction(),    vtxMass,
		  vtxPt,     vtx3dL,   vtxNtracks, vtx3DeL);
} 


float BJetRegression::mvaValue(int iNPV,double iCorr,const pat::Jet &iJet,float iVtxPt,float iVtxMass,float iVtx3DVal,float iVtxNtracks,float iVtx3DeL) { 
  float vtxPt      = iVtxPt;
  float vtxMass    = iVtxMass;
  float vtx3dL     = iVtx3DVal;
  float vtxNtracks = iVtxNtracks;
  float vtx3DeL    = iVtx3DeL;
  return mvaValue(iNPV,iJet.pt(),iJet.eta(),iCorr*iJet.mass(),
		  JetTools::leadPt(iJet),JetTools::leptons(iJet,1),JetTools::leptons(iJet,0),JetTools::leptons(iJet,2),
		  iJet.neutralEmEnergyFraction(),iJet.neutralHadronEnergyFraction(),    vtxMass,
		  vtxPt,     vtx3dL,   vtxNtracks, vtx3DeL);
} 
