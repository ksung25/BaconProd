#ifndef BACONPROD_UTILS_BJETREGRESSION_HH
#define BACONPROD_UTILS_BJETREGRESSION_HH

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <string>

// forward class declarations
namespace TMVA {
  class Reader;
}

namespace baconhep {

  class BJetRegression  {
    public:
      BJetRegression();
      ~BJetRegression();
      
    void initialize(const std::string iMethodTag, const std::string iPtWeightFile);
    bool isInitialized() const {return fIsInitialized;}
    float mvaValue(const float nvtx,      const float jetPt,    const float jetEta, const float jetMass,
		   const float leadTrack, const float lepPtRel, const float lepPt,  const float lepDR,
		   const float nHEF,      const float nEmEF,    const float vtxMass,
		   const float vtxPt,     const float vtx3dL,   const float vtxNtrk, const float vtx3deL,
		   const bool printDebug=false);    

    float mvaValue(int iNPV,double iCorr,const pat::Jet    &iJet);
    float mvaValue(int iNPV,double iCorr,const reco::PFJet &iJet,float iVtxPt,float iVtxMass,float iVtx3DVal,float iVtxNtracks,float iVtx3DeL);
    float mvaValue(int iNPV,double iCorr,const pat::Jet    &iJet,float iVtxPt,float iVtxMass,float iVtx3DVal,float iVtxNtracks,float iVtx3DeL);

  private:
    void initReader(TMVA::Reader *reader, const std::string filename);
    bool fIsInitialized;
    TMVA::Reader *fReader;
    std::string   fMethodTag;

    float _jetPt; 
    float _nvtx; 
    float _jetEta; 
    float _jetMass;
    float _leadTrackPt;
    float _leptonPtRel;
    float _leptonPt;
    float _leptonDR;
    float _nHEF;
    float _nEmEF;
    float _vtxPt;
    float _vtxMass;
    float _vtx3dL;
    float _vtxNtrk;
    float _vtx3deL;      
  };
}
#endif
