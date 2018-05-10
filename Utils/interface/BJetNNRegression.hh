#ifndef BACONPROD_UTILS_BJETNNREGRESSION_HH
#define BACONPROD_UTILS_BJETNNREGRESSION_HH

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <string>

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

namespace baconhep {


  class BJetNNRegression  {
    public:
      BJetNNRegression();
      ~BJetNNRegression();
      
      void initialize(const std::string iWeightFile,float mean,float std);
      void SetNNVectorVar();
      std::pair<float,float> EvaluateNN();

      float Jet_pt ;
      float Jet_eta ;
      float rho ;
      float Jet_mt ;
      float Jet_leadTrackPt ;
      float Jet_leptonPtRel ;
      float Jet_leptonDeltaR ;
      float Jet_neHEF ;
      float Jet_neEmEF ;
      float Jet_vtxPt ;
      float Jet_vtxMass ;
      float Jet_vtx3dL ;
      float Jet_vtxNtrk ;
      float Jet_vtx3deL ;
      float Jet_numDaughters_pt03 ;
      float Jet_energyRing_dR0_em_Jet_e ;
      float Jet_energyRing_dR1_em_Jet_e ;
      float Jet_energyRing_dR2_em_Jet_e ;
      float Jet_energyRing_dR3_em_Jet_e ;
      float Jet_energyRing_dR4_em_Jet_e ;
      float Jet_energyRing_dR0_neut_Jet_e ;
      float Jet_energyRing_dR1_neut_Jet_e ;
      float Jet_energyRing_dR2_neut_Jet_e ;
      float Jet_energyRing_dR3_neut_Jet_e ;
      float Jet_energyRing_dR4_neut_Jet_e ;
      float Jet_energyRing_dR0_ch_Jet_e ;
      float Jet_energyRing_dR1_ch_Jet_e ;
      float Jet_energyRing_dR2_ch_Jet_e ;
      float Jet_energyRing_dR3_ch_Jet_e ;
      float Jet_energyRing_dR4_ch_Jet_e ;
      float Jet_energyRing_dR0_mu_Jet_e ;
      float Jet_energyRing_dR1_mu_Jet_e ;
      float Jet_energyRing_dR2_mu_Jet_e ;
      float Jet_energyRing_dR3_mu_Jet_e ;
      float Jet_energyRing_dR4_mu_Jet_e ;
      float Jet_chHEF;//implement from here
      float Jet_chEmEF;
      float Jet_leptonPtRelInv;
      int isEle;
      int isMu;
      int isOther;
      float Jet_mass;
      float Jet_withPtd;

    private:
      tensorflow::Session* session;
      tensorflow::GraphDef* graphDef;
      std::vector<float> NNvectorVar_; 
      //add vector of mva for eache jet

      float y_mean, y_std;

      //float y_mean = 1.0454729795455933;
      //float y_std  = 0.31628304719924927;
      //std::string fBjetNNRegressionPBFile_ = "BaconProd/Utils/data/breg_training_2017.pb";
  };
}
#endif
