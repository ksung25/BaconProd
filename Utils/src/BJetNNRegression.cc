#include "BaconProd/Utils/interface/BJetNNRegression.hh"
#include "BaconProd/Utils/interface/JetTools.hh"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <iostream>
#include <cmath>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
BJetNNRegression::BJetNNRegression()
{
  NNvectorVar_.clear();
}

//--------------------------------------------------------------------------------------------------
BJetNNRegression::~BJetNNRegression() {
  tensorflow::closeSession(session);
  delete graphDef;
}

void BJetNNRegression::initialize(const std::string iWeightFile,float mean,float std){
  graphDef= tensorflow::loadGraphDef(iWeightFile.c_str());
  session = tensorflow::createSession(graphDef);
  y_mean = mean;
  y_std = std;
}

void BJetNNRegression::SetNNVectorVar(){

    NNvectorVar_.clear();
    NNvectorVar_.push_back(Jet_pt) ;//0
    NNvectorVar_.push_back(Jet_eta) ;
    NNvectorVar_.push_back(rho) ;
    NNvectorVar_.push_back(Jet_mt) ;
    NNvectorVar_.push_back(Jet_leadTrackPt) ;
    NNvectorVar_.push_back(Jet_leptonPtRel) ;//5
    NNvectorVar_.push_back(Jet_leptonDeltaR) ;
    NNvectorVar_.push_back(Jet_neHEF) ;
    NNvectorVar_.push_back(Jet_neEmEF) ;
    NNvectorVar_.push_back(Jet_vtxPt) ;
    NNvectorVar_.push_back(Jet_vtxMass) ;//10
    NNvectorVar_.push_back(Jet_vtx3dL) ;
    NNvectorVar_.push_back(Jet_vtxNtrk) ;
    NNvectorVar_.push_back(Jet_vtx3deL) ;
    NNvectorVar_.push_back(Jet_numDaughters_pt03) ;//this variable has changed order, in bdt it was last, check why
    NNvectorVar_.push_back(Jet_energyRing_dR0_em_Jet_e) ;//15
    NNvectorVar_.push_back(Jet_energyRing_dR1_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_neut_Jet_e) ;//20
    NNvectorVar_.push_back(Jet_energyRing_dR1_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_ch_Jet_e) ;//25
    NNvectorVar_.push_back(Jet_energyRing_dR1_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_mu_Jet_e) ;//30
    NNvectorVar_.push_back(Jet_energyRing_dR1_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_chHEF);//35
    NNvectorVar_.push_back(Jet_chEmEF);
    NNvectorVar_.push_back(Jet_leptonPtRelInv);
    NNvectorVar_.push_back(isEle);
    NNvectorVar_.push_back(isMu);
    NNvectorVar_.push_back(isOther);//40
    NNvectorVar_.push_back(Jet_mass);
    NNvectorVar_.push_back(Jet_withPtd);

}

std::pair<float,float> BJetNNRegression::EvaluateNN(){
    tensorflow::Tensor input(tensorflow::DT_FLOAT, {1,(unsigned int)NNvectorVar_.size()});//was {1,35} but get size mismatch, CHECK
    for (unsigned int i = 0; i < NNvectorVar_.size(); i++){
        //std::cout<<"i:"<<i<<" x:"<<NNvectorVar_[i]<<std::endl;
        input.matrix<float>()(0,i) =  float(NNvectorVar_[i]);
    }
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session, { { "ffwd_inp:0",input } }, { "ffwd_out/BiasAdd:0" }, &outputs);


    std::pair<float,float> corr_reg;
    corr_reg.first = y_mean+(outputs[0].matrix<float>()(0, 0)*y_std);
    corr_reg.second = 0.5*(outputs[0].matrix<float>()(0, 2)-outputs[0].matrix<float>()(0, 1))*y_std;
    //std::cout<<"\tReg: "<<corr_reg.first<<" = "<<y_mean<<"+("<<outputs[0].matrix<float>()(0, 0)<<"*"<<y_std<<") ,  "<<corr_reg.second<<" = 0.5*("<<outputs[0].matrix<float>()(0, 2)<<"-"<<outputs[0].matrix<float>()(0, 1)<<")*"<<y_std<<std::endl;

    return corr_reg;

}//end EvaluateNN
