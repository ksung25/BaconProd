#include "BaconProd/Utils/interface/BoostedBtaggingMVACalculator.hh"
#include "TMVA/Reader.h"
#include <iostream>
#include <cmath>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
BoostedBtaggingMVACalculator::BoostedBtaggingMVACalculator():
  fIsInitialized(false),
  fReader(0),
  fMethodTag(""),
{}

//--------------------------------------------------------------------------------------------------
BoostedBtaggingMVACalculator::~BoostedBtaggingMVACalculator() {
  delete fReader;
  fIsInitialized = false;
}

//--------------------------------------------------------------------------------------------------
void BoostedBtaggingMVACalculator::initialize(
                                      const std::string MethodTag, const std::string WeightFile)
  fMethodTag  = MethodTag;
  
  if(WeightFile.length()>0) {
    if(fReader !=0) delete fReader;
    fReader = new TMVA::Reader();

    fReader->AddSpectator("massPruned",&_massPruned);
    fReader->AddSpectator("flavour",&_flavour);
    fReader->AddSpectator("nbHadrons",&_nbHadrons);
    fReader->AddSpectator("ptPruned",&_ptPruned);
    fReader->AddSpectator("etaPruned",&_etaPruned);
    fReader->AddVariable("SubJet_csv",&_SubJet_csv);
    fReader->AddVariable("z_ratio",&_z_ratio);
    fReader->AddVariable("trackSipdSig_3",&_trackSip3dSig_3);
    fReader->AddVariable("trackSipdSig_2",&_trackSip3dSig_2);
    fReader->AddVariable("trackSipdSig_1",&_trackSip3dSig_1);
    fReader->AddVariable("trackSipdSig_0",&_trackSip3dSig_0);
    fReader->AddVariable("trackSipdSig_1_0",&_tau2_trackSip3dSig_0);
    fReader->AddVariable("trackSipdSig_0_0",&_tau1_trackSip3dSig_0);
    fReader->AddVariable("trackSipdSig_1_1",&_tau2_trackSip3dSig_1);
    fReader->AddVariable("trackSipdSig_0_1",&_tau1_trackSip3dSig_1);
    fReader->AddVariable("trackSip2dSigAboveCharm_0",&_trackSip2dSigAboveCharm);
    fReader->AddVariable("trackSip2dSigAboveBottom_0",&_trackSip2dSigAboveBottom_0);
    fReader->AddVariable("trackSip2dSigAboveBottom_1",&_trackSip2dSigAboveBottom_1);
    fReader->AddVariable("tau1_trackEtaRel_0",&_tau2_trackEtaRel_0);
    fReader->AddVariable("tau1_trackEtaRel_1",&_tau2_trackEtaRel_1);
    fReader->AddVariable("tau1_trackEtaRel_2",&_tau2_trackEtaRel_2);
    fReader->AddVariable("tau0_trackEtaRel_0",&_tau1_trackEtaRel_0);
    fReader->AddVariable("tau0_trackEtaRel_1",&_tau1_trackEtaRel_1);
    fReader->AddVariable("tau0_trackEtaRel_2",&_tau1_trackEtaRel_2);
    fReader->AddVariable("tau_vertexMass_0",&_tau1_vertexMass);
    fReader->AddVariable("tau_vertexEnergyRatio_0",&_tau1_vertexEnergyRatio);
    fReader->AddVariable("tau_vertexDeltaR_0",&_tau1_vertexDeltaR);
    fReader->AddVariable("tau_flightDistance2dSig_0",&_tau1_flightDistance2dSig);
    fReader->AddVariable("tau_vertexMass_1",&_tau2_vertexMass);
    fReader->AddVariable("tau_vertexEnergyRatio_1",&_tau2_vertexEnergyRatio);
    fReader->AddVariable("tau_flightDistance2dSig_1",&_tau2_flightDistance2dSig);
    fReader->AddVariable("jetNTracks",&_jetNTracks);
    fReader->AddVariable("nSV",&_jetNSecondaryVertices);



    fReader->BookMVA(fMethodTag, WeightFile);
  }


fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
float BoostedBtaggingMVACalculator::mvaValue(
		const float massPruned, const float flavour, const int nbHadrons, const float ptPruned, const float etaPruned,
		const float SubJet_csv,const float z_ratio, const float trackSipdSig_3, const float trackSipdSig_2, const float trackSipdSig_1,
		const float trackSipdSig_0, const float trackSipdSig_1_0, const float trackSipdSig_0_0, const float trackSipdSig_1_1,
		const float trackSipdSig_0_1, const float trackSip2dSigAboveCharm_0, const float trackSip2dSigAboveBottom_0,
		const float trackSip2dSigAboveBottom_1, const float tau0_trackEtaRel_0, const float tau0_trackEtaRel_1, const float tau0_trackEtaRel_2,
		const float tau1_trackEtaRel_0, const float tau1_trackEtaRel_1, const float tau1_trackEtaRel_2, const float tau_vertexMass_0,
		const float tau_vertexEnergyRatio_0, const float tau_vertexDeltaR_0, const float tau_flightDistance2dSig_0, const float tau_vertexMass_1,
		const float tau_vertexEnergyRatio_1, const float tau_flightDistance2dSig_1, const int jetNTracks, const int nSV,
		const bool printDebug)
{
	_massPruned=massPruned;
	_flavour=flavour;
	_nbHadrons=nbHadrons;
	_ptPruned=ptPruned;
	_etaPruned=etaPruned;
	_z_ratio=z_ratio;
	_trackSipdSig_3=trackSip3dSig_3;
	_trackSipdSig_2=trackSip3dSig_2;
	_trackSipdSig_1=trackSip3dSig_1;
	_trackSipdSig_0=trackSip3dSig_0;
	_trackSipdSig_1_0=tau2_trackSip3dSig_0;
	_trackSipdSig_0_0=tau1_trackSip3dSig_0;
	_trackSipdSig_1_1=tau2_trackSip3dSig_1;
	_trackSipdSig_0_1=tau1_trackSip3dSig_1;
	_trackSip2dSigAboveCharm_0=trackSip2dSigAboveCharm;
	_trackSip2dSigAboveBottom_0=trackSip2dSigAboveBottom_0;
	_trackSip2dSigAboveBottom_1=trackSip2dSigAboveBottom_1;
	_tau1_trackEtaRel_0=tau2_trackEtaRel_0;
	_tau1_trackEtaRel_1=tau2_trackEtaRel_1;
	_tau1_trackEtaRel_2=tau2_trackEtaRel_2;
	_tau0_trackEtaRel_0=tau1_trackEtaRel_0;
	_tau0_trackEtaRel_1=tau1_trackEtaRel_1;
	_tau0_trackEtaRel_2=tau1_trackEtaRel_2;
	_tau_vertexMass_0=tau1_vertexMass;
	_tau_vertexEnergyRatio_0=tau1_vertexEnergyRatio;
	_tau_vertexDeltaR_0=tau1_vertexDeltaR;
	_tau_flightDistance2dSig_0=tau1_flightDistance2dSig;
	_tau_vertexMass_1=tau2_vertexMass;
	_tau_vertexEnergyRatio_1=tau2_vertexEnergyRatio;
	_tau_flightDistance2dSig_1=tau2_flightDistance2dSig;
	_jetNTracks=jetNTracks;
	_nSV=jetNSecondaryVertices;



	double val = -2;
	val = (fReader  ? fReader->EvaluateMVA(fMethodTag)   : -2); 

	if(printDebug) {
		std::cout << "[BoostedBtaggingMVACalculator]" << std::endl;
		/*   std::cout << "Inputs: nvtx= " << _nvtx;
		     std::cout << "  jetPt= " << _jetPt << "  jetEta= " << _jetEta << "  jetPhi= " << _jetPhi;
		     std::cout << "  |d0|= " << _d0 << "  |dZ|= " << _dZ;
		     std::cout << "  beta= " << _beta << "  betaStar= " << _betaStar;
		     std::cout << "  nCharged= " << _nCharged << "  nNeutrals= " << nNeutrals;
		     std::cout << "  dRMean= " << _dRMean << "  dR2Mean= " << _dR2Mean << "  ptD= " << _ptD;
		     std::cout << "  frac01= " << _frac01 << "  frac02= " << _frac02 << "  frac03= " << _frac03 << "  frac04= " << _frac04 << "  frac05= " << _frac05;
		     */  std::cout << std::endl;
		std::cout << " > MVA value = " << val << std::endl;
	}

	return val;
}
