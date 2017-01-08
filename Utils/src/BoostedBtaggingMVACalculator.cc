#include "BaconProd/Utils/interface/BoostedBtaggingMVACalculator.hh"
#include "TMVA/Reader.h"
#include <iostream>
#include <cmath>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
BoostedBtaggingMVACalculator::BoostedBtaggingMVACalculator():
  fIsInitialized(false),
  fReader(0),
  fMethodTag("")
{}

//--------------------------------------------------------------------------------------------------
BoostedBtaggingMVACalculator::~BoostedBtaggingMVACalculator() {
  delete fReader;
  fIsInitialized = false;
}

//--------------------------------------------------------------------------------------------------
void BoostedBtaggingMVACalculator::initialize(const std::string MethodTag, const std::string WeightFile,bool iUseSubJet)
{
   fMethodTag  = MethodTag;
  
  if(WeightFile.length()>0) {
    delete fReader;
    fReader = new TMVA::Reader();

    
    if(iUseSubJet) fReader->AddVariable("SubJet_csv",&_SubJet_csv);
    fReader->AddVariable("z_ratio",&_z_ratio);
    fReader->AddVariable("trackSipdSig_3",&_trackSipdSig_3);
    fReader->AddVariable("trackSipdSig_2",&_trackSipdSig_2);
    fReader->AddVariable("trackSipdSig_1",&_trackSipdSig_1);
    fReader->AddVariable("trackSipdSig_0",&_trackSipdSig_0);
    fReader->AddVariable("trackSipdSig_1_0",&_trackSipdSig_1_0);
    fReader->AddVariable("trackSipdSig_0_0",&_trackSipdSig_0_0);
    fReader->AddVariable("trackSipdSig_1_1",&_trackSipdSig_1_1);
    fReader->AddVariable("trackSipdSig_0_1",&_trackSipdSig_0_1);
    fReader->AddVariable("trackSip2dSigAboveCharm_0",&_trackSip2dSigAboveCharm_0);
    fReader->AddVariable("trackSip2dSigAboveBottom_0",&_trackSip2dSigAboveBottom_0);
    fReader->AddVariable("trackSip2dSigAboveBottom_1",&_trackSip2dSigAboveBottom_1);
    fReader->AddVariable("tau0_trackEtaRel_0",&_tau0_trackEtaRel_0);
    fReader->AddVariable("tau0_trackEtaRel_1",&_tau0_trackEtaRel_1);
    fReader->AddVariable("tau0_trackEtaRel_2",&_tau0_trackEtaRel_2);
    fReader->AddVariable("tau1_trackEtaRel_0",&_tau1_trackEtaRel_0);
    fReader->AddVariable("tau1_trackEtaRel_1",&_tau1_trackEtaRel_1);
    fReader->AddVariable("tau1_trackEtaRel_2",&_tau1_trackEtaRel_2);
    fReader->AddVariable("tau_vertexMass_0",&_tau_vertexMass_0);
    fReader->AddVariable("tau_vertexEnergyRatio_0",&_tau_vertexEnergyRatio_0);
    fReader->AddVariable("tau_vertexDeltaR_0",&_tau_vertexDeltaR_0);
    fReader->AddVariable("tau_flightDistance2dSig_0",&_tau_flightDistance2dSig_0);
    fReader->AddVariable("tau_vertexMass_1",&_tau_vertexMass_1);
    fReader->AddVariable("tau_vertexEnergyRatio_1",&_tau_vertexEnergyRatio_1);
    fReader->AddVariable("tau_flightDistance2dSig_1",&_tau_flightDistance2dSig_1);
    fReader->AddVariable("jetNTracks",&_jetNTracks);
    fReader->AddVariable("nSV",&_nSV);

    fReader->AddSpectator("massPruned",&_massPruned);
    fReader->AddSpectator("flavour",&_flavour);
    fReader->AddSpectator("nbHadrons",&_nbHadrons);
    fReader->AddSpectator("ptPruned",&_ptPruned);
    fReader->AddSpectator("etaPruned",&_etaPruned);




    fReader->BookMVA(fMethodTag, WeightFile);
  }


fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
float BoostedBtaggingMVACalculator::mvaValue(
		const float massPruned, const float flavour, const float nbHadrons, const float ptPruned, const float etaPruned,
		const float SubJet_csv,const float z_ratio, const float trackSipdSig_3, const float trackSipdSig_2, const float trackSipdSig_1,
		const float trackSipdSig_0, const float trackSipdSig_1_0, const float trackSipdSig_0_0, const float trackSipdSig_1_1,
		const float trackSipdSig_0_1, const float trackSip2dSigAboveCharm_0, const float trackSip2dSigAboveBottom_0,
		const float trackSip2dSigAboveBottom_1, const float tau0_trackEtaRel_0, const float tau0_trackEtaRel_1, const float tau0_trackEtaRel_2,
		const float tau1_trackEtaRel_0, const float tau1_trackEtaRel_1, const float tau1_trackEtaRel_2, const float tau_vertexMass_0,
		const float tau_vertexEnergyRatio_0, const float tau_vertexDeltaR_0, const float tau_flightDistance2dSig_0, const float tau_vertexMass_1,
		const float tau_vertexEnergyRatio_1, const float tau_flightDistance2dSig_1, const float jetNTracks, const float nSV,
		const bool printDebug)
{


	_massPruned=massPruned;
	_flavour=flavour;
	_nbHadrons=nbHadrons;
	_ptPruned=ptPruned;
	_etaPruned=etaPruned;
	_z_ratio=z_ratio;
	_trackSipdSig_3=trackSipdSig_3;
	_trackSipdSig_2=trackSipdSig_2;
	_trackSipdSig_1=trackSipdSig_1;
	_trackSipdSig_0=trackSipdSig_0;
	_trackSipdSig_1_0=trackSipdSig_1_0;
	_trackSipdSig_0_0=trackSipdSig_0_0;
	_trackSipdSig_1_1=trackSipdSig_1_1;
	_trackSipdSig_0_1=trackSipdSig_0_1;
	_trackSip2dSigAboveCharm_0=trackSip2dSigAboveCharm_0;
	_trackSip2dSigAboveBottom_0=trackSip2dSigAboveBottom_0;
	_trackSip2dSigAboveBottom_1=trackSip2dSigAboveBottom_1;
	_tau1_trackEtaRel_0=tau1_trackEtaRel_0;
	_tau1_trackEtaRel_1=tau1_trackEtaRel_1;
	_tau1_trackEtaRel_2=tau1_trackEtaRel_2;
	_tau0_trackEtaRel_0=tau0_trackEtaRel_0;
	_tau0_trackEtaRel_1=tau0_trackEtaRel_1;
	_tau0_trackEtaRel_2=tau0_trackEtaRel_2;
	_tau_vertexMass_0=tau_vertexMass_0;
	_tau_vertexEnergyRatio_0=tau_vertexEnergyRatio_0;
	_tau_vertexDeltaR_0=tau_vertexDeltaR_0;
	_tau_flightDistance2dSig_0=tau_flightDistance2dSig_0;
	_tau_vertexMass_1=tau_vertexMass_1;
	_tau_vertexEnergyRatio_1=tau_vertexEnergyRatio_1;
	_tau_flightDistance2dSig_1=tau_flightDistance2dSig_1;
	_jetNTracks=jetNTracks;
	_nSV=nSV;
	_SubJet_csv=SubJet_csv;



        double val = -2;
        val = (fReader  ? fReader->EvaluateMVA(fMethodTag)   : -2);

        if(printDebug) {
                std::cout << "[BoostedBtaggingMVACalculator]" << std::endl;
                   std::cout << "Inputs: nvtx= " << _nSV;
                     std::cout << "  jetNTracks= " << _jetNTracks << "  tau_flightDistance2dSig_1= " << _tau_flightDistance2dSig_1 << "  tau_vertexEnergyRatio_1= " << _tau_vertexEnergyRatio_1;
                     std::cout << "  trackSipdSig_3= " << _trackSipdSig_3 << "  trackSipdSig_2= " << _trackSipdSig_2<< "  trackSipdSig_1=  "<<_trackSipdSig_1<<"   trackSipdSig_0   "<<_trackSipdSig_0;
                     std::cout << "  trackSipdSig_1_0= "<<_trackSipdSig_1_0 <<" trackSipdSig_1_1 "<< _trackSipdSig_1_1<< " trackSipdSig_0_0= "<<_trackSipdSig_0_0<<" trackSipdSig_0_1 "<<_trackSipdSig_0_1;
                     std::cout << "  trackSip2dSigAboveCharm_0= " << _trackSip2dSigAboveCharm_0 << "  z_ratio= " << _z_ratio;
                     std::cout << "  trackSip2dSigAboveBottom_1= " << _trackSip2dSigAboveBottom_1 << "  trackSip2dSigAboveBottom_0= " << trackSip2dSigAboveBottom_0;
                     std::cout << "  tau0_trackEtaRel_2= " << _tau0_trackEtaRel_2 << "  tau0_trackEtaRel_1= " << _tau0_trackEtaRel_1 << "  tau0_trackEtaRel_0= " << _tau0_trackEtaRel_0;
                     std::cout << "  tau1_trackEtaRel_2= " << _tau1_trackEtaRel_2 << "  tau1_trackEtaRel_1= " << _tau1_trackEtaRel_1 << "  tau1_trackEtaRel_0= " << _tau1_trackEtaRel_0;
                     std::cout << "  tau_vertexMass_1= " << _tau_vertexMass_1 << "  tau_flightDistance2dSig_0= " << _tau_flightDistance2dSig_0 << "  tau_vertexDeltaR_0= " << _tau_vertexDeltaR_0 << "  tau_vertexEnergyRatio_0= " << _tau_vertexEnergyRatio_0 << "  tau_vertexMass_0= " << _tau_vertexMass_0;
                     std::cout << std::endl;
		     std::cout << " > MVA value = " << val << std::endl;
	}

	return val;
}
