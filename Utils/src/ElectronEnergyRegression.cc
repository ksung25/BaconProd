#include "BaconProd/Utils/interface/ElectronEnergyRegression.hh"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "Cintex/Cintex.h"
#include <TFile.h>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
ElectronEnergyRegression::ElectronEnergyRegression() : 
  fIsInitialized(kFALSE),
  fVersionType(kNoTrkVar),
  forestCorrection_eb(0), 
  forestCorrection_ee(0), 
  forestUncertainty_eb(0), 
  forestUncertainty_ee(0) {
}

//--------------------------------------------------------------------------------------------------
ElectronEnergyRegression::~ElectronEnergyRegression()
{
  if(forestCorrection_eb)  delete forestCorrection_eb;
  if(forestCorrection_ee)  delete forestCorrection_ee;
  if(forestUncertainty_eb) delete forestUncertainty_eb;
  if(forestUncertainty_ee) delete forestUncertainty_ee;
}

//--------------------------------------------------------------------------------------------------
void ElectronEnergyRegression::initialize(const char *regressionFilename, RegressionType type) {

  ROOT::Cintex::Cintex::Enable();  // Needed to read non-TObject classes from a ROOT file!
  
  // Load "forests"
  TFile file(regressionFilename);
  forestCorrection_eb  = (GBRForest*)file.Get("EBCorrection");  assert(forestCorrection_eb);
  forestCorrection_ee  = (GBRForest*)file.Get("EECorrection");  assert(forestCorrection_ee);
  forestUncertainty_eb = (GBRForest*)file.Get("EBUncertainty"); assert(forestUncertainty_eb);
  forestUncertainty_ee = (GBRForest*)file.Get("EEUncertainty"); assert(forestUncertainty_ee);  
  
  // Updating type and marking as initialized
  fVersionType   = type;
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> ElectronEnergyRegression::evaluate(const reco::GsfElectron *ele, double rho, int nvertices, 
                                                            const edm::EventSetup &iSetup, EcalClusterLazyTools &lazyTools,
							    bool printDebug)
{  
  assert(ele);
  
  // Checking if instance has been initialized
  assert(fIsInitialized);
  
  const reco::SuperClusterRef sc   = ele->superCluster();
  const reco::CaloClusterPtr& seed = sc->seed();  
  
  // compute shower shape variables for the supercluster seed
  std::vector<float> vCov = lazyTools.localCovariances(*seed);
  
  // compute ECAL local quantites for the supercluster seed
  float etacry, phicry, thetatilt, phitilt;
  int ieta, iphi;
  EcalClusterLocal local;
  if(seed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel) {
    local.localCoordsEB(*seed,iSetup,etacry,phicry,ieta,iphi,thetatilt,phitilt); 
  } else {
    local.localCoordsEE(*seed,iSetup,etacry,phicry,ieta,iphi,thetatilt,phitilt);
  }
  
  //
  // Base variables
  //
  double SCRawEnergy      = sc->rawEnergy();                        // SuperCluster raw energy
  double scEta            = sc->eta();                              // SuperCluster eta
  double scPhi            = sc->phi();                              // SuperCluster phi
  double R9               = lazyTools.e3x3(*seed)/sc->rawEnergy();  // SuperCluster R9
  double etawidth         = sc->etaWidth();                         // SuperCluster width in eta
  double phiwidth         = sc->phiWidth();                         // SuperCluster width in phi
  double NClusters        = sc->clustersSize();                     // number of basic clusters in the SuperCluster
  double HoE              = ele->hcalOverEcal();                    // electron H/E
  double EtaSeed          = seed->eta();                            // seed cluster eta
  double PhiSeed          = seed->phi();                            // seed cluster phi
  double ESeed            = seed->energy();                         // seed cluster energy
  double E3x3Seed         = lazyTools.e3x3(*seed);                  // seed cluster E3x3
  double E5x5Seed         = lazyTools.e5x5(*seed);                  // seed cluster E5x5
  double see              = isnan(vCov[0]) ? 0 : sqrt(vCov[0]);     // SuperCluster sigma_ieta_ieta
  double spp              = isnan(vCov[2]) ? 0 : sqrt(vCov[2]);     // SuperCluster sigma_iphi_iphi
  double sep              = isnan(vCov[1]) ? 0 : vCov[1]/see/spp;   // SuperCluster cov(ieta,iphi)/(sigma_ieta_ieta * sigma_iphi_iphi)
  double EMaxSeed         = lazyTools.eMax(*seed);                  // seed cluster EMax
  double E2ndSeed         = lazyTools.e2nd(*seed);                  // seed cluster E2nd
  double ETopSeed         = lazyTools.eTop(*seed);                  // seed cluster ETop
  double EBottomSeed      = lazyTools.eBottom(*seed);               // seed cluster EBottom
  double ELeftSeed        = lazyTools.eLeft(*seed);                 // seed cluster ELeft
  double ERightSeed       = lazyTools.eRight(*seed);                // seed cluster ERight
  double E2x5MaxSeed      = lazyTools.e2x5Max(*seed);               // seed cluster E2x5Max
  double E2x5TopSeed      = lazyTools.e2x5Top(*seed);               // seed cluster E2x5Top
  double E2x5BottomSeed   = lazyTools.e2x5Bottom(*seed);            // seed cluster E2x5Bottom
  double E2x5LeftSeed     = lazyTools.e2x5Left(*seed);              // seed cluster E2x5Left
  double E2x5RightSeed    = lazyTools.e2x5Right(*seed);             // seed cluster E2x5Right
  double IEtaSeed         = ieta;                                   // seed cluster IEta
  double IPhiSeed         = iphi;                                   // seed cluster IPhi
  double EtaCrySeed       = etacry;                                 // eta of highest energy crystal in seed cluster
  double PhiCrySeed       = phicry;                                 // phi of highest energy crystal in seed cluster
  double PreShowerOverRaw = sc->preshowerEnergy()/sc->rawEnergy();  // (SuperCluster PreShower energy) / (SuperCluster raw energy)  
  int	 IsEcalDriven     = ele->ecalDrivenSeed();                  // IsEcalDriven electron?
  
  //
  // Track variables
  //
  double pt                     = ele->pt();                                              // electron pT
  double GsfTrackPIn            = ele->trackMomentumAtVtx().R();                          // electron GSF track momentum at PCA to beam line
  double fbrem                  = fmax(ele->fbrem(), -1.0);                               // electron fbrem
  double Charge                 = ele->charge();                                          // electron charge
  double EoP                    = fmin(ele->eSuperClusterOverP(), 20.);                   // E/p  
  double TrackMomentumError     = ele->trackMomentumError();                              // electron track momentum error
  double EcalEnergyError        = ele->correctedEcalEnergyError();                        // electron ECAL energy error
  int	 Classification         = ele->classification();                                  // electron classification  
  double detaIn                 = fmin(fabs(ele->deltaEtaSuperClusterTrackAtVtx()),0.6);  // Delta-eta between SC and track position at PCA to beam line
  double dphiIn                 = ele->deltaPhiSuperClusterTrackAtVtx();                  // Delta-phi between SC and track position at PCA to beam line
  double detaCalo               = ele->deltaEtaSeedClusterTrackAtCalo();                  // Delta-eta between SC and track position at calorimeter extracted from outermost state
  double dphiCalo               = ele->deltaPhiSeedClusterTrackAtCalo();                  // Delta-phi between SC and track position at calorimeter extracted from outermost state
  double GsfTrackChiSqr         = ele->gsfTrack()->normalizedChi2();                      // GSF track fit normalized chisquare (chi2/ndf)
  double ElectronEnergyOverPout = fmin(ele->eEleClusterOverPout(), 20.);                  // (electron cluster energy)/(track momentum at calorimeter)
  
  // number of layers with measurements in CTF track
  const reco::TrackRef kfTrackRef = ele->closestCtfTrackRef();
  double KFTrackNLayers = kfTrackRef.isNonnull() ? kfTrackRef->hitPattern().trackerLayersWithMeasurement() : -1;
    
  //
  // Sub-cluster variables
  //
  double ESubs=0;  // Energy sum of all basic clusters except seed
  const reco::CaloCluster *sub1=0, *sub2=0, *sub3=0;
  for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus!=sc->clustersEnd(); ++itClus) {
    const reco::CaloCluster *clus = (*itClus).get();
    if(clus == seed.get()) continue;
    ESubs+=clus->energy();
    if(!sub1 || clus->energy() > sub1->energy()) {
      sub3 = sub2;
      sub2 = sub1;
      sub1 = clus;
      
    } else if(!sub2 || clus->energy() > sub2->energy()) {
      sub3 = sub2;
      sub2 = clus;
    
    } else if(!sub3 || clus->energy() > sub3->energy()) {
      sub3 = clus;
    }
  }

  // preshower clusters
  unsigned int nPsClus=0;
  double EPshwSubs=0;
  const reco::CaloCluster *pshwsub1=0, *pshwsub2=0, *pshwsub3=0;
  for(reco::CaloCluster_iterator itPsClus = sc->preshowerClustersBegin(); itPsClus!=sc->preshowerClustersEnd(); ++itPsClus) {
    const reco::CaloCluster* psclus = (*itPsClus).get();
    nPsClus++;
    EPshwSubs+=psclus->energy();
    if(!pshwsub1 || psclus->energy() > pshwsub1->energy()) {
      pshwsub3 = pshwsub2;
      pshwsub2 = pshwsub1;
      pshwsub1 = psclus;
      
    } else if(!pshwsub2 || psclus->energy() > pshwsub2->energy()) {
      pshwsub3 = pshwsub2;
      pshwsub2 = psclus;
    
    } else if(!pshwsub3 || psclus->energy() > pshwsub3->energy()) {
      pshwsub3 = psclus;
    }
  }
  
  double isEtaGap      = ele->isEBEtaGap();                     // is in EB eta module gap
  double isPhiGap      = ele->isEBPhiGap();                     // is in EB phi module gap
  double isDeeGap      = ele->isEEDeeGap();                     // is in EE Dee gap
  double ESub1         = sub1 ? sub1->energy()        : 0.;     // 1st sub-cluster energy
  double EtaSub1       = sub1 ? sub1->eta()           : 999.;   // 1st sub-cluster eta
  double PhiSub1       = sub1 ? sub1->phi()           : 999.;   // 1st sub-cluster phi
  double EMaxSub1      = sub1 ? lazyTools.eMax(*sub1) : 0.;     // 1st sub-cluster EMax
  double E3x3Sub1      = sub1 ? lazyTools.e3x3(*sub1) : 0.;     // 1st sub-cluster E3x3
  double ESub2         = sub2 ? sub2->energy()        : 0.;     // 2nd sub-cluster energy
  double EtaSub2       = sub2 ? sub2->eta()           : 999.;   // 2nd sub-cluster eta
  double PhiSub2       = sub2 ? sub2->phi()           : 999.;   // 2nd sub-cluster phi
  double EMaxSub2      = sub2 ? lazyTools.eMax(*sub2) : 0.;     // 2nd sub-cluster EMax
  double E3x3Sub2      = sub2 ? lazyTools.e3x3(*sub2) : 0.;     // 2nd sub-cluster E3x3
  double ESub3         = sub3 ? sub3->energy()        : 0.;     // 3rd sub-cluster energy
  double EtaSub3       = sub3 ? sub3->eta()           : 999.;   // 3rd sub-cluster eta
  double PhiSub3       = sub3 ? sub3->phi()           : 999.;   // 3rd sub-cluster phi
  double EMaxSub3      = sub3 ? lazyTools.eMax(*sub3) : 0.;     // 3rd sub-cluster EMax
  double E3x3Sub3      = sub3 ? lazyTools.e3x3(*sub3) : 0.;     // 3rd sub-cluster E3x3
  double NPshwClusters = nPsClus;                               // number of preshower clusters
  double EPshwSub1     = pshwsub1 ? pshwsub1->energy() : 0.;    // 1st preshower cluster energy
  double EtaPshwSub1   = pshwsub1 ? pshwsub1->eta()    : 999.;  // 1st preshower cluster eta
  double PhiPshwSub1   = pshwsub1 ? pshwsub1->phi()    : 999.;  // 1st preshower cluster phi
  double EPshwSub2     = pshwsub2 ? pshwsub2->energy() : 0.;    // 2nd preshower cluster energy
  double EtaPshwSub2   = pshwsub2 ? pshwsub2->eta()    : 999.;  // 2nd preshower cluster eta
  double PhiPshwSub2   = pshwsub2 ? pshwsub2->phi()    : 999.;  // 2nd preshower cluster phi
  double EPshwSub3     = pshwsub3 ? pshwsub3->energy() : 0.;    // 3rd preshower cluster energy
  double EtaPshwSub3   = pshwsub3 ? pshwsub3->eta()    : 999.;  // 3rd preshower cluster eta
  double PhiPshwSub3   = pshwsub3 ? pshwsub3->phi()    : 999.;  // 3rd preshower cluster phi
  
    
  bool isBarrel = false;
  uint nvars = 0;
  float *vals = 0; 
  
  if(fVersionType == kNoTrkVar) {  ///// Regression without tracker variables /////    
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 38 : 31;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IEtaSeed;
      vals[31] = IPhiSeed;
      vals[32] = ((int) IEtaSeed)%5;
      vals[33] = ((int) IPhiSeed)%2;
      vals[34] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[35] = ((int) IPhiSeed)%20;
      vals[36] = EtaCrySeed;
      vals[37] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = PreShowerOverRaw;
    }
    
  } else if(fVersionType == kNoTrkVarV1) {  ///// Regression without tracker variables (V1) /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 39 : 32;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = IEtaSeed;
      vals[32] = IPhiSeed;
      vals[33] = ((int) IEtaSeed)%5;
      vals[34] = ((int) IPhiSeed)%2;
      vals[35] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[36] = ((int) IPhiSeed)%20;
      vals[37] = EtaCrySeed;
      vals[38] = PhiCrySeed;    
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = PreShowerOverRaw;    
    }
  
  } else if(fVersionType == kWithTrkVar) {  ///// Regression with tracker variables /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 43 : 36;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = pt;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = IEtaSeed;
      vals[36] = IPhiSeed;
      vals[37] = ((int) IEtaSeed)%5;
      vals[38] = ((int) IPhiSeed)%2;
      vals[39] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[40] = ((int) IPhiSeed)%20;
      vals[41] = EtaCrySeed;
      vals[42] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = pt;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = PreShowerOverRaw;    
    }
  
  } else if(fVersionType == kWithTrkVarV1) {  ///// Regression with tracker variables (V1) /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 46 : 39;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = IEtaSeed;
      vals[39] = IPhiSeed;
      vals[40] = ((int) IEtaSeed)%5;
      vals[41] = ((int) IPhiSeed)%2;
      vals[42] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[43] = ((int) IPhiSeed)%20;
      vals[44] = EtaCrySeed;
      vals[45] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = PreShowerOverRaw;
    }
  
  } else if(fVersionType == kWithTrkVarV2) {  ///// Regression with tracker variables (V2) /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 53 : 46;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = detaIn;
      vals[39] = dphiIn;
      vals[40] = detaCalo;
      vals[41] = dphiCalo;
      vals[42] = GsfTrackChiSqr;
      vals[43] = KFTrackNLayers;
      vals[44] = ElectronEnergyOverPout;
      vals[45] = IEtaSeed;
      vals[46] = IPhiSeed;
      vals[47] = ((int) IEtaSeed)%5;
      vals[48] = ((int) IPhiSeed)%2;
      vals[49] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[50] = ((int) IPhiSeed)%20;
      vals[51] = EtaCrySeed;
      vals[52] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = detaIn;
      vals[39] = dphiIn;
      vals[40] = detaCalo;
      vals[41] = dphiCalo;
      vals[42] = GsfTrackChiSqr;
      vals[43] = KFTrackNLayers;
      vals[44] = ElectronEnergyOverPout;
      vals[45] = PreShowerOverRaw;    
    } 
  
  } else if(fVersionType == kWithSubCluVar) {  ///// Regression with subcluster variables and without track variables /////
    isBarrel = ele->isEB();
    nvars = isBarrel ? 61 : 65;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = rho;
      vals[1]  = nvertices;
      vals[2]  = IsEcalDriven;
      vals[3]  = isEtaGap;
      vals[4]  = isPhiGap;
      vals[5]  = isDeeGap;
      vals[6]  = SCRawEnergy;
      vals[7]  = scEta;
      vals[8]  = scPhi;
      vals[9]  = R9;
      vals[10] = etawidth;
      vals[11] = phiwidth;
      vals[12] = NClusters;
      vals[13] = HoE;
      vals[14] = EtaSeed - scEta;
      vals[15] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[16] = ESeed/SCRawEnergy;
      vals[17] = E3x3Seed/ESeed;
      vals[18] = E5x5Seed/SCRawEnergy;
      vals[19] = E5x5Seed/ESeed;
      vals[20] = EMaxSeed/ESeed;
      vals[21] = E2ndSeed/ESeed;
      vals[22] = ETopSeed/ESeed;
      vals[23] = EBottomSeed/ESeed;
      vals[24] = ELeftSeed/ESeed;
      vals[25] = ERightSeed/ESeed;
      vals[26] = E2x5MaxSeed/ESeed;
      vals[27] = E2x5TopSeed/ESeed;
      vals[28] = E2x5BottomSeed/ESeed;
      vals[29] = E2x5LeftSeed/ESeed;
      vals[30] = E2x5RightSeed/ESeed;
      vals[31] = see;
      vals[32] = spp;
      vals[33] = sep;
      vals[34] = phiwidth/etawidth;
      vals[35] = (ELeftSeed+ERightSeed==0. ? 0. : (ELeftSeed-ERightSeed)/(ELeftSeed+ERightSeed));
      vals[36] = (ETopSeed+EBottomSeed==0. ? 0. : (ETopSeed-EBottomSeed)/(ETopSeed+EBottomSeed));
      vals[37] = ESubs/SCRawEnergy;
      vals[38] = ESub1/SCRawEnergy;
      vals[39] = (NClusters<=1 ? 999. : EtaSub1-EtaSeed);
      vals[40] = (NClusters<=1 ? 999. : atan2(sin(PhiSub1-PhiSeed),cos(PhiSub1-PhiSeed)));
      vals[41] = (NClusters<=1 ? 0.   : EMaxSub1/ESub1);
      vals[42] = (NClusters<=1 ? 0.   : E3x3Sub1/ESub1);
      vals[43] = ESub2/SCRawEnergy;
      vals[44] = (NClusters<=2 ? 999. : EtaSub2-EtaSeed);
      vals[45] = (NClusters<=2 ? 999. : atan2(sin(PhiSub2-PhiSeed),cos(PhiSub2-PhiSeed)));
      vals[46] = (NClusters<=2 ? 0.   : EMaxSub2/ESub2);
      vals[47] = (NClusters<=2 ? 0.   : E3x3Sub2/ESub2);
      vals[48] = ESub3/SCRawEnergy;
      vals[49] = (NClusters<=3 ? 999. : EtaSub3-EtaSeed);
      vals[50] = (NClusters<=3 ? 999. : atan2(sin(PhiSub3-PhiSeed),cos(PhiSub3-PhiSeed)));
      vals[51] = (NClusters<=3 ? 0.   : EMaxSub3/ESub3);
      vals[52] = (NClusters<=3 ? 0.   : E3x3Sub3/ESub3);
      vals[53] = IEtaSeed;
      vals[54] = ((int) IEtaSeed)%5;
      vals[55] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[56] = IPhiSeed;
      vals[57] = ((int) IPhiSeed)%2;
      vals[58] = ((int) IPhiSeed)%20;
      vals[59] = EtaCrySeed;
      vals[60] = PhiCrySeed;
    
    } else {
      vals[0]  = rho;
      vals[1]  = nvertices;
      vals[2]  = IsEcalDriven;
      vals[3]  = isEtaGap;
      vals[4]  = isPhiGap;
      vals[5]  = isDeeGap;
      vals[6]  = SCRawEnergy;
      vals[7]  = scEta;
      vals[8]  = scPhi;
      vals[9]  = R9;
      vals[10] = etawidth;
      vals[11] = phiwidth;
      vals[12] = NClusters;
      vals[13] = HoE;
      vals[14] = EtaSeed - scEta;
      vals[15] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[16] = ESeed/SCRawEnergy;
      vals[17] = E3x3Seed/ESeed;
      vals[18] = E5x5Seed/SCRawEnergy;
      vals[19] = E5x5Seed/ESeed;
      vals[20] = EMaxSeed/ESeed;
      vals[21] = E2ndSeed/ESeed;
      vals[22] = ETopSeed/ESeed;
      vals[23] = EBottomSeed/ESeed;
      vals[24] = ELeftSeed/ESeed;
      vals[25] = ERightSeed/ESeed;
      vals[26] = E2x5MaxSeed/ESeed;
      vals[27] = E2x5TopSeed/ESeed;
      vals[28] = E2x5BottomSeed/ESeed;
      vals[29] = E2x5LeftSeed/ESeed;
      vals[30] = E2x5RightSeed/ESeed;
      vals[31] = see;
      vals[32] = spp;
      vals[33] = sep;
      vals[34] = phiwidth/etawidth;
      vals[35] = (ELeftSeed+ERightSeed==0. ? 0. : (ELeftSeed-ERightSeed)/(ELeftSeed+ERightSeed));
      vals[36] = (ETopSeed+EBottomSeed==0. ? 0. : (ETopSeed-EBottomSeed)/(ETopSeed+EBottomSeed));
      vals[37] = ESubs/SCRawEnergy;
      vals[38] = ESub1/SCRawEnergy;
      vals[39] = (NClusters<=1 ? 999. : EtaSub1-EtaSeed);
      vals[40] = (NClusters<=1 ? 999. : atan2(sin(PhiSub1-PhiSeed),cos(PhiSub1-PhiSeed)));
      vals[41] = (NClusters<=1 ? 0.   : EMaxSub1/ESub1);
      vals[42] = (NClusters<=1 ? 0.   : E3x3Sub1/ESub1);
      vals[43] = ESub2/SCRawEnergy;
      vals[44] = (NClusters<=2 ? 999. : EtaSub2-EtaSeed);
      vals[45] = (NClusters<=2 ? 999. : atan2(sin(PhiSub2-PhiSeed),cos(PhiSub2-PhiSeed)));
      vals[46] = (NClusters<=2 ? 0.   : EMaxSub2/ESub2);
      vals[47] = (NClusters<=2 ? 0.   : E3x3Sub2/ESub2);
      vals[48] = ESub3/SCRawEnergy;
      vals[49] = (NClusters<=3 ? 999. : EtaSub3-EtaSeed);
      vals[50] = (NClusters<=3 ? 999. : atan2(sin(PhiSub3-PhiSeed),cos(PhiSub3-PhiSeed)));
      vals[51] = (NClusters<=3 ? 0.   : EMaxSub3/ESub3);
      vals[52] = (NClusters<=3 ? 0.   : E3x3Sub3/ESub3);
      vals[53] = PreShowerOverRaw;
      vals[54] = NPshwClusters;
      vals[55] = EPshwSubs/SCRawEnergy;
      vals[56] = EPshwSub1/SCRawEnergy;
      vals[57] = (NPshwClusters==0 ? 999. : EtaPshwSub1-EtaSeed);
      vals[58] = (NPshwClusters==0 ? 999. : atan2(sin(PhiPshwSub1-PhiSeed),cos(PhiPshwSub1-PhiSeed)));
      vals[59] = EPshwSub2/SCRawEnergy;
      vals[60] = (NPshwClusters<=1 ? 999. : EtaPshwSub2-EtaSeed);
      vals[61] = (NPshwClusters<=1 ? 999. : atan2(sin(PhiPshwSub2-PhiSeed),cos(PhiPshwSub2-PhiSeed)));
      vals[62] = EPshwSub3/SCRawEnergy;
      vals[63] = (NPshwClusters<=2 ? 999. : EtaPshwSub3-EtaSeed);
      vals[64] = (NPshwClusters<=2 ? 999. : atan2(sin(PhiPshwSub3-PhiSeed),cos(PhiPshwSub3-PhiSeed)));    
    }
  
  } else {
    std::cout << "Error: Regression VersionType " << fVersionType << " is not supported!" << std::endl;
    assert(0);    
  }

  // evaluate regression
  double value       = isBarrel ? (SCRawEnergy * forestCorrection_eb->GetResponse(vals))
                                : (SCRawEnergy * (1+PreShowerOverRaw) * forestCorrection_ee->GetResponse(vals));
  double uncertainty = isBarrel ? (SCRawEnergy * forestUncertainty_eb->GetResponse(vals))
                                : (SCRawEnergy * (1+PreShowerOverRaw) * forestUncertainty_ee->GetResponse(vals));
  Int_t BinIndex     = isBarrel ? 0 : 1;

  // print debug
  if(printDebug) {    
    std::cout << "[ElectronEnergyRegression]" << std::endl;
    std::cout << (isBarrel ? "Barrel :" : "Endcap :") << std::endl;
    for(uint v=0; v < nvars; ++v) std::cout << v << "=" << vals[v] << ", ";
    std::cout << std::endl;

    std::cout << "BinIndex : " << BinIndex << "\n";
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy = " << value << std::endl;
    std::cout << "regression energy uncertainty = " << uncertainty << std::endl;
  }
  
  // Cleaning up and returning
  delete[] vals;
  return std::pair<double,double>(value,uncertainty);
}
