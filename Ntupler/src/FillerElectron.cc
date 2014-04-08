#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <utility>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerElectron::FillerElectron(const edm::ParameterSet &iConfig):
  fMinPt       (iConfig.getUntrackedParameter<double>("minPt",7)),
  fEleName     (iConfig.getUntrackedParameter<std::string>("edmName","gsfElectrons")),
  fPFCandName  (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName   (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fConvName    (iConfig.getUntrackedParameter<std::string>("edmConversionName","allConversions")),
  fRhoName     (iConfig.getUntrackedParameter<std::string>("edmRhoForEnergyRegression","kt6PFJets")),
  fEBSCName    (iConfig.getUntrackedParameter<std::string>("edmEBSuperClusterName","correctedHybridSuperClusters")),
  fEESCName    (iConfig.getUntrackedParameter<std::string>("edmEESuperClusterName","correctedMulti5x5SuperClustersWithPreshower")),
  fEBRecHitName(iConfig.getUntrackedParameter<std::string>("edmEBRecHitName","reducedEcalRecHitsEB")),
  fEERecHitName(iConfig.getUntrackedParameter<std::string>("edmEERecHitName","reducedEcalRecHitsEE")),
  fDoEleCorr   (iConfig.getUntrackedParameter<bool>("doHZZ4lCorr",true))
{
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  
  std::vector<std::string> empty_vstring;  
  std::vector<std::string> eleIDFilenames;
  eleIDFilenames = iConfig.getUntrackedParameter< std::vector<std::string> >("eleIDFilenames",empty_vstring);
  for(unsigned int ifile=0; ifile<eleIDFilenames.size(); ifile++) {
    eleIDFilenames[ifile] = cmssw_base_src + eleIDFilenames[ifile];
  }  
  
  fEleIDMVA.initialize("BDT", EGammaMvaEleEstimator::kNonTrig, true, eleIDFilenames);
  
  if(fDoEleCorr) {
    bool isData = iConfig.getUntrackedParameter<bool>("isData",false);
    bool doRand = iConfig.getUntrackedParameter<bool>("doRandForMC",false);
    fEleCorr.initialize((cmssw_base_src + iConfig.getUntrackedParameter<std::string>("eleEnergyRegWeights","")).c_str(),
  		        baconhep::ElectronEnergyRegression::kWithSubCluVar,
        	        isData ? baconhep::ElectronEnergySmearingScaling::k22Jan2013ReReco : baconhep::ElectronEnergySmearingScaling::kSummer12_LegacyPaper,
        	        2,
        	        (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("scalesCorr","")).c_str(),
		        (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("smearsCorrType1","")).c_str(),
        	        (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("smearsCorrType2","")).c_str(),
		        (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("smearsCorrType3","")).c_str(),
		        (cmssw_base_src + iConfig.getUntrackedParameter<std::string>("linearityCorr","")).c_str(),
        	        !isData && doRand);

  }
}

//--------------------------------------------------------------------------------------------------
FillerElectron::~FillerElectron(){}

//--------------------------------------------------------------------------------------------------
void FillerElectron::fill(TClonesArray *array,	    
	                  const edm::Event &iEvent, const edm::EventSetup &iSetup,      
	                  const reco::Vertex &pv, const int nvtx,
			  const std::vector<const reco::PFCandidate*> &pfNoPU,
			  const std::vector<TriggerRecord> &triggerRecords,
			  const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
  assert(!fDoEleCorr || fEleCorr.isInitialized());
  assert(fEleIDMVA.isInitialized());
  
  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> hEleProduct;
  iEvent.getByLabel(fEleName,hEleProduct);
  assert(hEleProduct.isValid());
  const reco::GsfElectronCollection *eleCol = hEleProduct.product();

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByLabel(fPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
  
  // Get track collection
  edm::Handle<reco::TrackCollection> hTrackProduct;
  iEvent.getByLabel(fTrackName,hTrackProduct);
  assert(hTrackProduct.isValid());
  const reco::TrackCollection *trackCol = hTrackProduct.product();
  
  // Get conversions collection
  edm::Handle<reco::ConversionCollection> hConvProduct;
  iEvent.getByLabel(fConvName,hConvProduct);
  assert(hConvProduct.isValid());

  // event energy density
  edm::Handle<double> hRho;
  edm::InputTag rhoTag(fRhoName,"rho","RECO");
  iEvent.getByLabel(rhoTag,hRho);
  
  // Get SuperCluster collections
  edm::Handle<reco::SuperClusterCollection> hEBSCProduct;
  iEvent.getByLabel(fEBSCName,hEBSCProduct);
  assert(hEBSCProduct.isValid());
  const reco::SuperClusterCollection *ebSCCol = hEBSCProduct.product();
  
  edm::Handle<reco::SuperClusterCollection> hEESCProduct;
  iEvent.getByLabel(fEESCName,hEESCProduct);
  assert(hEESCProduct.isValid());
  const reco::SuperClusterCollection *eeSCCol = hEESCProduct.product();
  
    
  const double ELE_MASS = 0.000511;

  
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();
  
  edm::InputTag ebRecHitTag(fEBRecHitName);
  edm::InputTag eeRecHitTag(fEERecHitName);
  EcalClusterLazyTools lazyTools(iEvent, iSetup, ebRecHitTag, eeRecHitTag);
  
  for(reco::GsfElectronCollection::const_iterator itEle = eleCol->begin(); itEle!=eleCol->end(); ++itEle) {
    
    const reco::GsfTrackRef gsfTrack = itEle->gsfTrack();
    const reco::SuperClusterRef sc   = itEle->superCluster();
    
    // electron pT cut
    std::pair<double,double> result = fDoEleCorr ? fEleCorr.evaluate(&(*itEle), *hRho, nvtx, iEvent.id().run(), iSetup, lazyTools, false)
                                                 : std::pair<double,double>(itEle->p(), 0);
    TLorentzVector elevec;
    elevec.SetPtEtaPhiM(itEle->pt(), itEle->eta(), itEle->phi(), ELE_MASS);
    double ptCorr = result.first*TMath::Sin(elevec.Theta());
    if(itEle->pt() < fMinPt && ptCorr < fMinPt) continue;
    
    // construct object and place in array    
    TClonesArray &rElectronArr = *array;
    assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
    const int index = rElectronArr.GetEntries();  
    new(rElectronArr[index]) baconhep::TElectron();
    baconhep::TElectron *pElectron = (baconhep::TElectron*)rElectronArr[index];
    
    //
    // Kinematics
    //==============================    
    pElectron->pt         = itEle->pt();
    pElectron->eta        = itEle->eta();
    pElectron->phi        = itEle->phi();
    pElectron->ptHZZ4l    = ptCorr;
    pElectron->ptErrHZZ4l = result.second*TMath::Sin(elevec.Theta());;
    pElectron->q          = itEle->charge();
    pElectron->ecalEnergy = itEle->correctedEcalEnergy();
    pElectron->scEt       = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pElectron->scEta      = sc->eta();
    pElectron->scPhi      = sc->phi();
    pElectron->scEtHZZ4l  = (sc->energy())*(sc->position().Rho())/(sc->position().R())*(pElectron->ptHZZ4l)/(pElectron->pt);
    pElectron->r9         = lazyTools.e3x3(*(sc->seed())) / sc->rawEnergy();
    
    pElectron->pfPt  = 0;
    pElectron->pfEta = 0;
    pElectron->pfPhi = 0;
    for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
      if( (itEle->track().isNonnull() && itPF->trackRef().isNonnull() && itEle->track() == itPF->trackRef()) ||
          (itPF->gsfTrackRef().isNonnull() && itPF->gsfTrackRef() == gsfTrack) ) {
        
	pElectron->pfPt  = itPF->pt();
	pElectron->pfEta = itPF->eta();
	pElectron->pfPhi = itPF->phi();
      }
    }
    
    //
    // Isolation
    //==============================
    pElectron->trkIso03  = itEle->dr03TkSumPt();
    pElectron->ecalIso03 = itEle->dr03EcalRecHitSumEt();
    pElectron->hcalIso03 = itEle->dr03HcalTowerSumEt();
    
    computeIso(*itEle, 0.3, pfNoPU,
               pElectron->chHadIso03,
               pElectron->gammaIso03,
               pElectron->neuHadIso03);   
    
    computeIso(*itEle, 0.4, pfNoPU,
               pElectron->chHadIso04,
               pElectron->gammaIso04,
               pElectron->neuHadIso04);   
    
    //
    // Impact Parameter
    //==============================
    if(gsfTrack.isNonnull()) { 
      pElectron->d0 = -gsfTrack->dxy(pv.position());
      pElectron->dz =  gsfTrack->dz(pv.position());
      const reco::TransientTrack &tt = transientTrackBuilder->build(gsfTrack);
      const double gsfsign = (pElectron->d0 >= 0) ? 1. : -1.;
      const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
      pElectron->sip3d = ip3d.first ? gsfsign*ip3d.second.value() / ip3d.second.error() : -999.;
    }    
    //
    // Identification
    //==============================
    pElectron->sieie  = itEle->sigmaIetaIeta();
    pElectron->hovere = itEle->hcalOverEcal();
    pElectron->eoverp = itEle->eSuperClusterOverP();
    pElectron->fbrem  = itEle->fbrem();
    pElectron->dEtaIn = itEle->deltaEtaSuperClusterTrackAtVtx();
    pElectron->dPhiIn = itEle->deltaPhiSuperClusterTrackAtVtx();
    
    pElectron->mva = evalEleIDMVA(*itEle, lazyTools);
    
    pElectron->classification = itEle->classification();
    
    pElectron->isConv = ConversionTools::hasMatchedConversion(*itEle, hConvProduct, pv.position(), true, 2.0, 1e-6, 0);
    
    if(gsfTrack.isNonnull()) pElectron->nMissingHits = gsfTrack->trackerExpectedHitsInner().numberOfHits();
    
    pElectron->typeBits=0;
    if(itEle->ecalDrivenSeed())    pElectron->typeBits |= baconhep::kEcalDriven;
    if(itEle->trackerDrivenSeed()) pElectron->typeBits |= baconhep::kTrackerDriven;
    
    pElectron->fiducialBits=0;
    if(itEle->isEB())        pElectron->fiducialBits |= kIsEB;
    if(itEle->isEE())        pElectron->fiducialBits |= kIsEE;
    if(itEle->isGap())       pElectron->fiducialBits |= kIsGap;
    if(itEle->isEBEEGap())   pElectron->fiducialBits |= kIsEBEEGap;
    if(itEle->isEBGap())     pElectron->fiducialBits |= kIsEBGap;
    if(itEle->isEBEtaGap())  pElectron->fiducialBits |= kIsEBEtaGap;
    if(itEle->isEBPhiGap())  pElectron->fiducialBits |= kIsEBPhiGap;
    if(itEle->isEEGap())     pElectron->fiducialBits |= kIsEEGap;
    if(itEle->isEEDeeGap())  pElectron->fiducialBits |= kIsEEDeeGap;
    if(itEle->isEERingGap()) pElectron->fiducialBits |= kIsEERingGap;
    
    pElectron->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator iS = ebSCCol->begin(); iS!=ebSCCol->end(); ++iS) {
      scIndex++;
      if(itEle->superCluster().get() == &(*iS)) {
        pElectron->scID = scIndex;
	break;
      }
    }
    for(reco::SuperClusterCollection::const_iterator iS = eeSCCol->begin(); iS!=eeSCCol->end(); ++iS) {
      scIndex++;
      if(itEle->superCluster().get() == &(*iS)) {
        pElectron->scID = scIndex;
	break;
      }
    }
    
    pElectron->trkID = -1;
    if(itEle->closestTrack().isNonnull()) {
      int trkIndex = -1;    
      for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
        trkIndex++;
        if(itEle->closestTrack().get() == &(*itTrk)) {
          pElectron->trkID = trkIndex;
	  break;
        }
      }
    }
    
    pElectron->hltMatchBits = TriggerTools::matchHLT(pElectron->eta, pElectron->phi, triggerRecords, triggerEvent);
  }
}

//--------------------------------------------------------------------------------------------------
void FillerElectron::computeIso(const reco::GsfElectron &ele, const double extRadius,
                                const std::vector<const reco::PFCandidate*> &pfNoPU,
                                float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double intRadiusChHad  = 0.015;
  const double intRadiusGamma  = 0.08;
  const double intRadiusNeuHad = 0.;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<pfNoPU.size(); ipf++) {
    const reco::PFCandidate *pfcand = pfNoPU.at(ipf);
    
    if(pfcand->particleId() == reco::PFCandidate::gamma &&
       pfcand->mva_nothing_gamma()>0.99 && 
       ele.gsfTrack().isNonnull() && ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits()>0 &&
       pfcand->superClusterRef().isNonnull() && ele.superCluster().isNonnull() && pfcand->superClusterRef() == ele.superCluster())
     continue;
    
    if     (pfcand->particleId() == reco::PFCandidate::h)     { intRadius = ele.isEB() ? 0 : intRadiusChHad; }
    else if(pfcand->particleId() == reco::PFCandidate::gamma) { intRadius = ele.isEB() ? 0 : intRadiusGamma;  }
    else if(pfcand->particleId() == reco::PFCandidate::h0)    { intRadius = ele.isEB() ? 0 : intRadiusNeuHad; }
    
    double dr = reco::deltaR(pfcand->eta(), pfcand->phi(), ele.eta(), ele.phi());
    if(dr>=extRadius || dr<intRadius) continue;

    if     (pfcand->particleId() == reco::PFCandidate::h)     { chHadIso  += pfcand->pt(); }
    else if(pfcand->particleId() == reco::PFCandidate::gamma) { gammaIso  += pfcand->pt(); }
    else if(pfcand->particleId() == reco::PFCandidate::h0)    { neuHadIso += pfcand->pt(); }
  }
  
  out_chHadIso  = chHadIso;
  out_gammaIso  = gammaIso;
  out_neuHadIso = neuHadIso;
}

//--------------------------------------------------------------------------------------------------
double FillerElectron::evalEleIDMVA(const reco::GsfElectron &ele, EcalClusterLazyTools& lazyTools)
{
  assert(fEleIDMVA.isInitialized());
  
  const reco::TrackRef kfTrackRef   = ele.closestCtfTrackRef();
  const reco::SuperClusterRef scRef = ele.superCluster();
  std::vector<float> vCov = lazyTools.localCovariances(*(scRef->seed()));
  
  double _fbrem          = (ele.fbrem() < -1.) ? -1. : ele.fbrem();
  double _kftrk_chisq    = kfTrackRef.isNonnull() ? TMath::Min(kfTrackRef->normalizedChi2(),10.) : 0;
  double _kftrk_nhits    = kfTrackRef.isNonnull() ? kfTrackRef->hitPattern().trackerLayersWithMeasurement() : -1; 
  double _gsftrk_chisq   = 1; 
  if(ele.gsfTrack().isNonnull()) _gsftrk_chisq  = TMath::Min(ele.gsfTrack()->normalizedChi2(),200.); 
  double _deta           = TMath::Min(fabs(ele.deltaEtaSuperClusterTrackAtVtx()),0.06);
  double _dphi           = TMath::Min(fabs(ele.deltaPhiSuperClusterTrackAtVtx()),0.6);
  double _detacalo       = TMath::Min(fabs(ele.deltaEtaSeedClusterTrackAtCalo()),0.2);
  double _sigieie        = ele.sigmaIetaIeta();
  double _sigiphiiphi    = isnan(vCov[2]) ? 0 : sqrt(vCov[2]);
  double _etawidth       = scRef->etaWidth(); 
  double _phiwidth       = scRef->phiWidth();
  double _e1x5e5x5       = (ele.e5x5() != 0) ? 1 - ele.e1x5()/ele.e5x5() : -1;
  double _r9             = TMath::Min(lazyTools.e3x3(*(scRef->seed()))/scRef->rawEnergy(), 5.);
  double _h_o_e          = ele.hcalOverEcal();
  double _e_o_p          = (ele.eSuperClusterOverP()  > 20.) ? 20. : ele.eSuperClusterOverP();
  double _eeleclu_o_pout = (ele.eEleClusterOverPout() > 20.) ? 20. : ele.eEleClusterOverPout();
  double _IoEmIoP        = (double)(1./ele.correctedEcalEnergy()) - (double)(1./ele.p());
  double _epreoraw       = (double)scRef->preshowerEnergy() / scRef->rawEnergy();

  return fEleIDMVA.mvaValue(_fbrem,
                            _kftrk_chisq,
                            _kftrk_nhits,
                            _gsftrk_chisq,
                            _deta,
                            _dphi,
                            _detacalo,
                            _sigieie,
                            _sigiphiiphi,
                            _etawidth,
                            _phiwidth,
                            _e1x5e5x5,
                            _r9,
                            _h_o_e,
                            _e_o_p,
                            _IoEmIoP,
                            _eeleclu_o_pout,
                            _epreoraw,
                            scRef->eta(), 
                            ele.pt(),
                            false);
}
