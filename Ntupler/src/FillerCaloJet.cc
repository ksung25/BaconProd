#include "BaconProd/Ntupler/interface/FillerCaloJet.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TCaloJet.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace baconhep;


//--------------------------------------------------------------------------------------------------
FillerCaloJet::FillerCaloJet(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fMinPt              (iConfig.getUntrackedParameter<double>("minPt",20)),
  fUseGen             (iConfig.getUntrackedParameter<bool>("doGenJet",true)),
  fRhoName            (iConfig.getUntrackedParameter<std::string>("edmRhoName","fixedGridRhoFastjetAll")),
  fJetName            (iConfig.getUntrackedParameter<std::string>("jetName","ak4CaloJets")),
  fGenJetName         (iConfig.getUntrackedParameter<std::string>("genJetName","AK4GenJetsCHS")),
  fJetFlavorName      (iConfig.getUntrackedParameter<std::string>("jetFlavorName","AK4byValAlgoCHS")),
  fConeSize           (iConfig.getUntrackedParameter<double>("coneSize",0.4)),
  fJetCorr            (0),
  fJetUnc             (0)
{
  std::vector<std::string> empty_vstring;
  initJetCorr(iConfig.getUntrackedParameter< std::vector<std::string> >("jecFiles",empty_vstring),
	      iConfig.getUntrackedParameter< std::vector<std::string> >("jecUncFiles",empty_vstring));

  fTokJetName         = iC.consumes<reco::CaloJetCollection>(fJetName);
  fTokGenJetName      = iC.consumes<reco::GenJetCollection> (fGenJetName);
  fTokRhoTag          = iC.consumes<double>                 (fRhoName);
}

//--------------------------------------------------------------------------------------------------
FillerCaloJet::~FillerCaloJet()
{
  delete fJetCorr;
  delete fJetUnc;
}

//--------------------------------------------------------------------------------------------------
void FillerCaloJet::initJetCorr(const std::vector<std::string> &jecFiles,
				const std::vector<std::string> &jecUncFiles)
{
  assert(jecFiles.size()>0);
  assert(jecUncFiles.size()>0);

  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
 
  std::vector<JetCorrectorParameters> corrParams;
  for(unsigned int icorr=0; icorr<jecFiles.size(); icorr++) {
    corrParams.push_back(JetCorrectorParameters(cmssw_base_src + jecFiles[icorr]));
  }
  fJetCorr = new FactorizedJetCorrector(corrParams);
  
//  JetCorrectorParameters param(cmssw_base_src + jecUncFiles[0]);
//  fJetUnc = new JetCorrectionUncertainty(param);
}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerCaloJet::fill(TClonesArray *array,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, 
		     const std::vector<TriggerRecord> &triggerRecords,
		     const trigger::TriggerEvent &triggerEvent) 
{
  assert(array);
  
  // Get jet collection
  edm::Handle<reco::CaloJetCollection> hCaloJetProduct;
  iEvent.getByToken(fTokJetName,hCaloJetProduct);
  assert(hCaloJetProduct.isValid());
  const reco::CaloJetCollection *jetCol = hCaloJetProduct.product();
  
  // Get gen jet collection
  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJetCol = 0;
  if(fUseGen) { 
    iEvent.getByToken(fTokGenJetName,hGenJetProduct);
    assert(hGenJetProduct.isValid());
    genJetCol = hGenJetProduct.product();
  }

  // Get Jet Flavor Match
  //edm::Handle<reco::JetFlavourInfoMatchingCollection> hJetFlavourMatch;
  //if(fUseGen) {
  //  iEvent.getByToken(fTokJetFlavorName, hJetFlavourMatch);
  //  assert(hJetFlavourMatch.isValid());
  //}
  
  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid()); 
 
  TClonesArray &rArray      = *array;
  for(reco::CaloJetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {
    const double ptRaw = itJet->pt();

    // input to jet corrections
    double jetcorr = 1;
    if(fabs(itJet->eta()) < 5.191) {
      fJetCorr->setJetPt(ptRaw);
      fJetCorr->setJetEta(itJet->eta());
      fJetCorr->setJetPhi(itJet->phi());
      fJetCorr->setJetE(itJet->energy());
      fJetCorr->setRho(*hRho);
      fJetCorr->setJetA(itJet->jetArea());
      fJetCorr->setJetEMF(-99.0);
      jetcorr = fJetCorr->getCorrection();
    }

    // jet pT cut (BOTH raw AND corrected pT must exceed threshold)
    if(ptRaw*jetcorr < fMinPt || ptRaw < fMinPt) continue;
    //    fJetUnc->setJetPt ( ptRaw  );
    //    fJetUnc->setJetEta( itJet->eta() );
    //    double jetunc = fJetUnc->getUncertainty(true);
    
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TCaloJet();
    baconhep::TCaloJet *pJet = (baconhep::TCaloJet*)rArray[index];
 
    //
    // Kinematics
    //==============================    
    pJet->pt    = ptRaw * jetcorr;
    pJet->eta   = itJet->eta();
    pJet->phi   = itJet->phi();
    pJet->mass  = itJet->mass() * jetcorr;
    pJet->ptRaw = ptRaw;
    pJet->area  = itJet->jetArea();
    //    pJet->unc   = jetunc;
    

    // Identification
    //==============================
    pJet->nParticles = itJet->nConstituents();
    // Basic Noise Variables
    pJet->neuEmEBFrac = itJet->emEnergyInEB() / itJet->energy();
    pJet->neuEmEEFrac = itJet->emEnergyInEE() / itJet->energy();
    pJet->neuEmHFFrac = itJet->emEnergyInHF() / itJet->energy();
    pJet->neuHadHBFrac = itJet->hadEnergyInHB() / itJet->energy();
    pJet->neuHadHEFrac = itJet->hadEnergyInHE() / itJet->energy();
    pJet->neuHadHFFrac = itJet->hadEnergyInHF() / itJet->energy();

    //
    // Generator matching
    //==============================
    const reco::GenJet * matchGenJet = 0; 
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    if(matchGenJet != 0) { 
      //reco::CaloJetRef jetRef(hCaloJetProduct, itJet - jetCol->begin());
      //reco::JetBaseRef jetBaseRef(jetRef);
      //pJet->partonFlavor = (*hJetFlavourMatch)[jetBaseRef].getPartonFlavour();
      //pJet->hadronFlavor = (*hJetFlavourMatch)[jetBaseRef].getHadronFlavour();
      pJet->genpt        = matchGenJet->pt();
      pJet->geneta       = matchGenJet->eta();
      pJet->genphi       = matchGenJet->phi();
      pJet->genm         = matchGenJet->mass();
    }
    
    pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, triggerEvent);
  } 
}
const reco::BasicJet* FillerCaloJet::match( const reco::CaloJet *iJet,const reco::BasicJetCollection *jets ) { 
  int lId = -1;
  double dRmin = 999.;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::BasicJet* jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if(dR > fConeSize) continue;
    if ( dR < dRmin ) {
      dRmin = dR;
      lId = i;
    }
  }
  const reco::BasicJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}
const reco::GenJet* FillerCaloJet::match( const reco::CaloJet *iJet,const reco::GenJetCollection *jets ) { 
  int lId = -1;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::GenJet *jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if ( dR < 0.25 ) {
      lId = i;
      break;
    }
  }
  const reco::GenJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}
