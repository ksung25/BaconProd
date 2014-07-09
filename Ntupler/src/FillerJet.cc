#include "BaconProd/Ntupler/interface/FillerJet.hh"
#include "BaconProd/Utils/interface/EnergyCorrelator.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/JetTools.hh"
#include "BaconProd/Utils/interface/CMSTopTagger.hh"
#include "BaconProd/Utils/interface/HEPTopTaggerWrapper.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace baconhep;


//--------------------------------------------------------------------------------------------------
FillerJet::FillerJet(const edm::ParameterSet &iConfig, const double coneSize, const std::string prefix,std::string postfix):
  fMinPt              (iConfig.getUntrackedParameter<double>("minPt",20)),
  fUseGen             (iConfig.getUntrackedParameter<bool>("doGenJet",true)),
  fPVName             (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fRhoName            (iConfig.getUntrackedParameter<std::string>("edmRhoName","kt6PFJets")),
  fJetName            (prefix + iConfig.getUntrackedParameter<std::string>("jetName","ak5PFJets")+postfix),
  fGenJetName         (prefix + iConfig.getUntrackedParameter<std::string>("genJetName","AKT5GenJets")),
  fJetFlavorName      (prefix + iConfig.getUntrackedParameter<std::string>("jetFlavorName","AKT5byValAlgo")+postfix),
  fJetFlavorPhysName  (prefix + iConfig.getUntrackedParameter<std::string>("jetFlavorPhysName","AKT5byValPhys")+postfix),
  fPruneJetName       (prefix + iConfig.getUntrackedParameter<std::string>("pruneJetName","ca8PrunedJets")+postfix),
  fSubJetName         (prefix + iConfig.getUntrackedParameter<std::string>("subJetName","ca8PrunedJets__SubJets")+postfix),
  fCSVbtagName        (prefix + iConfig.getUntrackedParameter<std::string>("csvBTagName","myCombinedSecondaryVertexBJetTags")+postfix),
  fCSVbtagSubJetName  (prefix + iConfig.getUntrackedParameter<std::string>("csvBTagSubJetName","myCombinedSecondaryVertexBJetTagsSubJets")+postfix),
  fJettinessName      (prefix + iConfig.getUntrackedParameter<std::string>("jettiness","Njettiness")+postfix),
  fQGLikelihood       (prefix + iConfig.getUntrackedParameter<std::string>("qgLikelihood","QGLikelihood")+postfix),
  fQGLikelihoodSubJets(prefix + iConfig.getUntrackedParameter<std::string>("qgLikelihoodSubjet","QGLikelihood")+postfix),
  fConeSize           (coneSize),
  fComputeFullJetInfo (iConfig.getUntrackedParameter<bool>("doComputeFullJetInfo",true)),  
  fJetDef             (0),
  fGenJetDef          (0),
  fCAJetDef           (0),
  fActiveArea         (0),
  fAreaDefinition     (0),
  fClustering         (0),
  fPruner1            (0),
  fPruner2            (0),
  fFilter1            (0),
  fFilter2            (0),
  fSoftDrop1          (0),
  fSoftDrop2          (0),
  fSoftDrop3          (0),
  fTrimmer1           (0),
  fTrimmer2           (0),
  fTrimmer3           (0),
  fTrimmer4           (0), 
  fCMSTopTagger       (0),
  fHEPTopTagger       (0),
  fJetCorr            (0),
  fJetCorrForID       (0),
  fJetUnc             (0)
{
  fCAJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);

  int activeAreaRepeats = 1;
  double ghostArea      = 0.01;
  double ghostEtaMax    = 7.0;
  fActiveArea           = new fastjet::ActiveAreaSpec (ghostEtaMax,activeAreaRepeats,ghostArea);  
  fAreaDefinition       = new fastjet::AreaDefinition (fastjet::active_area_explicit_ghosts, *fActiveArea );
  //Only 2 subjets?
  fPruner1 = new fastjet::Pruner( fastjet::cambridge_algorithm, 0.1, 0.5); //CMS Default
  fPruner2 = new fastjet::Pruner( fastjet::cambridge_algorithm, 0.1, 0.2); //CMS Default
   
  fFilter1 = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorNHardest(3)));
  fFilter2 = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3)));


  fSoftDrop1 = new fastjet::contrib::SoftDropTagger( 0. ,0.1,1.0);
  fSoftDrop2 = new fastjet::contrib::SoftDropTagger( 2. ,0.1,1.0);
  fSoftDrop3 = new fastjet::contrib::SoftDropTagger(-1. ,0.1,1.0);

  fTrimmer1  = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2),  fastjet::SelectorPtFractionMin(0.05)));
  fTrimmer2  = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2),  fastjet::SelectorPtFractionMin(0.03)));
  fTrimmer3  = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.1),  fastjet::SelectorPtFractionMin(0.03)));
  fTrimmer4  = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.05), fastjet::SelectorPtFractionMin(0.03)));
  
  fCMSTopTagger = new fastjet::CMSTopTagger();//0.05,0.8,0.19);
  fHEPTopTagger = new fastjet::HEPTopTagger(0.8,30.,false);

  
  std::vector<std::string> empty_vstring;
  initJetCorr(iConfig.getUntrackedParameter< std::vector<std::string> >("jecFiles",empty_vstring),
              iConfig.getUntrackedParameter< std::vector<std::string> >("jecUncFiles",empty_vstring),
	      iConfig.getUntrackedParameter< std::vector<std::string> >("jecFilesForID",empty_vstring),
	      (postfix.size() > 0) );

  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";

  std::vector<std::string> puIDFiles = iConfig.getUntrackedParameter< std::vector<std::string> >("jetPUIDFiles",empty_vstring);
  assert(puIDFiles.size()==2);
  std::string lowPtWeightFile  = (puIDFiles[0].length()>0) ? (cmssw_base_src + puIDFiles[0]) : "";
  std::string highPtWeightFile = (puIDFiles[1].length()>0) ? (cmssw_base_src + puIDFiles[1]) : "";
  if(postfix.size() > 0)  { 
    TString lTmp = highPtWeightFile;
    lTmp.ReplaceAll("53X","53X_chs");
    highPtWeightFile = lTmp.Data();
  }
  
  fJetPUIDMVACalc.initialize(baconhep::JetPUIDMVACalculator::k53,
                             "BDT",lowPtWeightFile,
			     "BDT",highPtWeightFile);

  fRand = new TRandom2();
}

//--------------------------------------------------------------------------------------------------
FillerJet::~FillerJet()
{
  delete fActiveArea;
  delete fAreaDefinition;
  delete fPruner1;
  delete fPruner2;
  delete fSoftDrop1;
  delete fSoftDrop2; 
  delete fSoftDrop3; 
  delete fTrimmer1;
  delete fTrimmer2;
  delete fTrimmer3;
  delete fTrimmer4;  
  delete fCMSTopTagger;
  delete fHEPTopTagger;
  
  delete fJetCorr;
  delete fJetCorrForID;
  delete fJetUnc;
}

//--------------------------------------------------------------------------------------------------
void FillerJet::initJetCorr(const std::vector<std::string> &jecFiles,
                            const std::vector<std::string> &jecUncFiles,
			    const std::vector<std::string> &jecFilesForID,
			    bool iCHS)
{
  assert(jecFiles.size()>0);
  assert(jecUncFiles.size()>0);
  assert(jecFilesForID.size()>0);
  
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  
  std::vector<JetCorrectorParameters> corrParams;
  for(unsigned int icorr=0; icorr<jecFiles.size(); icorr++) {
    TString lTmp = jecFiles[icorr];
    lTmp.ReplaceAll("PF","PFchs");
    corrParams.push_back(JetCorrectorParameters( (cmssw_base_src + std::string(lTmp.Data())).c_str() ));
  }
  fJetCorr = new FactorizedJetCorrector(corrParams);
  
  JetCorrectorParameters param(cmssw_base_src + jecUncFiles[0]);
  fJetUnc = new JetCorrectionUncertainty(param);
  
  std::vector<JetCorrectorParameters> corrParamsForID;
  for(unsigned int icorr=0; icorr<jecFilesForID.size(); icorr++) {
    corrParamsForID.push_back(JetCorrectorParameters( (cmssw_base_src + jecFilesForID[icorr]).c_str() ));
  }
  fJetCorrForID = new FactorizedJetCorrector(corrParamsForID);
}

//--------------------------------------------------------------------------------------------------
double FillerJet::correction(fastjet::PseudoJet &iJet,double iRho) { 
  fJetCorr->setJetEta(iJet.eta());
  fJetCorr->setJetPt (iJet.pt());
  fJetCorr->setJetPhi(iJet.phi());
  fJetCorr->setJetE  (iJet.e());
  fJetCorr->setRho   (iRho);
  fJetCorr->setJetA  (iJet.area());
  fJetCorr->setJetEMF(-99.0);     
  return fJetCorr->getCorrection();
}

//--------------------------------------------------------------------------------------------------
void FillerJet::fill(TClonesArray *array, TClonesArray *iExtraArray,TClonesArray *iTopArray,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, 
		     const reco::Vertex	&pv,
		     const std::vector<TriggerRecord> &triggerRecords,
		     const trigger::TriggerEvent &triggerEvent) 
{
  assert(array);
  assert(!fComputeFullJetInfo || iExtraArray);
  assert(fJetPUIDMVACalc.isInitialized());
  fRand->SetSeed(iEvent.id().event());
  
  // Get jet collection
  edm::Handle<reco::PFJetCollection> hJetProduct;
  iEvent.getByLabel(fJetName,hJetProduct);
  assert(hJetProduct.isValid());
  const reco::PFJetCollection *jetCol = hJetProduct.product();
  
  // Get gen jet collection
  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJetCol = 0;
  if(fUseGen) { 
    iEvent.getByLabel(fGenJetName,hGenJetProduct);
    assert(hGenJetProduct.isValid());
    genJetCol = hGenJetProduct.product();
  }

  // Get Jet Flavor Match
  edm::Handle<reco::JetFlavourMatchingCollection> jetFlavourMatch;
  if(fUseGen) iEvent.getByLabel(fJetFlavorName, jetFlavourMatch);
  edm::Handle<reco::JetFlavourMatchingCollection> jetFlavourMatchPhys;
  if(fUseGen) iEvent.getByLabel(fJetFlavorPhysName, jetFlavourMatchPhys);

  // Get pruned jet collection
  edm::Handle<reco::BasicJetCollection> hPruneJetProduct;
  iEvent.getByLabel(fPruneJetName,hPruneJetProduct);
  assert(hPruneJetProduct.isValid());
  const reco::BasicJetCollection *pruneJetCol = hPruneJetProduct.product();

  // Get pruned sub jet collection
  edm::Handle<reco::PFJetCollection> hSubJetProduct;
  edm::InputTag subJetTag(fSubJetName,"SubJets");
  iEvent.getByLabel(subJetTag,hSubJetProduct);
  assert(hSubJetProduct.isValid());
  //const reco::PFJetCollection *subJetCol = hSubJetProduct.product();
  
  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  edm::InputTag rhoTag(fRhoName,"rho","RECO");
  iEvent.getByLabel(rhoTag,hRho);
  assert(hRho.isValid()); 
 
  // Get b-jet tagger
  edm::Handle<reco::JetTagCollection> hCSVbtags;
  iEvent.getByLabel(fCSVbtagName, hCSVbtags);
  assert(hCSVbtags.isValid());

  // Get b sub-jets 
  edm::Handle<reco::JetTagCollection> hCSVbtagsSubJets;
  iEvent.getByLabel(fCSVbtagSubJetName, hCSVbtagsSubJets);
  assert(hCSVbtagsSubJets.isValid());
  reco::JetTagCollection hCSVbtagSubJets = *(hCSVbtagsSubJets.product());
  
  // Get N-subjettiness moments
  edm::Handle<edm::ValueMap<float> > hTau1;
  iEvent.getByLabel(fJettinessName,"tau1",hTau1);                                                                                                                                                      
  assert(hTau1.isValid());
  edm::Handle<edm::ValueMap<float> > hTau2;
  iEvent.getByLabel(fJettinessName,"tau2",hTau2);
  assert(hTau2.isValid());
  edm::Handle<edm::ValueMap<float> > hTau3;
  iEvent.getByLabel(fJettinessName,"tau3",hTau3); 
  assert(hTau3.isValid());
  edm::Handle<edm::ValueMap<float> > hTau4;
  iEvent.getByLabel(fJettinessName,"tau3",hTau4); 
  assert(hTau4.isValid());

  //Get Quark Gluon Likelihood
  edm::Handle<edm::ValueMap<float> > hQGLikelihood; 
  iEvent.getByLabel(fQGLikelihood,"qgLikelihood",hQGLikelihood); 
  assert(hQGLikelihood.isValid());

  //Get Quark Gluon Likelihood on subjets
  edm::Handle<edm::ValueMap<float> > hQGLikelihoodSubJets;
  iEvent.getByLabel(fQGLikelihoodSubJets,"qgLikelihood",hQGLikelihoodSubJets);
  assert(hQGLikelihoodSubJets.isValid());
  
  int pId = 0; 
  TClonesArray &rArray      = *array;
  TClonesArray &rExtraArray = *iExtraArray;
  TClonesArray &rTopArray   = *iTopArray;
  for(reco::PFJetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {
    const double ptRaw = itJet->pt();
    pId++;
    // input to jet corrections
    fJetCorr->setJetPt(ptRaw);
    fJetCorr->setJetEta(itJet->eta());
    fJetCorr->setJetPhi(itJet->phi());
    fJetCorr->setJetE(itJet->energy());
    fJetCorr->setRho(*hRho);
    fJetCorr->setJetA(itJet->jetArea());
    fJetCorr->setJetEMF(-99.0);
    double jetcorr = fJetCorr->getCorrection();

    // jet pT cut (both raw and corrected pT must exceed threshold)
    if(ptRaw*jetcorr < fMinPt || ptRaw < fMinPt) continue;
    
    // jet eta cut to prevent exception in computing JEC uncertainty
    if(fabs(itJet->eta()) >= 5.4) continue;

    fJetUnc->setJetPt ( ptRaw  );
    fJetUnc->setJetEta( itJet->eta() );
    double jetunc = fJetUnc->getUncertainty(true);

    bool passLoose = JetTools::passPFLooseID(*itJet);
    
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TJet();
    baconhep::TJet    *pJet = (baconhep::TJet*)rArray[index];
 
    //
    // Kinematics
    //==============================    
    pJet->pt    = ptRaw * jetcorr;
    pJet->eta   = itJet->eta();
    pJet->phi   = itJet->phi();
    pJet->mass  = itJet->mass() * jetcorr;
    pJet->ptRaw = ptRaw;
    pJet->area  = itJet->jetArea();
    pJet->unc   = jetunc;
    
    //
    // Impact Parameter
    //==============================
    pJet->d0 = JetTools::jetD0(*itJet, pv);
    pJet->dz = JetTools::jetDz(*itJet, pv);
    //
    // Identification
    //==============================
    reco::PFJetRef jetRef(hJetProduct, itJet - jetCol->begin());
    reco::JetBaseRef jetBaseRef(jetRef);
    pJet->csv        = (*(hCSVbtags.product()))     [jetBaseRef];
    pJet->qgid       = (*(hQGLikelihood.product())) [jetBaseRef];
    pJet->tau1       = (*(hTau1.product())) [jetBaseRef];
    pJet->tau2       = (*(hTau2.product())) [jetBaseRef];
    pJet->tau3       = (*(hTau3.product())) [jetBaseRef];
    pJet->tau4       = (*(hTau4.product())) [jetBaseRef];
    const reco::BasicJet* matchJet = match(&(*itJet),pruneJetCol);
    double *lQG  = JetTools::subJetBTag(*itJet,hCSVbtagSubJets                                   ,fConeSize );
    double *lCSV = JetTools::subJetQG  (*itJet,hSubJetProduct,(*(hQGLikelihoodSubJets.product())),fConeSize);
    pJet->qg1        = lQG[0];
    pJet->qg2        = lQG[1];
    pJet->csv1       = lCSV[0];
    pJet->csv2       = lCSV[1];
    if(matchJet) pJet->prunedm = matchJet->mass();
    pJet->nCharged   = itJet->chargedMultiplicity();
    pJet->nNeutrals  = itJet->neutralMultiplicity();
    pJet->nParticles = itJet->getPFConstituents().size();
    pJet->beta       = JetTools::beta(*itJet, pv);
    pJet->betaStar   = JetTools::betaStar(*itJet, pv, pvCol);
    pJet->dR2Mean    = JetTools::dR2Mean(*itJet);
    pJet->ptD        = JetTools::jetWidth(*itJet);
    pJet->q          = JetTools::jetCharge(*itJet);
    TVector2 lPull = JetTools::jetPull(*itJet);   //Color Flow observables
    pJet->pullPt     = lPull.Mod();
    pJet->pullEta    = lPull.X();
    pJet->pullPhi    = lPull.Y();
    pJet->pullAngle  = JetTools::jetPullAngle(*itJet,hSubJetProduct,fConeSize);
    pJet->mva = -2;
    if(passLoose) {
      double dRMean = JetTools::dRMean(*itJet);
      double frac01 = JetTools::frac(*itJet,0.1);
      double frac02 = JetTools::frac(*itJet,0.2);
      double frac03 = JetTools::frac(*itJet,0.3);
      double frac04 = JetTools::frac(*itJet,0.4);
      double frac05 = JetTools::frac(*itJet,0.5);
      
      fJetCorrForID->setJetPt(ptRaw);
      fJetCorrForID->setJetEta(itJet->eta());
      fJetCorrForID->setJetPhi(itJet->phi());
      fJetCorrForID->setJetE(itJet->energy());
      fJetCorrForID->setRho(*hRho);
      fJetCorrForID->setJetA(itJet->jetArea());
      fJetCorrForID->setJetEMF(-99.0);
      double jetcorrForID = fJetCorrForID->getCorrection();
      
      pJet->mva = fJetPUIDMVACalc.mvaValue((float)pvCol->size(), ptRaw*jetcorrForID, itJet->eta(), itJet->phi(),
			                   pJet->d0, pJet->dz, pJet->beta, pJet->betaStar, itJet->chargedMultiplicity(), itJet->neutralMultiplicity(),
			                   dRMean, pJet->dR2Mean, pJet->ptD, frac01, frac02, frac03, frac04, frac05);
    }
    
    // Basic Noise Variables
    pJet->chEmFrac   = itJet->chargedEmEnergy() / itJet->energy();
    pJet->neuEmFrac  = itJet->neutralEmEnergy() / itJet->energy();
    pJet->chHadFrac  = itJet->chargedHadronEnergy() / itJet->energy();
    pJet->neuHadFrac = itJet->neutralHadronEnergy() / itJet->energy();
    //
    // Generator matching
    //==============================
    const reco::GenJet * matchGenJet   = 0; 
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    if(matchGenJet != 0) { 
      pJet->mcFlavor               = (*jetFlavourMatch)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour();
      pJet->mcFlavorPhys           = (*jetFlavourMatchPhys)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour();
      pJet->genpt                  = matchGenJet->pt();
      pJet->geneta                 = matchGenJet->eta();
      pJet->genphi                 = matchGenJet->phi();
      pJet->genm                   = matchGenJet->mass();
    }
    pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, triggerEvent);

    ////Add Extras
    baconhep::TAddJet *pAddJet = 0; 
    if(fComputeFullJetInfo && itJet->pt() > 100 ) {
      assert(rExtraArray.GetEntries() < rExtraArray.GetSize());
      const int extraIndex = rExtraArray.GetEntries();
      new(rExtraArray[extraIndex]) baconhep::TAddJet();
      pAddJet = (baconhep::TAddJet*)rExtraArray[extraIndex];
      pAddJet->index = index;
    }
    if(fComputeFullJetInfo && itJet->pt() > 100)                        addJet(pAddJet,*itJet,*(hRho.product()));

    baconhep::TTopJet *pTopJet = 0; 
    if(fComputeFullJetInfo && itJet->pt() > 150.) {
      assert(rTopArray.GetEntries() < rTopArray.GetSize());
      const int topIndex = rTopArray.GetEntries();
      new(rTopArray[topIndex]) baconhep::TTopJet();
      pTopJet = (baconhep::TTopJet*)rTopArray[topIndex];
      pTopJet->index = index;
    }
    if(fComputeFullJetInfo && itJet->pt() >  150.) topJet(pTopJet,*itJet,*(hRho.product()));
  } 
}

//--------------------------------------------------------------------------------------------------
void FillerJet::addJet(baconhep::TAddJet *pPFJet,const reco::PFJet &itJet,double iRho) { 
  std::vector<reco::PFCandidatePtr> pfConstituents = itJet.getPFConstituents(); 
  std::vector<fastjet::PseudoJet>  lClusterParticles;
  for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {
    reco::PFCandidatePtr pfcand = pfConstituents[ic];
    fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());
    lClusterParticles.push_back(pPart);
  }
  fClustering = new fastjet::ClusterSequenceArea(lClusterParticles, *fCAJetDef, *fAreaDefinition);
  fastjet::PseudoJet iJet = CACluster(iJet,*fClustering);
  fastjet::PseudoJet pP1Jet = (*fPruner1)( iJet);
  double pCorr        = correction(pP1Jet,iRho);
  pPFJet->pt_p1       = pP1Jet.pt()*pCorr;
  pPFJet->ptraw_p1    = pP1Jet.pt();
  pPFJet->eta_p1      = pP1Jet.eta();
  pPFJet->phi_p1      = pP1Jet.phi();
  pPFJet->mass_p1     = pP1Jet.m()*pCorr;
  pPFJet->area_p1     = pP1Jet.area();

  fastjet::PseudoJet pP2Jet = (*fPruner2)( iJet);
  pCorr               = correction(pP2Jet,iRho);
  pPFJet->pt_p2       = pP2Jet.pt()*pCorr;
  pPFJet->ptraw_p2    = pP2Jet.pt();
  pPFJet->eta_p2      = pP2Jet.eta();
  pPFJet->phi_p2      = pP2Jet.phi();
  pPFJet->mass_p2     = pP2Jet.m()*pCorr;
  pPFJet->area_p2     = pP2Jet.area();

  fastjet::PseudoJet pT1Jet = (*fTrimmer1)( iJet);
  pCorr               = correction(pT1Jet,iRho);
  pPFJet->pt_t1       = pT1Jet.pt()*pCorr;
  pPFJet->ptraw_t1    = pT1Jet.pt();
  pPFJet->eta_t1      = pT1Jet.eta();
  pPFJet->phi_t1      = pT1Jet.phi();
  pPFJet->mass_t1     = pT1Jet.m()*pCorr;
  pPFJet->area_t1     = pT1Jet.area();
  
  fastjet::PseudoJet pT2Jet = (*fTrimmer2)( iJet);
  pCorr               = correction(pT2Jet,iRho);
  pPFJet->pt_t2       = pT2Jet.pt()*pCorr;
  pPFJet->ptraw_t2    = pT2Jet.pt();
  pPFJet->eta_t2      = pT2Jet.eta();
  pPFJet->phi_t2      = pT2Jet.phi();
  pPFJet->mass_t2     = pT2Jet.m()*pCorr;
  pPFJet->area_t2     = pT2Jet.area();
  
  //fastjet::PseudoJet pT3Jet = (*fTrimmer3)( iJet);
  //pCorr               = correction(pT3Jet,iRho);
  //pPFJet->pt_t3       = pT3Jet.pt()*pCorr;
  //pPFJet->ptraw_t3    = pT3Jet.pt();
  //pPFJet->eta_t3      = pT3Jet.eta();
  //pPFJet->phi_t3      = pT3Jet.phi();
  //pPFJet->mass_t3     = pT3Jet.m()*pCorr;
  //pPFJet->area_t3     = pT3Jet.area();
  
  fastjet::PseudoJet pT4Jet = (*fTrimmer4)( iJet);
  pCorr               = correction(pT4Jet,iRho);
  pPFJet->pt_t4       = pT4Jet.pt()*pCorr;
  pPFJet->ptraw_t4    = pT4Jet.pt();
  pPFJet->eta_t4      = pT4Jet.eta();
  pPFJet->phi_t4      = pT4Jet.phi();
  pPFJet->mass_t4     = pT4Jet.m()*pCorr;
  pPFJet->area_t4     = pT4Jet.area();

  /*
  fastjet::PseudoJet pF1Jet = (*fFilter1)( iJet);
  pCorr               = correction(pF1Jet,iRho);
  pPFJet->pt_f1       = pF1Jet.pt()*pCorr;
  pPFJet->ptraw_f1    = pF1Jet.pt();
  pPFJet->eta_f1      = pF1Jet.eta();
  pPFJet->phi_f1      = pF1Jet.phi();
  pPFJet->mass_f1     = pF1Jet.m()*pCorr;
  pPFJet->area_f1     = pF1Jet.area();
  
  fastjet::PseudoJet pF2Jet = (*fFilter2)( iJet);
  pCorr               = correction(pF2Jet,iRho);
  pPFJet->pt_f2       = pF2Jet.pt()*pCorr;
  pPFJet->ptraw_f2    = pF2Jet.pt();
  pPFJet->eta_f2      = pF2Jet.eta();
  pPFJet->phi_f2      = pF2Jet.phi();
  pPFJet->mass_f2     = pF2Jet.m()*pCorr;
  pPFJet->area_f2     = pF2Jet.area();
  */
                                            
  fastjet::PseudoJet pM1Jet = (*fSoftDrop1)( iJet);
  pCorr               = correction(pM1Jet,iRho);
  pPFJet->pt_m1       = pM1Jet.pt()*pCorr;
  pPFJet->ptraw_m1    = pM1Jet.pt(); 
  pPFJet->eta_m1      = pM1Jet.eta(); 
  pPFJet->phi_m1      = pM1Jet.phi(); 
  pPFJet->mass_m1     = pM1Jet.m()*pCorr; 
  pPFJet->area_m1     = pM1Jet.area();                                                                                                                                                                       
  fastjet::PseudoJet pM2Jet = (*fSoftDrop2)( iJet);                                                                                                                                       
  pCorr               = correction(pM2Jet,iRho);
  pPFJet->pt_m2       = pM2Jet.pt()*pCorr;
  pPFJet->ptraw_m2    = pM2Jet.pt(); 
  pPFJet->eta_m2      = pM2Jet.eta(); 
  pPFJet->phi_m2      = pM2Jet.phi(); 
  pPFJet->mass_m2     = pM2Jet.m()*pCorr; 
  pPFJet->area_m2     = pM2Jet.area();      

  fastjet::PseudoJet pM3Jet = (*fSoftDrop3)( iJet);                                                                                                                                                 
  pCorr               = correction(pM3Jet,iRho);
  pPFJet->pt_m3       = pM3Jet.pt()*pCorr;
  pPFJet->ptraw_m3    = pM3Jet.pt(); 
  pPFJet->eta_m3      = pM3Jet.eta(); 
  pPFJet->phi_m3      = pM3Jet.phi(); 
  pPFJet->mass_m3     = pM3Jet.m()*pCorr; 
  pPFJet->area_m3     = pM3Jet.area();      
    
  //Jet Shape Correlation observables
  fastjet::JetDefinition lCJet_def    (fastjet::cambridge_algorithm, 2.0);
  fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
  std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
  fastjet::EnergyCorrelatorRatio C2beta0 (2,0. ,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorRatio C2beta02(2,0.2,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorRatio C2beta05(2,0.5,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorRatio C2beta10(2,1.0,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorRatio C2beta20(2,2.0,fastjet::EnergyCorrelator::pt_R);
  pPFJet->c2_0    = C2beta0 (inclusive_jets[0]);
  pPFJet->c2_0P2  = C2beta02(inclusive_jets[0]);
  pPFJet->c2_0P5  = C2beta05(inclusive_jets[0]);
  pPFJet->c2_1P0  = C2beta10(inclusive_jets[0]);
  pPFJet->c2_2P0  = C2beta20(inclusive_jets[0]);
  //QJets
  pPFJet->qjet    = 0; 
  if(itJet.pt() > 100) pPFJet->qjet  = JetTools::qJetVolatility(lClusterParticles,25.,fRand->Rndm());
  //Subjet q/g
  std::vector<fastjet::PseudoJet>  lSubJets = pP1Jet.pieces();//fClustering->exclusive_subjets_up_to(iJet, 2.);
  if(lSubJets.size() > 0) pPFJet->sj1_npart = float(lSubJets[0].constituents().size());
  if(lSubJets.size() > 0) pPFJet->sj1_ptd   = JetTools::jetWidth(lSubJets[0],6); 
  if(lSubJets.size() > 0) pPFJet->sj1_maxW  = JetTools::jetWidth(lSubJets[0],1); 
  if(lSubJets.size() > 0) pPFJet->sj1_minW  = JetTools::jetWidth(lSubJets[0],2); 
  if(lSubJets.size() > 1) pPFJet->sj2_npart = float(lSubJets[1].constituents().size());
  if(lSubJets.size() > 1) pPFJet->sj2_ptd   = JetTools::jetWidth(lSubJets[1],6); 
  if(lSubJets.size() > 1) pPFJet->sj2_maxW  = JetTools::jetWidth(lSubJets[1],1); 
  if(lSubJets.size() > 1) pPFJet->sj2_minW  = JetTools::jetWidth(lSubJets[1],2); 
  delete fClustering;
}
void FillerJet::topJet(TTopJet *pPFJet,const reco::PFJet &itJet,double iRho) {
  std::vector<reco::PFCandidatePtr> pfConstituents = itJet.getPFConstituents(); 
  std::vector<fastjet::PseudoJet>  lClusterParticles;
  for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {
    reco::PFCandidatePtr pfcand = pfConstituents[ic];
    fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());
    lClusterParticles.push_back(pPart);
  }
  fClustering = new fastjet::ClusterSequenceArea(lClusterParticles, *fCAJetDef, *fAreaDefinition);
  fastjet::PseudoJet iJet = CACluster(iJet,*fClustering);
  //Top Tagging
  //fastjet::JetDefinition::Plugin *plugin = new fastjet::SISConePlugin(0.6, 0.75);
  //fastjet::JetDefinition NsubJetDef(plugin);
  //fNSUBTagger   = std::auto_ptr<fastjet::RestFrameNSubjettinessTagger>(new fastjet::RestFrameNSubjettinessTagger(NsubJetDef));

  fastjet::PseudoJet cmsTopJet = fCMSTopTagger->result(iJet);
  bool lCheckTop = (cmsTopJet.structure_non_const_ptr() == 0);
  if(!lCheckTop) { 
    std::vector<fastjet::PseudoJet> lPieces =  cmsTopJet.pieces();
    fastjet::PseudoJet pW1jet    = lPieces[0];//((fastjet::CMSTopTaggerStructure*) cmsTopJet.structure_non_const_ptr())->W1();
    fastjet::PseudoJet pW2jet    = lPieces[1];//((fastjet::CMSTopTaggerStructure*) cmsTopJet.structure_non_const_ptr())->W2();
    fastjet::PseudoJet pNWjet    = lPieces[2];//((fastjet::CMSTopTaggerStructure*) cmsTopJet.structure_non_const_ptr())->non_W();
    if(cmsTopJet.pieces().size() > 3) std::cout << cmsTopJet.pt() << " ===> Missing Pt " << pW1jet.pt() << " - " << pW2jet.pt() << " - " << pNWjet.pt()  << " -- " << lPieces[3].pt() << std::endl;
    double pCorr        = correction(cmsTopJet,iRho);
    pPFJet->pt_cms      = float(cmsTopJet.pt())*float(pCorr);
    pPFJet->ptraw_cms   = cmsTopJet.pt();
    pPFJet->eta_cms     = cmsTopJet.eta(); 
    pPFJet->phi_cms     = cmsTopJet.phi();
    pPFJet->mass_cms    = cmsTopJet.m()*pCorr; 
    pPFJet->area_cms    = cmsTopJet.area();        
    pCorr               = correction(pW1jet,iRho);
    pPFJet->pt_cms1     = pW1jet.pt()*pCorr;
    pPFJet->ptraw_cms1  = pW1jet.pt(); 
    pPFJet->eta_cms1    = pW1jet.eta(); 
    pPFJet->phi_cms1    = pW1jet.phi(); 
    pPFJet->mass_cms1   = pW1jet.m()*pCorr; 
    pPFJet->area_cms1   = pW1jet.area();        
    pCorr               = correction(pW2jet,iRho);
    pPFJet->pt_cms2     = pW2jet.pt()*pCorr;
    pPFJet->ptraw_cms2  = pW2jet.pt(); 
    pPFJet->eta_cms2    = pW2jet.eta(); 
    pPFJet->phi_cms2    = pW2jet.phi(); 
    pPFJet->mass_cms2   = pW2jet.m()*pCorr; 
    pPFJet->area_cms2   = pW2jet.area();        
    pCorr               = correction(pNWjet,iRho);
    pPFJet->pt_cms3     = pNWjet.pt()*pCorr;
    pPFJet->ptraw_cms3  = pNWjet.pt(); 
    pPFJet->eta_cms3    = pNWjet.eta(); 
    pPFJet->phi_cms3    = pNWjet.phi(); 
    pPFJet->mass_cms3   = pNWjet.m()*pCorr; 
    pPFJet->area_cms3   = pNWjet.area();        
  }
  fastjet::PseudoJet hepTopJet = fHEPTopTagger->result(iJet);
  bool lCheckTop2 = (hepTopJet.structure_non_const_ptr() == 0);
  if(!lCheckTop2) { 
    std::vector<fastjet::PseudoJet> lPieces = hepTopJet.pieces();
    fastjet::PseudoJet pHW1jet    = ((fastjet::HEPTopTaggerStructure*) hepTopJet.structure_non_const_ptr())->W1();
    fastjet::PseudoJet pHW2jet    = ((fastjet::HEPTopTaggerStructure*) hepTopJet.structure_non_const_ptr())->W2();
    fastjet::PseudoJet pHNWjet    = ((fastjet::HEPTopTaggerStructure*) hepTopJet.structure_non_const_ptr())->non_W();
   
    double pCorr              = 1.;//correction(hepTopJet,iRho);
    pPFJet->pt_htt     = hepTopJet.pt()*pCorr;
    pPFJet->ptraw_htt  = hepTopJet.pt(); 
    pPFJet->eta_htt    = hepTopJet.eta(); 
    pPFJet->phi_htt    = hepTopJet.phi(); 
    pPFJet->mass_htt   = hepTopJet.m()*pCorr; 
    //pPFJet->area_htt   = hepTopJet.area();        
    pCorr              = 1.;//correction(pHW1jet,iRho);
    pPFJet->pt_htt1    = pHW1jet.pt()*pCorr;
    pPFJet->ptraw_htt1 = pHW1jet.pt(); 
    pPFJet->eta_htt1   = pHW1jet.eta(); 
    pPFJet->phi_htt1   = pHW1jet.phi(); 
    pPFJet->mass_htt1  = pHW1jet.m()*pCorr; 
    //pPFJet->area_htt1  = pHW1jet.area();        
   
    pCorr               = 1.;///correction(pHW2jet,iRho);
    pPFJet->pt_htt2     = pHW2jet.pt()*pCorr;
    pPFJet->ptraw_htt2  = pHW2jet.pt(); 
    pPFJet->eta_htt2    = pHW2jet.eta(); 
    pPFJet->phi_htt2    = pHW2jet.phi(); 
    pPFJet->mass_htt2   = pHW2jet.m()*pCorr; 
    //pPFJet->area_htt2   = pHW2jet.area();        
   	
    pCorr               = 1.;//correction(pHNWjet,iRho);
    pPFJet->pt_htt3     = pHNWjet.pt()*pCorr;
    pPFJet->ptraw_htt3  = pHNWjet.pt(); 
    pPFJet->eta_htt3    = pHNWjet.eta(); 
    pPFJet->phi_htt3    = pHNWjet.phi(); 
    pPFJet->mass_htt3   = pHNWjet.m()*pCorr; 
    //pPFJet->area_htt3   = pHNWjet.area();        
  }
  delete fClustering;
}

//--------------------------------------------------------------------------------------------------
fastjet::PseudoJet FillerJet::CACluster   (fastjet::PseudoJet &iJet, fastjet::ClusterSequenceArea &iCAClustering) { 
  std::vector<fastjet::PseudoJet>  lOutJets = sorted_by_pt(iCAClustering.inclusive_jets(0.0));
  return lOutJets[0];
}

//--------------------------------------------------------------------------------------------------
const reco::BasicJet* FillerJet::match( const reco::PFJet *iJet,const reco::BasicJetCollection *jets ) { 
  int lId = -1;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::BasicJet* jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if ( dR < 0.25 ) {
      lId = i;
      break;
    }
  }
  const reco::BasicJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

//--------------------------------------------------------------------------------------------------
const reco::GenJet* FillerJet::match( const reco::PFJet *iJet,const reco::GenJetCollection *jets ) { 
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
