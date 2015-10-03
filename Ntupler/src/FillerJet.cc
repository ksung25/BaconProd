#include "BaconProd/Ntupler/interface/FillerJet.hh"
#include "BaconProd/Utils/interface/EnergyCorrelator.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/JetTools.hh"
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
FillerJet::FillerJet(const edm::ParameterSet &iConfig):
  fMinPt              (iConfig.getUntrackedParameter<double>("minPt",20)),
  fConeSize           (iConfig.getUntrackedParameter<double>("coneSize",0.5)),
  fUseGen             (iConfig.getUntrackedParameter<bool>("doGenJet",true)),
  fPVName             (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fRhoName            (iConfig.getUntrackedParameter<std::string>("edmRhoName","kt6PFJets")),
  fJetName            (iConfig.getUntrackedParameter<std::string>("jetName","ak5PFJets")),
  fGenJetName         (iConfig.getUntrackedParameter<std::string>("genJetName","AKT5GenJets")),
  fJetFlavorName      (iConfig.getUntrackedParameter<std::string>("jetFlavorName","AKT5byValAlgo")),
  fJetFlavorPhysName  (iConfig.getUntrackedParameter<std::string>("jetFlavorPhysName","AKT5byValPhys")),
  fPruneJetName       (iConfig.getUntrackedParameter<std::string>("pruneJetName","ca8PrunedJets")),
  fSubJetName         (iConfig.getUntrackedParameter<std::string>("subJetName","ca8PrunedJets__SubJets")),
  fCSVbtagName        (iConfig.getUntrackedParameter<std::string>("csvBTagName","myCombinedSecondaryVertexBJetTags")),
  fCSVbtagSubJetName  (iConfig.getUntrackedParameter<std::string>("csvBTagSubJetName","myCombinedSecondaryVertexBJetTagsSubJets")),
  fJettinessName      (iConfig.getUntrackedParameter<std::string>("jettiness","Njettiness")),
  fQGLikelihood       (iConfig.getUntrackedParameter<std::string>("qgLikelihood","QGLikelihood")),
  fQGLikelihoodSubJets(iConfig.getUntrackedParameter<std::string>("qgLikelihoodSubjet","QGLikelihood")),
  fTopTagType         (iConfig.getUntrackedParameter<std::string>("topTagType","CMS")),
  fComputeFullJetInfo (iConfig.getUntrackedParameter<bool>("doComputeFullJetInfo",true)),  
  fJetDef             (0),
  fGenJetDef          (0),
  fCAJetDef           (0),
  fActiveArea         (0),
  fAreaDefinition     (0),
  fClustering         (0),
  fPruner             (0),
  fSoftDrop0          (0),
  fSoftDrop1          (0),
  fTrimmer            (0), 
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

  fPruner = new fastjet::Pruner( fastjet::cambridge_algorithm, 0.1, 0.5); //CMS Default

  fSoftDrop0 = new fastjet::contrib::SoftDropTagger(0., 0.10, 1.0);
  fSoftDrop1 = new fastjet::contrib::SoftDropTagger(1., 0.15, 1.0);

  fTrimmer  = new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2),  fastjet::SelectorPtFractionMin(0.05)));
  
  fCMSTopTagger = new fastjet::CMSTopTagger();//0.05,0.8,0.19);
  fHEPTopTagger = new fastjet::HEPTopTagger(0.8,30.,false);

  
  std::vector<std::string> empty_vstring;
  initJetCorr(iConfig.getUntrackedParameter< std::vector<std::string> >("jecFiles",empty_vstring),
              iConfig.getUntrackedParameter< std::vector<std::string> >("jecUncFiles",empty_vstring),
	      iConfig.getUntrackedParameter< std::vector<std::string> >("jecFilesForID",empty_vstring));

  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";

  std::vector<std::string> puIDFiles = iConfig.getUntrackedParameter< std::vector<std::string> >("jetPUIDFiles",empty_vstring);
  assert(puIDFiles.size()==2);
  std::string lowPtWeightFile  = (puIDFiles[0].length()>0) ? (cmssw_base_src + puIDFiles[0]) : "";
  std::string highPtWeightFile = (puIDFiles[1].length()>0) ? (cmssw_base_src + puIDFiles[1]) : "";
  
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
  delete fPruner;
  delete fSoftDrop0;
  delete fSoftDrop1;
  delete fTrimmer;
  delete fCMSTopTagger;
  delete fHEPTopTagger;
  
  delete fJetCorr;
  delete fJetCorrForID;
  delete fJetUnc;
}

//--------------------------------------------------------------------------------------------------
void FillerJet::initJetCorr(const std::vector<std::string> &jecFiles,
                            const std::vector<std::string> &jecUncFiles,
			    const std::vector<std::string> &jecFilesForID)
{
  assert(jecFiles.size()>0);
  assert(jecUncFiles.size()>0);
  assert(jecFilesForID.size()>0);
  
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  
  std::vector<JetCorrectorParameters> corrParams;
  for(unsigned int icorr=0; icorr<jecFiles.size(); icorr++) {
    corrParams.push_back(JetCorrectorParameters( (cmssw_base_src + jecFiles[icorr]).c_str() ));
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
double FillerJet::correction(fastjet::PseudoJet &iJet, double iRho) { 
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
void FillerJet::fill(TClonesArray *array, TClonesArray *iExtraArray,
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
  
  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  edm::InputTag rhoTag(fRhoName,"rho");
  iEvent.getByLabel(rhoTag,hRho);
  assert(hRho.isValid()); 
 
  // Get b-jet tagger
  edm::Handle<reco::JetTagCollection> hCSVbtags;
  iEvent.getByLabel(fCSVbtagName, hCSVbtags);
  assert(hCSVbtags.isValid());

  //Get Quark Gluon Likelihood
  edm::Handle<edm::ValueMap<float> > hQGLikelihood; 
  iEvent.getByLabel(fQGLikelihood,"qgLikelihood",hQGLikelihood); 
  assert(hQGLikelihood.isValid());


  TClonesArray &rArray      = *array;
  TClonesArray &rExtraArray = *iExtraArray;
  for(reco::PFJetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {
    const double ptRaw = itJet->pt();

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
    baconhep::TJet *pJet = (baconhep::TJet*)rArray[index];
 
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
    pJet->csv  = (*(hCSVbtags.product()))    [jetBaseRef];
    pJet->qgid = (*(hQGLikelihood.product()))[jetBaseRef];
    
    pJet->nCharged   = itJet->chargedMultiplicity();
    pJet->nNeutrals  = itJet->neutralMultiplicity();
    pJet->nParticles = itJet->getPFConstituents().size();
    pJet->beta       = JetTools::beta(*itJet, pv);
    pJet->betaStar   = JetTools::betaStar(*itJet, pv, pvCol);
    pJet->dR2Mean    = JetTools::dR2Mean(*itJet);
    pJet->ptD        = JetTools::jetWidth(*itJet);
    pJet->q          = JetTools::jetCharge(*itJet,false);

    TVector2 lPull    = JetTools::jetPull(*itJet,0);
    pJet->pullY       = lPull.X();
    pJet->pullPhi     = lPull.Y();
    TVector2 lChPull  = JetTools::jetPull(*itJet,1);
    pJet->chPullY     = lChPull.X();
    pJet->chPullPhi   = lChPull.Y();
    TVector2 lNeuPull = JetTools::jetPull(*itJet,2);
    pJet->neuPullY    = lNeuPull.X();
    pJet->neuPullPhi  = lNeuPull.Y();

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
    const reco::GenJet *matchGenJet = 0; 
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    if(matchGenJet != 0) { 
      pJet->mcFlavor     = (*jetFlavourMatch)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour();
      pJet->mcFlavorPhys = (*jetFlavourMatchPhys)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour();
      pJet->genpt        = matchGenJet->pt();
      pJet->geneta       = matchGenJet->eta();
      pJet->genphi       = matchGenJet->phi();
      pJet->genm         = matchGenJet->mass();
    }
    pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, triggerEvent);

    ////Add Extras
    baconhep::TAddJet *pAddJet = 0; 
    if(fComputeFullJetInfo) {
      assert(rExtraArray.GetEntries() < rExtraArray.GetSize());
      const int extraIndex = rExtraArray.GetEntries();
      new(rExtraArray[extraIndex]) baconhep::TAddJet();
      pAddJet = (baconhep::TAddJet*)rExtraArray[extraIndex];
      pAddJet->index = index;
      addJet(pAddJet, iEvent, *itJet, jetBaseRef);
    }
  } 
}

//--------------------------------------------------------------------------------------------------
void FillerJet::addJet(baconhep::TAddJet *pAddJet, const edm::Event &iEvent, 
                       const reco::PFJet &itJet, const reco::JetBaseRef &jetBaseRef)
{ 
  // Get pruned jet collection
//  edm::Handle<reco::BasicJetCollection> hPruneJetProduct;
//  iEvent.getByLabel(fPruneJetName,hPruneJetProduct);
//  assert(hPruneJetProduct.isValid());
//  const reco::BasicJetCollection *pruneJetCol = hPruneJetProduct.product();

  // Get pruned sub jet collection
  edm::Handle<reco::PFJetCollection> hSubJetProduct;
  edm::InputTag subJetTag(fSubJetName,"SubJets");
  iEvent.getByLabel(subJetTag,hSubJetProduct);
  assert(hSubJetProduct.isValid());
  const reco::PFJetCollection *subJetCol = hSubJetProduct.product();

  // Get b sub-jets 
  edm::Handle<reco::JetTagCollection> hCSVbtagsSubJets;
  iEvent.getByLabel(fCSVbtagSubJetName, hCSVbtagsSubJets);
  assert(hCSVbtagsSubJets.isValid());

  //Get Quark Gluon Likelihood on subjets
  edm::Handle<edm::ValueMap<float> > hQGLikelihoodSubJets;
  iEvent.getByLabel(fQGLikelihoodSubJets,"qgLikelihood",hQGLikelihoodSubJets);
  assert(hQGLikelihoodSubJets.isValid());
  
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

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  edm::InputTag rhoTag(fRhoName,"rho");
  iEvent.getByLabel(rhoTag,hRho);
  assert(hRho.isValid());

  pAddJet->pullAngle = JetTools::jetPullAngle(itJet,hSubJetProduct,fConeSize);
  pAddJet->tau1 = (*(hTau1.product()))[jetBaseRef];
  pAddJet->tau2 = (*(hTau2.product()))[jetBaseRef];
  pAddJet->tau3 = (*(hTau3.product()))[jetBaseRef];


  std::vector<reco::PFCandidatePtr> pfConstituents = itJet.getPFConstituents(); 
  std::vector<fastjet::PseudoJet> lClusterParticles;
  for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {
    reco::PFCandidatePtr pfcand = pfConstituents[ic];
    fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());
    lClusterParticles.push_back(pPart);
  }
  fClustering = new fastjet::ClusterSequenceArea(lClusterParticles, *fCAJetDef, *fAreaDefinition);
  fastjet::PseudoJet iJet = CACluster(iJet,*fClustering);

  double pCorr=1;

  // Pruning
  fastjet::PseudoJet pP1Jet = (*fPruner)(iJet);
  pCorr = 1;//correction(pP1Jet,*hRho);
  pAddJet->mass_prun = pP1Jet.m()*pCorr;

  // Trimming
  fastjet::PseudoJet pT1Jet = (*fTrimmer)(iJet);
  pCorr = 1;//correction(pT1Jet,*hRho);
  pAddJet->mass_trim = pT1Jet.m()*pCorr;

  // Soft drop                                          
  fastjet::PseudoJet pM0Jet = (*fSoftDrop0)(iJet);
  pCorr = 1;//correction(pM0Jet,*hRho);
  pAddJet->mass_sd0 = pM0Jet.m()*pCorr;
  
  fastjet::PseudoJet pM1Jet = (*fSoftDrop1)(iJet);
  pCorr = 1;//correction(pM1Jet,*hRho);
  pAddJet->mass_sd1 = pM1Jet.m()*pCorr;																					       
    
  // Jet Shape Correlation observables
  fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 2.0);
  fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
  std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
  fastjet::EnergyCorrelatorDoubleRatio C2beta0 (2,0. ,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorDoubleRatio C2beta02(2,0.2,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorDoubleRatio C2beta05(2,0.5,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorDoubleRatio C2beta10(2,1.0,fastjet::EnergyCorrelator::pt_R);
  fastjet::EnergyCorrelatorDoubleRatio C2beta20(2,2.0,fastjet::EnergyCorrelator::pt_R);
  pAddJet->c2_0   = C2beta0 (inclusive_jets[0]);
  pAddJet->c2_0P2 = C2beta02(inclusive_jets[0]);
  pAddJet->c2_0P5 = C2beta05(inclusive_jets[0]);
  pAddJet->c2_1P0 = C2beta10(inclusive_jets[0]);
  pAddJet->c2_2P0 = C2beta20(inclusive_jets[0]);

  // Q-Jets
  pAddJet->qjet = 0; 
  if(itJet.pt() > 100) pAddJet->qjet = JetTools::qJetVolatility(lClusterParticles,25.,fRand->Rndm());  // (!) why pT > 100 cut? computation time?

  //
  // Subjets
  //

  // find/sort up to 4 hardest subjets
  const reco::PFJet *subjet1=0, *subjet2=0, *subjet3=0, *subjet4=0;
  double csv1=-2, csv2=-2, csv3=-2, csv4=-2;
  double qgid1=-2, qgid2=-2, qgid3=-2, qgid4=-2;
  double q1=-100, q2=-100, q3=-100, q4=-100;
  for(reco::PFJetCollection::const_iterator itSubJet = subJetCol->begin(); itSubJet!=subJetCol->end(); ++itSubJet) {
    if(reco::deltaR(itJet.eta(),itJet.phi(),itSubJet->eta(),itSubJet->phi())>fConeSize) continue;  // (!) get associated subjets by dR...is there a better way???

    reco::PFJetRef subjetRef(hSubJetProduct, itSubJet - subJetCol->begin());
    reco::JetBaseRef subjetBaseRef(subjetRef);

    if(!subjet1 || itSubJet->pt() > subjet1->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;
      
      subjet3 = subjet2;
      csv3    = csv2;
      qgid3   = qgid2;
      q3      = q2;
      
      subjet2 = subjet1;
      csv2    = csv1;
      qgid2   = qgid1;
      q2      = q1;
      
      subjet1 = &(*itSubJet);      
      csv1    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];      
      qgid1   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q1      = JetTools::jetCharge(*itSubJet,false);
      
    } else if(!subjet2 || itSubJet->pt() > subjet2->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;

      subjet3 = subjet2;
      csv3    = csv2;
      qgid3   = qgid2;
      q3      = q2;

      subjet2 = &(*itSubJet);      
      csv2    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef]; 
      qgid2   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
      q2      = JetTools::jetCharge(*itSubJet,false);
      
    } else if(!subjet3 || itSubJet->pt() > subjet3->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;

      subjet3 = &(*itSubJet);
      csv3    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      qgid3   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
      q3      = JetTools::jetCharge(*itSubJet,false);
      
    } else if(!subjet4 || itSubJet->pt() > subjet4->pt()) {
      subjet4 = &(*itSubJet);
      csv4    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      qgid4   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q4      = JetTools::jetCharge(*itSubJet,false);
    }
  }
  if(subjet1) {
    pAddJet->sj1_pt   = subjet1->pt();
    pAddJet->sj1_eta  = subjet1->eta();
    pAddJet->sj1_phi  = subjet1->phi();
    pAddJet->sj1_m    = subjet1->mass();
    pAddJet->sj1_csv  = csv1;
    pAddJet->sj1_qgid = qgid1;
    pAddJet->sj1_q    = q1;
  }
  if(subjet2) {
    pAddJet->sj2_pt   = subjet2->pt();
    pAddJet->sj2_eta  = subjet2->eta();
    pAddJet->sj2_phi  = subjet2->phi();
    pAddJet->sj2_m    = subjet2->mass();
    pAddJet->sj2_csv  = csv2;
    pAddJet->sj2_qgid = qgid2;
    pAddJet->sj2_q    = q2;
  }
  if(subjet3) {
    pAddJet->sj3_pt   = subjet3->pt();
    pAddJet->sj3_eta  = subjet3->eta();
    pAddJet->sj3_phi  = subjet3->phi();
    pAddJet->sj3_m    = subjet3->mass();
    pAddJet->sj3_csv  = csv3;
    pAddJet->sj3_qgid = qgid3;
    pAddJet->sj3_q    = q3;
  }
  if(subjet4) {
    pAddJet->sj4_pt   = subjet4->pt();
    pAddJet->sj4_eta  = subjet4->eta();
    pAddJet->sj4_phi  = subjet4->phi();
    pAddJet->sj4_m    = subjet4->mass();
    pAddJet->sj4_csv  = csv4;
    pAddJet->sj4_qgid = qgid4;
    pAddJet->sj4_q    = q4;
  } 
  
  //
  // Top Tagging
  //
  if(fTopTagType.compare("CMS")==0) {  // CMS Top tagger (R=0.8)
    fastjet::PseudoJet cmsTopJet = fCMSTopTagger->result(iJet);
    bool lCheckTop = (cmsTopJet.structure_non_const_ptr() == 0);
    if(!lCheckTop) {
      pAddJet->topTagType |= kCMSTT;

      double pCorr= 1;//correction(cmsTopJet,*hRho);

      // order the subjets in the following order:
      //  - hardest of the W subjets
      //  - softest of the W subjets
      //  - hardest of the remaining subjets
      //  - softest of the remaining subjets (if any)
      std::vector<fastjet::PseudoJet> lPieces = cmsTopJet.pieces();
      fastjet::PseudoJet pW1jet = lPieces[0];
      fastjet::PseudoJet pW2jet = lPieces[1];
      fastjet::PseudoJet pNWjet = lPieces[2];
      
      TLorentzVector vSubjet1; vSubjet1.SetPtEtaPhiM(pW1jet.pt()*pCorr, pW1jet.eta(), pW1jet.phi(), pW1jet.m()*pCorr);
      TLorentzVector vSubjet2; vSubjet2.SetPtEtaPhiM(pW2jet.pt()*pCorr, pW2jet.eta(), pW2jet.phi(), pW2jet.m()*pCorr);
      TLorentzVector vSubjet3; vSubjet3.SetPtEtaPhiM(pNWjet.pt()*pCorr, pNWjet.eta(), pNWjet.phi(), pNWjet.m()*pCorr);
      double m12 = (vSubjet1+vSubjet2).M();
      double m23 = (vSubjet2+vSubjet3).M();
      double m31 = (vSubjet3+vSubjet1).M();
      pAddJet->top_m_min = 0;
      if     (m12<m23 && m12<m31) { pAddJet->top_m_min = m12; }
      else if(m23<m12 && m23<m31) { pAddJet->top_m_min = m23; }
      else if(m31<m12 && m31<m23) { pAddJet->top_m_min = m31; }
      
      pAddJet->top_m_123 = (vSubjet1 + vSubjet2 + vSubjet3).M();
    }
    
  } else if(fTopTagType.compare("HEP")==0) {  // HEP Top tagger (R=1.5)
    fastjet::PseudoJet hepTopJet = fHEPTopTagger->result(iJet);
    bool lCheckTop = (hepTopJet.structure_non_const_ptr() == 0);
    if(!lCheckTop) { 
      pAddJet->topTagType |= kHEPTT;
   
      double pCorr = 1;//correction(hepTopJet,*hRho);

      // Only 3 subjets in HEP top tagger
      fastjet::PseudoJet pW1jet = ((fastjet::HEPTopTaggerStructure*) hepTopJet.structure_non_const_ptr())->W1();
      fastjet::PseudoJet pW2jet = ((fastjet::HEPTopTaggerStructure*) hepTopJet.structure_non_const_ptr())->W2();
      fastjet::PseudoJet pNWjet = ((fastjet::HEPTopTaggerStructure*) hepTopJet.structure_non_const_ptr())->non_W();
      
      TLorentzVector vSubjet1; vSubjet1.SetPtEtaPhiM(pW1jet.pt()*pCorr, pW1jet.eta(), pW1jet.phi(), pW1jet.m()*pCorr);
      TLorentzVector vSubjet2; vSubjet2.SetPtEtaPhiM(pW2jet.pt()*pCorr, pW2jet.eta(), pW2jet.phi(), pW2jet.m()*pCorr);
      TLorentzVector vSubjet3; vSubjet3.SetPtEtaPhiM(pNWjet.pt()*pCorr, pNWjet.eta(), pNWjet.phi(), pNWjet.m()*pCorr);
      double m12 = (vSubjet1+vSubjet2).M();
      double m23 = (vSubjet2+vSubjet3).M();
      double m31 = (vSubjet3+vSubjet1).M();
      pAddJet->top_m_min = 0;
      if     (m12<m23 && m12<m31) { pAddJet->top_m_min = m12; }
      else if(m23<m12 && m23<m31) { pAddJet->top_m_min = m23; }
      else if(m31<m12 && m31<m23) { pAddJet->top_m_min = m31; }
      
      pAddJet->top_m_123 = (vSubjet1 + vSubjet2 + vSubjet3).M();   
    }
  }
  
  delete fClustering;
}

//--------------------------------------------------------------------------------------------------
fastjet::PseudoJet FillerJet::CACluster(fastjet::PseudoJet &iJet, fastjet::ClusterSequenceArea &iCAClustering) { 
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(iCAClustering.inclusive_jets(0.0));
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
