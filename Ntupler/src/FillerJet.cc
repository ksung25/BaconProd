#include "BaconProd/Ntupler/interface/FillerJet.hh"
//#include "BaconProd/Utils/interface/EnergyCorrelator.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/JetTools.hh"
#include "BaconProd/Utils/interface/BJetNNRegression.hh"
//#include "BaconProd/Utils/interface/RecursiveSoftDrop.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TSVtx.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"

//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
//#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include <fastjet/JetDefinition.hh>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace baconhep;


//--------------------------------------------------------------------------------------------------
FillerJet::FillerJet(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt              (iConfig.getUntrackedParameter<double>("minPt",20)),
  fUseGen             (iConfig.getUntrackedParameter<bool>("doGenJet",true)),
  fApplyJEC           (iConfig.getUntrackedParameter<bool>("applyJEC",true)),
  fPVName             (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fRhoName            (iConfig.getUntrackedParameter<std::string>("edmRhoName","fixedGridRhoFastjetAll")),
  fJetName            (iConfig.getUntrackedParameter<std::string>("jetName","ak4PFJetsCHS")),
  fJECName            (iConfig.getUntrackedParameter<std::string>("jecName","")),
  fJECUName           (iConfig.getUntrackedParameter<std::string>("jecUncName","")),
  fGenJetName         (iConfig.getUntrackedParameter<std::string>("genJetName","slimmedGenJets")),
  fJetFlavorName      (iConfig.getUntrackedParameter<std::string>("jetFlavorName","AK4byValAlgoCHS")),
  fPrunedJetName      (iConfig.getUntrackedParameter<std::string>("prunedJetName","AK4caPFJetsPrunedCHS")),
  fTrimmedJetName     (iConfig.getUntrackedParameter<std::string>("trimmedJetName","AK4caPFJetsTrimmedCHS")),
  fSoftDropJetName    (iConfig.getUntrackedParameter<std::string>("softdropJetName","AK4caPFJetsSoftDropCHS")),
  fSubJetName         (iConfig.getUntrackedParameter<std::string>("subJetName","AK4caPFJetsTrimmedCHS__SubJets")),
  fCVLctagName        (iConfig.getUntrackedParameter<std::string>("cvlcTagName","AK4PFCombinedCvsLJetTagsCHS")),
  fSVName             (iConfig.getUntrackedParameter<std::string>("secVertices","slimmedSecondaryVertices")),
  fCVBctagName        (iConfig.getUntrackedParameter<std::string>("cvbcTagName","AK4PFCombinedCvsBJetTagsCHS")),
  fMVAbtagName        (iConfig.getUntrackedParameter<std::string>("mvaBTagName","AK4PFCombinedMVAV2BJetTagsCHS")),
  fCSVbtagName        (iConfig.getUntrackedParameter<std::string>("csvBTagName","combinedInclusiveSecondaryVertexV2BJetTags")),
  fCSVbtagSubJetName  (iConfig.getUntrackedParameter<std::string>("csvBTagSubJetName","AK4CombinedInclusiveSecondaryVertexV2BJetTagsSJCHS")),
  fCSVDoubleBtagName  (iConfig.getUntrackedParameter<std::string>("csvDoubleBTagName","AK8PFBoostedDoubleSecondaryVertexBJetTagsCHS")),
  fDeepCSVBtagName    (iConfig.getUntrackedParameter<std::string>("deepCSVBTagName","AK4PFDeepCSVJetTagsCHS")),
  fDeepCSVBtagNameb   (iConfig.getUntrackedParameter<std::string>("deepCSVBTagNameb","AK4PFDeepCSVJetTagsCHS:probb")),
  fDeepCSVBtagNamec   (iConfig.getUntrackedParameter<std::string>("deepCSVBTagNamec","AK4PFDeepCSVJetTagsCHS:probc")),
  fDeepCSVBtagNamel   (iConfig.getUntrackedParameter<std::string>("deepCSVBTagNamel","AK4PFDeepCSVJetTagsCHS:probudsg")),
  fDeepCSVBtagNamebb  (iConfig.getUntrackedParameter<std::string>("deepCSVBTagNamebb","AK4PFDeepCSVJetTagsCHS:probbb")),
  fDeepCMVABtagName   (iConfig.getUntrackedParameter<std::string>("deepCMVABTagName","AK4PFDeepCMVAJetTagsCHS")),
  fDeepCMVABtagNameb  (iConfig.getUntrackedParameter<std::string>("deepCMVABTagNameb","AK4PFDeepCMVAJetTagsCHS:probb")),
  fDeepCMVABtagNamec  (iConfig.getUntrackedParameter<std::string>("deepCMVABTagNamec","AK4PFDeepCMVAJetTagsCHS:probc")),
  fDeepCMVABtagNamel  (iConfig.getUntrackedParameter<std::string>("deepCMVABTagNamel","AK4PFDeepCMVAJetTagsCHS:probudsg")),
  fDeepCMVABtagNamebb (iConfig.getUntrackedParameter<std::string>("deepCMVABTagNamebb","AK4PFDeepCMVAJetTagsCHS:probbb")),
  //fSVTagInfoName      (iConfig.getUntrackedParameter<std::string>("svTagInfoName","AK4PFSecondaryVertexTagInfosCHS")),
  fBoostedDoubleSVTagInfoName (iConfig.getUntrackedParameter<std::string>("boostedDoubleSVTagInfoName","AK8PFBoostedDoubleSVTagInfosCHS")),
  fDeepDoubleBvLtagName (iConfig.getUntrackedParameter<std::string>("deepDoubleBvLTagName","AK8PFBoostedDeepDoubleBvLJetTagsCHS:probH")),
  fDeepDoubleBvLNoMassSculptPentagName (iConfig.getUntrackedParameter<std::string>("deepDoubleBvLNoMassSculptPenTagName","AK8PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsCHS:probH")),
  fDeepDoubleCvLtagName (iConfig.getUntrackedParameter<std::string>("deepDoubleCvLTagName","AK8PFBoostedDeepDoubleCvLJetTagsCHS:probH")),
  fDeepDoubleCvLNoMassSculptPentagName (iConfig.getUntrackedParameter<std::string>("deepDoubleCvLNoMassSculptPenTagName","AK8PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsCHS:probH")),
  fDeepDoubleCvBtagName (iConfig.getUntrackedParameter<std::string>("deepDoubleCvBTagName","AK8PFBoostedDeepDoubleCvBJetTagsCHS:probH")),
  fDeepDoubleCvBNoMassSculptPentagName (iConfig.getUntrackedParameter<std::string>("deepDoubleCvBNoMassSculptPenTagName","AK8PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsCHS:probH")),
  fBRegNNFile         (iConfig.getUntrackedParameter<std::string>("BRegNNFileName","BaconProd/Utils/data/breg_training_2017.pb")),
  fBRegNNMean         (iConfig.getUntrackedParameter<double>("BRegNNMean",1.0610932111740112)),
  fBRegNNStd          (iConfig.getUntrackedParameter<double>("BRegNNStd",0.39077115058898926)),
  //fMuonName           (iConfig.getUntrackedParameter<std::string>("edmMuonName","muons")),
  //fEleName            (iConfig.getUntrackedParameter<std::string>("edmElectronName","gedGsfElectrons")),
  //fsoftPFMuonTagInfoName    (iConfig.getUntrackedParameter<std::string>("softPFMuonTagInfoName","AK4PFSoftPFMuonsTagInfosCHS")),
  //fsoftPFElectronTagInfoName(iConfig.getUntrackedParameter<std::string>("softPFElectronTagInfoName","AK4PFSoftPFElectronsTagInfosCHS")),
  fJettinessName      (iConfig.getUntrackedParameter<std::string>("jettiness","AK4NjettinessCHS")),
  fQGLikelihood       (iConfig.getUntrackedParameter<std::string>("qgLikelihood","QGLikelihood")),
  fQGLikelihoodSubJets(iConfig.getUntrackedParameter<std::string>("qgLikelihoodSubjet","QGLikelihood")),
  fTopTaggerName      (iConfig.getUntrackedParameter<std::string>("topTaggerName","")),
  //fShowerDecoConf     (iConfig.getUntrackedParameter<std::string>("showerDecoConf","")),
  fConeSize           (iConfig.getUntrackedParameter<double>("coneSize",0.4)),
  fComputeFullJetInfo (iConfig.getUntrackedParameter<bool>("doComputeFullJetInfo",false)),  
  fAddPFCand          (iConfig.getUntrackedParameter<bool>("addPFCand",true)),
  fComputeSVInfo      (iConfig.getUntrackedParameter<bool>("doComputeSVInfo",false)),  
  fUseTO              (iConfig.getUntrackedParameter<bool>("useTriggerObject",false)),
  //fShowerDeco         (0),
  fJetCorr            (0),
  fJetUnc             (0),
  fUseAOD             (useAOD)
  //fRecursiveSoftDrop1 (0),
  //fRecursiveSoftDrop2 (0)
{
    std::vector<std::string> empty_vstring;
    /* ===> Switching to DB
    initJetCorr(iConfig.getUntrackedParameter< std::vector<std::string> >("jecFiles",empty_vstring),
                iConfig.getUntrackedParameter< std::vector<std::string> >("jecUncFiles",empty_vstring));
    */
    std::string cmssw_base_src = getenv("CMSSW_BASE");
    cmssw_base_src += "/src/";


   std::string empty_string; 
   std::string BoostedBtaggingFiles = iConfig.getUntrackedParameter< std::string >("jetBoostedBtaggingFiles",empty_string);
   fWeightFile  =  (cmssw_base_src + "BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml");//BoostedBtaggingFiles) ;
   if(fConeSize < 1.0)   fWeightFile  =  (cmssw_base_src + "BaconProd/Utils/data/BoostedDoubleSV_AK8_BDT_v4.weights.xml");
   std::cout<<"Double B-tag: "<< fWeightFile <<std::endl;

   fBRegFile  =  (cmssw_base_src + "BaconProd/Utils/data/ttbar-G25-500k-13d-300t.weights.xml");
   std::cout<<"BJet Regression "<< fBRegFile <<std::endl;
   initBReg();
   
   fBRegNNFile  =  (cmssw_base_src + fBRegNNFile);
   std::cout<<"BJet NN Regression "<< fBRegNNFile <<std::endl;
   initBRegNN();
   
   if(fUseAOD) {
      std::vector<std::string> puIDFiles = iConfig.getUntrackedParameter< std::vector<std::string> >("jetPUIDFiles",empty_vstring);
      assert(puIDFiles.size()==2);
      fLowPtWeightFile  = (puIDFiles[0].length()>0) ? (cmssw_base_src + puIDFiles[0]) : "";
      fHighPtWeightFile = (puIDFiles[1].length()>0) ? (cmssw_base_src + puIDFiles[1]) : "";
      //initPUJetId();
    }

  fRand = new TRandom2();
  //if(fShowerDecoConf.size() > 0) { 
  //  fShowerDeco = new ShowerDeco(cmssw_base_src+fShowerDecoConf);
  //}
  fTokRhoTag        = iC.consumes<double>               (fRhoName);
  if(fUseAOD)  fTokJetName       = iC.consumes<reco::PFJetCollection>(fJetName);
  if(!fUseAOD) fTokPatJetName    = iC.consumes<pat::JetCollection>   (fJetName);
  fTokJECName       = iC.consumes<reco::JetCorrector>(fJECName);
  fTokGenJetName    = iC.consumes<reco::GenJetCollection>(fGenJetName);
  fTokJetFlavorName = iC.consumes<reco::JetFlavourInfoMatchingCollection>(fJetFlavorName);
  fTokPVName        = iC.consumes<reco::VertexCollection>(fPVName);
  fTokCSVbtagName   = iC.consumes<reco::JetTagCollection>(fCSVbtagName);
  fTokMVAbtagName   = iC.consumes<reco::JetTagCollection>(fMVAbtagName);
  fTokCVBctagName   = iC.consumes<reco::JetTagCollection>(fCVBctagName);
  fTokCVLctagName   = iC.consumes<reco::JetTagCollection>(fCVLctagName);
  fTokSVName        = iC.consumes<reco::VertexCompositePtrCandidateCollection>(fSVName);
  //fTokSVTagInfoCollection = iC.consumes<std::vector<reco::CandIPTagInfo> >(fSVTagInfoName);
  //if(fUseAOD)  
  //fTokMuonName       = iC.consumes<reco::MuonCollection>       (fMuonName);
  //if(!fUseAOD) 
  //fTokPatMuonName    = iC.consumes<pat::MuonCollection>(iConfig.getUntrackedParameter<std::string>("edmMuonName","muons"));
  //if(fUseAOD)  
  //fTokEleName        = iC.consumes<reco::GsfElectronCollection>(iConfig.getUntrackedParameter<std::string>("edmElectronName","gedGsfElectrons"));
  //if(!fUseAOD) 
  //fTokPatEleName     = iC.consumes<pat::ElectronCollection>(iConfig.getUntrackedParameter<std::string>("edmElectronName","gedGsfElectrons"));
  
  edm::InputTag lQGLikelihood(fQGLikelihood,"qgLikelihood");
  edm::InputTag lQGLAxis2    (fQGLikelihood,"axis2");
  edm::InputTag lQGLPtD      (fQGLikelihood,"ptD");
  edm::InputTag lQGLMult     (fQGLikelihood,"mult");
  fTokQGLikelihood        = iC.consumes<edm::ValueMap<float> >   (lQGLikelihood);
  fTokQGLAxis2            = iC.consumes<edm::ValueMap<float> >   (lQGLAxis2);
  fTokQGLPtD              = iC.consumes<edm::ValueMap<float> >   (lQGLPtD);
  fTokQGLMult             = iC.consumes<edm::ValueMap<int> >     (lQGLMult);
  fTokDeepCSVBtagNameb    = iC.consumes<reco::JetTagCollection>  (fDeepCSVBtagNameb);
  fTokDeepCSVBtagNamec    = iC.consumes<reco::JetTagCollection>  (fDeepCSVBtagNamec);
  fTokDeepCSVBtagNamel    = iC.consumes<reco::JetTagCollection>  (fDeepCSVBtagNamel);
  fTokDeepCSVBtagNamebb   = iC.consumes<reco::JetTagCollection>  (fDeepCSVBtagNamebb);
  fTokDeepCMVABtagNameb   = iC.consumes<reco::JetTagCollection>  (fDeepCMVABtagNameb);
  fTokDeepCMVABtagNamec   = iC.consumes<reco::JetTagCollection>  (fDeepCMVABtagNamec);
  fTokDeepCMVABtagNamel   = iC.consumes<reco::JetTagCollection>  (fDeepCMVABtagNamel);
  fTokDeepCMVABtagNamebb  = iC.consumes<reco::JetTagCollection>  (fDeepCMVABtagNamebb);
  if(fComputeFullJetInfo) { 
    fTokPrunedJetName     = iC.consumes<reco::BasicJetCollection>(fPrunedJetName);
    fTokTrimmedJetName    = iC.consumes<reco::BasicJetCollection>(fTrimmedJetName);
    fTokSoftDropJetName   = iC.consumes<reco::BasicJetCollection>(fSoftDropJetName);
    fTokCSVbtagSubJetName = iC.consumes<reco::JetTagCollection>  (fCSVbtagSubJetName);
    fTokCSVDoubleBtagName = iC.consumes<reco::JetTagCollection>  (fCSVDoubleBtagName);
    fTokBoostedDoubleSVTagInfo = iC.consumes<reco::BoostedDoubleSVTagInfoCollection> (fBoostedDoubleSVTagInfoName);
    fTokDeepDoubleBvLtagName = iC.consumes<reco::JetTagCollection>  (fDeepDoubleBvLtagName);
    fTokDeepDoubleBvLNoMassSculptPentagName = iC.consumes<reco::JetTagCollection>  (fDeepDoubleBvLNoMassSculptPentagName);
    fTokDeepDoubleCvLtagName = iC.consumes<reco::JetTagCollection>  (fDeepDoubleCvLtagName);
    fTokDeepDoubleCvLNoMassSculptPentagName = iC.consumes<reco::JetTagCollection>  (fDeepDoubleCvLNoMassSculptPentagName);
    fTokDeepDoubleCvBtagName = iC.consumes<reco::JetTagCollection>  (fDeepDoubleCvBtagName);
    fTokDeepDoubleCvBNoMassSculptPentagName = iC.consumes<reco::JetTagCollection>  (fDeepDoubleCvBNoMassSculptPentagName);
    //fToksoftPFMuonTagInfo     = iC.consumes<reco::CandSoftLeptonTagInfoCollection>     (fsoftPFMuonTagInfoName);
    //fToksoftPFElectronTagInfo = iC.consumes<reco::CandSoftLeptonTagInfoCollection>     (fsoftPFElectronTagInfoName);
    edm::InputTag lTau1(fJettinessName,"tau1");
    edm::InputTag lTau2(fJettinessName,"tau2");
    edm::InputTag lTau3(fJettinessName,"tau3");
    edm::InputTag lTau4(fJettinessName,"tau4");
    edm::InputTag lQGSubJets(fQGLikelihoodSubJets,"qgLikelihood");
    edm::InputTag lSubJets  (fSubJetName,"SubJets");
    edm::InputTag lTopTagSubJet(fTopTaggerName,"caTopSubJets"); //"cmsTopTagPFJetsCHS"
    edm::InputTag lTopTag      (fTopTaggerName);
    fTokTau1Name           = iC.consumes<edm::ValueMap<float> >   (lTau1);
    fTokTau2Name           = iC.consumes<edm::ValueMap<float> >   (lTau2);
    fTokTau3Name           = iC.consumes<edm::ValueMap<float> >   (lTau3);
    fTokTau4Name           = iC.consumes<edm::ValueMap<float> >   (lTau4);
    fTokQGLSubJets         = iC.consumes<edm::ValueMap<float> >   (lQGSubJets);
    fTokSubJets            = iC.consumes<reco::PFJetCollection>   (lSubJets);
    fTokCMSTTJetProduct    = iC.consumes<reco::BasicJetCollection>(lTopTag);
    fTokCMSTTSubJetProduct = iC.consumes<reco::PFJetCollection>   (lTopTagSubJet);
    fECF = new EnergyCorrelations();
    //fRecursiveSoftDrop1 = new fastjet::RecursiveSoftDrop( 1. ,0.05,fConeSize,-1); // beta = 1, zcut=0.05, n = Inf
    //fRecursiveSoftDrop2 = new fastjet::RecursiveSoftDrop( 2. ,0.1,fConeSize,-1); // beta = 2, zcut=0.1, n = Inf
  }
}

//--------------------------------------------------------------------------------------------------
FillerJet::~FillerJet()
{
  delete fJetCorr;
  delete fJetUnc;
  //delete fRecursiveSoftDrop1;
  //delete fRecursiveSoftDrop2;
  //delete fShowerDeco;
  //delete fECF;
  //delete fRand;
}
void FillerJet::initPUJetId() { 
  if(!fUseAOD) return;
  std::cout << "===> Re-initializing" << std::endl;
  fJetPUIDMVACalc.initialize(baconhep::JetPUIDMVACalculator::k53,
			     "BDT",fLowPtWeightFile,
			     "BDT",fHighPtWeightFile);
  initBoostedBtaggingJetId();
  initBReg();
}
void FillerJet::initBReg(){
  fBReg.initialize("BDT",fBRegFile);
}
void FillerJet::initBRegNN(){
  fBRegNN.initialize(fBRegNNFile,fBRegNNMean,fBRegNNStd);
}
void FillerJet::initBoostedBtaggingJetId(){
  fJetBoostedBtaggingMVACalc.initialize("BDT",fWeightFile,fWeightFile.find("Subjet") < std::string::npos);

}
//--------------------------------------------------------------------------------------------------
void FillerJet::initJetCorr(const std::vector<std::string> &jecFiles,
                            const std::vector<std::string> &jecUncFiles)
{
  assert(jecFiles.size()>0);
  assert(jecUncFiles.size()>0);

  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
 
  std::vector<JetCorrectorParameters> corrParams;
  for(unsigned int icorr=0; icorr<jecFiles.size(); icorr++) {
    std::cout << "JEC===> " << (cmssw_base_src + jecFiles[icorr]) << std::endl;
    corrParams.push_back(JetCorrectorParameters(cmssw_base_src + jecFiles[icorr]));
  }
  fJetCorr = new FactorizedJetCorrector(corrParams);
  std::cout << "JEC U===> " << (cmssw_base_src + jecUncFiles[0]) << std::endl;
  JetCorrectorParameters param(cmssw_base_src + jecUncFiles[0]);
  fJetUnc = new JetCorrectionUncertainty(param);
}
Measurement1D FillerJet::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
  VertexDistanceXY dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

Measurement1D FillerJet::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
  VertexDistance3D dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerJet::fill(TClonesArray *array, TClonesArray *iExtraArray,TClonesArray *iSVArray,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, 
		     const reco::Vertex	&pv,
		     int iNPV,
		     const TClonesArray *iPFArr,
		     const std::vector<TriggerRecord> &triggerRecords,
		     const trigger::TriggerEvent *triggerEvent,
                     const pat::TriggerObjectStandAloneCollection *patTriggerObjects)
{
  assert(array);
  assert(!fComputeFullJetInfo || iExtraArray);
  //if(fUseAOD) { assert(fJetPUIDMVACalc.isInitialized()); }
  //assert(fJetBoostedBtaggingMVACalc.isInitialized()); 
  fRand->SetSeed(iEvent.id().event());
 
  // Get jet collection
  edm::Handle<reco::PFJetCollection> hJetProduct;
  iEvent.getByToken(fTokJetName,hJetProduct);
  assert(hJetProduct.isValid());
  const reco::PFJetCollection *jetCol = hJetProduct.product();

  // Get JEC collection
  edm::Handle<reco::JetCorrector> hJetCorrector;
  iEvent.getByToken(fTokJECName,hJetCorrector);
  assert(hJetCorrector.isValid());
  const reco::JetCorrector *jetCorr = hJetCorrector.product();

  // Get JEC Unc collection
  edm::ESHandle<JetCorrectorParametersCollection> hJetCorrectorParsCol;
  iSetup.get<JetCorrectionsRecord>().get(fJECUName,hJetCorrectorParsCol);
  assert(hJetCorrectorParsCol.isValid());
  const JetCorrectorParameters    jetPars = (*hJetCorrectorParsCol)["Uncertainty"];
  delete fJetUnc;
  fJetUnc  = new JetCorrectionUncertainty(jetPars);

  // Get gen jet collection
  //std::cout << fGenJetName << std::endl;
  //std::cout << "fTokGenJetName " << fJetName << fGenJetName << fConeSize << std::endl;

  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJetCol = 0;
  if(fUseGen) { 
    iEvent.getByToken(fTokGenJetName,hGenJetProduct);
    assert(hGenJetProduct.isValid());
    genJetCol = hGenJetProduct.product();
  }
  // Get Jet Flavor Match
  edm::Handle<reco::JetFlavourInfoMatchingCollection> hJetFlavourMatch;
  if(fUseGen) {
    iEvent.getByToken(fTokJetFlavorName, hJetFlavourMatch);
    assert(hJetFlavourMatch.isValid());
  }

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(fTokPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid()); 
 
  // Get b-jet tagger
  edm::Handle<reco::JetTagCollection> hCSVbtags;
  iEvent.getByToken(fTokCSVbtagName, hCSVbtags);
  assert(hCSVbtags.isValid());

  // Get b-jet MVA tagger
  edm::Handle<reco::JetTagCollection> hMVAbtags;
  iEvent.getByToken(fTokMVAbtagName, hMVAbtags);
  assert(hMVAbtags.isValid());

  // Get c-jet vs B tagger
  edm::Handle<reco::JetTagCollection> hCVBctags;
  iEvent.getByToken(fTokCVBctagName, hCVBctags);
  assert(hCVBctags.isValid());

  // Get c-jet vs Light tagger
  edm::Handle<reco::JetTagCollection> hCVLctags;
  iEvent.getByToken(fTokCVLctagName, hCVLctags);
  assert(hCVLctags.isValid());

  //Get Quark Gluon Likelihood
  edm::Handle<edm::ValueMap<float> > hQGLikelihood; 
  iEvent.getByToken(fTokQGLikelihood,hQGLikelihood); 
  assert(hQGLikelihood.isValid());

  edm::Handle<edm::ValueMap<float> > hQGLaxis2;
  iEvent.getByToken(fTokQGLAxis2,hQGLaxis2);
  assert(hQGLaxis2.isValid());

  edm::Handle<edm::ValueMap<float> > hQGLptD;
  iEvent.getByToken(fTokQGLPtD,hQGLptD);
  assert(hQGLptD.isValid());

  edm::Handle<edm::ValueMap<int> > hQGLmult;
  iEvent.getByToken(fTokQGLMult,hQGLmult);
  assert(hQGLmult.isValid());

  // Get DeepCSV tagger
  edm::Handle<reco::JetTagCollection> hDeepCSVBtagb;
  iEvent.getByToken(fTokDeepCSVBtagNameb, hDeepCSVBtagb);
  assert(hDeepCSVBtagb.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCSVBtagc;
  iEvent.getByToken(fTokDeepCSVBtagNamec, hDeepCSVBtagc);
  assert(hDeepCSVBtagc.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCSVBtagl;
  iEvent.getByToken(fTokDeepCSVBtagNamel, hDeepCSVBtagl);
  assert(hDeepCSVBtagl.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCSVBtagbb;
  iEvent.getByToken(fTokDeepCSVBtagNamebb, hDeepCSVBtagbb);
  assert(hDeepCSVBtagbb.isValid());

  // Get DeepCMVA tagger
  edm::Handle<reco::JetTagCollection> hDeepCMVABtagb;
  iEvent.getByToken(fTokDeepCMVABtagNameb, hDeepCMVABtagb);
  assert(hDeepCMVABtagb.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCMVABtagc;
  iEvent.getByToken(fTokDeepCMVABtagNamec, hDeepCMVABtagc);
  assert(hDeepCMVABtagc.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCMVABtagl;
  iEvent.getByToken(fTokDeepCMVABtagNamel, hDeepCMVABtagl);
  assert(hDeepCMVABtagl.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCMVABtagbb;
  iEvent.getByToken(fTokDeepCMVABtagNamebb, hDeepCMVABtagbb);
  assert(hDeepCMVABtagbb.isValid());

  //edm::Handle<std::vector<reco::CandIPTagInfo> > hSVTagInfoCollection;
  //iEvent.getByToken(fTokSVTagInfoCollection, hSVTagInfoCollection);
  //assert(hSVTagInfoCollection.isValid());
  //const std::vector<reco::CandIPTagInfo> *trackIPTagInfos = hSVTagInfoCollection.product();

  VertexDistance3D vdist;
  TClonesArray &rArray      = *array;
  TClonesArray &rExtraArray = *iExtraArray;
  for(reco::PFJetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {
    const double ptRaw = itJet->pt();
    if(ptRaw < fMinPt) continue;
    // input to jet corrections
    double jetcorr = 1;
    if(fabs(itJet->eta()) < 5.191 && fApplyJEC) {
      /* Switching to DB 
      fJetCorr->setJetPt(ptRaw);
      fJetCorr->setJetEta(itJet->eta());
      fJetCorr->setJetPhi(itJet->phi());
      fJetCorr->setJetE(itJet->energy());
      fJetCorr->setRho(*hRho);
      fJetCorr->setJetA(itJet->jetArea());
      fJetCorr->setJetEMF(-99.0);
      jetcorr = fJetCorr->getCorrection();
      */
      jetcorr = jetCorr->correction(*itJet);
    }

    // jet pT cut (BOTH raw AND corrected pT must exceed threshold)
    if(ptRaw*jetcorr < fMinPt || ptRaw < fMinPt) continue;
    fJetUnc->setJetPt ( ptRaw*jetcorr  );
    fJetUnc->setJetEta( itJet->eta() );
    double jetunc = fJetUnc->getUncertainty(true);
    
    //bool passLoose = JetTools::passPFLooseID(*itJet);
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
    // Impact Parameter and leptons
    //==============================
    pJet->d0 = JetTools::jetD0(*itJet, pv);
    pJet->dz = JetTools::jetDz(*itJet, pv);
    pJet->leadPt = JetTools::leadPt(*itJet);
    pJet->lepPt  = JetTools::leptons(*itJet,0);
    pJet->lepDR  = JetTools::leptons(*itJet,2);
    
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
    iEvent.getByToken(fTokSVName, secVertices);
    const reco::VertexCompositePtrCandidateCollection svtx=*secVertices;
    float vtxPt=0,vtxMass=0,vtx3DVal=0,vtx3DeVal=0,vtxNTrks=0,maxFoundSignificance=0;
    for (const reco::VertexCompositePtrCandidate &sv : svtx) {
      GlobalVector flightDir(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(),sv.vertex().z() - pv.z());
      GlobalVector jetDir(itJet->px(),itJet->py(),itJet->pz());
      if( Geom::deltaR2( flightDir, jetDir ) < 0.09 ){
	Measurement1D dl= vdist.distance(pv,VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
	if(dl.significance() > maxFoundSignificance){
	  maxFoundSignificance=dl.significance();
	  vtxPt    = sv.pt();
	  vtxMass  = sv.p4().M();
	  vtx3DVal = dl.value();
	  vtx3DeVal = dl.error();
	  vtxNTrks  = sv.numberOfSourceCandidatePtrs();
	}
      }
    }
    pJet->vtx3DSig = maxFoundSignificance;
    pJet->vtxMass  = vtxMass;
    pJet->vtxPt    = vtxPt;
    pJet->vtxNtk   = vtxNTrks;
    pJet->vtx3DVal = vtx3DVal;
    pJet->ptreg    = fBReg.mvaValue(iNPV,jetcorr,*itJet,vtxPt,vtxMass,vtx3DVal,vtxNTrks,vtx3DeVal);
      
    //
    // Bjet NN Regression
    //
    fBRegNN.Jet_pt = ptRaw;
    fBRegNN.Jet_eta = itJet->eta();
    fBRegNN.rho = *hRho;
    fBRegNN.Jet_mt = sqrt(itJet->energy()*itJet->energy()-itJet->pz()*itJet->pz());
    fBRegNN.Jet_leadTrackPt = JetTools::leadTrkPt(*itJet);
    fBRegNN.Jet_leptonDeltaR = JetTools::leptons(*itJet,7);
    if (fBRegNN.Jet_leptonDeltaR<1e-4) {
      fBRegNN.Jet_leptonPtRel = 0.;
      fBRegNN.Jet_leptonPtRelInv = 0.;
    }
    else {
      fBRegNN.Jet_leptonPtRel = JetTools::leptons(*itJet,8);
      fBRegNN.Jet_leptonPtRelInv = JetTools::leptons(*itJet,9);
    }
    fBRegNN.Jet_neHEF = itJet->neutralHadronEnergyFraction();
    fBRegNN.Jet_neEmEF = itJet->neutralEmEnergyFraction();
    fBRegNN.Jet_vtxPt = vtxPt;
    fBRegNN.Jet_vtxMass = vtxMass;
    fBRegNN.Jet_vtx3dL = vtx3DVal;
    fBRegNN.Jet_vtxNtrk = vtxNTrks;
    fBRegNN.Jet_vtx3deL = vtx3DeVal;
    fBRegNN.Jet_numDaughters_pt03 = JetTools::nDaughters(*itJet,0.03);
    std::vector<float> chvec, emvec, nevec, muvec;
    JetTools::energyRings(*itJet,chvec,emvec,nevec,muvec);
    fBRegNN.Jet_energyRing_dR0_em_Jet_e = emvec[0]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_em_Jet_e = emvec[1]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_em_Jet_e = emvec[2]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_em_Jet_e = emvec[3]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_em_Jet_e = emvec[4]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR0_neut_Jet_e = nevec[0]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_neut_Jet_e = nevec[1]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_neut_Jet_e = nevec[2]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_neut_Jet_e = nevec[3]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_neut_Jet_e = nevec[4]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR0_ch_Jet_e = chvec[0]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_ch_Jet_e = chvec[1]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_ch_Jet_e = chvec[2]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_ch_Jet_e = chvec[3]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_ch_Jet_e = chvec[4]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR0_mu_Jet_e = muvec[0]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_mu_Jet_e = muvec[1]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_mu_Jet_e = muvec[2]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_mu_Jet_e = muvec[3]/(itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_mu_Jet_e = muvec[4]/(itJet->energy()) ;
    fBRegNN.Jet_chHEF = itJet->chargedHadronEnergyFraction();
    fBRegNN.Jet_chEmEF = itJet->chargedEmEnergyFraction();
    fBRegNN.isEle = 0;
    fBRegNN.isMu = 0;
    fBRegNN.isOther = 1;
    int softLepID = abs(JetTools::leptons(*itJet,4));
    if (softLepID==13) fBRegNN.isMu = 1;
    else if (softLepID==11) fBRegNN.isEle = 1;
    fBRegNN.Jet_mass = itJet->mass();
    fBRegNN.Jet_withPtd = JetTools::ptD(*itJet);


    fBRegNN.SetNNVectorVar();
    std::pair<float,float> bjetnnout = fBRegNN.EvaluateNN();
    //std::cout<<"NN BJet Correction = "<<bjetnnout.first<<std::endl;
    pJet->bjetcorr = bjetnnout.first;
    pJet->bjetres  = bjetnnout.second;

    //
    // Identification
    //==============================
    reco::PFJetRef jetRef(hJetProduct, itJet - jetCol->begin());
    reco::JetBaseRef jetBaseRef(jetRef);
    //unsigned int idx = itJet - jetCol->begin();
    //edm::RefToBase<reco::Jet> jetRef_ = jetCol->refAt(idx);  
    pJet->csv  = (*(hCSVbtags.product()))[jetBaseRef];
    pJet->bmva = (*(hMVAbtags.product()))[jetBaseRef];
    pJet->cvb  = (*(hCVBctags.product()))[jetBaseRef];
    pJet->cvl  = (*(hCVLctags.product()))[jetBaseRef];
    pJet->qgid  = (*(hQGLikelihood.product()))[jetBaseRef];
    pJet->axis2 = (*(hQGLaxis2.product()))[jetBaseRef];
    pJet->ptD   = (*(hQGLptD.product()))[jetBaseRef];
    pJet->mult  = (*(hQGLmult.product()))[jetBaseRef];
    pJet->q = JetTools::jetCharge(*itJet);
    
    pJet->deepcsvb = (*(hDeepCSVBtagb.product()))[jetBaseRef];
    pJet->deepcsvc = (*(hDeepCSVBtagc.product()))[jetBaseRef];
    pJet->deepcsvl = (*(hDeepCSVBtagl.product()))[jetBaseRef];
    pJet->deepcsvbb = (*(hDeepCSVBtagbb.product()))[jetBaseRef];

    pJet->deepcmvab = (*(hDeepCMVABtagb.product()))[jetBaseRef];
    pJet->deepcmvac = (*(hDeepCMVABtagc.product()))[jetBaseRef];
    pJet->deepcmval = (*(hDeepCMVABtagl.product()))[jetBaseRef];
    pJet->deepcmvabb = (*(hDeepCMVABtagbb.product()))[jetBaseRef];

    pJet->beta     = JetTools::beta(*itJet, pv);
    pJet->betaStar = JetTools::betaStar(*itJet, pv, pvCol);
    pJet->dR2Mean  = JetTools::dR2Mean(*itJet);
    pJet->mva = -2;
    /*
    if(passLoose) {
      double dRMean = JetTools::dRMean(*itJet);
      double frac01 = JetTools::frac(*itJet,0.1);
      double frac02 = JetTools::frac(*itJet,0.2);
      double frac03 = JetTools::frac(*itJet,0.3);
      double frac04 = JetTools::frac(*itJet,0.4);
      double frac05 = JetTools::frac(*itJet,0.5);
      pJet->mva = fJetPUIDMVACalc.mvaValue((float)pvCol->size(), ptRaw*jetcorr, itJet->eta(), itJet->phi(),
                                           pJet->d0, pJet->dz, pJet->beta, pJet->betaStar, itJet->chargedMultiplicity(), itJet->neutralMultiplicity(),
					   dRMean, pJet->dR2Mean, pJet->ptD, frac01, frac02, frac03, frac04, frac05);
    }
    */
    TVector2 lPull = JetTools::jetPull(*itJet,0);
    pJet->pullY      = lPull.X();
    pJet->pullPhi    = lPull.Y();
    TVector2 lChPull = JetTools::jetPull(*itJet,1);
    pJet->chPullY    = lChPull.X();
    pJet->chPullPhi  = lChPull.Y();
    TVector2 lNeuPull = JetTools::jetPull(*itJet,2);
    pJet->neuPullY   = lNeuPull.X();
    pJet->neuPullPhi = lNeuPull.Y();

    // Basic Noise Variables
    pJet->chEmFrac   = itJet->chargedEmEnergy() / itJet->energy();
    pJet->neuEmFrac  = itJet->neutralEmEnergy() / itJet->energy();
    pJet->chHadFrac  = itJet->chargedHadronEnergy() / itJet->energy();
    pJet->neuHadFrac = itJet->neutralHadronEnergy() / itJet->energy();
    pJet->muonFrac   = itJet->muonEnergy() / itJet->energy();
        
    //
    // Generator matching
    //==============================
    const reco::GenJet * matchGenJet = 0; 
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    int lNB =0,lNC = 0;
    if(matchGenJet != 0) { 
      pJet->partonFlavor = (*hJetFlavourMatch)[jetBaseRef].getPartonFlavour();
      pJet->hadronFlavor = (*hJetFlavourMatch)[jetBaseRef].getHadronFlavour();
      pJet->genpt        = matchGenJet->pt();
      pJet->geneta       = matchGenJet->eta();
      pJet->genphi       = matchGenJet->phi();
      pJet->genm         = matchGenJet->mass();
      //float pMsd = 0; float pe2sdb1 =0; float pe3_v2_sdb1 =0;
      //if(pJet->genpt > 100) softdrop(&(*matchGenJet),pMsd,pe2sdb1,pe3_v2_sdb1);
      //pJet->genmsd       = pMsd;
      //pJet->gene2sdb1    = pe2sdb1;
      //pJet->gene3_v2_sdb1= pe3_v2_sdb1;
      lNB = (*hJetFlavourMatch)[jetBaseRef].getbHadrons().size();
      lNC = (*hJetFlavourMatch)[jetBaseRef].getcHadrons().size();
    }
    pJet->vtxFlavInfo = lNB*1000 + lNC*100;
    pJet->nCharged   = itJet->chargedMultiplicity();
    pJet->nNeutrals  = itJet->neutralMultiplicity();
    pJet->nParticles = itJet->nConstituents ();
    if(fUseTO) { 
      if(triggerEvent      != 0) {pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, *triggerEvent); } 
      else                       {pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, *patTriggerObjects); }
    }
     if(fAddPFCand) { 
      pJet->pfCands.clear();
      std::vector<reco::CandidatePtr> pfConstituents = itJet->getJetConstituents();                                                                                                   
      for(unsigned int i0 = 0; i0 < pfConstituents.size(); i0++) { 
	reco::CandidatePtr pfcand = pfConstituents[i0]; 
	for(       int i1 = 0; i1 < iPFArr->GetEntriesFast();  i1++) { 
	  baconhep::TPFPart *pPF = (baconhep::TPFPart*)(*iPFArr)[i1];
	  if(pfcand->pdgId() != pPF->pfType)  continue;
	  if(fabs(pfcand->pt() - pPF->pt)  > 0.1)  continue;
	  if(fabs(pfcand->eta()- pPF->eta) > 0.01) continue;
	  if(fabs(reco::deltaPhi(pfcand->phi(),pPF->phi)) > 0.01) continue;
	  pJet->pfCands.push_back(i1);
	}
      }
    }
    ////Add Extras
    baconhep::TAddJet *pAddJet = 0; 
    if(fComputeFullJetInfo) {
      assert(rExtraArray.GetEntries() < rExtraArray.GetSize());
      const int extraIndex = rExtraArray.GetEntries();
      new(rExtraArray[extraIndex]) baconhep::TAddJet();
      pAddJet = (baconhep::TAddJet*)rExtraArray[extraIndex];
      pAddJet->index = index;
      addJet(pAddJet, iSVArray, pv, iEvent, *itJet, jetBaseRef);
      
    }
  } 
}

// === filler for MINIAOD ===
void FillerJet::fill(TClonesArray *array, TClonesArray *iExtraArray,TClonesArray *iSVArray,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup,
                     const reco::Vertex &pv,
		     int iNPV,
		     const TClonesArray *iPFArr,
                     const std::vector<TriggerRecord> &triggerRecords,
                     const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);
  assert(!fComputeFullJetInfo || iExtraArray);
  fRand->SetSeed(iEvent.id().event());

  // Get jet collection
  edm::Handle<pat::JetCollection> hJetProduct;
  iEvent.getByToken(fTokPatJetName,hJetProduct);
  assert(hJetProduct.isValid());
  const pat::JetCollection *jetCol = hJetProduct.product();

  // Get JEC collection
  edm::Handle<reco::JetCorrector> hJetCorrector;
  iEvent.getByToken(fTokJECName,hJetCorrector);
  assert(hJetCorrector.isValid());
  const reco::JetCorrector *jetCorr = hJetCorrector.product();

  // Get JEC Unc collection
  edm::ESHandle<JetCorrectorParametersCollection> hJetCorrectorParsCol;
  iSetup.get<JetCorrectionsRecord>().get(fJECUName,hJetCorrectorParsCol);
  assert(hJetCorrectorParsCol.isValid());
  const JetCorrectorParameters    jetPars = (*hJetCorrectorParsCol)["Uncertainty"];
  delete fJetUnc;
  fJetUnc  = new JetCorrectionUncertainty(jetPars);

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(fTokPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();
  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid()); 

  // Get DeepCSV tagger
  /*
  edm::Handle<reco::JetTagCollection> hDeepCMVABtagb;
  iEvent.getByToken(fTokDeepCMVABtagNameb, hDeepCMVABtagb);
  assert(hDeepCMVABtagb.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCMVABtagc;
  iEvent.getByToken(fTokDeepCSVBtagNamec, hDeepCSVBtagc);
  assert(hDeepCSVBtagc.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCMVABtagl;
  iEvent.getByToken(fTokDeepCSVBtagNamel, hDeepCSVBtagl);
  assert(hDeepCSVBtagl.isValid());

  edm::Handle<reco::JetTagCollection> hDeepCMVABtagbb;
  iEvent.getByToken(fTokDeepCSVBtagNamebb, hDeepCSVBtagbb);
  assert(hDeepCSVBtagbb.isValid());
  */
  // Get Quark Gluon Likelihood
  edm::Handle<edm::ValueMap<float> > hQGLikelihood;
  edm::Handle<edm::ValueMap<float> > hQGLaxis2;
  edm::Handle<edm::ValueMap<float> > hQGLptD;
  edm::Handle<edm::ValueMap<int> >   hQGLmult;
  if(fQGLikelihood.size() > 0) { 
    iEvent.getByToken(fTokQGLikelihood,hQGLikelihood);
    assert(hQGLikelihood.isValid());

    iEvent.getByToken(fTokQGLAxis2,hQGLaxis2);
    assert(hQGLaxis2.isValid());
    
    iEvent.getByToken(fTokQGLPtD,hQGLptD);
    assert(hQGLptD.isValid());
    
    iEvent.getByToken(fTokQGLMult,hQGLmult);
    assert(hQGLmult.isValid());
  }
  // Get gen jet collection
  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJetCol = 0;
  //std::cout << "fTokGenJetName " << fJetName << fGenJetName << std::endl;

  if(fUseGen) {
    iEvent.getByToken(fTokGenJetName,hGenJetProduct);
    assert(hGenJetProduct.isValid());
    genJetCol = hGenJetProduct.product();
  }
  TClonesArray &rArray      = *array;
  TClonesArray &rExtraArray = *iExtraArray;
  VertexDistance3D vdist;
  for(pat::JetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {

    float uncorr = itJet->jecFactor("Uncorrected");
    double ptRaw = itJet->pt()*itJet->jecFactor("Uncorrected");
    if(ptRaw < fMinPt) continue;
    // input to jet corrections
    double jetcorr = 1;
    if(fabs(itJet->eta()) < 5.191 && fApplyJEC) {
      /* Switching to DB 
      fJetCorr->setJetPt(ptRaw);
      fJetCorr->setJetEta(itJet->eta());
      fJetCorr->setJetPhi(itJet->phi());
      fJetCorr->setJetE(itJet->energy());
      fJetCorr->setRho(*hRho);
      fJetCorr->setJetA(itJet->jetArea());
      fJetCorr->setJetEMF(-99.0);
      jetcorr = fJetCorr->getCorrection();
      */
      jetcorr = jetCorr->correction(*itJet);
    }	
    // jet pT cut (BOTH raw AND corrected pT must exceed threshold)
    if(ptRaw*jetcorr < fMinPt || ptRaw < fMinPt) continue;
    fJetUnc->setJetPt ( ptRaw*jetcorr  );
    fJetUnc->setJetEta( itJet->eta() );
    double jetunc = fJetUnc->getUncertainty(true);

    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TJet();
    baconhep::TJet *pJet = (baconhep::TJet*)rArray[index];

    //
    // Kinematics
    //==============================
    pJet->pt    = itJet->pt();
    if(fApplyJEC) pJet->pt    = ptRaw*jetcorr;
    pJet->eta   = itJet->eta();
    pJet->phi   = itJet->phi();
    pJet->mass  = itJet->mass();
    //if(fApplyJEC) pJet->mass    = itJet->mass()*jetcorr;
    pJet->ptRaw = ptRaw;
    pJet->area  = itJet->jetArea();
    pJet->unc   = jetunc;

    //
    // Impact Parameter
    //==============================
    pJet->d0 = JetTools::jetD0(*itJet, pv);
    pJet->dz = JetTools::jetDz(*itJet, pv);
    pJet->leadPt = JetTools::leadPt(*itJet);
    pJet->lepPt  = JetTools::leptons(*itJet,0);
    pJet->lepDR  = JetTools::leptons(*itJet,2);
    
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
    iEvent.getByToken(fTokSVName, secVertices);
    const reco::VertexCompositePtrCandidateCollection svtx=*secVertices;
    float vtxPt=0,vtxMass=0,vtx3DVal=0,vtx3DeVal=0,vtxNTrks=0,maxFoundSignificance=0;
    for (const reco::VertexCompositePtrCandidate &sv : svtx) {
      GlobalVector flightDir(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(),sv.vertex().z() - pv.z());
      GlobalVector jetDir(itJet->px(),itJet->py(),itJet->pz());
      if( Geom::deltaR2( flightDir, jetDir ) < 0.09 ){
	Measurement1D dl= vdist.distance(pv,VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
	if(dl.significance() > maxFoundSignificance){
	  maxFoundSignificance=dl.significance();
	  vtxPt    = sv.pt();
	  vtxMass  = sv.p4().M();
	  vtx3DVal = dl.value();
	  vtx3DeVal = dl.error();
	  vtxNTrks  = sv.numberOfSourceCandidatePtrs();
	}
      }
    }
    pJet->vtx3DSig = maxFoundSignificance;
    pJet->vtxMass  = vtxMass;
    pJet->vtxPt    = vtxPt;
    pJet->vtxNtk   = vtxNTrks;
    pJet->vtx3DVal = vtx3DVal;
    pJet->ptreg    = fBReg.mvaValue(iNPV,jetcorr,*itJet,vtxPt,vtxMass,vtx3DVal,vtxNTrks,vtx3DeVal);
      
    //
    // Bjet NN Regression
    //
    fBRegNN.Jet_pt = ptRaw;
    fBRegNN.Jet_eta = itJet->eta();
    fBRegNN.rho = *hRho;
    fBRegNN.Jet_mt = uncorr*sqrt(itJet->energy()*itJet->energy()-itJet->pz()*itJet->pz());
    fBRegNN.Jet_leadTrackPt = JetTools::leadTrkPt(*itJet);
    fBRegNN.Jet_leptonDeltaR = JetTools::leptons(*itJet,7);
    if (fBRegNN.Jet_leptonDeltaR<1e-4) {
      fBRegNN.Jet_leptonPtRel = 0.;
      fBRegNN.Jet_leptonPtRelInv = 0.;
    }
    else {
      fBRegNN.Jet_leptonPtRel = JetTools::leptons(*itJet,8);
      fBRegNN.Jet_leptonPtRelInv = uncorr*JetTools::leptons(*itJet,9);
    }
    fBRegNN.Jet_neHEF = itJet->neutralHadronEnergyFraction();
    fBRegNN.Jet_neEmEF = itJet->neutralEmEnergyFraction();
    fBRegNN.Jet_vtxPt = vtxPt;
    fBRegNN.Jet_vtxMass = vtxMass;
    fBRegNN.Jet_vtx3dL = vtx3DVal;
    fBRegNN.Jet_vtxNtrk = vtxNTrks;
    fBRegNN.Jet_vtx3deL = vtx3DeVal;
    fBRegNN.Jet_numDaughters_pt03 = JetTools::nDaughters(*itJet,0.03);
    std::vector<float> chvec, emvec, nevec, muvec;
    JetTools::energyRings(*itJet,chvec,emvec,nevec,muvec);
    fBRegNN.Jet_energyRing_dR0_em_Jet_e = emvec[0]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_em_Jet_e = emvec[1]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_em_Jet_e = emvec[2]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_em_Jet_e = emvec[3]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_em_Jet_e = emvec[4]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR0_neut_Jet_e = nevec[0]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_neut_Jet_e = nevec[1]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_neut_Jet_e = nevec[2]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_neut_Jet_e = nevec[3]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_neut_Jet_e = nevec[4]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR0_ch_Jet_e = chvec[0]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_ch_Jet_e = chvec[1]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_ch_Jet_e = chvec[2]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_ch_Jet_e = chvec[3]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_ch_Jet_e = chvec[4]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR0_mu_Jet_e = muvec[0]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR1_mu_Jet_e = muvec[1]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR2_mu_Jet_e = muvec[2]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR3_mu_Jet_e = muvec[3]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_energyRing_dR4_mu_Jet_e = muvec[4]/(uncorr*itJet->energy()) ;
    fBRegNN.Jet_chHEF = itJet->chargedHadronEnergyFraction();
    fBRegNN.Jet_chEmEF = itJet->chargedEmEnergyFraction();
    fBRegNN.isEle = 0;
    fBRegNN.isMu = 0;
    fBRegNN.isOther = 1;
    int softLepID = abs(JetTools::leptons(*itJet,4));
    if (softLepID==13) fBRegNN.isMu = 1;
    else if (softLepID==11) fBRegNN.isEle = 1;
    fBRegNN.Jet_mass = itJet->mass()*uncorr;
    fBRegNN.Jet_withPtd = JetTools::ptD(*itJet);


    fBRegNN.SetNNVectorVar();
    std::pair<float,float> bjetnnout = fBRegNN.EvaluateNN();
    //std::cout<<"NN BJet Correction = "<<bjetnnout.first<<std::endl;
    pJet->bjetcorr = bjetnnout.first;
    pJet->bjetres  = bjetnnout.second;
    //
    // Identification
    //==============================
    pJet->csv  = itJet->bDiscriminator(fCSVbtagName);
    pJet->bmva = itJet->bDiscriminator(fMVAbtagName);
    pJet->cvb  = itJet->bDiscriminator(fCVBctagName);
    pJet->cvl  = itJet->bDiscriminator(fCVLctagName);

    // Deep-CSV
    pJet->deepcsvb = itJet->bDiscriminator((fDeepCSVBtagName+":probb").c_str());
    pJet->deepcsvc = itJet->bDiscriminator((fDeepCSVBtagName+":probc").c_str());
    pJet->deepcsvl = itJet->bDiscriminator((fDeepCSVBtagName+":probudsg").c_str());
    pJet->deepcsvbb = itJet->bDiscriminator((fDeepCSVBtagName+":probbb").c_str());

    pJet->deepcmvab = itJet->bDiscriminator((fDeepCMVABtagName+":probb").c_str());
    pJet->deepcmvac = itJet->bDiscriminator((fDeepCMVABtagName+":probc").c_str());
    pJet->deepcmval = itJet->bDiscriminator((fDeepCMVABtagName+":probudsg").c_str());
    pJet->deepcmvabb = itJet->bDiscriminator((fDeepCMVABtagName+":probbb").c_str());
    edm::RefToBase<pat::Jet> jetBaseRef( edm::Ref<pat::JetCollection>(hJetProduct, itJet - jetCol->begin()) );
    //std::cout << "===> deep cmva" << pJet->deepcmvab << " -- " << pJet->deepcmvac << " -- " <<  pJet->deepcmval << " -- " << pJet->deepcmvabb << " -- " << pJet->deepcsvb  << std::endl;
    if(fQGLikelihood.size() > 0) { 
      pJet->qgid  = (*hQGLikelihood)[jetBaseRef];
      pJet->axis2 = (*hQGLaxis2)[jetBaseRef];
      pJet->ptD   = (*hQGLptD)[jetBaseRef];
      pJet->mult  = (*hQGLmult)[jetBaseRef];
    }
    pJet->q = JetTools::jetCharge(*itJet);

    pJet->beta     = JetTools::beta(*itJet, pv);
    pJet->betaStar = JetTools::betaStar(*itJet, pv, pvCol);
    pJet->dR2Mean  = JetTools::dR2Mean(*itJet);
    pJet->mva = itJet->userFloat("pileupJetId:fullDiscriminant");

    TVector2 lPull    = JetTools::jetPull(*itJet,0);
    pJet->pullY       = lPull.X();
    pJet->pullPhi     = lPull.Y();
    TVector2 lChPull  = JetTools::jetPull(*itJet,1);
    pJet->chPullY     = lChPull.X();
    pJet->chPullPhi   = lChPull.Y();
    TVector2 lNeuPull = JetTools::jetPull(*itJet,2);
    pJet->neuPullY    = lNeuPull.X();
    pJet->neuPullPhi  = lNeuPull.Y();

    // Basic Noise Variables
    pJet->chEmFrac   = itJet->chargedEmEnergyFraction();
    pJet->neuEmFrac  = itJet->neutralEmEnergyFraction();
    pJet->chHadFrac  = itJet->chargedHadronEnergyFraction();
    pJet->neuHadFrac = itJet->neutralHadronEnergyFraction();
    pJet->muonFrac   = itJet->muonEnergyFraction();

    //
    // Generator matching
    //==============================
    const reco::GenJet * matchGenJet = 0;
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    if(matchGenJet != 0) {
      pJet->partonFlavor = itJet->partonFlavour();
      pJet->hadronFlavor = itJet->hadronFlavour();
      pJet->genpt        = matchGenJet->pt();
      pJet->geneta       = matchGenJet->eta();
      pJet->genphi       = matchGenJet->phi();
      pJet->genm         = matchGenJet->mass();
      //float pMsd = 0; float pe2sdb1 =0; float pe3_v2_sdb1 =0;
      //if(pJet->genpt > 10) softdrop(&(*matchGenJet),pMsd,pe2sdb1,pe3_v2_sdb1);
      //pJet->genmsd       = pMsd;
      //pJet->gene2sdb1    = pe2sdb1;
      //pJet->gene3_v2_sdb1= pe3_v2_sdb1;
    }
    int lNB = itJet->jetFlavourInfo().getbHadrons().size();
    int lNC = itJet->jetFlavourInfo().getcHadrons().size();
    pJet->vtxFlavInfo = lNB*1000 + lNC*100;
    pJet->nCharged   = itJet->chargedMultiplicity();
    pJet->nNeutrals  = itJet->neutralMultiplicity();
    pJet->nParticles = itJet->numberOfDaughters();
    if(fUseTO) pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, triggerObjects);

    if(fAddPFCand) { 
      pJet->pfCands.clear();
      std::vector<reco::CandidatePtr> pfConstituents = itJet->getJetConstituents();
      for(unsigned int i0 = 0; i0 < pfConstituents.size(); i0++) { 
	reco::CandidatePtr pfcand = pfConstituents[i0];    
	for(int i1 = 0; i1 < iPFArr->GetEntries(); i1++) { 
	  baconhep::TPFPart *pPF = (baconhep::TPFPart*)(*iPFArr)[i1];
	  if(fabs(pfcand->pt() - pPF->pt)  > 0.1) continue;
	  if(fabs(pfcand->eta()- pPF->eta) > 0.01) continue;
	  if(reco::deltaPhi(pfcand->phi(),pPF->phi) > 0.01) continue;
	  if(pfcand->pdgId() != pPF->pfType) continue;
	  pJet->pfCands.push_back(i1);
	}
      }
    }
    ////Add Extras
    baconhep::TAddJet *pAddJet = 0;
    if(fComputeFullJetInfo) {
      assert(rExtraArray.GetEntries() < rExtraArray.GetSize());
      const int extraIndex = rExtraArray.GetEntries();
      new(rExtraArray[extraIndex]) baconhep::TAddJet();
      pAddJet = (baconhep::TAddJet*)rExtraArray[extraIndex];
      pAddJet->index = index;
      addJet(pAddJet, iSVArray, pv, iEvent, *itJet);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void FillerJet::addJet(baconhep::TAddJet *pAddJet, TClonesArray *iSVArr,const reco::Vertex &pv, const edm::Event &iEvent, 
                       const reco::PFJet &itJet, const reco::JetBaseRef &jetBaseRef)
{ 
  // Get pruned jet collection
  edm::Handle<reco::BasicJetCollection> hPrunedJetProduct;
  iEvent.getByToken(fTokPrunedJetName,hPrunedJetProduct);
  assert(hPrunedJetProduct.isValid());
  const reco::BasicJetCollection *prunedJetCol = hPrunedJetProduct.product();

  // Get trimmed jet collection
  edm::Handle<reco::BasicJetCollection> hTrimmedJetProduct;
  iEvent.getByToken(fTokTrimmedJetName,hTrimmedJetProduct);
  assert(hTrimmedJetProduct.isValid());
  const reco::BasicJetCollection *trimmedJetCol = hTrimmedJetProduct.product();

  // Get soft drop jet collection
  edm::Handle<reco::BasicJetCollection> hSoftDropJetProduct;
  iEvent.getByToken(fTokSoftDropJetName,hSoftDropJetProduct);
  assert(hSoftDropJetProduct.isValid());
  const reco::BasicJetCollection *softdropJetCol = hSoftDropJetProduct.product();

  // Get sub-jet collection
  edm::Handle<reco::PFJetCollection> hSubJetProduct;
  iEvent.getByToken(fTokSubJets,hSubJetProduct);
  assert(hSubJetProduct.isValid());
  const reco::PFJetCollection *subJetCol = hSubJetProduct.product();

  // Get b sub-jets 
  edm::Handle<reco::JetTagCollection> hCSVbtagsSubJets;
  iEvent.getByToken(fTokCSVbtagSubJetName, hCSVbtagsSubJets);
  assert(hCSVbtagsSubJets.isValid());

  // Get double b tag
  edm::Handle<reco::JetTagCollection> hCSVDoubleBtag;
  iEvent.getByToken(fTokCSVDoubleBtagName, hCSVDoubleBtag);
  assert(hCSVDoubleBtag.isValid());

  // Get deep double b tag
  edm::Handle<reco::JetTagCollection> hDeepDoubleBvLtag;
  iEvent.getByToken(fTokDeepDoubleBvLtagName, hDeepDoubleBvLtag);
  assert(hDeepDoubleBvLtag.isValid());
  edm::Handle<reco::JetTagCollection> hDeepDoubleBvLNoMassSculptPentag;
  iEvent.getByToken(fTokDeepDoubleBvLNoMassSculptPentagName, hDeepDoubleBvLNoMassSculptPentag);
  assert(hDeepDoubleBvLNoMassSculptPentag.isValid());

  // Get deep double c tag
  edm::Handle<reco::JetTagCollection> hDeepDoubleCvLtag;
  iEvent.getByToken(fTokDeepDoubleCvLtagName, hDeepDoubleCvLtag);
  assert(hDeepDoubleCvLtag.isValid());
  edm::Handle<reco::JetTagCollection> hDeepDoubleCvLNoMassSculptPentag;
  iEvent.getByToken(fTokDeepDoubleCvLNoMassSculptPentagName, hDeepDoubleCvLNoMassSculptPentag);
  assert(hDeepDoubleCvLNoMassSculptPentag.isValid());

  // Get deep double cvb tag
  edm::Handle<reco::JetTagCollection> hDeepDoubleCvBtag;
  iEvent.getByToken(fTokDeepDoubleCvBtagName, hDeepDoubleCvBtag);
  assert(hDeepDoubleCvBtag.isValid());
  edm::Handle<reco::JetTagCollection> hDeepDoubleCvBNoMassSculptPentag;
  iEvent.getByToken(fTokDeepDoubleCvBNoMassSculptPentagName, hDeepDoubleCvBNoMassSculptPentag);
  assert(hDeepDoubleCvBNoMassSculptPentag.isValid());

  //Get Quark Gluon Likelihood on subjets
  edm::Handle<edm::ValueMap<float> > hQGLikelihoodSubJets;
  iEvent.getByToken(fTokQGLSubJets,hQGLikelihoodSubJets);
  assert(hQGLikelihoodSubJets.isValid());
  
  // Get N-subjettiness moments
  edm::Handle<edm::ValueMap<float> > hTau1;
  iEvent.getByToken(fTokTau1Name,hTau1);                                                                                                                                                      
  assert(hTau1.isValid());
  edm::Handle<edm::ValueMap<float> > hTau2;
  iEvent.getByToken(fTokTau2Name,hTau2);
  assert(hTau2.isValid());
  edm::Handle<edm::ValueMap<float> > hTau3;
  iEvent.getByToken(fTokTau3Name,hTau3); 
  assert(hTau3.isValid());
  edm::Handle<edm::ValueMap<float> > hTau4;
  iEvent.getByToken(fTokTau4Name,hTau4); 
  assert(hTau4.isValid());

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid());

  pAddJet->pullAngle = JetTools::jetPullAngle(itJet,hSubJetProduct,fConeSize);
  pAddJet->tau1 = (*(hTau1.product()))[jetBaseRef];
  pAddJet->tau2 = (*(hTau2.product()))[jetBaseRef];
  pAddJet->tau3 = (*(hTau3.product()))[jetBaseRef];
  pAddJet->tau4 = (*(hTau4.product()))[jetBaseRef];
  pAddJet->doublecsv = (*(hCSVDoubleBtag.product()))[jetBaseRef];
  pAddJet->deepdoubleb = (*(hDeepDoubleBvLtag.product()))[jetBaseRef];
  pAddJet->deepdoubleb_nomasssculptpen = (*(hDeepDoubleBvLNoMassSculptPentag.product()))[jetBaseRef];
  pAddJet->deepdoublec = (*(hDeepDoubleCvLtag.product()))[jetBaseRef];
  pAddJet->deepdoublec_nomasssculptpen = (*(hDeepDoubleCvLNoMassSculptPentag.product()))[jetBaseRef];
  pAddJet->deepdoublecvb = (*(hDeepDoubleCvBtag.product()))[jetBaseRef];
  pAddJet->deepdoublecvb_nomasssculptpen = (*(hDeepDoubleCvBNoMassSculptPentag.product()))[jetBaseRef];

  //if(fShowerDeco != 0) { 
  std::vector<reco::CandidatePtr> pfConstituents = itJet.getJetConstituents();                                                                                                                     
  std::vector<fastjet::PseudoJet>   lClusterParticles;                                                                                                                                     
  for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {                                                                                                                                         
    reco::CandidatePtr pfcand = pfConstituents[ic];                                                                                                                                  
    fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());                                                                                                      
    lClusterParticles.emplace_back(pPart);                                                                                                                                                       
  }                                                                           
  std::sort(lClusterParticles.begin(),lClusterParticles.end(),JetTools::orderPseudoJet);
  //if(fShowerDeco != 0) pAddJet->topchi2 = fShowerDeco->chi(itJet.pt(),lClusterParticles);
  fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 2.0);
  fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
  std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
  /*
  fastjet::contrib::EnergyCorrelatorDoubleRatio C2beta0 (2,0. ,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio C2beta02(2,0.2,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio C2beta05(2,0.5,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio C2beta10(2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio C2beta20(2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
  pAddJet->c2_0   = C2beta0 (inclusive_jets[0]);
  pAddJet->c2_0P2 = C2beta02(inclusive_jets[0]);
  pAddJet->c2_0P5 = C2beta05(inclusive_jets[0]);
  pAddJet->c2_1P0 = C2beta10(inclusive_jets[0]);
  pAddJet->c2_2P0 = C2beta20(inclusive_jets[0]);
  */
  //fastjet::contrib::EnergyCorrelator::Measure measure = fastjet::contrib::EnergyCorrelator::pt_R;
  double beta=1;
  int nFilter = TMath::Min(100,(int)lClusterParticles.size());
  std::vector<fastjet::PseudoJet> lFiltered(lClusterParticles.begin(),lClusterParticles.begin()+nFilter);
  fECF->calcECFN(beta,lFiltered,true);
  pAddJet->e2_b1      = float(fECF->manager->ecfns["2_2"]);
  pAddJet->e3_b1      = float(fECF->manager->ecfns["3_3"]);
  pAddJet->e3_v1_b1   = float(fECF->manager->ecfns["3_1"]);
  pAddJet->e3_v2_b1   = float(fECF->manager->ecfns["3_2"]);
  pAddJet->e4_v1_b1   = float(fECF->manager->ecfns["4_1"]);
  pAddJet->e4_v2_b1   = float(fECF->manager->ecfns["4_2"]);
  /*
  beta=2;
  fECF->calcECFN(beta,lFiltered);
  pAddJet->e2_b2      = float(fECF->manager->ecfns["2_2"]);
  pAddJet->e3_b2      = float(fECF->manager->ecfns["3_3"]);
  pAddJet->e3_v1_b2   = float(fECF->manager->ecfns["3_1"]);
  pAddJet->e3_v2_b2   = float(fECF->manager->ecfns["3_2"]);
  pAddJet->e4_v1_b2   = float(fECF->manager->ecfns["4_1"]);
  pAddJet->e4_v2_b2   = float(fECF->manager->ecfns["4_2"]);
  beta=4;
  fECF->calcECFN(beta,lFiltered); 
  pAddJet->e2_sdb4      = float(fECF->manager->ecfns["2_2"]);
  pAddJet->e3_sdb4      = float(fECF->manager->ecfns["3_3"]);
  pAddJet->e3_v1_sdb4   = float(fECF->manager->ecfns["3_1"]);
  pAddJet->e3_v2_sdb4   = float(fECF->manager->ecfns["3_2"]);
  pAddJet->e4_v1_sdb4   = float(fECF->manager->ecfns["4_1"]);
  pAddJet->e4_v2_sdb4   = float(fECF->manager->ecfns["4_2"]);
  beta=0.5;
  fECF->calcECFN(beta,lFiltered);
  pAddJet->e2_sdb05      = float(fECF->manager->ecfns["2_2"]);
  pAddJet->e3_sdb05      = float(fECF->manager->ecfns["3_3"]);
  pAddJet->e3_v1_sdb05   = float(fECF->manager->ecfns["3_1"]);
  pAddJet->e3_v2_sdb05   = float(fECF->manager->ecfns["3_2"]);
  pAddJet->e4_v1_sdb05   = float(fECF->manager->ecfns["4_1"]);
  pAddJet->e4_v2_sdb05   = float(fECF->manager->ecfns["4_2"]);
  */
  double pCorr=1;
  const reco::BasicJet* matchJet = 0;
  // Pruning
  matchJet = match(&itJet,prunedJetCol);
  if(matchJet) {
    //pCorr = correction(pP1Jet,*hRho);
    pAddJet->mass_prun  = matchJet->mass()*pCorr;
  }

  // Trimming
  matchJet = match(&itJet,trimmedJetCol);
  if(matchJet) {
    //pCorr = correction(pP1Jet,*hRho);
    pAddJet->mass_trim  = matchJet->mass()*pCorr;
  }

  // Soft drop
  matchJet = match(&itJet,softdropJetCol);
  if(matchJet) {
    //pCorr = correction(pP1Jet,*hRho);
    pAddJet->mass_sd0  = matchJet->mass()*pCorr;
    fastjet::contrib::SoftDrop sd(0.,0.1,2.0);
    fastjet::PseudoJet sd_jet = sd(inclusive_jets[0]);
    //fastjet::contrib::EnergyCorrelator::Measure measure = fastjet::contrib::EnergyCorrelator::pt_R;
    double beta=1;
    std::vector<fastjet::PseudoJet> lSDClusterParticles = sd_jet.constituents();
    std::sort(lSDClusterParticles.begin(),lSDClusterParticles.end(),JetTools::orderPseudoJet);
    nFilter = TMath::Min(100,(int)lSDClusterParticles.size());
    std::vector<fastjet::PseudoJet> lSDFiltered(lSDClusterParticles.begin(),lSDClusterParticles.begin()+nFilter);

    fECF->calcECFN(beta,lSDFiltered,true);
    pAddJet->e2_sdb1      = float(fECF->manager->ecfns["2_2"]);
    pAddJet->e3_sdb1      = float(fECF->manager->ecfns["3_3"]);
    pAddJet->e3_v1_sdb1   = float(fECF->manager->ecfns["3_1"]);
    pAddJet->e3_v2_sdb1   = float(fECF->manager->ecfns["3_2"]);
    pAddJet->e4_v1_sdb1   = float(fECF->manager->ecfns["4_1"]);
    pAddJet->e4_v2_sdb1   = float(fECF->manager->ecfns["4_2"]);
    /*
    beta=2;
    fECF->calcECFN(beta,lSDFiltered);
    pAddJet->e2_sdb2      = float(fECF->manager->ecfns["2_2"]);
    pAddJet->e3_sdb2      = float(fECF->manager->ecfns["3_3"]);
    pAddJet->e3_v1_sdb2   = float(fECF->manager->ecfns["3_1"]);
    pAddJet->e3_v2_sdb2   = float(fECF->manager->ecfns["3_2"]);
    pAddJet->e4_v1_sdb2   = float(fECF->manager->ecfns["4_1"]);
    pAddJet->e4_v2_sdb2   = float(fECF->manager->ecfns["4_2"]);
    */
  }
  
  // Recursive Soft drop
  /*
  fastjet::PseudoJet pRSM1Jet = (*fRecursiveSoftDrop1) (inclusive_jets[0]);
  fastjet::PseudoJet pRSM2Jet = (*fRecursiveSoftDrop2) (inclusive_jets[0]);
  pAddJet->mass_rsd0  = pRSM1Jet.m()*pCorr;
  pAddJet->mass_rsd1  = pRSM2Jet.m()*pCorr;
  */
  /*
  // Q-Jets
  pAddJet->qjet = 0;
  if(itJet.pt() > 100) pAddJet->qjet = JetTools::qJetVolatility(lClusterParticles,25.,fRand->Rndm());  // (!) why pT > 100 cut? computation time?
  //
  // Subjets
  //
  */
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
      q1      = JetTools::jetCharge(*itSubJet);
      
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
      q2      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet3 || itSubJet->pt() > subjet3->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;

      subjet3 = &(*itSubJet);
      csv3    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      qgid3   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
      q3      = JetTools::jetCharge(*itSubJet);
      
    } else if(!subjet4 || itSubJet->pt() > subjet4->pt()) {
      subjet4 = &(*itSubJet);
      csv4    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      qgid4   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q4      = JetTools::jetCharge(*itSubJet);
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
  //CA15 double-b with subjet
  //
  //
  edm::Handle<reco::BoostedDoubleSVTagInfoCollection> hBoostedDoubleSVTagInfo;
  iEvent.getByToken(fTokBoostedDoubleSVTagInfo,hBoostedDoubleSVTagInfo);  // 
  assert(hBoostedDoubleSVTagInfo.isValid());
  //

  //match to jet 
  bool matched ;
  reco::BoostedDoubleSVTagInfoCollection::const_iterator matchTI = hBoostedDoubleSVTagInfo->end();
  for( reco::BoostedDoubleSVTagInfoCollection::const_iterator itTI = hBoostedDoubleSVTagInfo->begin(); itTI != hBoostedDoubleSVTagInfo->end(); ++itTI ) {
      matched =false;
      const reco::JetBaseRef jetTI = itTI->jet();
      if( jetTI->px() ==  jetBaseRef->px()  && jetTI->pz() ==  jetBaseRef->pz() ) {
        matchTI = itTI;
        matched = true;
        break;
      }
  }
  if( matchTI != hBoostedDoubleSVTagInfo->end() && matched) {
    const reco::TaggingVariableList vars = matchTI->taggingVariables();
    float SubJet_csv__ =  std::min( pAddJet->sj2_csv , pAddJet->sj1_csv) ;
    if (SubJet_csv__ < -1 || SubJet_csv__ >1.0) SubJet_csv__ =-1;
    float z_ratio__ = vars.get(reco::btau::z_ratio);
    float trackSipdSig_3__ = vars.get(reco::btau::trackSip3dSig_3);
    float trackSipdSig_2__ = vars.get(reco::btau::trackSip3dSig_2);
    float trackSipdSig_1__ = vars.get(reco::btau::trackSip3dSig_1);
    float trackSipdSig_0__ = vars.get(reco::btau::trackSip3dSig_0);
    float trackSipdSig_1_0__ = vars.get(reco::btau::tau2_trackSip3dSig_0);
    float trackSipdSig_0_0__ = vars.get(reco::btau::tau1_trackSip3dSig_0);
    float trackSipdSig_1_1__ = vars.get(reco::btau::tau2_trackSip3dSig_1);
    float trackSipdSig_0_1__ = vars.get(reco::btau::tau1_trackSip3dSig_1);
    float trackSip2dSigAboveCharm_0__ = vars.get(reco::btau::trackSip2dSigAboveCharm);
    float trackSip2dSigAboveBottom_0__ = vars.get(reco::btau::trackSip2dSigAboveBottom_0);
    float trackSip2dSigAboveBottom_1__ = vars.get(reco::btau::trackSip2dSigAboveBottom_1);
    float tau1_trackEtaRel_0__ = vars.get(reco::btau::tau2_trackEtaRel_0);
    float tau1_trackEtaRel_1__ = vars.get(reco::btau::tau2_trackEtaRel_1);
    float tau1_trackEtaRel_2__ = vars.get(reco::btau::tau2_trackEtaRel_2);
    float tau0_trackEtaRel_0__ = vars.get(reco::btau::tau1_trackEtaRel_0);
    float tau0_trackEtaRel_1__ = vars.get(reco::btau::tau1_trackEtaRel_1);
    float tau0_trackEtaRel_2__ = vars.get(reco::btau::tau1_trackEtaRel_2);
    float tau_vertexMass_0__ = vars.get(reco::btau::tau1_vertexMass);
    float tau_vertexEnergyRatio_0__ = vars.get(reco::btau::tau1_vertexEnergyRatio);
    float tau_vertexDeltaR_0__ = vars.get(reco::btau::tau1_vertexDeltaR);
    float tau_flightDistance2dSig_0__ = vars.get(reco::btau::tau1_flightDistance2dSig);
    float tau_vertexMass_1__ = vars.get(reco::btau::tau2_vertexMass);
    float tau_vertexEnergyRatio_1__ = vars.get(reco::btau::tau2_vertexEnergyRatio);
    float tau_flightDistance2dSig_1__ = vars.get(reco::btau::tau2_flightDistance2dSig);
    float jetNTracks__ = vars.get(reco::btau::jetNTracks);
    float nSV__ = vars.get(reco::btau::jetNSecondaryVertices);
    float massPruned__ =pAddJet->mass_prun;
    float flavour__ = -1;//itJet.partonFlavor();   // they're spectator variables
    float nbHadrons__ = -1;//itJet.hadronFlavor(); // 
    float ptPruned__ = itJet.pt();
    float etaPruned__ =itJet.eta();
    
    pAddJet->Double_sub = fJetBoostedBtaggingMVACalc.mvaValue(massPruned__, flavour__, nbHadrons__, ptPruned__, etaPruned__,SubJet_csv__,z_ratio__,trackSipdSig_3__,trackSipdSig_2__,trackSipdSig_1__,trackSipdSig_0__,trackSipdSig_1_0__,trackSipdSig_0_0__,trackSipdSig_1_1__,trackSipdSig_0_1__,trackSip2dSigAboveCharm_0__,trackSip2dSigAboveBottom_0__,trackSip2dSigAboveBottom_1__,tau0_trackEtaRel_0__,tau0_trackEtaRel_1__,tau0_trackEtaRel_2__,tau1_trackEtaRel_0__,tau1_trackEtaRel_1__,tau1_trackEtaRel_2__,tau_vertexMass_0__,tau_vertexEnergyRatio_0__,tau_vertexDeltaR_0__,tau_flightDistance2dSig_0__,tau_vertexMass_1__,tau_vertexEnergyRatio_1__,tau_flightDistance2dSig_1__,jetNTracks__,nSV__,false);
  }
  else std::cout<< "   not found matched double-b tag info  "<<std::endl;	
  //Secondary Vertices
  
  //sec vtx                    
  pAddJet->svtx.clear();
  if(fComputeSVInfo) { 
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
    iEvent.getByToken(fTokSVName, secVertices);
    assert(secVertices.isValid());
    const reco::VertexCompositePtrCandidateCollection svtx=*secVertices;
    pAddJet->svtx.clear();
    TClonesArray &rSVArray = *iSVArr;
    for (const reco::VertexCompositePtrCandidate &sv : svtx) {
      if (reco::deltaR(sv,itJet)>0.7) { continue; }
      bool pFind = false;
      for(int i0 = 0; i0 < iSVArr->GetEntries(); i0++) { 
	baconhep::TSVtx* pLSV = (baconhep::TSVtx*) iSVArr->At(i0);
	if(fabs(sv.pt()-pLSV->pt) < 0.01 && fabs(pLSV->eta-sv.eta()) < 0.01 && fabs(pLSV->phi-sv.phi()) < 0.01) {pFind=true; break;}
      }
      if(pFind) continue;
      assert(rSVArray.GetEntries() < rSVArray.GetSize());
      const int svIndex = rSVArray.GetEntries();
      new(rSVArray[svIndex]) baconhep::TSVtx();
      baconhep::TSVtx* pSV = (baconhep::TSVtx*)rSVArray[svIndex];
      pAddJet->svtx.push_back(svIndex);
      pSV->pt          = sv.pt();
      pSV->eta         = sv.eta();
      pSV->phi         = sv.phi();
      pSV->mass        = sv.mass();
      pSV->etarel      = fabs(sv.eta()-itJet.eta());
      pSV->phirel      = fabs(reco::deltaPhi(sv.phi(),itJet.phi()));
      pSV->sv_deltaR   = fabs(reco::deltaR(sv,itJet));
      pSV->sv_ntracks  = sv.numberOfDaughters();
      pSV->sv_chi2     = sv.vertexChi2();
      pSV->sv_ndf      = sv.vertexNdof();
      pSV->sv_normchi2 = pSV->sv_chi2/pSV->sv_ndf; // really?!!!
      pSV->sv_dxy      = vertexDxy(sv,pv).value();
      pSV->sv_dxyerr   = vertexDxy(sv,pv).error();
      pSV->sv_dxysig   = pSV->sv_dxy/pSV->sv_dxyerr; //really?!!!
      pSV->sv_d3d      = vertexD3d(sv,pv).value();
      pSV->sv_d3derr   = vertexD3d(sv,pv).error();
      pSV->sv_d3dsig   = vertexD3d(sv,pv).value()/vertexD3d(sv,pv).error();
      pSV->sv_enratio  = sv.energy()/itJet.energy();
    }
  }

  // Lepton in Jet

  /*
  edm::Handle<pat::MuonCollection> hMuonProduct;
  iEvent.getByToken(fTokPatMuonName,hMuonProduct);
  assert(hMuonProduct.isValid());
  const pat::MuonCollection *muonCol = hMuonProduct.product();

  edm::Handle<reco::CandSoftLeptonTagInfoCollection> hsoftPFMuonTagInfo;
  iEvent.getByToken(fToksoftPFMuonTagInfo,hsoftPFMuonTagInfo);
  assert(hsoftPFMuonTagInfo.isValid());

  edm::Handle<pat::ElectronCollection> hEleProduct;
  iEvent.getByToken(fTokPatEleName,hEleProduct);
  assert(hEleProduct.isValid());
  const pat::ElectronCollection *eleCol = hEleProduct.product();

  edm::Handle<reco::CandSoftLeptonTagInfoCollection> hsoftPFElectronTagInfo;
  iEvent.getByToken(fToksoftPFElectronTagInfo, hsoftPFElectronTagInfo);
  assert(hsoftPFElectronTagInfo.isValid());

  // match to jet and find highest pT matched lepton
  bool matchedMu,matchedEle;
  int indexTM(-1);
  float lepPt(-100), lepEta(-100), lepPhi(-100);
  float lepRPt(-100), lepREta(-100), lepRPhi(-100);
  float lepId(0), lepRId(0);
  reco::CandSoftLeptonTagInfoCollection::const_iterator matchSM = hsoftPFMuonTagInfo->end();
  for( reco::CandSoftLeptonTagInfoCollection::const_iterator itTM = hsoftPFMuonTagInfo->begin(); itTM != hsoftPFMuonTagInfo->end(); ++itTM ) {
    indexTM +=1;
    matchedMu =false;
    const reco::JetBaseRef jetTM = itTM->jet();
    if( jetTM->px() ==  jetBaseRef->px()  && jetTM->pz() ==  jetBaseRef->pz() ) {
      matchSM = itTM;
      matchedMu = true;
      break;
    }
  }

  reco::CandSoftLeptonTagInfoCollection::const_iterator matchSE = hsoftPFElectronTagInfo->end();
  for( reco::CandSoftLeptonTagInfoCollection::const_iterator itTE = hsoftPFElectronTagInfo->begin(); itTE != hsoftPFElectronTagInfo->end(); ++itTE ) {
    matchedEle =false;
    const reco::JetBaseRef jetTE = itTE->jet();
    if( jetTE->px() ==  jetBaseRef->px()  && jetTE->pz() ==  jetBaseRef->pz() ) {
      matchSE = itTE;
      matchedEle = true;
      break;
    }
  }

  
  if( matchSM != hsoftPFMuonTagInfo->end() && matchedMu) {
    for (size_t PFmu = 0; PFmu < (size_t)matchSM->leptons(); ++PFmu) {
      if(matchSM->lepton(PFmu)->pt() > lepPt) { 
	lepPt = matchSM->lepton(PFmu)->pt(); 
	lepEta = matchSM->lepton(PFmu)->eta();
        lepPhi = matchSM->lepton(PFmu)->phi();
	lepId = 13; }
    }
  }
  else std::cout<< "   not found matched soft muon tag info  "<<std::endl;

  if( matchSE != hsoftPFElectronTagInfo->end() && matchedEle) {
    for (size_t PFele = 0; PFele < (size_t)matchSE->leptons(); ++PFele) {
      if(matchSE->lepton(PFele)->pt() > lepPt) { 
	lepPt = matchSE->lepton(PFele)->pt(); 
	lepEta = matchSE->lepton(PFele)->eta();
        lepPhi = matchSE->lepton(PFele)->phi();
	lepId = 11; 
      }
    }
  }
  else std::cout<< "   not found matched soft electron tag info  "<<std::endl;

  for(pat::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {
    if(itMu->pt() > lepRPt && deltaR( itJet.eta(), itJet.phi(), itMu->eta(), itMu->phi()) < fConeSize){
      lepRPt = itMu->pt(); 
      lepREta = itMu->eta();
      lepRPhi = itMu->phi();
      lepRId = 13;
    }
  }
  for(pat::ElectronCollection::const_iterator itEle = eleCol->begin(); itEle!=eleCol->end(); ++itEle) {
    if(itEle->pt() > lepRPt  && deltaR( itJet.eta(), itJet.phi(), itEle->eta(), itEle->phi()) < fConeSize) {
      lepRPt = itEle->pt(); 
      lepREta = itEle->eta();
      lepRPhi = itEle->phi();
      lepRId = 11;
    }
  }
  */

  float lepCPt(-100), lepCEta(-100), lepCPhi(-100);
  float lepCId(0);

  if(JetTools::leptons(itJet,3)> 0 && JetTools::leptons(itJet,7)<fConeSize) {
    lepCPt = JetTools::leptons(itJet,3);
    lepCEta = JetTools::leptons(itJet,5);
    lepCPhi = JetTools::leptons(itJet,6);
    lepCId = JetTools::leptons(itJet,4);
  }
  
  pAddJet->lepCPt = lepCPt;
  pAddJet->lepCEta = lepCEta;
  pAddJet->lepCPhi = lepCPhi;
  pAddJet->lepCId = lepCId;

  // LSF
  std::vector<fastjet::PseudoJet> vSubCInc; pAddJet->lsfCInc = JetTools::lsf(lClusterParticles, vSubCInc, lepCPt, lepCEta, lepCPhi, lepCId, 0.2, 2);
  std::vector<fastjet::PseudoJet> vSubC_2;  pAddJet->lsfC_2 = JetTools::lsf(lClusterParticles, vSubC_2, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 2);
  std::vector<fastjet::PseudoJet> vSubC_3;  pAddJet->lsfC_3 = JetTools::lsf(lClusterParticles, vSubC_3, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 3);
  std::vector<fastjet::PseudoJet> vSubC_4;  pAddJet->lsfC_4 = JetTools::lsf(lClusterParticles, vSubC_4, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 4);

  if(vSubC_3.size() > 0) {
    pAddJet->lsfC_3_sj1_pt = vSubC_3[0].pt();
    pAddJet->lsfC_3_sj1_eta = vSubC_3[0].eta();
    pAddJet->lsfC_3_sj1_phi = vSubC_3[0].phi();
    pAddJet->lsfC_3_sj1_m = vSubC_3[0].m();
  }
  if(vSubC_3.size() > 1) {
    pAddJet->lsfC_3_sj2_pt = vSubC_3[1].pt();
    pAddJet->lsfC_3_sj2_eta = vSubC_3[1].eta();
    pAddJet->lsfC_3_sj2_phi = vSubC_3[1].phi();
    pAddJet->lsfC_3_sj2_m = vSubC_3[1].m();
  }
  if(vSubC_3.size() > 2) {
    pAddJet->lsfC_3_sj3_pt = vSubC_3[2].pt();
    pAddJet->lsfC_3_sj3_eta = vSubC_3[2].eta();
    pAddJet->lsfC_3_sj3_phi = vSubC_3[2].phi();
    pAddJet->lsfC_3_sj3_m = vSubC_3[2].m();
  }

  pAddJet->lmdCInc = JetTools::lsf(lClusterParticles, vSubCInc, lepCPt, lepCEta, lepCPhi, lepCId, 0.2, 2, 1);
  pAddJet->lmdC_2 = JetTools::lsf(lClusterParticles, vSubC_2, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 2, 1);
  pAddJet->lmdC_3 = JetTools::lsf(lClusterParticles, vSubC_3, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 3, 1);
  pAddJet->lmdC_4 = JetTools::lsf(lClusterParticles, vSubC_4, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 4, 1);

  //
  // Top Tagging
  //
  pAddJet->topTagType=0;

  if(fTopTaggerName.compare("CMS")==0) {  // CMS Top Tagger

    edm::Handle<reco::BasicJetCollection> hCMSTTJetProduct;
    iEvent.getByToken(fTokCMSTTJetProduct,hCMSTTJetProduct);  // (!) hard-code
    assert(hCMSTTJetProduct.isValid());
    const reco::BasicJetCollection *cmsttJetCol = hCMSTTJetProduct.product();

    edm::Handle<reco::PFJetCollection> hCMSTTSubJetProduct;
    iEvent.getByToken(fTokCMSTTSubJetProduct,hCMSTTSubJetProduct);  // (!) hard-code
    assert(hCMSTTSubJetProduct.isValid());
    const reco::PFJetCollection *cmsttSubJetCol = hCMSTTSubJetProduct.product();

    matchJet = match(&itJet,cmsttJetCol);
    if(matchJet) {
      pAddJet->topTagType |= kCMSTT;

      unsigned int nsubjets=0;
      const reco::PFJet *sub1=0, *sub2=0, *sub3=0;
      for(reco::PFJetCollection::const_iterator itSub = cmsttSubJetCol->begin(); itSub!=cmsttSubJetCol->end(); ++itSub) {
        if(reco::deltaR(itJet.eta(),itJet.phi(),itSub->eta(),itSub->phi())>fConeSize) continue;  // (!) get associated subjets by dR...is there a better way???

        nsubjets++;

        if(!sub1 || itSub->pt() > sub1->pt()) {
          sub3 = sub2;
          sub2 = sub1;
          sub1 = &(*itSub);

        } else if(!sub2 || itSub->pt() > sub2->pt()) {
          sub3 = sub2;
          sub2 = &(*itSub);

        } else if(!sub3 || itSub->pt() > sub3->pt()) {
          sub3 = &(*itSub);
        }
      }
      pAddJet->top_n_subjets = nsubjets;

      TLorentzVector vSub1; if(sub1) { vSub1.SetPtEtaPhiM(sub1->pt(), sub1->eta(), sub1->phi(), sub1->mass()); }
      TLorentzVector vSub2; if(sub2) { vSub2.SetPtEtaPhiM(sub2->pt(), sub2->eta(), sub2->phi(), sub2->mass()); }
      TLorentzVector vSub3; if(sub3) { vSub3.SetPtEtaPhiM(sub3->pt(), sub3->eta(), sub3->phi(), sub3->mass()); }
      double m12 = (sub1 && sub2) ? (vSub1+vSub2).M() : 0;
      double m23 = (sub2 && sub3) ? (vSub2+vSub3).M() : 0;
      double m31 = (sub3 && sub1) ? (vSub3+vSub1).M() : 0;
      pAddJet->top_m_min=0;
      if     (m12 < m23 && m12 < m31) { pAddJet->top_m_min = m12; }
      else if(m23 < m12 && m23 < m31) { pAddJet->top_m_min = m23; }
      else if(m31 < m12 && m31 < m23) { pAddJet->top_m_min = m31; }

      pAddJet->top_m_123 = 0;
      pAddJet->top_fRec  = 0;
    }
/*
  } else if(fTopTaggerName.compare("HEP")==0) {  // HEP Top Tagger

    edm::Handle<reco::HTTTopJetTagInfoCollection> hHTTTagInfos;
    //iEvent.getByLabel("CA15HTTCHSAOD",hHTTTagInfos);  // (!) hard-code
    assert(hHTTTagInfos.isValid());
    const reco::HTTTopJetTagInfoCollection *httTagInfoCol = hHTTTagInfos.product();

    double dRmin = fConeSize;
    for(reco::HTTTopJetTagInfoCollection::const_iterator itTagInfo = httTagInfoCol->begin(); itTagInfo!=httTagInfoCol->end(); ++itTagInfo) {
      double dR = reco::deltaR(itJet.eta(), itJet.phi(),
                               itTagInfo->properties().fjEta, itTagInfo->properties().fjPhi);

      if(dR < dRmin) {
        dRmin = dR;
        pAddJet->topTagType |= kHEPTT;
        pAddJet->top_n_subjets = 3;
        pAddJet->top_m_min     = 0;
        pAddJet->top_m_123     = itTagInfo->properties().topMass;
        pAddJet->top_fRec      = itTagInfo->properties().fRec;
      }
    }*/
  }
}

void FillerJet::addJet(baconhep::TAddJet *pAddJet, TClonesArray *iSVArr, const reco::Vertex &pv, const edm::Event &iEvent, const pat::Jet &itJet) {
  //pAddJet->pullAngle = JetTools::jetPullAngle(itJet,hSubJetProduct,fConeSize);
  pAddJet->tau1 = itJet.userFloat(fJettinessName + std::string(":tau1"));
  pAddJet->tau2 = itJet.userFloat(fJettinessName + std::string(":tau2"));
  pAddJet->tau3 = itJet.userFloat(fJettinessName + std::string(":tau3"));
  pAddJet->tau4 = -1;  //(!) Not in MINIAOD

  // Pruning
  pAddJet->mass_prun = itJet.userFloat(fPrunedJetName + std::string("Mass"));

  // Trimming
  pAddJet->mass_trim = itJet.userFloat(fTrimmedJetName + std::string("Mass"));

  // Soft drop
  pAddJet->mass_sd0 = itJet.userFloat(fSoftDropJetName + std::string("Mass"));


  //Bosted b tagging for CA15

  //reco::BoostedDoubleSVTagInfo const *bdsvTagInfo = itJet.tagInfoBoostedDoubleSV();//dynamic_cast<reco::BoostedDoubleSVTagInfo const *>(itJet.tagInfo("pfBoostedDoubleSVCA15"));
  reco::BoostedDoubleSVTagInfo const *bdsvTagInfo = dynamic_cast<reco::BoostedDoubleSVTagInfo const *>(itJet.tagInfo("pfBoostedDoubleSVCA15"));
  const reco::TaggingVariableList vars = bdsvTagInfo->taggingVariables();
  
  float SubJet_csv_ =  std::min( pAddJet->sj2_csv , pAddJet->sj1_csv) ;
  float z_ratio_ = vars.get(reco::btau::z_ratio);
  float trackSipdSig_3_ = vars.get(reco::btau::trackSip3dSig_3);
  float trackSipdSig_2_ = vars.get(reco::btau::trackSip3dSig_2);
  float trackSipdSig_1_ = vars.get(reco::btau::trackSip3dSig_1);
  float trackSipdSig_0_ = vars.get(reco::btau::trackSip3dSig_0);
  float trackSipdSig_1_0_ = vars.get(reco::btau::tau2_trackSip3dSig_0);
  float trackSipdSig_0_0_ = vars.get(reco::btau::tau1_trackSip3dSig_0);
  float trackSipdSig_1_1_ = vars.get(reco::btau::tau2_trackSip3dSig_1);
  float trackSipdSig_0_1_ = vars.get(reco::btau::tau1_trackSip3dSig_1);
  float trackSip2dSigAboveCharm_0_ = vars.get(reco::btau::trackSip2dSigAboveCharm);
  float trackSip2dSigAboveBottom_0_ = vars.get(reco::btau::trackSip2dSigAboveBottom_0);
  float trackSip2dSigAboveBottom_1_ = vars.get(reco::btau::trackSip2dSigAboveBottom_1);
  float tau1_trackEtaRel_0_ = vars.get(reco::btau::tau2_trackEtaRel_0);
  float tau1_trackEtaRel_1_ = vars.get(reco::btau::tau2_trackEtaRel_1);
  float tau1_trackEtaRel_2_ = vars.get(reco::btau::tau2_trackEtaRel_2);
  float tau0_trackEtaRel_0_ = vars.get(reco::btau::tau1_trackEtaRel_0);
  float tau0_trackEtaRel_1_ = vars.get(reco::btau::tau1_trackEtaRel_1);
  float tau0_trackEtaRel_2_ = vars.get(reco::btau::tau1_trackEtaRel_2);
  float tau_vertexMass_0_ = vars.get(reco::btau::tau1_vertexMass);
  float tau_vertexEnergyRatio_0_ = vars.get(reco::btau::tau1_vertexEnergyRatio);
  float tau_vertexDeltaR_0_ = vars.get(reco::btau::tau1_vertexDeltaR);
  float tau_flightDistance2dSig_0_ = vars.get(reco::btau::tau1_flightDistance2dSig);
  float tau_vertexMass_1_ = vars.get(reco::btau::tau2_vertexMass);
  float tau_vertexEnergyRatio_1_ = vars.get(reco::btau::tau2_vertexEnergyRatio);
  float tau_flightDistance2dSig_1_ = vars.get(reco::btau::tau2_flightDistance2dSig);
  float jetNTracks_ = vars.get(reco::btau::jetNTracks);
  float nSV_ = vars.get(reco::btau::jetNSecondaryVertices);
  float massPruned_ =pAddJet->mass_prun;
  float flavour_ = -1;//itJet.partonFlavor();   // they're spectator variables
  float nbHadrons_ = -1;//itJet.hadronFlavor(); // 
  float ptPruned_ =itJet.pt();
  float etaPruned_ =itJet.eta();
    
  pAddJet->Double_sub = fJetBoostedBtaggingMVACalc.mvaValue(massPruned_, flavour_, nbHadrons_, ptPruned_, etaPruned_,SubJet_csv_,z_ratio_,trackSipdSig_3_,trackSipdSig_2_,trackSipdSig_1_,trackSipdSig_0_,trackSipdSig_1_0_,trackSipdSig_0_0_,trackSipdSig_1_1_,trackSipdSig_0_1_,trackSip2dSigAboveCharm_0_,trackSip2dSigAboveBottom_0_,trackSip2dSigAboveBottom_1_,tau0_trackEtaRel_0_,tau0_trackEtaRel_1_,tau0_trackEtaRel_2_,tau1_trackEtaRel_0_,tau1_trackEtaRel_1_,tau1_trackEtaRel_2_,tau_vertexMass_0_,tau_vertexEnergyRatio_0_,tau_vertexDeltaR_0_,tau_flightDistance2dSig_0_,tau_vertexMass_1_,tau_vertexEnergyRatio_1_,tau_flightDistance2dSig_1_,jetNTracks_,nSV_, true);


  if(fComputeSVInfo) { 
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
    iEvent.getByToken(fTokSVName, secVertices);
    assert(secVertices.isValid());
    const reco::VertexCompositePtrCandidateCollection svtx=*secVertices;
    pAddJet->svtx.clear();
    TClonesArray &rSVArray = *iSVArr;
    for (const reco::VertexCompositePtrCandidate &sv : svtx) {
      if (reco::deltaR(sv,itJet)>0.7) { continue; }
      bool pFind = false;
      for(int i0 = 0; i0 < iSVArr->GetEntries(); i0++) { 
	baconhep::TSVtx* pLSV = (baconhep::TSVtx*) iSVArr->At(i0);
	if(fabs(sv.pt()-pLSV->pt) < 0.01 && fabs(pLSV->eta-sv.eta()) < 0.01 && fabs(pLSV->phi-sv.phi()) < 0.01) {pFind=true; break;}
      }
      if(pFind) continue;
      assert(rSVArray.GetEntries() < rSVArray.GetSize());
      const int svIndex = rSVArray.GetEntries();
      new(rSVArray[svIndex]) baconhep::TSVtx();
      baconhep::TSVtx* pSV = (baconhep::TSVtx*)rSVArray[svIndex];
      pAddJet->svtx.push_back(svIndex);
      pSV->pt          = sv.pt();
      pSV->eta         = sv.eta();
      pSV->phi         = sv.phi();
      pSV->mass        = sv.mass();
      pSV->etarel      = fabs(sv.eta()-itJet.eta());
      pSV->phirel      = fabs(reco::deltaPhi(sv.phi(),itJet.phi()));
      pSV->sv_deltaR   = fabs(reco::deltaR(sv,itJet));
      pSV->sv_ntracks  = sv.numberOfDaughters();
      pSV->sv_chi2     = sv.vertexChi2();
      pSV->sv_ndf      = sv.vertexNdof();
      pSV->sv_normchi2 = pSV->sv_chi2/pSV->sv_ndf; // really?!!!
      pSV->sv_dxy      = vertexDxy(sv,pv).value();
      pSV->sv_dxyerr   = vertexDxy(sv,pv).error();
      pSV->sv_dxysig   = pSV->sv_dxy/pSV->sv_dxyerr; //really?!!!
      pSV->sv_d3d      = vertexD3d(sv,pv).value();
      pSV->sv_d3derr   = vertexD3d(sv,pv).error();
      pSV->sv_d3dsig   = vertexD3d(sv,pv).value()/vertexD3d(sv,pv).error();
      pSV->sv_enratio  = sv.energy()/itJet.energy();
    }
  }

  // Jet Shape Correlation observables
//  fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 2.0);
//  fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
//  std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta0 (2,0. ,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta02(2,0.2,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta05(2,0.5,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta10(2,1.0,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta20(2,2.0,fastjet::EnergyCorrelator::pt_R);
//  pAddJet->c2_0   = C2beta0 (inclusive_jets[0]);
//  pAddJet->c2_0P2 = C2beta02(inclusive_jets[0]);
//  pAddJet->c2_0P5 = C2beta05(inclusive_jets[0]);
//  pAddJet->c2_1P0 = C2beta10(inclusive_jets[0]);
//  pAddJet->c2_2P0 = C2beta20(inclusive_jets[0]);

  // Q-Jets
//  pAddJet->qjet = 0;
//  if(itJet->pt() > 100) pAddJet->qjet = JetTools::qJetVolatility(lClusterParticles,25.,fRand->Rndm());  // (!) why pT > 100 cut? computation time?

  //
  // Subjets
  //
/*
  // find/sort up to 4 hardest subjets
  const pat::Jet *subjet1=0, *subjet2=0, *subjet3=0, *subjet4=0;
  double csv1=-2, csv2=-2, csv3=-2, csv4=-2;
  double qgid1=-2, qgid2=-2, qgid3=-2, qgid4=-2;
  double q1=-100, q2=-100, q3=-100, q4=-100;
  pat::JetPtrCollection const &subjetCol = itJet.subjets(fSubJetName);
  for(edm::Ptr<pat::Jet> const &itSubJet : subjetCol) {

    if(!subjet1 || itSubJet->pt() > subjet1->pt()) {
      subjet4 = subjet3;
      //csv4    = csv3;
      //qgid4   = qgid3;
      q4      = q3;

      subjet3 = subjet2;
      //csv3    = csv2;
      //qgid3   = qgid2;
      q3      = q2;

      subjet2 = subjet1;
      //csv2    = csv1;
      //qgid2   = qgid1;
      q2      = q1;

      subjet1 = itSubJet.get();
      //csv1    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid1   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q1      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet2 || itSubJet->pt() > subjet2->pt()) {
      subjet4 = subjet3;
      //csv4    = csv3;
      //qgid4   = qgid3;
      q4      = q3;

      subjet3 = subjet2;
      //csv3    = csv2;
      //qgid3   = qgid2;
      q3      = q2;

      subjet2 = itSubJet.get();
      //csv2    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid2   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q2      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet3 || itSubJet->pt() > subjet3->pt()) {
      subjet4 = subjet3;
      //csv4    = csv3;
      //qgid4   = qgid3;
      q4      = q3;

      subjet3 = itSubJet.get();
      //csv3    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid3   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q3      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet4 || itSubJet->pt() > subjet4->pt()) {
      subjet4 = itSubJet.get();
      //csv4    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid4   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q4      = JetTools::jetCharge(*itSubJet);
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
*/
  //
  // Top Tagging
  //
  if(fTopTaggerName.compare("CMS")==0) {  // CMS Top Tagger
    reco::CATopJetTagInfo const *tagInfo = dynamic_cast<reco::CATopJetTagInfo const *>(itJet.tagInfo("caTop"));  // (!) hard-code
    if(tagInfo) {
      pAddJet->topTagType |= kCMSTT;
      pAddJet->top_n_subjets = tagInfo->properties().nSubJets;
      pAddJet->top_m_min     = tagInfo->properties().minMass;
      pAddJet->top_m_123     = 0;
      pAddJet->top_fRec      = 0;
    }
/*  } else if(fTopTaggerName.compare("HEP")==0) {  // HEP Top Tagger
    reco::HTTTopJetTagInfo const *tagInfo = dynamic_cast<reco::HTTTopJetTagInfo const *>(itJet.tagInfo("CA15HTTCHSMINIAOD"));  // (!) hard-code
    if(tagInfo) {
      pAddJet->topTagType |= kHEPTT;
      pAddJet->top_n_subjets = 3;
      pAddJet->top_m_min     = 0;
      pAddJet->top_m_123     = tagInfo->properties().topMass;
      pAddJet->top_fRec      = tagInfo->properties().fRec;
    }*/
  }
}

//--------------------------------------------------------------------------------------------------
const reco::PFJet* FillerJet::matchPF( const reco::PFJet *iJet,const reco::PFJetCollection *jets ) { 
  int lId = -1;
  double dRmin = 999.;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::PFJet* jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if(dR > fConeSize) continue;
    if ( dR < dRmin ) {
      dRmin = dR;
      lId = i;
    }
  }
  const reco::PFJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

const reco::BasicJet* FillerJet::match( const reco::PFJet *iJet,const reco::BasicJetCollection *jets ) { 
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

const reco::BasicJet* FillerJet::match(const pat::Jet *iJet, const reco::BasicJetCollection *jets) {
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

const reco::GenJet* FillerJet::match(const pat::Jet *iJet, const reco::GenJetCollection *jets) {
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
  if(lId != -1) {
    lJet = &((*jets)[lId]);
  }
  return lJet;
}

void FillerJet::softdrop(const reco::GenJet *iJet,float &iMsd,float &ie2,float &ie3) {
  std::vector<fastjet::PseudoJet>  lClusterParticles;
  for (unsigned i = 0;  i < iJet->numberOfDaughters (); i++) {
    const reco::Candidate * daughter = iJet->daughter( i );
    fastjet::PseudoJet   pPart(daughter->px(),daughter->py(),daughter->pz(),daughter->energy());
    lClusterParticles.emplace_back(pPart);
  }
  fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 0.8);
  fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
  std::vector<fastjet::PseudoJet>  lOutJets = lCClust_seq.inclusive_jets(0.0);
  if(lOutJets.size() == 0) {
    std::cout << "no jets " << std::endl;
    return;
  }
  double beta=1;
  fastjet::contrib::SoftDrop SD(0.,0.1,0.8);
  fastjet::PseudoJet SD_jet = SD(lOutJets[0]);
  iMsd = SD_jet.m();
  //std::cout << iMsd << std::endl;
  std::vector<fastjet::PseudoJet> lSDClusterParticles = SD_jet.constituents();
  std::sort(lSDClusterParticles.begin(),lSDClusterParticles.end(),JetTools::orderPseudoJet);
  int nFilter = TMath::Min(100,(int)lSDClusterParticles.size());
  std::vector<fastjet::PseudoJet> lSDFilter(lSDClusterParticles.begin(),lSDClusterParticles.begin()+nFilter);
  fECFn->calcECFN(beta,lSDFilter);
  ie2 = float(fECFn->manager->ecfns["2_2"]);
  ie3 = float(fECFn->manager->ecfns["3_2"]);
}
