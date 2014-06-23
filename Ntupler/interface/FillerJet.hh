#ifndef BACONPROD_NTUPLER_FILLERJET_HH
#define BACONPROD_NTUPLER_FILLERJET_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/JetPUIDMVACalculator.hh"
#include "BaconProd/Utils/interface/SoftDrop.hh"
#include "BaconProd/Utils/interface/CMSTopTagger.hh"
#include "BaconProd/Utils/interface/HEPTopTaggerWrapper.h"
#include "BaconAna/DataFormats/interface/TTopJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "TRandom2.h"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
class FactorizedJetCorrector;
class JetCorrectionUncertainty;
namespace trigger {
  class TriggerEvent;
}

namespace baconhep
{
  class FillerJet
  {
    public:
       FillerJet(const edm::ParameterSet &iConfig, const double coneSize=0.5, const std::string prefix="",std::string postfix="");		
      ~FillerJet();
      
      
      void fill(TClonesArray                     *array,           // output array to be filled
		TClonesArray                     *iExtraArray,     // Extra Array to be filled
		TClonesArray                     *iTopArray,       // Top Jet Array to be filled
                const edm::Event                 &iEvent,          // event info
		const edm::EventSetup            &iSetup,          // event setup info
	        const reco::Vertex		 &pv,	           // event primary vertex
		const std::vector<TriggerRecord> &triggerRecords,  // list of trigger names and objects
		const trigger::TriggerEvent      &triggerEvent);   // event trigger objects
            
    protected:
      void initJetCorr(const std::vector<std::string> &jecFiles, 
                       const std::vector<std::string> &jecUncFiles,
		       const std::vector<std::string> &jecFilesForID,
		       bool iCHS=false);
      
      double correction(fastjet::PseudoJet &iJet,double iRho);      
      void   addJet(TAddJet *pPFJet,const reco::PFJet &itJet,double iRho);
      void   topJet(TTopJet *pPFJet,const reco::PFJet &itJet,double iRho);
      fastjet::PseudoJet CACluster(fastjet::PseudoJet &iJet, fastjet::ClusterSequenceArea &iCAClustering); 
      //float              getTau( fastjet::PseudoJet &iJet,int iN, float iKappa );
      const reco::BasicJet*    match( const reco::PFJet *jet,const reco::BasicJetCollection  *jets );
      const reco::GenJet*      match( const reco::PFJet *jet,const reco::GenJetCollection    *jets );
      
      // Jet cuts
      double fMinPt;
      bool   fUseGen;
      
      // EDM object collection names
      std::string fPVName;
      std::string fRhoName;
      std::string fJetName;
      std::string fGenJetName;
      std::string fJetFlavorName;
      std::string fJetFlavorPhysName;
      std::string fPruneJetName;
      std::string fSubJetName;
      std::string fCSVbtagName;
      std::string fCSVbtagSubJetName;
      std::string fJettinessName;
      std::string fQGLikelihood;
      std::string fQGLikelihoodSubJets;
      double      fConeSize;
      bool        fComputeFullJetInfo;
      
      // Jet ID MVA
      JetPUIDMVACalculator fJetPUIDMVACalc;


      fastjet::JetDefinition*       fJetDef;
      fastjet::JetDefinition*       fGenJetDef;
      fastjet::JetDefinition*       fCAJetDef;
    
      fastjet::ActiveAreaSpec*      fActiveArea;
      fastjet::AreaDefinition*      fAreaDefinition;
      fastjet::ClusterSequenceArea* fClustering;
      
      fastjet::Pruner* fPruner1;
      fastjet::Pruner* fPruner2;

      fastjet::Filter* fFilter1;
      fastjet::Filter* fFilter2;

      fastjet::contrib::SoftDropTagger *fSoftDrop1;
      fastjet::contrib::SoftDropTagger *fSoftDrop2;
      fastjet::contrib::SoftDropTagger *fSoftDrop3;

      fastjet::Filter* fTrimmer1;
      fastjet::Filter* fTrimmer2;
      fastjet::Filter* fTrimmer3;
      fastjet::Filter* fTrimmer4;

      fastjet::CMSTopTagger* fCMSTopTagger;
      fastjet::HEPTopTagger* fHEPTopTagger;

      TRandom2*        fRand;

      bool passPFLooseID();
      
      // jet correctors
      FactorizedJetCorrector   *fJetCorr, *fJetCorrForID;
      JetCorrectionUncertainty *fJetUnc;
  };
}
#endif
