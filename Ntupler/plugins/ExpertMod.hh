#include "FWCore/Framework/interface/MakerMacros.h"    // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"   // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"     // EDAnalyzer class
#include <string>                                      // string class

// forward class declarations
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TFile;
class TH1D;
class TTree;
class TClonesArray;
namespace edm {
  class TriggerResults;
  class TriggerNames;
}
namespace baconhep {
  class TEventInfo;
  class TGenEventInfo;
  class TGenWeight;
  class TTrigger;
  class FillerEventInfo;
  class FillerGenInfo;
  class FillerVertex;
  class FillerPF;
}

//
class ExpertMod : public edm::EDAnalyzer {
  public:
    explicit ExpertMod(const edm::ParameterSet &iConfig);
    ~ExpertMod();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

    virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void endRun  (const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void beginLuminosityBlock(const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
    virtual void endLuminosityBlock  (const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);

    // specify trigger paths of interest
    void setTriggers();
    
    // initialization from HLT menu; needs to be called on every change in HLT menu
    void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);


    //--------------------------------------------------------------------------------------------------
    //  data members
    //==================================================================================================   
    bool fSkipOnHLTFail;
    
    // variables to handle triggers
    edm::ParameterSetID fTriggerNamesID;
    edm::InputTag	fHLTTag;
    edm::InputTag       fHLTObjTag;
    std::string         fHLTFile;
        
    // bacon fillers
    baconhep::FillerEventInfo *fFillerEvtInfo;
    baconhep::FillerGenInfo   *fFillerGenInfo;
    baconhep::FillerVertex    *fFillerPV;
    baconhep::FillerPF        *fFillerPF;
 
    baconhep::TTrigger        *fTrigger;
    
    bool fIsActiveEvtInfo;
    bool fIsActiveGenInfo;
    bool fIsActivePV;
    bool fIsActivePF;
    
    // Objects and arrays for output file
    std::string              fOutputName;
    TFile                   *fOutputFile;
    TH1D                    *fTotalEvents;
    TTree                   *fEventTree;
    baconhep::TEventInfo    *fEvtInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    baconhep::TGenWeight    *fGenWeight;
    TClonesArray            *fGenParArr;
    TClonesArray            *fPFParArr;
    TClonesArray	    *fPVArr;
};
