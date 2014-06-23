#ifndef BACONPROD_NTUPLER_FILLERTAU_HH
#define BACONPROD_NTUPLER_FILLERTAU_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
//#include "BaconProd/Utils/interface/TauIsoMVACalculator.hh"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/Common/interface/Handle.h"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class MyTauDiscHandle
  {
    public:
      MyTauDiscHandle(const std::string inName="", const unsigned long inFlag=0):name(inName),flag(inFlag){}
      ~MyTauDiscHandle(){}
    
      edm::Handle<reco::PFTauDiscriminator> handle;  // EDM handle
      std::string name;                              // EDM name
      unsigned long flag;                            // bacon flag bit value

      float value(reco::PFTauRef tauRef) {
  	return handle.isValid() ? (*handle)[tauRef] : 0;
      }
  };  
  
  class FillerTau
  {
    public:
      FillerTau(const edm::ParameterSet &iConfig);
      ~FillerTau();
            
      void fill(TClonesArray                     *array,           // output array to be filled
                const edm::Event                 &iEvent,          // event info
		const edm::EventSetup            &iSetup,          // event setup info
		const reco::Vertex               &pv,              // event primary vertex
		const std::vector<TriggerRecord> &triggerRecords,  // list of trigger names and objects
		const trigger::TriggerEvent      &triggerEvent);   // event trigger objects
      
      
      // Tau cuts
      double fMinPt;
            
      // EDM object collect names
      std::string fTauName;
      std::string fRhoName;
      
      // Tau isolation MVAs
      //TauIsoMVACalculator fRingIso;
      //TauIsoMVACalculator fRingIso2;
  
    protected:
      std::vector<MyTauDiscHandle*> fMyTauDiscHandles;
  };
}
#endif
