#ifndef BACONPROD_NTUPLER_FILLERRH_HH
#define BACONPROD_NTUPLER_FILLERRH_HH

#include <string>

// forward class declarations
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

class TClonesArray;


namespace baconhep
{
  class FillerRH
  {
    public:
       FillerRH(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerRH();
      
       void fill(TClonesArray       *array,    // output array to be filled
		 const edm::Event   &iEvent,const edm::EventSetup &iSetup);  // event info
   
    protected:
      // EDM object collection names
      edm::InputTag fRHName;
      const CaloSubdetectorGeometry *fHCAL;
      edm::EDGetTokenT<HBHERecHitCollection> fTokRH;
  };
}
#endif
