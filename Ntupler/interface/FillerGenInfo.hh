#ifndef BACONPROD_NTUPLER_FILLERGENINFO_HH
#define BACONPROD_NTUPLER_FILLERGENINFO_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

class TClonesArray;


namespace baconhep
{
  class TGenEventInfo;  // foward declaration

  class FillerGenInfo
  {
    public:
       FillerGenInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerGenInfo();
      
      void fill(TGenEventInfo    *genEvtInfo,     // output object to be filled
                TClonesArray     *particlesArr,   // output array of particles to be filled
                TClonesArray     *weightsArr,     // output array of LHE weights to be filled
		const edm::Event &iEvent,         // EDM event info
                const float       iXS);           // LHE cross section

    protected:        
      // EDM object collection names
      std::string fGenEvtInfoName;
      std::string fLHEEvtInfoName;
      std::string fGenParName;
      bool fFillAll;
      edm::EDGetTokenT<GenEventInfoProduct>         fTokGenEvent;
      edm::EDGetTokenT<reco::GenParticleCollection> fTokGenPar;
      edm::EDGetTokenT<LHEEventProduct>             fTokLHEEventInfo;
      bool fFillLHEWeights;
  };
}
#endif
