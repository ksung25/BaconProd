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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

class TClonesArray;


namespace baconhep
{
  class TGenEventInfo;  // foward declaration

  class simPrimaryVertex {
  public:
    simPrimaryVertex(double x1,double y1,double z1):x(x1),y(y1),z(z1),ptsq(0),nGenTrk(0){};
    double x,y,z;
    HepMC::FourVector ptot;
    //HepLorentzVector ptot;
    double ptsq;
    int nGenTrk;
    std::vector<int> finalstateParticles;
    std::vector<int> simTrackIndex;
    std::vector<int> genVertex;
    const reco::Vertex *recVtx;
  };

  class FillerGenInfo
  {
    public:
       FillerGenInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerGenInfo();
      
      void fill(TGenEventInfo    *genEvtInfo,     // output object to be filled
                TClonesArray     *particlesArr,   // output array of particles to be filled
                TClonesArray     *vtxArr,         // gen vertex
                TClonesArray     *weightsArr,     // output array of LHE weights to be filled
		const edm::Event &iEvent,         // EDM event info
                const float       iXS);           // LHE cross section

    protected:        
      // EDM object collection names
      std::string fGenEvtInfoName;
      std::string fLHEEvtInfoName;
      std::string fGenParName;
      std::string fPackGenParName;
      std::string fHepMCProduct;
      bool fFillAll;
      edm::EDGetTokenT<GenEventInfoProduct>         fTokGenEvent;
      edm::EDGetTokenT<reco::GenParticleCollection> fTokGenPar;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> fTokGenPackPar;
      edm::EDGetTokenT<LHEEventProduct>             fTokLHEEventInfo;
      bool fFillLHEWeights;
      int flavor(int iId);
  };
}
#endif
