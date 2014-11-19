#ifndef BACONPROD_NTUPLER_FILLERGENJET_HH
#define BACONPROD_NTUPLER_FILLERGENJET_HH

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TClonesArray;

namespace baconhep
{
  class FillerGenJets
  {
    public:
       FillerGenJets(const edm::ParameterSet &iConfig);
      ~FillerGenJets();
       void    fill(TClonesArray *array,const edm::Event &iEvent);
    protected:
       double* genCone(const reco::GenJet *iJet,const reco::GenParticleCollection &iGenParticles,double iDRMin,double iDRMax,int iType);
       int     flavor (const reco::GenJet *iJet,const reco::GenParticleCollection &iGenParticles);            
      
      // EDM object collection names
      std::string fGenParName;
      std::string fGenJetName;
  };
}
#endif
