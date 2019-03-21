#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerGenInfo::FillerGenInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fGenEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmGenEventInfoName","generator")),
  fLHEEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmLHEEventInfoName","externalLHEProducer")),
  fGenParName    (iConfig.getUntrackedParameter<std::string>("edmGenParticlesName","genParticles")),
  fFillAll       (iConfig.getUntrackedParameter<bool>("fillAllGen",false)),
  fFillLHEWeights(iConfig.getUntrackedParameter<bool>("fillLHEWeights",false))
{
  fTokGenEvent     = iC.mayConsume<GenEventInfoProduct>        ( fGenEvtInfoName);
  fTokGenPar       = iC.mayConsume<reco::GenParticleCollection>( fGenParName    );
  fTokLHEEventInfo = iC.mayConsume<LHEEventProduct>            ( fLHEEvtInfoName );
}

//--------------------------------------------------------------------------------------------------
FillerGenInfo::~FillerGenInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerGenInfo::fill(TGenEventInfo *genEvtInfo, TClonesArray *particlesArr, TClonesArray *weightsArr,
                         const edm::Event &iEvent, const float iXS)
{
  assert(particlesArr);
  assert(!fFillLHEWeights || weightsArr);

  // Get generator event information
  edm::Handle<GenEventInfoProduct> hGenEvtInfoProduct;
  iEvent.getByToken(fTokGenEvent,hGenEvtInfoProduct);
  assert(hGenEvtInfoProduct.isValid());

  const gen::PdfInfo *pdfInfo = (hGenEvtInfoProduct->hasPDF()) ? hGenEvtInfoProduct->pdf() : 0;
  genEvtInfo->id_1     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.first  : 0;
  genEvtInfo->id_2     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.second : 0;
  genEvtInfo->x_1      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.first   : 0;
  genEvtInfo->x_2      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.second  : 0;
  genEvtInfo->scalePDF = hGenEvtInfoProduct->qScale();
  genEvtInfo->xs       = iXS;
  genEvtInfo->weight   = hGenEvtInfoProduct->weight();

  // Get LHE event information
  if(fFillLHEWeights) {
    edm::Handle<LHEEventProduct> hLHEEvtInfoProduct;
    iEvent.getByToken(fTokLHEEventInfo,hLHEEvtInfoProduct);
    assert(hLHEEvtInfoProduct.isValid());
    TClonesArray &rWeightsArray = *weightsArr;
    unsigned int lMax = hLHEEvtInfoProduct->weights().size();
    if(lMax > 110) lMax = 110.;
    for(unsigned int iw=0; iw< lMax; iw++) {
      // construct object and place in array
      assert(rWeightsArray.GetEntries() < rWeightsArray.GetSize());
      const int index = rWeightsArray.GetEntries();
      new(rWeightsArray[index]) baconhep::TLHEWeight();
      baconhep::TLHEWeight *pWeight = (baconhep::TLHEWeight*)rWeightsArray[index];
      pWeight->weight = hLHEEvtInfoProduct->weights().at(iw).wgt;
      std::string pId = hLHEEvtInfoProduct->weights().at(iw).id;
      int id = -1;
      try {id = atoi(pId.c_str());} catch(int e) { std::cout << " ===> Error converting LHE to int" << std::endl;}
      pWeight->id     = id;
    }
  }

  // Get generator particles collection
  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByToken(fTokGenPar,hGenParProduct);
  assert(hGenParProduct.isValid());  
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  

  // loop over GEN particles
  std::vector<edm::Ptr<reco::GenParticle>> lMothers;
  TClonesArray &rGenParArray = *particlesArr;
  for (reco::GenParticleCollection::const_iterator itGenP = genParticles.begin(); itGenP!=genParticles.end(); ++itGenP) {

    // if not storing all gen particles, then do selective storing
    // Note: assuming Pythia8 status codes
    // Reference: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
    if(!fFillAll) {
      bool skip=true;

      if(itGenP->status()>20 && itGenP->status()<30)           { skip = false; }  // keep particles from hard scatter process
      if(abs(itGenP->pdgId())>= 5 && abs(itGenP->pdgId())<= 8) { skip = false; }  // keep b, t, b', t'
      if(abs(itGenP->pdgId())>=11 && abs(itGenP->pdgId())<=18) { skip = false; }  // keep leptons
      if(abs(itGenP->pdgId())>=23 && abs(itGenP->pdgId())<=39) { skip = false; }  // keep bosons except photons and gluons
      if(abs(itGenP->pdgId())>10000)                           { skip = false; }  // keep exotic particles

      // photons (e.g. for FSR/ISR) and u,d,c,s-quarks (e.g. hadronic boson decays) coming from a previously stored particle
      if(itGenP->pdgId()==22 || (abs(itGenP->pdgId())>0 && abs(itGenP->pdgId())<5)) {
        if(itGenP->numberOfMothers() > 0) {
          edm::Ptr<reco::GenParticle> lMomPtr = edm::refToPtr(itGenP->motherRef());
          for(unsigned int im=0; im < lMothers.size(); ++im) {
            if(lMothers[im] == lMomPtr) {
              skip = false;
              break;
            }
          }
        }
      }

      if(skip) continue;
    }

    // construct object and place in array
    assert(rGenParArray.GetEntries() < rGenParArray.GetSize());
    const int index = rGenParArray.GetEntries();
    new(rGenParArray[index]) baconhep::TGenParticle();
    baconhep::TGenParticle *pGenPart = (baconhep::TGenParticle*)rGenParArray[index];
    pGenPart->pdgId  = itGenP->pdgId();
    pGenPart->status = itGenP->status();
    pGenPart->pt     = itGenP->pt();
    pGenPart->eta    = itGenP->eta();
    pGenPart->phi    = itGenP->phi();
    pGenPart->y      = itGenP->rapidity();
    pGenPart->mass   = itGenP->mass();
    pGenPart->isPromptFinalState = itGenP->isPromptFinalState();
    pGenPart->isDirectPromptTauDecayProductFinalState = itGenP->isDirectPromptTauDecayProductFinalState();
    pGenPart->isHardProcess = itGenP->isHardProcess();
    pGenPart->fromHardProcessDecayed = itGenP->fromHardProcessDecayed();
    pGenPart->fromHardProcessFinalState = itGenP->fromHardProcessFinalState();
    pGenPart->isPromptDecayed = itGenP->isPromptDecayed();
    if(itGenP->numberOfMothers() >  0 ) {
      int lId = -2;
      edm::Ptr<reco::GenParticle> lMomPtr = edm::refToPtr(itGenP->motherRef()); 
      for(unsigned int im=0; im < lMothers.size(); ++im) { 
	if(lMothers[im] == lMomPtr) {
	  lId = im;  
	  break;
	}
      }
      pGenPart->parent = lId;
    }
    if(itGenP->numberOfMothers() == 0) { pGenPart->parent = -1; }
    edm::Ptr<reco::GenParticle> thePtr(hGenParProduct, itGenP - genParticles.begin());
    lMothers.emplace_back(thePtr);
  }
}
