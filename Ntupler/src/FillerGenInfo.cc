#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenWeight.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerGenInfo::FillerGenInfo(const edm::ParameterSet &iConfig):
  fGenEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmGenEventInfoName","generator")),
  fLHEEvtName    (iConfig.getUntrackedParameter<std::string>("edmLHEEventName"    ,"generator")),
  fGenParName    (iConfig.getUntrackedParameter<std::string>("edmGenParticlesName","genParticles")),
  fFillAll       (iConfig.getUntrackedParameter<bool>("fillAllGen",true))
{}

//--------------------------------------------------------------------------------------------------
FillerGenInfo::~FillerGenInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerGenInfo::fill(TGenEventInfo *genEvtInfo,TGenWeight *genWeight, TClonesArray *array,      
                         const edm::Event &iEvent,float &iXS)
{
  assert(array);
  edm::Handle<GenEventInfoProduct> hGenEvtInfoProduct;
  iEvent.getByLabel(fGenEvtInfoName,hGenEvtInfoProduct);
  assert(hGenEvtInfoProduct.isValid());

  const gen::PdfInfo *pdfInfo = (hGenEvtInfoProduct->hasPDF()) ? hGenEvtInfoProduct->pdf() : 0;
  genEvtInfo->id_1     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.first    : 0;
  genEvtInfo->id_2     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.second   : 0;
  genEvtInfo->x_1      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.first     : 0;
  genEvtInfo->x_2      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.second    : 0;
  genEvtInfo->scalePDF = hGenEvtInfoProduct->qScale();
  genEvtInfo->weight   = float(hGenEvtInfoProduct->weight());
  genEvtInfo->xs       = iXS;

  // Get LHE
  /*
  edm::Handle<LHEEventProduct> hLHEEvtProduct;
  iEvent.getByLabel(fLHEEvtName,hLHEEvtProduct);
  assert(hLHEEvtProduct.isValid());
  double baseWeight = (hLHEEvtProduct->weights()[0]).wgt;
  genWeight->scale11 = float(((hLHEEvtProduct->weights()[0]).wgt)/baseWeight);
  genWeight->scale12 = float(((hLHEEvtProduct->weights()[1]).wgt)/baseWeight);
  genWeight->scale10 = float(((hLHEEvtProduct->weights()[2]).wgt)/baseWeight);
  genWeight->scale21 = float(((hLHEEvtProduct->weights()[3]).wgt)/baseWeight);
  genWeight->scale22 = float(((hLHEEvtProduct->weights()[4]).wgt)/baseWeight);
  genWeight->scale20 = float(((hLHEEvtProduct->weights()[5]).wgt)/baseWeight);
  genWeight->scale01 = float(((hLHEEvtProduct->weights()[6]).wgt)/baseWeight);
  genWeight->scale02 = float(((hLHEEvtProduct->weights()[7]).wgt)/baseWeight);
  genWeight->scale00 = float(((hLHEEvtProduct->weights()[8]).wgt)/baseWeight);
  for(int i0 = 0; i0 < 100; i0++) genWeight->pdf[i0] = float(((hLHEEvtProduct->weights()[i0+9]).wgt)/baseWeight);
  */

  // Get generator particles collection
  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByLabel(fGenParName,hGenParProduct);
  assert(hGenParProduct.isValid());  
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  

  // loop over GEN particles
  std::vector<edm::Ptr<reco::GenParticle>> lMothers;
  TClonesArray &rArray = *array;
  for (reco::GenParticleCollection::const_iterator itGenP = genParticles.begin(); itGenP!=genParticles.end(); ++itGenP) {
    if((itGenP->status() == 1     || itGenP->status()      == 99)    &&   //Remove all Status 1 and 99 particles
       (fabs(itGenP->pdgId()) < 11 || fabs(itGenP->pdgId()) > 17)    &&   //Keep Letpons 
       !fFillAll) continue; //Require to not fill all particles
    if(!fFillAll && (itGenP->status() == 2 && itGenP->pdgId() == 21)) continue;
    if(!fFillAll && (itGenP->status() == 2 && itGenP->pdgId() == 22)) continue;
    if(!fFillAll && fabs(itGenP->pdgId()) > 50) continue;
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TGenParticle();
    baconhep::TGenParticle *pGenPart = (baconhep::TGenParticle*)rArray[index];
    pGenPart->pdgId  = itGenP->pdgId();
    pGenPart->status = itGenP->status();
    pGenPart->pt     = itGenP->pt();
    pGenPart->eta    = itGenP->eta();
    pGenPart->phi    = itGenP->phi();
    pGenPart->y      = itGenP->rapidity();
    pGenPart->mass   = itGenP->mass();
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
    if(itGenP->numberOfMothers() == 0 ) pGenPart->parent = -1;
    edm::Ptr<reco::GenParticle> thePtr(hGenParProduct, itGenP - genParticles.begin());
    lMothers.push_back(thePtr);
  }
}
