#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerGenInfo::FillerGenInfo(const edm::ParameterSet &iConfig):
  fGenEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmGenEventInfoName","generator")),
  fGenParName    (iConfig.getUntrackedParameter<std::string>("edmGenParticlesName","genParticles")),
  fFillAll       (iConfig.getUntrackedParameter<bool>("fillAllGen",true))
{}

//--------------------------------------------------------------------------------------------------
FillerGenInfo::~FillerGenInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerGenInfo::fill(TGenEventInfo *genEvtInfo, TClonesArray *array,      
                         const edm::Event &iEvent)
{
  assert(array);
  
  // Get generator event information
  edm::Handle<GenEventInfoProduct> hGenEvtInfoProduct;
  iEvent.getByLabel(fGenEvtInfoName,hGenEvtInfoProduct);
  assert(hGenEvtInfoProduct.isValid());

  const gen::PdfInfo *pdfInfo = (hGenEvtInfoProduct->hasPDF()) ? hGenEvtInfoProduct->pdf() : 0;
  genEvtInfo->id_1     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.first    : 0;
  genEvtInfo->id_2     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.second   : 0;
  genEvtInfo->x_1      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.first     : 0;
  genEvtInfo->x_2      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.second    : 0;
  genEvtInfo->scalePDF = hGenEvtInfoProduct->qScale();
  
  
  // Get generator particles collection
  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByLabel(fGenParName,hGenParProduct);
  assert(hGenParProduct.isValid());  
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  
  // loop over GEN particles
  std::vector<reco::GenParticle> lMothers;
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
      int pId = 0;
      reco::GenParticleRef lMom = itGenP->motherRef(); 
      for( std::vector<reco::GenParticle>::const_iterator itMomP = lMothers.begin(); itMomP != lMothers.end(); ++itMomP) { 
        // TODO: figure out how to match motherRef with GenParticle
	if(itMomP->pdgId() == lMom->pdgId() && itMomP->status() == lMom->status() &&
	   reco::deltaR(itMomP->eta(),itMomP->phi(),lMom->eta(),lMom->phi()) < 0.001) {
	  lId = pId;  
	  break;
	}
	pId++;
      }
      pGenPart->parent =  lId;
    }
    if(itGenP->numberOfMothers() == 0 ) pGenPart->parent = -1;
    lMothers.push_back(*itGenP);
  }
}
