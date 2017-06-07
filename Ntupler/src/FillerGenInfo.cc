#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenVtx.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include <TLorentzVector.h>
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerGenInfo::FillerGenInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fGenEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmGenEventInfoName","generator")),
  fLHEEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmLHEEventInfoName","externalLHEProducer")),
  fGenParName    (iConfig.getUntrackedParameter<std::string>("edmGenParticlesName","genParticles")),
  fPackGenParName(iConfig.getUntrackedParameter<std::string>("edmGenPackParticlesName","packedGenParticles")),
  fHepMCProduct  (iConfig.getUntrackedParameter<std::string>("edmHepMCProduct"    ,"genHepMCInfo")),
  fFillAll       (iConfig.getUntrackedParameter<bool>("fillAllGen",false)),
  fFillLHEWeights(iConfig.getUntrackedParameter<bool>("fillLHEWeights",false))
{
  fTokGenEvent     = iC.mayConsume<GenEventInfoProduct>        ( fGenEvtInfoName);
  fTokGenPar       = iC.mayConsume<reco::GenParticleCollection>( fGenParName    );
  fTokGenPackPar   = iC.mayConsume<pat::PackedGenParticleCollection>( fPackGenParName    );
  fTokLHEEventInfo = iC.mayConsume<LHEEventProduct>            ( fLHEEvtInfoName );
}

//--------------------------------------------------------------------------------------------------
FillerGenInfo::~FillerGenInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerGenInfo::fill(TGenEventInfo *genEvtInfo, TClonesArray *particlesArr, TClonesArray *vtxArr, TClonesArray *weightsArr,
                         const edm::Event &iEvent, const float iXS) {
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
  /*
  std::vector<simPrimaryVertex> simpv;
  edm::Handle<edm::HepMCProduct> hHepMCProduct;
  iEvent.getByToken(fHepMCProduct,hHepMCProduct);
  const HepMC::GenEvent* evt = hHepMCProduct->GetEvent();
  for( HepMC::GenEvent::vertex_const_iterator vitr= evt->vertices_begin(); vitr != evt->vertices_end(); ++vitr ) { 
    HepMC::FourVector pos = (*vitr)->position();
    bool hasMotherVertex=false;
    for ( HepMC::GenVertex::particle_iterator mother  = (*vitr)->particles_begin(HepMC::parents); mother != (*vitr)->particles_end(HepMC::parents); ++mother ) {
      HepMC::GenVertex * mv=(*mother)->production_vertex();
      if (!mv) continue;
      hasMotherVertex=true;
      break; //if verbose_, print all particles of gen vertices
    }
    if (hasMotherVertex) {continue;}
    // could be a new vertex, check  all primaries found so far to avoid multiple entries
    const double mm=0.1;
    simPrimaryVertex sv(pos.x()*mm,pos.y()*mm,pos.z()*mm);
    simPrimaryVertex *vp=NULL;  // will become non-NULL if a vertex is found and then point to it
    for (std::vector<simPrimaryVertex>::iterator v0=simpv.begin(); v0!=simpv.end(); ++v0) {
      if ( (fabs(sv.x-v0->x)<1e-5) && (fabs(sv.y-v0->y)<1e-5) && (fabs(sv.z-v0->z)<1e-5) ) {
	vp=&(*v0);
	break;
      }
    }
    if (!vp) simpv.push_back(sv);
    if (!vp) { 
      std::cout << "------------  Gen Vertex  : " << sv.x << " -- " << sv.y << " -- " << sv.z << " -- " << (*vitr)->id() << std::endl;
      for( const HepMC::GenVertex::particle_iterator pPart = (*vitr)->particles_begin(HepMC::parents);  pPart != (*vitr)->particles_end(HepMC::parents); pPart++) {   
	std::cout << " particles in : " << pPart->pdg_id() << " -- " << pPart->status() << " -- " << pPart->pt() << " -- " << pPart->eta() << " -- " << pPart->phi() << std::endl; 
      }
      for( const HepMC::GenVertex::particle_iterator pPart = (*vitr)->particles_begin(HepMC::children); pPart != (*vitr)->particles_end(HepMC::children); pPart++) {   
	std::cout << " particles out : " << pPart->pdg_id() << " -- " << pPart->status() << " -- " << pPart->pt() << " -- " << pPart->eta() << " -- " << pPart->phi() << std::endl; 
      }
    }
  }
  */

  // loop over GEN particles
  std::vector<edm::Ptr<reco::GenParticle> > lMothers;
  TClonesArray &rGenParArray = *particlesArr;
  TClonesArray &rVtxArray    = *vtxArr;
  std::cout <<" ===> size : " << genParticles.size() << std::endl;
  std::vector<std::pair<double,int> > lVtxs;
  std::vector<TLorentzVector>         lVecs;
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
    int lMFlav       = 0;
    if(itGenP->numberOfMothers() >  0 ) {
      int lId = -2;
      edm::Ptr<reco::GenParticle> lMomPtr = edm::refToPtr(itGenP->motherRef()); 
      for(unsigned int im=0; im < lMothers.size(); ++im) { 
	if(lMothers[im] == lMomPtr) {
	  lId = im;  
	  lMFlav = flavor(lMomPtr->pdgId());
	  break;
	}
      }
      pGenPart->parent = lId;
    }
    if(itGenP->numberOfMothers() == 0) { pGenPart->parent = -1; }
    int lVId = -2;
    for(unsigned int iv=0; iv < lVtxs.size(); ++iv) { 
      if(itGenP->vx() == lVtxs[iv].first) {
	lVId = iv;
	TLorentzVector pGenVec; pGenVec.SetPtEtaPhiM(itGenP->pt(),itGenP->eta(),itGenP->phi(),itGenP->mass());
	if(lVtxs[iv].second < lMFlav) lVtxs[iv].second = lMFlav;
	lVecs[iv] += pGenVec;
	break;
      }
    }
    if(lVId == -2) { 
      std::pair<double,int> pVtx(itGenP->vx(),lMFlav);
      lVtxs.push_back(pVtx);
      TLorentzVector pGenVec; pGenVec.SetPtEtaPhiM(itGenP->pt(),itGenP->eta(),itGenP->phi(),itGenP->mass());
      lVecs.push_back(pGenVec);
      const int pVIndex = rVtxArray.GetEntries();
      new(rVtxArray[pVIndex]) baconhep::TGenVtx();
      baconhep::TGenVtx *pGenVtx = (baconhep::TGenVtx*)rVtxArray[pVIndex];
      pGenVtx->vx = itGenP->vx();
      pGenVtx->vy = itGenP->vy();
      pGenVtx->vz = itGenP->vz();
      pGenVtx->pdgId = lMFlav;
      lVId = pVIndex;
    }
    pGenPart->vtxId = lVId;
    std::cout << rGenParArray.GetEntries()-1  << " -pdg--> " << itGenP->pdgId() << " -status- " << itGenP->status() << " -pt- " << itGenP->pt() << " -vx- " << itGenP->vx() << "  -vy- " << itGenP->vy() << " -vz- " << itGenP->vz() << " -parent- " << pGenPart->parent << " -- " << lVId << std::endl;
    edm::Ptr<reco::GenParticle> thePtr(hGenParProduct, itGenP - genParticles.begin());
    lMothers.emplace_back(thePtr);
  }
  double lVX0=0,lVY0,lVZ0=0;
  baconhep::TGenVtx *lGenVtx = (baconhep::TGenVtx*)rVtxArray[1];
  lVX0 =  lGenVtx->vx;
  lVY0 =  lGenVtx->vy;
  lVZ0 =  lGenVtx->vz;
  std::cout << "vertices: " << lVtxs.size() << std::endl;
  std::cout << " ---> " << lVtxs.size() << " -- " << lVecs.size() << " -- " << rVtxArray.GetEntries() << std::endl;
  for(int i0 = 0; i0 < rVtxArray.GetEntries(); i0++) { 
      baconhep::TGenVtx *pGenVtx = (baconhep::TGenVtx*)rVtxArray[i0];
      pGenVtx->pdgId = lVtxs[i0].second;
      pGenVtx->pt    = lVecs[i0].Pt();
      pGenVtx->eta   = lVecs[i0].Eta();
      pGenVtx->phi   = lVecs[i0].Phi();
      pGenVtx->mass  = lVecs[i0].M();
      pGenVtx->vx    = pGenVtx->vx-lVX0;
      pGenVtx->vy    = pGenVtx->vy-lVY0;
      pGenVtx->vz    = pGenVtx->vz-lVZ0;
      std::cout << "===> Vertex " << i0 << " -- " << sqrt((pGenVtx->vx)*(pGenVtx->vx)+(pGenVtx->vy)*(pGenVtx->vy))  << " -flav- " << lVtxs[i0].second << std::endl;
  }
  lGenVtx->pdgId = 0;
  edm::Handle<pat::PackedGenParticleCollection> hGenPackParProduct;
  iEvent.getByToken(fTokGenPackPar,hGenPackParProduct);
  assert(hGenPackParProduct.isValid());  
  const pat::PackedGenParticleCollection packedGenParticles = *(hGenPackParProduct.product());  
  //std::vector<edm::Ptr<pat::PackedGenParticle> > lPatMothers;
  std::cout <<" ===> packed size : " << packedGenParticles.size() << std::endl;
  for (pat::PackedGenParticleCollection::const_iterator itGenP = packedGenParticles.begin(); itGenP!= packedGenParticles.end(); ++itGenP) {
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
    pGenPart->vtxId  = 0;
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
    //std::cout << rGenParArray.GetEntries()-1  << " -pdg--> " << itGenP->pdgId() << " -status- " << itGenP->status() << " -pt- " << itGenP->pt() << " -vx- " << itGenP->vx() << "  -vy- " << itGenP->vy() << " -vz- " << itGenP->vz() << " -parent- " << pGenPart->parent << std::endl;
    edm::Ptr<reco::GenParticle> thePtr = (edm::Ptr<reco::GenParticle>) itGenP->sourceCandidatePtr(0);//edm::refToPtr(itGenP->masterRef()); 
    lMothers.emplace_back(thePtr);
  }
  for(int i0 = 0; i0 < rGenParArray.GetEntries(); i0++) { 
    baconhep::TGenParticle *pGenPart = (baconhep::TGenParticle*)rGenParArray[i0];
    baconhep::TGenVtx      *pVtx     = (baconhep::TGenVtx*)     rVtxArray[pGenPart->vtxId];
    pGenPart->d0      = sqrt(pVtx->vx*pVtx->vx+pVtx->vy*pVtx->vy);
    pGenPart->vtxFlav = pVtx->pdgId;
  }
}
int FillerGenInfo::flavor(int iId) { 
  int rid = abs(iId) % 10000;
  int flavor = 0;
  if (rid > 999) {
    flavor = rid / 1000;
  } else {
    int flavor1 = rid / 100;
    int flavor2 = (rid % 100) / 10;
    flavor = flavor1;
    if (flavor1 == flavor2 && flavor1 < 4) {
      flavor = 1;
    }
    if (iId == 130) flavor = 3;
  }
  return flavor;
}

