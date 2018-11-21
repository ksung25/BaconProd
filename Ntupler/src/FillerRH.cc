#include "BaconProd/Ntupler/interface/FillerRH.hh"
#include "BaconAna/DataFormats/interface/TRHPart.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include <TClonesArray.h>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerRH::FillerRH(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fRHName      (iConfig.getUntrackedParameter<edm::InputTag>("edmRecHitName"))
{
  fTokRH = iC.consumes<HBHERecHitCollection>(fRHName);
}

//--------------------------------------------------------------------------------------------------
FillerRH::~FillerRH(){}

//--------------------------------------------------------------------------------------------------
void FillerRH::fill(TClonesArray *array,const edm::Event &iEvent,const edm::EventSetup &iSetup)
		   
{
  assert(array);
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry *gTmp0 = geoHandle->getSubdetectorGeometry(DetId::Hcal,1);
  fHCAL = dynamic_cast< const CaloSubdetectorGeometry* > (gTmp0);

  const HBHERecHitCollection *recHitHCAL   = 0;
  edm::Handle<HBHERecHitCollection> hRecHitHCAL;
  iEvent.getByToken(fTokRH,hRecHitHCAL);
  recHitHCAL = hRecHitHCAL.product();

  TClonesArray &rArray = *array;
  int pId = 0; 
  for(HBHERecHitCollection::const_iterator itRH = recHitHCAL->begin(); itRH != recHitHCAL->end(); itRH++) {
    pId++;
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TRHPart();
    baconhep::TRHPart *pRH = (baconhep::TRHPart*) rArray[index];
    std::shared_ptr<const CaloCellGeometry> pCell = fHCAL->getGeometry(itRH->detid());
    TLorentzVector pVec; 
    double pR = sqrt(pCell->getPosition().x()*pCell->getPosition().x() + pCell->getPosition().y()*pCell->getPosition().y() + pCell->getPosition().z()*pCell->getPosition().z());
    pVec.SetPxPyPzE( pCell->getPosition().x()/pR * itRH->energy(),       pCell->getPosition().y()/pR * itRH->energy(),       pCell->getPosition().z()/pR * itRH->energy(),itRH->energy()); 
    // Kinematics
    //==============================    
    pRH->pt     = pVec.Pt();
    pRH->eta    = pVec.Eta();
    pRH->phi    = pVec.Phi();
    pRH->rho    = pCell->getPosition().perp();
    pRH->ieta   = itRH->id().ieta();
    pRH->iphi   = itRH->id().iphi();
    pRH->depth  = itRH->id().depth();
    //pRH->timefalling   = itRH->timeFalling();
    pRH->energy = itRH->energy();
    pRH->time   = itRH->time();
    pRH->x      = pCell->getPosition().x();
    pRH->y      = pCell->getPosition().y();
    pRH->z      = pCell->getPosition().z();
  } 
}
