#include "BaconProd/Ntupler/interface/FillerPF.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TClonesArray.h>
#include <TMath.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerPF::FillerPF(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fPFName      (iConfig.getUntrackedParameter<std::string>("edmName","particleFlow")),
  fPVName      (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fAddDepthTime(iConfig.getUntrackedParameter<bool>("doAddDepthTime",false))
{
  fTokPFName       = iC.consumes<reco::PFCandidateCollection>(fPFName);
  fTokPackCandName = iC.consumes<pat::PackedCandidateCollection>(fPFName);
  fTokPVName       = iC.consumes<reco::VertexCollection>     (fPVName);
  if(fAddDepthTime) { 
    std::string lTokPFRecHitECAL = "particleFlowRecHitECAL";
    std::string lTokPFRecHitHCAL = "particleFlowRecHitHCAL";
    std::string lTokPFRecHitHO   = "particleFlowRecHitHO";
    fTokPFRecHitECAL = iC.consumes<reco::PFRecHitCollection>(lTokPFRecHitECAL);
    fTokPFRecHitHCAL = iC.consumes<reco::PFRecHitCollection>(lTokPFRecHitHCAL);
    fTokPFRecHitHO   = iC.consumes<reco::PFRecHitCollection>(lTokPFRecHitHO);
  }
}

//--------------------------------------------------------------------------------------------------
FillerPF::~FillerPF(){}

//--------------------------------------------------------------------------------------------------
void FillerPF::fill(TClonesArray *array,TClonesArray *iVtxCol,
		    const edm::Event &iEvent) 
		   
{
  assert(array);
  // Get PF collection
  edm::Handle<reco::PFCandidateCollection> hPFProduct;
  iEvent.getByToken(fTokPFName,hPFProduct);
  assert(hPFProduct.isValid());
  const reco::PFCandidateCollection *PFCol = hPFProduct.product();

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(fTokPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  const reco::PFRecHitCollection *pfRecHitECAL = 0;
  const reco::PFRecHitCollection *pfRecHitHCAL = 0;
  const reco::PFRecHitCollection *pfRecHitHO   = 0;
  reco::PFRecHitCollection pfRecHitAll;
  if(fAddDepthTime) { 
    //Load all of the stupid PF Rec Hits
    edm::Handle<reco::PFRecHitCollection> hPFRecHitECAL;
    iEvent.getByToken(fTokPFRecHitECAL,hPFRecHitECAL);
    assert(hPFRecHitECAL.isValid());
    pfRecHitECAL = hPFRecHitECAL.product();
    
    edm::Handle<reco::PFRecHitCollection> hPFRecHitHCAL;
    iEvent.getByToken(fTokPFRecHitHCAL,hPFRecHitHCAL);
    assert(hPFRecHitHCAL.isValid());
    pfRecHitHCAL = hPFRecHitHCAL.product();
    
    edm::Handle<reco::PFRecHitCollection> hPFRecHitHO;
    iEvent.getByToken(fTokPFRecHitHO,hPFRecHitHO);
    assert(hPFRecHitHO.isValid());
    pfRecHitHO = hPFRecHitHO.product();
   
    for(unsigned int i0 = 0; i0 < pfRecHitECAL ->size(); i0++) pfRecHitAll.push_back((*pfRecHitECAL) [i0]);
    for(unsigned int i0 = 0; i0 < pfRecHitHCAL ->size(); i0++) pfRecHitAll.push_back((*pfRecHitHCAL) [i0]);
    for(unsigned int i0 = 0; i0 < pfRecHitHO   ->size(); i0++) pfRecHitAll.push_back((*pfRecHitHO)   [i0]);
  }
  /*
  edm::Handle<reco::PFRecHitCollection> hPFRecHitHFEM;
  //iEvent.getByLabel(edm::InputTag("particleFlowClusterHFEM"),hPFRecHitHFEM);
  assert(hPFRecHitHFEM.isValid());
  const reco::PFRecHitCollection *pfRecHitHFEM = hPFRecHitHFEM.product();

  edm::Handle<reco::PFRecHitCollection> hPFRecHitHFHAD;
  //iEvent.getByLabel(edm::InputTag("particleFlowClusterHFHAD"),hPFRecHitHFHAD);
  assert(hPFRecHitHFHAD.isValid());
  const reco::PFRecHitCollection *pfRecHitHFHAD = hPFRecHitHFHAD.product();
  */


  //for(unsigned int i0 = 0; i0 < pfRecHitHFEM ->size(); i0++) pfRecHitAll.push_back((*pfRecHitHFEM) [i0]);
  //for(unsigned int i0 = 0; i0 < pfRecHitHFHAD->size(); i0++) pfRecHitAll.push_back((*pfRecHitHFHAD)[i0]);
  //std::cout << "=====> Size " << pfRecHitAll.size() << " -- " << pfRecHitECAL->size() << " - " << pfRecHitHCAL->size() << " -- " << std::endl;//pfRecHitHFHAD->size()  << std::endl;

  TClonesArray &rArray = *array;
  int pId = 0; 
  for(reco::PFCandidateCollection::const_iterator itPF = PFCol->begin(); itPF!=PFCol->end(); itPF++) {
    if(itPF->pt() < 0.1) continue;
    pId++;
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TPFPart();
    baconhep::TPFPart *pPF = (baconhep::TPFPart*)rArray[index];

    //
    // Kinematics
    //==============================    
    pPF->pt  = itPF->pt();
    pPF->eta = itPF->eta();
    pPF->phi = itPF->phi();
    pPF->m   = itPF->mass();
    pPF->e   = itPF->energy();
    pPF->q   = itPF->charge();
    pPF->pfType = itPF->particleId();
    //pPF->deltaP = itPF->deltaP();
    pPF->ecalE  = itPF->ecalEnergy();
    pPF->hcalE  = itPF->hcalEnergy();
    //
    // Depth & Timing Info
    //==============================
    if(fAddDepthTime) { 
      pPF->time  = timeDeltaR (&(*itPF),pfRecHitAll);
      pPF->depth = depthDeltaR(&(*itPF),pfRecHitAll);
    }
    //pPF->time  = time (&(*itPF));
    //pPF->depth = depth(&(*itPF));
    //float pDepth =  depth(&(*itPF));
    //
    // TrackInfo
    //==============================
    const reco::TrackRef& pfTrack = itPF->trackRef();
    if(!pfTrack.isNonnull()) continue;
    int    ndof     = pfTrack->ndof();
    double chi2     = pfTrack->chi2();
    pPF->trkChi2 = TMath::Prob(chi2,ndof);
    const reco::Vertex *closestVtx = 0;
    bool lFirst = true;
    for(reco::VertexCollection::const_iterator iV = pvCol->begin(); iV!=pvCol->end(); ++iV) {
      if(lFirst) { 
	if      ( itPF->trackRef().isNonnull()    ) pPF->dz = itPF->trackRef()   ->dz (iV->position());
	else if ( itPF->gsfTrackRef().isNonnull() ) pPF->dz = itPF->gsfTrackRef()->dz (iV->position());
	if      ( itPF->trackRef().isNonnull()    ) pPF->d0 = itPF->trackRef()   ->dxy(iV->position());
	else if ( itPF->gsfTrackRef().isNonnull() ) pPF->d0 = itPF->gsfTrackRef()->dxy(iV->position());
	lFirst = false;
      }
      if(iV->trackWeight(itPF->trackRef())>0) {
	closestVtx  = &(*iV);
	break;
      }
    }
    if(closestVtx == 0) continue;
    pPF->vtxChi2 = closestVtx->trackWeight(itPF->trackRef());
    int lId = -1;
    for(int i0 = 0; i0 < iVtxCol->GetEntries(); i0++) { 
      baconhep::TVertex* pVertex = (TVertex*)(*iVtxCol)[i0];
      if(fabs(closestVtx->x() - pVertex->x) + 
	 fabs(closestVtx->y() - pVertex->y) + 
	 fabs(closestVtx->z() - pVertex->z) > 0.0001) continue;
      lId = i0;
      break;
    }
    pPF->vtxId = lId;
  } 
}
void FillerPF::fillMiniAOD(TClonesArray *array,TClonesArray *iVtxCol,
			   const edm::Event &iEvent) 		   
{
  assert(array);
  // Get PF collection
  edm::Handle<pat::PackedCandidateCollection> hPFProduct;
  iEvent.getByToken(fTokPackCandName,hPFProduct);
  assert(hPFProduct.isValid());
  const pat::PackedCandidateCollection *PFCol = hPFProduct.product();

  TClonesArray &rArray = *array;
  int pId = 0; 
  for(pat::PackedCandidateCollection::const_iterator itPF = PFCol->begin(); itPF!=PFCol->end(); itPF++) {
    pId++;
    // construct object and place in array
    //assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TPFPart();
    baconhep::TPFPart *pPF = (baconhep::TPFPart*)rArray[index];

    //
    // Kinematics
    //==============================    
    pPF->pt      = itPF->pt();
    pPF->eta     = itPF->eta();
    pPF->phi     = itPF->phi();
    pPF->m       = itPF->mass();
    pPF->e       = itPF->energy();
    pPF->q       = itPF->charge();
    pPF->pfType  = itPF->pdgId();
    pPF->pup     = itPF->puppiWeight();
    if (itPF->hasTrackDetails()) {
      pPF->vtxChi2 = itPF->vertexChi2();
      //pPF->pup     = itPF->puppiWeight();
      if(itPF->charge() == 0) continue;
      
      const reco::Track & pseudoTrack =  itPF->pseudoTrack();
      pPF->trkChi2 = pseudoTrack.normalizedChi2();
      pPF->dz      = itPF->dz();
      pPF->d0      = itPF->dxy();
      pPF->d0Err   = itPF->dxyError();
      reco::Track::CovarianceMatrix myCov = pseudoTrack.covariance ();
      pPF->dptdpt    = catchInfsAndBound(myCov[0][0],0,-1,1);
      pPF->detadeta  = catchInfsAndBound(myCov[1][1],0,-1,0.01);
      pPF->dphidphi  = catchInfsAndBound(myCov[2][2],0,-1,0.1);
      pPF->dxydxy    = catchInfsAndBound(myCov[3][3],7.,-1,7); 
      pPF->dzdz      = catchInfsAndBound(myCov[4][4],6.5,-1,6.5); 
      pPF->dxydz     = catchInfsAndBound(myCov[3][4],6.,-6,6); 
      pPF->dphidxy   = catchInfs(myCov[2][3],-0.03); 
      pPF->dlambdadz = catchInfs(myCov[1][4],-0.03);  
    }
    else continue;
   
  }
}
float FillerPF::depthDeltaR(const reco::PFCandidate *iPF,const reco::PFRecHitCollection &iPFCol,double iDR) { 
   //Get Calo Depth of PF Clusters in a cylinder using deltar R
  float lEta     = iPF->positionAtECALEntrance().eta();
  float lPhi     = iPF->positionAtECALEntrance().phi();
  //float lRhoE    = iPF->positionAtECALEntrance().rho();
  float lTotRho  = 0; 
  float lTotE    = 0;
  for(unsigned int i0 = 0; i0 < iPFCol.size(); i0++) { 
    //iPFCol[i0].calculatePositionREP();
    float pEta = iPFCol[i0].position().eta();
    float pPhi = iPFCol[i0].position().phi();
    if(reco::deltaR(pEta,pPhi,lEta,lPhi) > iDR) continue;
    //std::cout << " --> " << iPFCol[i0].energy() << " -- " << iPFCol[i0].position().rho() << " -- " << lRhoE << " -- Diff - " << iPFCol[i0].position().rho()-lRhoE << std::endl;
    //lTotRho += (iPFCol[i0].position().rho()-lRhoE)*iPFCol[i0].energy();
    lTotE   += iPFCol[i0].energy();
  }
  if(lTotE == 0) return 0;
  return lTotRho/lTotE;
}
float FillerPF::timeDeltaR(const reco::PFCandidate *iPF,const reco::PFRecHitCollection &iPFCol,double iDR) { 
  //highest Energy Time in a cylinder using delta R 
  float lEta     = iPF->positionAtECALEntrance().eta();
  float lPhi     = iPF->positionAtECALEntrance().phi();
  float lMaxT    = -999; 
  float lMaxE    = -999;
  for(unsigned int i0 = 0; i0 < iPFCol.size(); i0++) { 
    //iPFCol[i0].calculatePositionREP();
    float pEta = iPFCol[i0].position().eta();
    float pPhi = iPFCol[i0].position().phi();
    if(reco::deltaR(pEta,pPhi,lEta,lPhi) > iDR) continue;
    double pEnergy = iPFCol[i0].energy();
    if(lMaxE > pEnergy) continue;
    lMaxE    = pEnergy;
    lMaxT    = iPFCol[i0].time();
    //std::cout << "===> E " << lMaxE << " - " << lMaxT << std::endl;
  }
  return lMaxT;
}
float FillerPF::depth(const reco::PFCandidate *iPF) { 
  float lTotRho  = 0; 
  float lTotE    = 0;
  //Get Calo Depth of PF Clusters
  float lRhoE    = iPF->positionAtECALEntrance().rho();
  assert(!iPF->elementsInBlocks().empty() );
  //std::cout << "S:" << std::endl;
  for(unsigned int i0 = 0; i0 < iPF->elementsInBlocks().size(); i0++ ) { 
    reco::PFBlockRef blockRef = iPF->elementsInBlocks()[i0].first;
    if(blockRef.isNull()) continue;
    const reco::PFBlock& block = *blockRef;
    const edm::OwnVector<reco::PFBlockElement>& elements = block.elements();
    for(unsigned int iEle=0; iEle< elements.size(); iEle++) {
      // Find the tracks in the block
      reco::PFClusterRef pCluster = elements[iEle].clusterRef();
      if(pCluster.isNull()) continue;
      lTotRho += (pCluster->position().rho() - lRhoE)*pCluster->energy();
      lTotE   += pCluster->energy();
      //std::cout << " --> " << pCluster->energy()  << " -- " << pCluster->position().rho() << " -- " << lRhoE << " -- Diff - " << pCluster->position().rho()-lRhoE << " -- " << pCluster->layer() << std::endl;
    }
  }
  if(lTotE == 0) return 0;
  //std::cout << "E:" << std::endl;
  return lTotRho/lTotE;
}
float FillerPF::time(const reco::PFCandidate *iPF) { 
  //Get Max Time of the recHits
  float lMaxTime = -999;
  float lMaxE    = -999;
  assert(!iPF->elementsInBlocks().empty() );
  for(unsigned int i0 = 0; i0 < iPF->elementsInBlocks().size(); i0++) { 
    reco::PFBlockRef blockRef = iPF->elementsInBlocks()[i0].first;
    const reco::PFBlock& block = *blockRef;
    const edm::OwnVector<reco::PFBlockElement>& elements = block.elements();
    for(unsigned int iEle=0; iEle<elements.size(); iEle++) {
      // Find the tracks in the block
      reco::PFClusterRef pCluster = elements[iEle].clusterRef();
      if(pCluster.isNull()) continue;
      const std::vector< reco::PFRecHitFraction > pRecHits = pCluster->recHitFractions();
      for(unsigned  int iRec = 0; iRec < pRecHits.size(); iRec++) { 
	reco::PFRecHitRef pRHRef = pRecHits[iRec].recHitRef();
	double pEnergy = pRHRef->energy() * pRecHits[iRec].fraction();
	if(lMaxE > pEnergy) continue;
	lMaxE    = pEnergy;
	lMaxTime = pRHRef->time();
      }
    }
  }
  return lMaxTime;
}
const float& FillerPF::catchInfs(const float& in,const float& replace_value){
  if(in==in){
    if(std::isinf(in))
      return replace_value;
    else if(in < -1e32 || in > 1e32)
      return replace_value;
    return in;
  }
  return replace_value;
}
float FillerPF::catchInfsAndBound(const float& in,const float& replace_value, const float& lowerbound, const float& upperbound){
  float withoutinfs=catchInfs(in,replace_value);
  if(withoutinfs<lowerbound) return lowerbound;
  if(withoutinfs>upperbound) return upperbound;
  return withoutinfs;
}
