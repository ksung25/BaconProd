#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerVertex::FillerVertex(const edm::ParameterSet &iConfig):
  fPVName       (iConfig.getUntrackedParameter<std::string>("edmName","offlinePrimaryVertices")),
  fMinNTracksFit(0),
  fMinNdof      (4),
  fMaxAbsZ      (24),
  fMaxRho       (2)
{}

//--------------------------------------------------------------------------------------------------
FillerVertex::~FillerVertex(){}

//--------------------------------------------------------------------------------------------------
const reco::Vertex* FillerVertex::fill(TClonesArray *array, int &nvtx, const edm::Event &iEvent)
{
  assert(array);
  
  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  const reco::Vertex* pv = &(*pvCol->begin());
  nvtx = 0;
    
  for(reco::VertexCollection::const_iterator itVtx = pvCol->begin(); itVtx!=pvCol->end(); ++itVtx) {    
    if(itVtx->isFake())                          continue;
    if(itVtx->tracksSize()     < fMinNTracksFit) continue;
    if(itVtx->ndof()           < fMinNdof)       continue;
    if(fabs(itVtx->z())        > fMaxAbsZ)       continue;
    if(itVtx->position().Rho() > fMaxRho)        continue;

    // vertices are sorted by sum{pT^2}, so the first one passing cuts
    // is taken as the event primary vertex
    if(nvtx==0) {
      pv = &(*itVtx);
    }
    nvtx++;

    // construct object and place in array
    TClonesArray &rArray = *array;
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) TVertex();
    TVertex *pVertex = (TVertex*)rArray[index];
    
    pVertex->nTracksFit = itVtx->tracksSize();
    pVertex->ndof       = itVtx->ndof();
    pVertex->chi2       = itVtx->chi2();
    pVertex->x          = itVtx->x();
    pVertex->y          = itVtx->y();
    pVertex->z          = itVtx->z();
  }
  return pv;
}
