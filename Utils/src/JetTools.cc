#include "BaconProd/Utils/interface/JetTools.hh"
#include "BaconProd/Utils/interface/QjetsPlugin.h"
#include "BaconProd/Utils/interface/Qjets.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
double JetTools::beta(const reco::PFJet &jet, const reco::Vertex &pv, const double dzCut)
{
  double pt_jets=0, pt_jets_tot=0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    const reco::TrackRef       track  = pfcand->trackRef();
    if(track.isNull()) continue;
    
    pt_jets_tot += track->pt();
    
    if(fabs(track->dz(pv.position())) < dzCut) {
      pt_jets += track->pt();
    }
  }
  
  return (pt_jets_tot>0) ? pt_jets/pt_jets_tot : 0;
}

//--------------------------------------------------------------------------------------------------
double JetTools::betaStar(const reco::PFJet &jet, const reco::Vertex &pv, const reco::VertexCollection *pvCol, const double dzCut)
{
  double pileup=0, total=0;
  
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    const reco::TrackRef       track  = pfcand->trackRef();
    if(track.isNull()) continue;
    total += track->pt();
    
    double dzPV = fabs(track->dz(pv.position()));
    if(dzPV <= dzCut) continue;
    
    double dzMin = dzPV;
    for(reco::VertexCollection::const_iterator itVtx = pvCol->begin(); itVtx!=pvCol->end(); ++itVtx) {
      if(itVtx->ndof() < 4 || (pv.position() - itVtx->position()).R() < 0.02) continue;
      dzMin = TMath::Min(dzMin, fabs(track->dz(itVtx->position())));
    }
    if(dzMin < dzCut) pileup += track->pt();
  }
  if(total==0) total=1;
  
  return pileup/total;
}

//--------------------------------------------------------------------------------------------------
double JetTools::dRMean(const reco::PFJet &jet, const int pfType)
{
  double drmean=0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    if(pfType!=-1 && pfcand->particleId() != pfType) continue;
    
    double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand->eta(),pfcand->phi());    
    drmean += dr*(pfcand->pt())/(jet.pt());
  }
  
  return drmean;
}

//--------------------------------------------------------------------------------------------------
double JetTools::dR2Mean(const reco::PFJet &jet, const int pfType)
{
  double dr2mean=0;
  double sumpt2=0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    if(pfType!=-1 && pfcand->particleId() != pfType) continue;
    
    sumpt2 += pfcand->pt() * pfcand->pt();
    
    double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand->eta(),pfcand->phi());    
    dr2mean += dr*dr*(pfcand->pt() * pfcand->pt());
  }
  dr2mean/=sumpt2;
  
  return dr2mean;
}

//--------------------------------------------------------------------------------------------------
double JetTools::frac(const reco::PFJet &jet, const double dRMax, const int pfType)
{
  const double dRMin = dRMax - 0.1;
  
  double fraction = 0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    if(pfType!=-1 && pfcand->particleId() != pfType) continue;
    
    double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand->eta(),pfcand->phi());    
    if(dr > dRMax) continue;
    if(dr < dRMin) continue;
    
    fraction += pfcand->pt() / jet.pt();
  }
  
  return fraction;
}

//--------------------------------------------------------------------------------------------------
double JetTools::jetDz(const reco::PFJet &jet, const reco::Vertex &pv)
{
  // Assumes constituents are stored by descending pT
  double dz=-1000;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    const reco::TrackRef       track  = pfcand->trackRef();
    if(track.isNull()) continue;
    dz = track->dz(pv.position());
    break;
  }
  
  return dz;
}
//--------------------------------------------------------------------------------------------------
double JetTools::jetD0(const reco::PFJet &jet, const reco::Vertex &pv)
{
  // Assumes constituents are stored by descending pT
  double d0=-1000;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    const reco::TrackRef       track  = pfcand->trackRef();
    if(track.isNull()) continue;
    d0 = -track->dxy(pv.position());
    break;
  }
  
  return d0;
}

//--------------------------------------------------------------------------------------------------
double JetTools::jetWidth(const reco::PFJet &jet, const int varType, const int pfType)
{
  double ptD=0, sumPt=0, sumPt2=0;
  TMatrixDSym covMatrix(2); covMatrix=0.;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    if(pfType!=-1 && pfcand->particleId() != pfType) continue;
    
    double dEta = jet.eta() - pfcand->eta();
    double dPhi = reco::deltaPhi(jet.phi(), pfcand->phi());
    
    covMatrix(0,0) += pfcand->pt() * pfcand->pt() * dEta * dEta;
    covMatrix(0,1) += pfcand->pt() * pfcand->pt() * dEta * dPhi;
    covMatrix(1,1) += pfcand->pt() * pfcand->pt() * dPhi * dPhi;
    ptD            += pfcand->pt() * pfcand->pt();
    sumPt          += pfcand->pt();
    sumPt2         += pfcand->pt() * pfcand->pt();
  }
  covMatrix(0,0) /= sumPt2;
  covMatrix(0,1) /= sumPt2;
  covMatrix(1,1) /= sumPt2;
  covMatrix(1,0) = covMatrix(0,1);
  
  ptD /= sqrt(ptD);
  ptD /= sumPt;
  
  double etaW = sqrt(covMatrix(0,0));
  double phiW = sqrt(covMatrix(1,1));
  double jetW = 0.5*(etaW+phiW);
  
  TVectorD eigVals(2);
  eigVals = TMatrixDSymEigen(covMatrix).GetEigenValues();
  double majW = sqrt(fabs(eigVals(0)));
  double minW = sqrt(fabs(eigVals(1)));
  
  if     (varType==1) { return majW; }
  else if(varType==2) { return minW; }
  else if(varType==3) { return etaW; }
  else if(varType==4) { return phiW; }
  else if(varType==5) { return jetW; }  

  return ptD;
}
//--------------------------------------------------------------------------------------------------
double JetTools::jetWidth(fastjet::PseudoJet &jet, const int varType) { 
  double ptD=0, sumPt=0, sumPt2=0;
  TMatrixDSym covMatrix(2); covMatrix=0.;
  const unsigned int nPFCands = jet.constituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    fastjet::PseudoJet pfcand = jet.constituents()[ipf];
    
    double dEta = jet.eta() - pfcand.eta();
    double dPhi = reco::deltaPhi(jet.phi(), pfcand.phi());
    
    covMatrix(0,0) += pfcand.pt() * pfcand.pt() * dEta * dEta;
    covMatrix(0,1) += pfcand.pt() * pfcand.pt() * dEta * dPhi;
    covMatrix(1,1) += pfcand.pt() * pfcand.pt() * dPhi * dPhi;
    ptD            += pfcand.pt() * pfcand.pt();
    sumPt          += pfcand.pt();
    sumPt2         += pfcand.pt() * pfcand.pt();
  }
  covMatrix(0,0) /= sumPt2;
  covMatrix(0,1) /= sumPt2;
  covMatrix(1,1) /= sumPt2;
  covMatrix(1,0) = covMatrix(0,1);
  
  ptD /= sqrt(ptD);
  ptD /= sumPt;
  
  double etaW = sqrt(covMatrix(0,0));
  double phiW = sqrt(covMatrix(1,1));
  double jetW = 0.5*(etaW+phiW);
  
  TVectorD eigVals(2);
  eigVals = TMatrixDSymEigen(covMatrix).GetEigenValues();
  double majW = sqrt(fabs(eigVals(0)));
  double minW = sqrt(fabs(eigVals(1)));

  if     (varType==1) { return majW; }
  else if(varType==2) { return minW; }
  else if(varType==3) { return etaW; }
  else if(varType==4) { return phiW; }
  else if(varType==5) { return jetW; }
  return ptD;
}
//--------------------------------------------------------------------------------------------------
bool JetTools::passPFLooseID(const reco::PFJet &jet)
{
  if(jet.energy() == 0) return false;
  if(jet.neutralHadronEnergy() / jet.energy() > 0.99) return false;
  if(jet.neutralEmEnergy() / jet.energy()     > 0.99) return false;
  if(jet.getPFConstituents().size()           < 2)    return false;
  if(fabs(jet.eta())<2.4) {
    if(jet.chargedHadronEnergy() / jet.energy() <= 0)   return false;
    if(jet.chargedEmEnergy() / jet.energy()     > 0.99) return false;
    if(jet.chargedMultiplicity()                < 1)    return false;
  }
  
  return true;
}
//--------------------------------------------------------------------------------------------------
double JetTools::jetCharge(const reco::PFJet &jet, const bool iSquare)
{
  // Assumes constituents are stored by descending pT
  double charge=0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    const reco::TrackRef       track  = pfcand->trackRef();
    if(track.isNull()) continue;
    if(!iSquare) charge += pfcand->pt()*track->charge();
    if( iSquare) charge += pfcand->pt()*track->charge()*pfcand->pt();
  }
  if(iSquare) return charge/jet.pt()/jet.pt();
  else        return charge/jet.pt();
}
double JetTools::jetCharge(const reco::PFJet &jet, const double kappa)
{
  // Assumes constituents are stored by descending pT
  double charge=0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    const reco::TrackRef       track  = pfcand->trackRef();
    if(track.isNull()) continue;
    charge += pow(pfcand->pt(),kappa)*track->charge();
  }
  double denom = pow(jet.pt(),kappa);
  return charge/denom;
}
//--------------------------------------------------------------------------------------------------
double* JetTools::subJetBTag(const reco::PFJet &jet,reco::JetTagCollection &subJetVal,double iConeSize )
{
  double* vals = new double[2];
  double pt[2]; 
  int subjetIndex[2];
  vals[0] = vals[1] = -10; 
  pt[0]   = pt[1]   =   0; 
  subjetIndex[0] = subjetIndex[1] = -1;
  for (unsigned int i = 0; i != subJetVal.size(); ++i) {
    double pDR = reco::deltaR(subJetVal[i].first->eta(),subJetVal[i].first->phi(),jet.eta(),jet.phi());
    if(pDR  > iConeSize) continue;    
    double sjpt = subJetVal[i].first->pt();
    if(sjpt > pt[0]) {
      pt[1] = pt[0]; subjetIndex[1] = subjetIndex[0];
      pt[0] = sjpt;  subjetIndex[0] = i;
    } else if(sjpt > pt[1]) {
      pt[1] = sjpt;  subjetIndex[1] = i;
    }
  }
  if(subjetIndex[0]>=0) {
    vals[0] = subJetVal[subjetIndex[0]].second;
  }
  if(subjetIndex[1]>=0) {
    vals[1] = subJetVal[subjetIndex[1]].second;
  }  
  return vals;
}
//--------------------------------------------------------------------------------------------------
double* JetTools::subJetQG(const reco::PFJet &jet,edm::Handle<reco::PFJetCollection> &subJets,const edm::ValueMap<float> iQGLikelihood,double iConeSize)
{
  double* vals = new double[4];
  double pt[2];
  int subjetIndex[2];
  vals[0] = vals[1] = vals[2] = vals[3] = -10;
  pt[0]   = pt[1]   = 0;
  subjetIndex[0] = subjetIndex[1] = -1;
  for (unsigned int i = 0; i != subJets->size(); ++i) {
    double pDR = reco::deltaR((*subJets)[i].eta(),(*subJets)[i].phi(),jet.eta(),jet.phi());
    if(pDR  > iConeSize) continue;
    double sjpt = (*subJets)[i].pt();
    if(sjpt > pt[0]) {
      pt[1] = pt[0]; subjetIndex[1] = subjetIndex[0];
      pt[0] = sjpt;  subjetIndex[0] = i;
    } else if(sjpt > pt[1]) {
      pt[1] = sjpt;  subjetIndex[1] = i;
    }
  }
  if(subjetIndex[0]>=0) {
    edm::RefToBase<reco::Jet> jetRef0(edm::Ref<reco::PFJetCollection>(subJets,subjetIndex[0]));
    vals[0] = (iQGLikelihood)[jetRef0];
    vals[2] = jetCharge((*subJets)[subjetIndex[0]],true);
  }
  if(subjetIndex[1]>=0) {
    edm::RefToBase<reco::Jet> jetRef1(edm::Ref<reco::PFJetCollection>(subJets,subjetIndex[1]));
    vals[1] = (iQGLikelihood)[jetRef1];
    vals[3] = jetCharge((*subJets)[subjetIndex[1]],true);
  }
  return vals;  
}
//--------------------------------------------------------------------------------------------------
TVector2 JetTools::jetPull(const reco::PFJet &jet, const int type)
{
  double dYSum=0, dPhiSum=0;
  const unsigned int nPFCands = jet.getPFConstituents().size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
    
    double pt_i=0, y_i=0, phi_i=0;
    if(type==0) {  // PF jet pull
      pt_i  = pfcand->pt();
      y_i   = pfcand->rapidity();
      phi_i = pfcand->phi();
    
    } else if(type==1) {  // charged jet pull component
      if(pfcand->charge()!=0) {
        pt_i  = pfcand->pt();
        y_i   = pfcand->rapidity();
        phi_i = pfcand->phi();  
      }
          
    } else if(type==2) {  // neutral jet pull component
      if(pfcand->charge()==0) {
        pt_i  = pfcand->pt();
        y_i   = pfcand->rapidity();
        phi_i = pfcand->phi();  
      }
      
    } else {
      assert(0);
    }
    
    double dY     = y_i - jet.rapidity();
    double dPhi   = reco::deltaPhi(phi_i,jet.phi());
    double weight = pt_i*sqrt(dY*dY + dPhi*dPhi);
    dYSum   += weight*dY;
    dPhiSum += weight*dPhi;
  }
  
  return TVector2(dYSum/jet.pt(), dPhiSum/jet.pt());

}

//--------------------------------------------------------------------------------------------------
double JetTools::jetPullAngle(const reco::PFJet &jet ,edm::Handle<reco::PFJetCollection> &subJets,double iConeSize)
{
  const reco::PFJet *subjet0 = 0;
  const reco::PFJet *subjet1 = 0;
  for(unsigned int i0 = 0; i0 < subJets->size(); i0++) {
    const reco::PFJet *sj = &(*subJets)[i0];
    
    double pDR = reco::deltaR(sj->eta(),sj->phi(),jet.eta(),jet.phi());
    if(pDR > iConeSize) continue;
    
    if(!subjet0 || subjet0->pt() < sj->pt()) {
      subjet1 = subjet0;
      subjet0 = sj;
            
    } else if(!subjet1 || subjet1->pt() < sj->pt()) {
      subjet1 = sj;
    }
  }
  if(subjet0 == 0 || subjet1 == 0) return -20;
  
  // work in dy-dphi space of subjet0
  TVector2 lPull = jetPull(*subjet0);
  TVector2 lJet(subjet1->rapidity() - subjet0->rapidity(), reco::deltaPhi(subjet1->phi(), subjet0->phi()));
  double lThetaP = lPull.DeltaPhi(lJet);
  return lThetaP;
}
//--------------------------------------------------------------------------------------------------
float JetTools::findRMS( std::vector<float> &iQJetMass) {
  float total = 0.;
  float ctr = 0.;
  for (unsigned int i = 0; i < iQJetMass.size(); i++){
    total = total + iQJetMass[i];
    ctr++;
  }
  float mean =  total/ctr;
    
  float totalsquared = 0.;
  for (unsigned int i = 0; i < iQJetMass.size(); i++){
    totalsquared += (iQJetMass[i] - mean)*(iQJetMass[i] - mean) ;
  }
  float RMS = sqrt( totalsquared/ctr );
  return RMS;
}
//--------------------------------------------------------------------------------------------------
float JetTools::findMean( std::vector<float> &iQJetMass ){
  float total = 0.;
  float ctr = 0.;
  for (unsigned int i = 0; i < iQJetMass.size(); i++){
    total = total + iQJetMass[i];
    ctr++;
  }
  return total/ctr;
}
//-------------------------------------------------------------------------------------P-------------
double JetTools::qJetVolatility(std::vector <fastjet::PseudoJet> &iConstits, int iQJetsN, int iSeed){
  std::vector< float > qjetmasses;
  double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);          
  for(unsigned int ii = 0 ; ii < (unsigned int) iQJetsN ; ii++){
    QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
    qjet_plugin.SetRandSeed(iSeed+ii); // new feature in Qjets to set the random seed
    fastjet::JetDefinition qjet_def(&qjet_plugin);
    fastjet::ClusterSequence qjet_seq(iConstits, qjet_def);
    std::vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
    if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }
  }
  // find RMS of a vector
  float qjetsRMS = findRMS( qjetmasses );
  // find mean of a vector
  float qjetsMean = findMean( qjetmasses );
  float qjetsVolatility = qjetsRMS/qjetsMean;
  return qjetsVolatility;
}


