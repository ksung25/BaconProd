#include "BaconProd/Utils/interface/JetTools.hh"
#include "BaconProd/Utils/interface/QjetsPlugin.h"
#include "BaconProd/Utils/interface/Qjets.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      pt_jets += packCand->pt();
      if(fabs(packCand->dz(pv.position())) < dzCut) pt_jets_tot += packCand->pt();
    } else { 
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
      const reco::TrackRef       track  = pfcand->trackRef();
      if(track.isNull()) continue;
      pt_jets_tot += track->pt();
      if(fabs(track->dz(pv.position())) < dzCut) {
	pt_jets += track->pt();
      }
    }
  }
  return (pt_jets_tot>0) ? pt_jets/pt_jets_tot : 0;
}

double JetTools::beta(const pat::Jet &jet, const reco::Vertex &pv, const double dzCut)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  double pt_jets=0, pt_jets_tot=0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfcand.charge()==0) continue; 

    pt_jets_tot += pfcand.pt();

    if(fabs(pfcand.dz(pv.position())) < dzCut) {
      pt_jets += pfcand.pt();
    }
  }

  return (pt_jets_tot>0) ? pt_jets/pt_jets_tot : 0;
}

//--------------------------------------------------------------------------------------------------
double JetTools::betaStar(const reco::PFJet &jet, const reco::Vertex &pv, const reco::VertexCollection *pvCol, const double dzCut)
{
  double pileup=0, total=0;

  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) {
      total += packCand->pt();
      if(packCand->charge()==0) continue;
      total += packCand->pt();
      double dzPV = fabs(packCand->dz(pv.position()));
      if(dzPV <= dzCut) continue;
      double dzMin = dzPV;
      for(reco::VertexCollection::const_iterator itVtx = pvCol->begin(); itVtx!=pvCol->end(); ++itVtx) {
	if(itVtx->ndof() < 4 || (pv.position() - itVtx->position()).R() < 0.02) continue;
	dzMin = TMath::Min(dzMin, (double)fabs(packCand->dz(pv.position())));
      }
      if(dzMin < dzCut) pileup += packCand->pt();
    } else {  
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
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
  }
  return pileup/total;
}

double JetTools::betaStar(const pat::Jet &jet, const reco::Vertex &pv, const reco::VertexCollection *pvCol, const double dzCut)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  double pileup=0, total=0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfcand.charge()==0) continue;
    total += pfcand.pt();

    double dzPV = fabs(pfcand.dz(pv.position()));
    if(dzPV <= dzCut) continue;

    double dzMin = dzPV;
    for(reco::VertexCollection::const_iterator itVtx = pvCol->begin(); itVtx!=pvCol->end(); ++itVtx) {
      if(itVtx->ndof() < 4 || (pv.position() - itVtx->position()).R() < 0.02) continue;
      dzMin = TMath::Min(dzMin, (double)fabs(pfcand.dz(pv.position())));
    }
    if(dzMin < dzCut) pileup += pfcand.pt();
  }
  if(total==0) total=1;

  return pileup/total;
}

//--------------------------------------------------------------------------------------------------
double JetTools::dRMean(const reco::PFJet &jet, const int pfType)
{
  double drmean=0;
  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      if(pfType!=-1 && packCand->pdgId() != pdgid) continue;
      double dr = reco::deltaR(jet.eta(),jet.phi(),packCand->eta(),packCand->phi());    
      drmean += dr*(packCand->pt())/(jet.pt());
    } else { 
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
      if(pfType!=-1 && pfcand->particleId() != pfType) continue;
      double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand->eta(),pfcand->phi());    
      drmean += dr*(pfcand->pt())/(jet.pt());
    }
  }
  return drmean;
}

double JetTools::dRMean(const pat::Jet &jet, const int pfType)
{
  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  double drmean=0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfType!=-1 && pfcand.pdgId() != pdgid) continue;

    double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand.eta(),pfcand.phi());
    drmean += dr*(pfcand.pt())/(jet.pt());
  }

  return drmean;
}

//--------------------------------------------------------------------------------------------------
double JetTools::dR2Mean(const reco::PFJet &jet, const int pfType)
{
  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  double dr2mean=0;
  double sumpt2=0;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) {
      if(pfType!=-1 && packCand->pdgId() != pdgid) continue;
      sumpt2 += packCand->pt() * packCand->pt();
      double dr = reco::deltaR(jet.eta(),jet.phi(),packCand->eta(),packCand->phi());    
      dr2mean += dr*dr*(packCand->pt())/(jet.pt());
    } else { 
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
      if(pfType!=-1 && pfcand->particleId() != pfType) continue;
      sumpt2 += pfcand->pt() * pfcand->pt();
      double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand->eta(),pfcand->phi());    
      dr2mean += dr*dr*(pfcand->pt() * pfcand->pt());
    }
  }
  dr2mean/=sumpt2;

  return dr2mean;
}

double JetTools::dR2Mean(const pat::Jet &jet, const int pfType)
{
  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );
  double dr2mean=0;
  double sumpt2=0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfType!=-1 && pfcand.pdgId() != pdgid) continue;

    sumpt2 += pfcand.pt() * pfcand.pt();

    double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand.eta(),pfcand.phi());
    dr2mean += dr*dr*(pfcand.pt() * pfcand.pt());
  }
  dr2mean/=sumpt2;

  return dr2mean;
}

//--------------------------------------------------------------------------------------------------
double JetTools::frac(const reco::PFJet &jet, const double dRMax, const int pfType)
{

  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  const double dRMin = dRMax - 0.1;
  
  double fraction = 0;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      if(pfType!=-1 && packCand->pdgId() != pdgid) continue;
      double dr = reco::deltaR(jet.eta(),jet.phi(),packCand->eta(),packCand->phi());
      if(dr > dRMax) continue;
      if(dr < dRMin) continue;
      fraction += packCand->pt() / jet.pt();
    } else { 
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
      if(pfType!=-1 && pfcand->particleId() != pfType) continue;
      
      double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand->eta(),pfcand->phi());    
      if(dr > dRMax) continue;
      if(dr < dRMin) continue;
      fraction += pfcand->pt() / jet.pt();
    }
  }
  
  return fraction;
}

double JetTools::frac(const pat::Jet &jet, const double dRMax, const int pfType)
{
  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  const double dRMin = dRMax - 0.1;

  double fraction = 0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfType!=-1 && pfcand.pdgId() != pdgid) continue;

    double dr = reco::deltaR(jet.eta(),jet.phi(),pfcand.eta(),pfcand.phi());
    if(dr > dRMax) continue;
    if(dr < dRMin) continue;

    fraction += pfcand.pt() / jet.pt();
  }
  
  return fraction;
}

//--------------------------------------------------------------------------------------------------
double JetTools::jetDz(const reco::PFJet &jet, const reco::Vertex &pv)
{
  // Assumes constituents are stored by descending pT
  double dz=-1000;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      if(packCand->charge()==0) continue;
      dz = packCand->dz(pv.position());
      break;
    } else { 
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
      const reco::TrackRef       track  = pfcand->trackRef();
      if(track.isNull()) continue;
      dz = track->dz(pv.position());
      break;
    }
  }  
  return dz;
}

double JetTools::jetDz(const pat::Jet &jet, const reco::Vertex &pv)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  // Assumes constituents are stored by descending pT
  double dz=-1000;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfcand.charge()==0) continue;
    dz = pfcand.dz(pv.position());
    break;
  }

  return dz;
}

//--------------------------------------------------------------------------------------------------
double JetTools::jetD0(const reco::PFJet &jet, const reco::Vertex &pv)
{
  // Assumes constituents are stored by descending pT
  double d0=-1000;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      if(packCand->charge()==0) continue;
      d0 = -packCand->dxy(pv.position());
      break;
    } else { 
      const reco::PFCandidatePtr pfcand = jet.getPFConstituent(ipf);
      const reco::TrackRef       track  = pfcand->trackRef();
      if(track.isNull()) continue;
      d0 = -track->dxy(pv.position());
      break;
    }
  }
  return d0;
}
double JetTools::jetD0(const pat::Jet &jet, const reco::Vertex &pv)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  // Assumes constituents are stored by descending pT
  double d0=-1000;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfcand.charge()==0) continue;
    d0 = -pfcand.dxy(pv.position());
    break;
  }

  return d0;
}

//--------------------------------------------------------------------------------------------------
double JetTools::jetWidth(const reco::PFJet &jet, const int varType, const int pfType)
{

  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  double ptD=0, sumPt=0, sumPt2=0;
  TMatrixDSym covMatrix(2); covMatrix=0.;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      if(pfType!=-1 && packCand->pdgId() != pdgid) continue;
      
      double dEta = jet.eta() - packCand->eta();
      double dPhi = reco::deltaPhi(jet.phi(), packCand->phi());
      
      covMatrix(0,0) += packCand->pt() * packCand->pt() * dEta * dEta;
      covMatrix(0,1) += packCand->pt() * packCand->pt() * dEta * dPhi;
      covMatrix(1,1) += packCand->pt() * packCand->pt() * dPhi * dPhi;
      ptD            += packCand->pt() * packCand->pt();
      sumPt          += packCand->pt();
      sumPt2         += packCand->pt() * packCand->pt();
    } else { 
      const reco::PFCandidatePtr pfcand = jet.getPFConstituent(ipf);
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

double JetTools::jetWidth(const pat::Jet &jet, const int varType, const int pfType)
{
  int pdgid = 0;  // ParticleType not stored in MINIAOD
  if     (pfType==reco::PFCandidate::X)         { pdgid = 0;   }  // undefined (dummy code)
  else if(pfType==reco::PFCandidate::h)         { pdgid = 211; }  // charged hadron
  else if(pfType==reco::PFCandidate::e)         { pdgid = 11;  }  // electron
  else if(pfType==reco::PFCandidate::mu)        { pdgid = 13;  }  // muon
  else if(pfType==reco::PFCandidate::gamma)     { pdgid = 22;  }  // photon
  else if(pfType==reco::PFCandidate::h0)        { pdgid = 130; }  // neutral hadron
  else if(pfType==reco::PFCandidate::h_HF)      { pdgid = 1;   }  // HF tower identified as a hadron (dummy code)
  else if(pfType==reco::PFCandidate::egamma_HF) { pdgid = 2;   }  // HF tower identified as an EM particle (dummy code)

  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  double ptD=0, sumPt=0, sumPt2=0;
  TMatrixDSym covMatrix(2); covMatrix=0.;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfType!=-1 && pfcand.pdgId() != pdgid) continue;

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
  if(jet.nConstituents ()                     < 2)    return false;
  if(jet.muonEnergy() / jet.energy()          > 0.8)  return false;
  if(fabs(jet.eta())<2.4) {
    if(jet.chargedHadronEnergy() / jet.energy() <= 0)   return false;
    if(jet.chargedEmEnergy() / jet.energy()     > 0.99) return false;
    if(jet.chargedMultiplicity()                < 1)    return false;
  }
  
  return true;
}

bool JetTools::passPFLooseID(const pat::Jet &jet)
{
  if(jet.energy() == 0) return false;
  if(jet.neutralHadronEnergyFraction() > 0.99) return false;
  if(jet.neutralEmEnergyFraction()     > 0.99) return false;
  if(jet.getPFConstituents().size()    < 2)    return false;
  if(jet.muonEnergyFraction()          > 0.8)  return false;
  if(fabs(jet.eta())<2.4) {
    if(jet.chargedHadronEnergyFraction() <= 0)   return false;
    if(jet.chargedEmEnergyFraction()     > 0.99) return false;
    if(jet.chargedMultiplicity()         < 1)    return false;
  }

  return true;
}

//--------------------------------------------------------------------------------------------------
double JetTools::jetCharge(const reco::PFJet &jet, const double kappa)
{
  double charge=0;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    if(packCand != 0) { 
      if(packCand->charge()==0) continue;
      charge += pow(packCand->pt(),kappa)*packCand->charge();
    } else { 
      const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
      const reco::TrackRef       track  = pfcand->trackRef();
      if(track.isNull()) continue;
      charge += pow(pfcand->pt(),kappa)*track->charge();
    }
  }
  double denom = (kappa==1) ? jet.pt() : pow(jet.pt(),kappa);
  return charge/denom;
}

double JetTools::jetCharge(const pat::Jet &jet, const double kappa)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  double charge=0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);
    if(pfcand.charge()==0) continue;
    charge += pow(pfcand.pt(),kappa)*pfcand.charge();
  }
  double denom = (kappa==1) ? jet.pt() : pow(jet.pt(),kappa);
  return charge/denom;
}

//--------------------------------------------------------------------------------------------------
double* JetTools::subJetBTag(const reco::PFJet &jet,reco::JetTagCollection &subJetVal,double iConeSize )
{
  // Following CMS convention first two daughters are leading subjets
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
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    double tpt_i=0, ty_i=0, tphi_i=0;
    double pt_i=0,   y_i=0,  phi_i=0;
    int tq_i = 0;
    if(packCand != 0) { 
      tq_i   = packCand->charge();
      tpt_i  = packCand->pt();
      ty_i   = packCand->rapidity();
      tphi_i = packCand->phi();
    } else { 
      const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
      tq_i   = pfcand->charge();
      tpt_i  = pfcand->pt();
      ty_i   = pfcand->rapidity();
      tphi_i = pfcand->phi();
    }
    if(type==0) {  // PF jet pull
      pt_i    = tpt_i;
      y_i     = ty_i;
      phi_i   = tphi_i;
    } else if(type==1) {  // charged jet pull component
      if(tq_i!=0) {
        pt_i  = tpt_i;
        y_i   = ty_i;
        phi_i = tphi_i;
      }
    } else if(type==2) {  // neutral jet pull component
      if(tq_i==0) {
        pt_i  = tpt_i;
        y_i   = ty_i;
        phi_i = tphi_i;
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

TVector2 JetTools::jetPull(const pat::Jet &jet, const int type)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  double dYSum=0, dPhiSum=0;
  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);

    double pt_i=0, y_i=0, phi_i=0;
    if(type==0) {  // PF jet pull
      pt_i  = pfcand.pt();
      y_i   = pfcand.rapidity();
      phi_i = pfcand.phi();

    } else if(type==1) {  // charged jet pull component
      if(pfcand.charge()!=0) {
        pt_i  = pfcand.pt();
        y_i   = pfcand.rapidity();
        phi_i = pfcand.phi();
      }

    } else if(type==2) {  // neutral jet pull component
      if(pfcand.charge()==0) {
        pt_i  = pfcand.pt();
        y_i   = pfcand.rapidity();
        phi_i = pfcand.phi();
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
double JetTools::jetPullAngle(const reco::PFJet &jet, edm::Handle<reco::PFJetCollection> &subJets,double iConeSize)
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

//--------------------------------------------------------------------------------------------------
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

//--------------------------------------------------------------------------------------------------
void JetTools::calcQGLVars(const reco::PFJet &jet, float &out_axis2, float &out_ptD, int &out_mult)
{
  // based on: RecoJets/JetProducers/plugins/QGTagger.cc
  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int mult = 0;

  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
    const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
    float deta,dphi,partPt,charge=0;
    if(packCand != 0) { 
      charge = packCand->charge();
      deta   = packCand->eta() - jet.eta();
      dphi   = reco::deltaPhi(packCand->phi(), jet.phi());
      partPt = packCand->pt();
    } else { 
      const reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(ipf);
      const reco::TrackRef       track  = pfcand->trackRef();
      // multiplicity = all charged particles + neutrals above 1 GeV
      if(track.isNull()) charge = 0;
      deta   = pfcand->eta() - jet.eta();
      dphi   = reco::deltaPhi(pfcand->phi(), jet.phi());
      partPt = pfcand->pt();
    }
    if((charge == 0 && partPt > 1) || charge != 0) mult++;
    float weight = partPt*partPt;
    sum_weight   += weight;
    sum_pt       += partPt;
    sum_deta     += deta*weight;
    sum_dphi     += dphi*weight;
    sum_deta2    += deta*deta*weight;
    sum_detadphi += deta*dphi*weight;
    sum_dphi2    += dphi*dphi*weight;
  }

  //Calculate axis2 and ptD
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ave_deta  = sum_deta/sum_weight;
    ave_dphi  = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a         = ave_deta2 - ave_deta*ave_deta;  			
    b         = ave_dphi2 - ave_dphi*ave_dphi;  			
    c         = -(sum_detadphi/sum_weight - ave_deta*ave_dphi); 	       
  }
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  float axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
  float ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);

  out_axis2 = axis2;
  out_ptD   = ptD;
  out_mult  = mult;
}

void JetTools::calcQGLVars(const pat::Jet &jet, float &out_axis2, float &out_ptD, int &out_mult)
{
  std::vector<reco::Candidate const *> constituents;
  for(unsigned int ida=0; ida<jet.numberOfDaughters(); ida++) {
    reco::Candidate const *cand = jet.daughter(ida);
    if(cand->numberOfDaughters()==0) {  // AK8 jets in MINIAOD are SUBJETS
      constituents.push_back(cand);
    } else {
      for(unsigned int jda=0; jda<cand->numberOfDaughters(); jda++) {
        reco::Candidate const *cand2 = cand->daughter(jda);
        constituents.push_back(cand2);
      }
    }
  }
  std::sort(constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda) { return ida->pt() > jda->pt(); } );

  // based on: RecoJets/JetProducers/plugins/QGTagger.cc
  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int mult = 0;

  const unsigned int nPFCands = constituents.size();
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {
    const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ipf]);

    // multiplicity = all charged particles + neutrals above 1 GeV
    if(pfcand.charge()!=0) {
      if(pfcand.pt()>1) { mult++; }
    } else {
      mult++;
    }

    float deta   = pfcand.eta() - jet.eta();
    float dphi   = reco::deltaPhi(pfcand.phi(), jet.phi());
    float partPt = pfcand.pt();
    float weight = partPt*partPt;

    sum_weight   += weight;
    sum_pt       += partPt;
    sum_deta     += deta*weight;
    sum_dphi     += dphi*weight;
    sum_deta2    += deta*deta*weight;
    sum_detadphi += deta*dphi*weight;
    sum_dphi2    += dphi*dphi*weight;
  }

  //Calculate axis2 and ptD
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ave_deta  = sum_deta/sum_weight;
    ave_dphi  = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a         = ave_deta2 - ave_deta*ave_deta;
    b         = ave_dphi2 - ave_dphi*ave_dphi;
    c         = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);
  }
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  float axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
  float ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);

  out_axis2 = axis2;
  out_ptD   = ptD;
  out_mult  = mult;
}

double JetTools::angle_squared(const fastjet::PseudoJet& jet1, const fastjet::PseudoJet& jet2)  {
  return jet1.squared_distance(jet2);
}
double JetTools::e2_func(double beta, fastjet::contrib::EnergyCorrelator::Measure measurelist, fastjet::PseudoJet myjet ){
  fastjet::contrib::EnergyCorrelator ECF1(1,beta,measurelist);
  fastjet::contrib::EnergyCorrelator ECF2(2,beta,measurelist);
  double temp = ECF2(myjet)/(pow(ECF1(myjet),2));
  return temp;
}

double JetTools::e3_func(double beta, fastjet::contrib::EnergyCorrelator::Measure measurelist, fastjet::PseudoJet myjet ){
  fastjet::contrib::EnergyCorrelator ECF1(1,beta,measurelist);
  fastjet::contrib::EnergyCorrelator ECF3(3,beta,measurelist);
  double temp =  ECF3(myjet)/(pow(ECF1(myjet),3));
  return temp;
}
double JetTools::e3_vn_func(unsigned int n, double beta,fastjet::contrib::EnergyCorrelator::Measure measurelist, fastjet::PseudoJet myjet ) {
  // n is the number of angles here used in the function, so 1 for 1e3, 2 for 2e3 and 3 for the usual e3.
  int N_total = 3;
  double temp = 0.0;
  double angle1, angle2, angle3;
  fastjet::contrib::EnergyCorrelator ECF1(1,beta,measurelist);
  double EJ = ECF1(myjet);
  vector<fastjet::PseudoJet> jet_consts = myjet.constituents();
  for (unsigned int i = 0; i< jet_consts.size(); i++){
    for (unsigned int j = i + 1; j < jet_consts.size(); j++){
      for (unsigned int k = j + 1; k < jet_consts.size(); k++) {
	angle1 = pow(angle_squared(jet_consts[i], jet_consts[j]), beta/2.);
	angle2 = pow(angle_squared(jet_consts[i], jet_consts[k]), beta/2.);
	angle3 = pow(angle_squared(jet_consts[j], jet_consts[k]), beta/2.);
	// Let's sort the angles first, then figure out what to multiply
	double angle_list[] = {angle1, angle2, angle3};
	std::vector<double> angle_vector (angle_list, angle_list+N_total);
	std::sort (angle_vector.begin(), angle_vector.begin()+N_total);
	double angle = 1.0;
	for (unsigned int l=0 ; l<n; l++){angle = angle*angle_vector[l];}
	temp += ( jet_consts[i].pt() / EJ) * ( jet_consts[j].pt() / EJ)*( jet_consts[k].pt() / EJ) * angle;
      }
    }
  }
  return temp;
}
