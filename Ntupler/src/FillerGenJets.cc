#include "BaconProd/Ntupler/interface/FillerGenJets.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include <TClonesArray.h>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerGenJets::FillerGenJets(const edm::ParameterSet &iConfig):
  fGenParName    (iConfig.getUntrackedParameter<std::string>("edmGenParticlesName","genParticles")),
  fGenJetName    (iConfig.getUntrackedParameter<std::string>("genJetName","AK4GenJets"))
{}

//--------------------------------------------------------------------------------------------------
FillerGenJets::~FillerGenJets(){}

//--------------------------------------------------------------------------------------------------
void FillerGenJets::fill(TClonesArray *array,
                         const edm::Event &iEvent)
{
  assert(array);
  // Get generator jet collection
  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJets = 0;
  iEvent.getByLabel(fGenJetName,hGenJetProduct);
  assert(hGenJetProduct.isValid());
  genJets = hGenJetProduct.product();

  // Get Jet Flavor Match
  //edm::Handle<reco::JetFlavourMatchingCollection> jetFlavourMatch;
  //iEvent.getByLabel(fJetFlavorName, jetFlavourMatch);
  //edm::Handle<reco::JetFlavourMatchingCollection> jetFlavourMatchPhys;
  //iEvent.getByLabel(fJetFlavorPhysName, jetFlavourMatchPhys);

  // Get generator particles collection
  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByLabel(fGenParName,hGenParProduct);
  assert(hGenParProduct.isValid());  
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  
  // loop over GEN particles
  std::vector<edm::Ptr<reco::GenParticle>> lMothers;
  TClonesArray &rArray = *array;
  for (reco::GenJetCollection::const_iterator itGenJ = genJets->begin(); itGenJ!=genJets->end(); ++itGenJ) {
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TGenJet();
    baconhep::TGenJet *pGenJet = (baconhep::TGenJet*)rArray[index];
    pGenJet->pdgId   = flavor(&(*itGenJ),genParticles);
    pGenJet->pt      = itGenJ->pt();
    pGenJet->eta     = itGenJ->eta();
    pGenJet->phi     = itGenJ->phi();
    pGenJet->mass   = itGenJ->mass();
    double* lElectrons =  genCone(&(*itGenJ),genParticles,0,0.10,  -1);
    double* lMuons     =  genCone(&(*itGenJ),genParticles,0,0.10,  -2);
    double* lPhotons   =  genCone(&(*itGenJ),genParticles,0,0.10,   1);
    double* lTotal     =  genCone(&(*itGenJ),genParticles,0,0.10,   0);
    double* lIso03     =  genCone(&(*itGenJ),genParticles,0.1,0.3,  0);
    double* lIso04     =  genCone(&(*itGenJ),genParticles,0.1,0.4,  0);
    double* lIso05     =  genCone(&(*itGenJ),genParticles,0.1,0.5,  0);

    pGenJet->elept       =  lElectrons[1];
    pGenJet->eleeta      =  lElectrons[2];
    pGenJet->elephi      =  lElectrons[3];
    pGenJet->elem        =  lElectrons[4];

    pGenJet->mupt        =  lMuons[1];
    pGenJet->mueta       =  lMuons[2];
    pGenJet->muphi       =  lMuons[3];
    pGenJet->mum         =  lMuons[4];

    pGenJet->gapt        =  lPhotons[1];
    pGenJet->gaeta       =  lPhotons[2];
    pGenJet->gaphi       =  lPhotons[3];
    pGenJet->gam         =  lPhotons[4];

    pGenJet->totpt       =  lTotal[1];
    pGenJet->toteta      =  lTotal[2];
    pGenJet->totphi      =  lTotal[3];
    pGenJet->totm        =  lTotal[4];

    pGenJet->iso03       =  lIso03[0];
    pGenJet->iso04       =  lIso04[0];
    pGenJet->iso05       =  lIso05[0];

    delete lElectrons;
    delete lMuons    ; 
    delete lPhotons  ; 
    delete lTotal    ; 
    delete lIso03    ; 
    delete lIso04    ; 
    delete lIso05;
  }
}
double* FillerGenJets::genCone(const reco::GenJet *iJet,const reco::GenParticleCollection &iGenParticles,double iDRMin,double iDRMax,int iType) {
  double lIso = 0;
  TLorentzVector lVec; lVec.SetPtEtaPhiM(0,0,0,0);
  for (reco::GenParticleCollection::const_iterator itGenP = iGenParticles.begin(); itGenP!=iGenParticles.end(); ++itGenP) {
    if(itGenP->status() != 1) continue;
    double pDPhi1 = fabs(iJet->phi()-itGenP->phi()); if(pDPhi1 > 2.*TMath::Pi()-pDPhi1) pDPhi1 = 2.*TMath::Pi()-pDPhi1;
    double pDEta1 = fabs(iJet->eta()-itGenP->eta());
    double pDR    = sqrt(pDPhi1*pDPhi1 + pDEta1*pDEta1);
    if(pDR  < iDRMin  || pDR > iDRMax) continue;
    int pPdgId = itGenP->pdgId();
    if(fabs(pPdgId) < 17  && fabs(pPdgId) > 10 && pPdgId % 2 == 0 )  continue; //skip neutrinos
    if(itGenP->charge()   == 0  && iType <   0) continue;
    if(itGenP->charge()   != 0  && iType >   0) continue;
    if(fabs(pPdgId)       != 22 && iType ==  1) continue;
    if(fabs(pPdgId)       != 11 && iType == -1) continue;
    if(fabs(pPdgId)       != 13 && iType == -2) continue;
    TLorentzVector pVec; pVec.SetPtEtaPhiM(itGenP->pt(),itGenP->eta(),itGenP->phi(),itGenP->mass());
    lVec += pVec;
    lIso += itGenP->pt();
  }
  double *lReturn = new double[5];
  lReturn[0] = lIso;
  if(lVec.Px() > 0) {
    lReturn[1] = lVec.Pt();
    lReturn[2] = lVec.Eta();
    lReturn[3] = lVec.Phi();
    lReturn[4] = lVec.M();
  }
  return lReturn;
}
int FillerGenJets::flavor(const reco::GenJet *iJet,const reco::GenParticleCollection &iGenParticles) { 
  int    lId    = -1;
  double lPtMax = -1; 
  for (reco::GenParticleCollection::const_iterator itGenP = iGenParticles.begin(); itGenP!=iGenParticles.end(); ++itGenP) {
    double pDPhi1 = fabs(iJet->phi()-itGenP->phi()); if(pDPhi1 > 2.*TMath::Pi()-pDPhi1) pDPhi1 = 2.*TMath::Pi()-pDPhi1;
    double pDEta1 = fabs(iJet->eta()-itGenP->eta());
    double pDR    = sqrt(pDPhi1*pDPhi1 + pDEta1*pDEta1);
    if(pDR  > 0.25) continue;
    if(itGenP->pt() > lPtMax) lId    = itGenP->pdgId();
    if(itGenP->pt() > lPtMax) lPtMax = itGenP->pt(); 
  }  
  return lId;
}
