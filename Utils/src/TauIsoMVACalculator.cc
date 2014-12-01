#include "BaconProd/Utils/interface/TauIsoMVACalculator.hh"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "TMVA/Reader.h"
#include "Cintex/Cintex.h"
#include <TFile.h>

using namespace baconhep;



//--------------------------------------------------------------------------------------------------
TauIsoMVACalculator::TauIsoMVACalculator():
  fIsInitialized(false),
  fUseGBR       (false),
  fReader       (0),
  fGBRReader    (0)
{}

//--------------------------------------------------------------------------------------------------
TauIsoMVACalculator::~TauIsoMVACalculator()
{
  delete fReader;
  delete fGBRReader;
}

//--------------------------------------------------------------------------------------------------
void TauIsoMVACalculator::initialize(const std::string weightFile, const bool useGBR)
{
  fUseGBR = useGBR;
  
  if(useGBR) {    
    ROOT::Cintex::Cintex::Enable();  // Needed to read non-TObject classes from a ROOT file!
    TFile *forest = new TFile(weightFile.c_str(),"READ");
    fGBRReader = (GBRForest*)forest->Get("gbrfTauIso");
    forest->Close();
  
  } else {
    fReader = new TMVA::Reader("!Color:!Silent");
    fReader->AddVariable("tau_nchiso",        (float *)0);
    fReader->AddVariable("tau_ngiso",         (float *)0);
    fReader->AddVariable("tau_nneuiso",       (float *)0);
    fReader->AddVariable("ring_ch_0*tau_pt",  (float *)0);
    fReader->AddVariable("ring_ch_1*tau_pt",  (float *)0);
    fReader->AddVariable("ring_ch_2*tau_pt",  (float *)0);
    fReader->AddVariable("ring_ch_3*tau_pt",  (float *)0);
    fReader->AddVariable("ring_ch_4*tau_pt",  (float *)0);
    fReader->AddVariable("ring_g_0*tau_pt",   (float *)0);
    fReader->AddVariable("ring_g_1*tau_pt",   (float *)0);
    fReader->AddVariable("ring_g_2*tau_pt",   (float *)0);
    fReader->AddVariable("ring_g_3*tau_pt",   (float *)0);
    fReader->AddVariable("ring_g_4*tau_pt",   (float *)0);
    fReader->AddVariable("ring_neu_0*tau_pt", (float *)0);
    fReader->AddVariable("ring_neu_1*tau_pt", (float *)0);
    fReader->AddVariable("ring_neu_2*tau_pt", (float *)0);
    fReader->AddVariable("ring_neu_3*tau_pt", (float *)0);
    fReader->AddVariable("ring_neu_4*tau_pt", (float *)0);
    fReader->AddVariable("shape_ch_eta",      (float *)0);
    fReader->AddVariable("shape_ch_phi",      (float *)0);
    fReader->AddVariable("shape_ch_etaeta",   (float *)0);
    fReader->AddVariable("shape_ch_phiphi",   (float *)0);
    fReader->AddVariable("shape_ch_etaphi",   (float *)0);
    fReader->AddVariable("shape_g_eta",       (float *)0);
    fReader->AddVariable("shape_g_phi",       (float *)0);
    fReader->AddVariable("shape_g_etaeta",    (float *)0);
    fReader->AddVariable("shape_g_phiphi",    (float *)0);
    fReader->AddVariable("shape_g_etaphi",    (float *)0);
    fReader->AddVariable("shape_neu_eta",     (float *)0);
    fReader->AddVariable("shape_neu_phi",     (float *)0);
    fReader->AddVariable("shape_neu_etaeta",  (float *)0);
    fReader->AddVariable("shape_neu_phiphi",  (float *)0);
    fReader->AddVariable("shape_neu_etaphi",  (float *)0);
    fReader->AddVariable("rho",               (float *)0);
    
    fReader->AddSpectator("tau_pt",           (float *)0);
    fReader->AddSpectator("tau_eta",          (float *)0);
    fReader->AddSpectator("tau_iso",          (float *)0);
    fReader->AddSpectator("gen_pt",           (float *)0);
    fReader->AddSpectator("jet_pt",           (float *)0);
    fReader->AddSpectator("pv",               (float *)0);
    
    fReader->BookMVA("BDTG", weightFile);
  }
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
float TauIsoMVACalculator::mvaValue(const reco::PFTau &tau, const double rho)
{
  baconhep::TauIsoRings rings = computeIsoRings(tau);
  
  std::vector<float> mvaInput = rings.getVector();
  mvaInput.push_back(rho);
  if(!fUseGBR) {
    mvaInput.insert(mvaInput.end(), 6, 0);
    return fReader->EvaluateMVA(mvaInput,"BDTG");
  }
  
  return fGBRReader->GetClassifier(&mvaInput[0]);
}

//--------------------------------------------------------------------------------------------------
baconhep::TauIsoRings TauIsoMVACalculator::computeIsoRings(const reco::PFTau &tau) {
  std::vector<int>                  niso    (3);
  std::vector< std::vector<float> > rings   (3, std::vector<float>(5));
  std::vector< std::vector<float> > shapes  (3, std::vector<float>(5));
  std::vector<float>                isoPtSum(3);
  
  for(reco::PFCandidateRefVector::const_iterator itPF = tau.isolationPFCands().begin(); 
      itPF!=tau.isolationPFCands().end(); 
      ++itPF)
  {
    reco::PFCandidateRef pfcand = *itPF;
    
    float dEta = tau.eta() - pfcand->eta();
    float dPhi = reco::deltaPhi(tau.phi(), pfcand->phi());
    float dR   = reco::deltaR(tau.eta(),tau.phi(), pfcand->eta(),pfcand->phi());
    
    int pftype = 0;
    
    if     (pfcand->charge() != 0)                            { pftype = 0; }
    else if(pfcand->particleId() == reco::PFCandidate::gamma) { pftype = 1; }
    else                                                      { pftype = 2; }
    
    // number of isolation candidates by type
    niso[pftype]++;
    
    // isolation rings
    if     (dR < 0.1) { rings[pftype][0] += pfcand->pt(); }
    else if(dR < 0.2) { rings[pftype][1] += pfcand->pt(); } 
    else if(dR < 0.3) { rings[pftype][2] += pfcand->pt(); }
    else if(dR < 0.4) { rings[pftype][3] += pfcand->pt(); }
    else if(dR < 0.5) { rings[pftype][4] += pfcand->pt(); }
    
    // angle shape variables
    shapes[pftype][0] += pfcand->pt() * dEta;
    shapes[pftype][1] += pfcand->pt() * dPhi;
    shapes[pftype][2] += pfcand->pt() * dEta*dEta;
    shapes[pftype][3] += pfcand->pt() * dPhi*dPhi;
    shapes[pftype][4] += pfcand->pt() * dEta*dPhi;
    
    // overall pT sum
    isoPtSum[pftype] += pfcand->pt();
  }
  
  // mean and variance of angle variables weighted by pT
  for(unsigned int i=0; i<shapes.size(); i++) {
    for(unsigned int j=0; j<shapes[i].size(); j++) {
      shapes[i][j] = (isoPtSum[i]>0) ? fabs(shapes[i][j]/isoPtSum[i]) : 0;
    }
  }
  
  // Fill TauIsoRings object
  baconhep::TauIsoRings isoRings;
  isoRings.niso   = niso;
  isoRings.rings  = rings;
  isoRings.shapes = shapes;
  
  return isoRings;
}
