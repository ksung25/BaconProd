#ifndef BACONPROD_NTUPLER_FILLEREVENTINFO_HH
#define BACONPROD_NTUPLER_FILLEREVENTINFO_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class TEventInfo;  // foward declaration

  class FillerEventInfo
  {
    public:
       FillerEventInfo(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC);
      ~FillerEventInfo();
      
      void fill(TEventInfo         *evtInfo,       // output object to be filled
                const edm::Event   &iEvent,        // EDM event info
		const reco::Vertex &pv,            // event primary vertex
		const bool          hasGoodPV,     // flag for if PV passing cuts is found
		const TriggerBits   triggerBits);  // bits for corresponding fired triggers
	       
    protected:
      void computeTrackMET(const unsigned int ipv,
			   const reco::VertexCollection *pvCol,
			   const reco::PFCandidateCollection *pfCandCol,
			   float &out_met, float &out_metphi);
    
      void computeTrackMET(const pat::PackedCandidateCollection *pfCandCol,
                           float &out_met, float &out_metphi);
    
      // EDM object collection names
      std::string fPFCandName;
      std::string fPUInfoName;
      std::string fPVName;
      std::string fBSName;
      std::string fCaloMETName;
      std::string fMETName;
      edm::InputTag fPFMETName;
      std::string fPFMETCName;
      std::string fMVAMETName;
      edm::InputTag fPUPPETName;
      std::string fPUPPETCName;
      std::string fALPACAMETName;
      std::string fPALPACAMETName;
      std::string fRhoIsoName;
      std::string fRhoJetName;
      bool fUseFilters;
      bool fUseAOD;

    edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  fTokPUInfoName;
    edm::EDGetTokenT<reco::BeamSpot>              fTokBSName         ;
    edm::EDGetTokenT<reco::CaloMETCollection>     fTokCaloMETName    ;
    edm::EDGetTokenT<pat::METCollection>          fTokCaloMETPATName ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokPFMETName      ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokPFMETCName     ;
    edm::EDGetTokenT<pat::METCollection>          fTokPFMETPATName     ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokMVAMETName     ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokPUPPETName     ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokPUPPETCName    ;
    edm::EDGetTokenT<pat::METCollection>          fTokPUPPETPATName     ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokALPACAMETName  ;
    edm::EDGetTokenT<reco::PFMETCollection>       fTokPALPACAMETName ;
    edm::EDGetTokenT<reco::VertexCollection>      fTokPVName         ;
    edm::EDGetTokenT<reco::PFCandidateCollection> fTokPFCandName     ;
    edm::EDGetTokenT<pat::PackedCandidateCollection> fTokPackCandName;
    edm::EDGetTokenT<double>                      fTokRhoIso;
    edm::EDGetTokenT<double>                      fTokRhoJet;
    edm::EDGetTokenT<reco::BeamHaloSummary>       fTokBeamHaloSummary;
    edm::EDGetTokenT<bool>                        fTokHBHENoiseFilterResultProducer;
    edm::EDGetTokenT<bool>                        fTokHcalLaserEventFilter         ;
    edm::EDGetTokenT<bool>                        fTokEEBadScFilter                ;
    edm::EDGetTokenT<bool>                        fTokEcalDeadCellTriggerPrimitiveFilter;
    edm::EDGetTokenT<bool>                        fToktrackingFailureFilter        ;
    edm::EDGetTokenT<bool>                        fTokManystripClus53X             ;
    edm::EDGetTokenT<bool>                        fTokTooManyStripClus53X          ;
    edm::EDGetTokenT<bool>                        fToklogErrorTooManyClusters      ;
    edm::EDGetTokenT<bool>                        fTokBadChCand                    ;
    edm::EDGetTokenT<bool>                        fTokBadPFMuon                    ;
    edm::EDGetTokenT<edm::TriggerResults>         fTokMetFiltersTag                ;
    edm::EDGetTokenT<double>                      fTokPrefWeight                    ;
    edm::EDGetTokenT<double>                      fTokPrefWeightUp                  ;
    edm::EDGetTokenT<double>                      fTokPrefWeightDown                ;
  };
}
#endif
