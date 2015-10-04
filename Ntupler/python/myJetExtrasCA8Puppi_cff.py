import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
CA8GenJetsPuppi = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8)
  )

# Jet Flavour
CA8FlavorPuppi = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("CA8PFJetsPuppi"),
    groomedJets              = cms.InputTag("CA8caPFJetsSoftDropPuppi"),
    subjets                  = cms.InputTag("CA8caPFJetsSoftDropPuppi", "SubJets"),
    bHadrons                 = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
    cHadrons                 = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
    partons                  = cms.InputTag("selectedHadronsAndPartons","partons"),
    jetAlgorithm             = cms.string("CambridgeAachen"),
    rParam                   = cms.double(0.8),
    ghostRescaling           = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(True)
  )

# take the default AK4 PFJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
CA8PFJetsPuppi = ak4PFJets.clone(
    src          = cms.InputTag('puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8),
    jetPtMin     = cms.double(150)
  )

# Pruned
CA8caPFJetsPrunedPuppi = CA8PFJetsPuppi.clone(
    cms.PSet(nFilt = cms.int32(2),
             zcut = cms.double(0.1),
             rcut_factor = cms.double(0.5)),
    jetAlgorithm        = cms.string("CambridgeAachen"),
    usePruning          = cms.bool(True),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Trimmed
CA8caPFJetsTrimmedPuppi = CA8PFJetsPuppi.clone(
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
CA8caPFJetsSoftDropPuppi = CA8PFJetsPuppi.clone(
    useSoftDrop         = cms.bool(True),
    zcut                = cms.double(0.1),
    beta                = cms.double(0.0),
    R0                  = cms.double(0.8),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# b-tagging
from RecoBTag.Configuration.RecoBTag_cff import *
CA8PFImpactParameterTagInfosPuppi = pfImpactParameterTagInfos.clone(
    jets      = cms.InputTag('CA8PFJetsPuppi'),
    maxDeltaR = cms.double(0.8)
  )
CA8PFInclusiveSecondaryVertexFinderTagInfosPuppi = pfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("CA8PFImpactParameterTagInfosPuppi")
  )
CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("CA8PFImpactParameterTagInfosPuppi"), cms.InputTag("CA8PFInclusiveSecondaryVertexFinderTagInfosPuppi") )
  )

# subjet b-tagging
CA8PFImpactParameterTagInfosSJPuppi = pfImpactParameterTagInfos.clone(
    jets      = cms.InputTag('CA8caPFJetsSoftDropPuppi','SubJets'),
    maxDeltaR = cms.double(0.8)
  )
CA8PFInclusiveSecondaryVertexFinderTagInfosSJPuppi = pfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("CA8PFImpactParameterTagInfosSJPuppi")
  )
CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("CA8PFImpactParameterTagInfosSJPuppi"), cms.InputTag("CA8PFInclusiveSecondaryVertexFinderTagInfosSJPuppi") )
  )

# q/g discriminator
# Note: need to provide JECs (or corrected jet collection) to QGL calculator
#from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff import *
from BaconProd.Ntupler.myPUPPICorrections_cff                  import *
from RecoJets.JetProducers.QGTagger_cfi import *
CA8QGTaggerPuppi           = QGTagger.clone()
CA8QGTaggerPuppi.srcJets   = cms.InputTag('CA8PFJetsPuppi')
CA8QGTaggerPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
CA8QGTaggerPuppi.jec       = cms.InputTag("ak8PuppiL1FastL2L3Corrector") # NOTE: use "ca8PFPuppiL1FastL2L3Corrector" for MC / "ca8PFPuppiL1FastL2L3ResidualCorrector" for Data

CA8QGTaggerSubJetsPuppi           = CA8QGTaggerPuppi.clone()
CA8QGTaggerSubJetsPuppi.srcJets   = cms.InputTag('CA8caPFJetsSoftDropPuppi','SubJets')
CA8QGTaggerSubJetsPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
CA8QGTaggerSubJetsPuppi.jec       = cms.InputTag("ak8PuppiL1FastL2L3Corrector") # NOTE: use "ca8PFPuppiL1FastL2L3Corrector" for MC / "ca8PFPuppiL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
CA8NjettinessPuppi = Njettiness.clone(
    src   = cms.InputTag('CA8PFJetsPuppi'),
    cone  = cms.double(0.8),
    Njets = cms.vuint32(1,2,3,4)
  )

CA8genjetsequencePuppi = cms.Sequence(
  CA8GenJetsPuppi* 
  CA8FlavorPuppi
)

#
# Define sequences
#
CA8jetsequencePuppi = cms.Sequence(
    CA8PFJetsPuppi*
    CA8caPFJetsPrunedPuppi*
    CA8caPFJetsTrimmedPuppi*
    CA8caPFJetsSoftDropPuppi*
    CA8PFImpactParameterTagInfosPuppi*
    CA8PFInclusiveSecondaryVertexFinderTagInfosPuppi*
    CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi*
    CA8PFImpactParameterTagInfosSJPuppi*
    CA8PFInclusiveSecondaryVertexFinderTagInfosSJPuppi*
    CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi*
    ak8PuppiL1FastL2L3Chain*
    CA8QGTaggerPuppi*
    CA8QGTaggerSubJetsPuppi*                
    CA8NjettinessPuppi*
    CA8FlavorPuppi
    )

CA8jetsequencePuppiData = cms.Sequence(
    CA8PFJetsPuppi*
    CA8caPFJetsPrunedPuppi*
    CA8caPFJetsTrimmedPuppi*
    CA8caPFJetsSoftDropPuppi*
    CA8PFImpactParameterTagInfosPuppi*
    CA8PFInclusiveSecondaryVertexFinderTagInfosPuppi*
    CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi*
    CA8PFImpactParameterTagInfosSJPuppi*
    CA8PFInclusiveSecondaryVertexFinderTagInfosSJPuppi*
    CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi*
    #ak8PuppiL1FastL2L3ResidualChain*
    ak8PuppiL1FastL2L3Chain*
    CA8QGTaggerPuppi*
    CA8QGTaggerSubJetsPuppi*                
    CA8NjettinessPuppi
  )

def setMiniAODCA8Puppi(process) :
    process.CA8PFImpactParameterTagInfosPuppi.primaryVertex                    = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA8PFImpactParameterTagInfosPuppi.candidates                       = cms.InputTag("packedPFCandidates")
    process.CA8PFInclusiveSecondaryVertexFinderTagInfosPuppi.extSVCollection   = cms.InputTag('slimmedSecondaryVertices')
    process.CA8PFImpactParameterTagInfosSJPuppi.primaryVertex                  = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA8PFImpactParameterTagInfosSJPuppi.candidates                     = cms.InputTag("packedPFCandidates")
    process.CA8PFInclusiveSecondaryVertexFinderTagInfosSJPuppi.extSVCollection = cms.InputTag('slimmedSecondaryVertices')
