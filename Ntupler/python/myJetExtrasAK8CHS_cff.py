import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
AK8GenJetsCHS = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.8)
  )

# Jet Flavour
AK8FlavorCHS = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("AK8PFJetsCHS"),
    groomedJets              = cms.InputTag("AK8caPFJetsSoftDropCHS"),
    subjets                  = cms.InputTag("AK8caPFJetsSoftDropCHS", "SubJets"),
    bHadrons                 = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
    cHadrons                 = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
    partons                  = cms.InputTag("selectedHadronsAndPartons","partons"),
    jetAlgorithm             = cms.string("AntiKt"),
    rParam                   = cms.double(0.8),
    ghostRescaling           = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(True)
  )

# take the default AK4 PFJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
AK8PFJetsCHS = ak4PFJets.clone(
    src          = cms.InputTag('pfNoPileUpJME'),
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.8),
    jetPtMin     = cms.double(150)
  )

# Pruned
AK8caPFJetsPrunedCHS = AK8PFJetsCHS.clone(
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
AK8caPFJetsTrimmedCHS = AK8PFJetsCHS.clone(
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
AK8caPFJetsSoftDropCHS = AK8PFJetsCHS.clone(
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
AK8PFImpactParameterTagInfosCHS = pfImpactParameterTagInfos.clone(
    jets      = cms.InputTag('AK8PFJetsCHS'),
    maxDeltaR = cms.double(0.8)
  )
AK8PFInclusiveSecondaryVertexFinderTagInfosCHS = pfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("AK8PFImpactParameterTagInfosCHS")
  )
AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("AK8PFImpactParameterTagInfosCHS"), cms.InputTag("AK8PFInclusiveSecondaryVertexFinderTagInfosCHS") )
  )

# subjet b-tagging
AK8PFImpactParameterTagInfosSJCHS = pfImpactParameterTagInfos.clone(
    jets      = cms.InputTag('AK8caPFJetsSoftDropCHS','SubJets'),
    maxDeltaR = cms.double(0.8)
  )
AK8PFInclusiveSecondaryVertexFinderTagInfosSJCHS = pfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("AK8PFImpactParameterTagInfosSJCHS")
  )
AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("AK8PFImpactParameterTagInfosSJCHS"), cms.InputTag("AK8PFInclusiveSecondaryVertexFinderTagInfosSJCHS") )
  )

# q/g discriminator
# Note: need to provide JECs (or corrected jet collection) to QGL calculator
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff import *
ak8PFCHSL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFCHSL1FastjetCorrector','ak8PFCHSL2RelativeCorrector','ak8PFCHSL3AbsoluteCorrector')
  )
ak8PFCHSL1FastL2L3CorrectorChain = cms.Sequence(
  ak8PFCHSL1FastjetCorrector * ak8PFCHSL2RelativeCorrector * ak8PFCHSL3AbsoluteCorrector * ak8PFCHSL1FastL2L3Corrector
)

ak8PFCHSL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFCHSL1FastjetCorrector','ak8PFCHSL2RelativeCorrector','ak8PFCHSL3AbsoluteCorrector','ak8PFCHSResidualCorrector')
  )
ak8PFCHSL1FastL2L3ResidualCorrectorChain = cms.Sequence(
  ak8PFCHSL1FastjetCorrector * ak8PFCHSL2RelativeCorrector * ak8PFCHSL3AbsoluteCorrector * ak8PFCHSResidualCorrector * ak8PFCHSL1FastL2L3ResidualCorrector
)

from RecoJets.JetProducers.QGTagger_cfi import *
AK8QGTaggerCHS           = QGTagger.clone()
AK8QGTaggerCHS.srcJets   = cms.InputTag('AK8PFJetsCHS')
AK8QGTaggerCHS.jetsLabel = cms.string('QGL_AK4PFchs')
AK8QGTaggerCHS.jec       = cms.InputTag("ak8PFCHSL1FastL2L3Corrector") # NOTE: use "ak8PFCHSL1FastL2L3Corrector" for MC / "ak8PFCHSL1FastL2L3ResidualCorrector" for Data

AK8QGTaggerSubJetsCHS           = AK8QGTaggerCHS.clone()
AK8QGTaggerSubJetsCHS.srcJets   = cms.InputTag('AK8caPFJetsSoftDropCHS','SubJets')
AK8QGTaggerSubJetsCHS.jetsLabel = cms.string('QGL_AK4PFchs')
AK8QGTaggerSubJetsCHS.jec       = cms.InputTag("ak8PFCHSL1FastL2L3Corrector") # NOTE: use "ak8PFCHSL1FastL2L3Corrector" for MC / "ak8PFCHSL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
AK8NjettinessCHS = Njettiness.clone(
    src   = cms.InputTag('AK8PFJetsCHS'),
    cone  = cms.double(0.8),
    Njets = cms.vuint32(1,2,3,4)
  )

#
# Define sequences
#
AK8genjetsequenceCHS = cms.Sequence(
  AK8GenJetsCHS*
  AK8FlavorCHS
)

AK8jetsequenceCHS = cms.Sequence(
    AK8PFJetsCHS*
    AK8caPFJetsPrunedCHS*
    AK8caPFJetsTrimmedCHS*
    AK8caPFJetsSoftDropCHS*
    AK8PFImpactParameterTagInfosCHS*
    AK8PFInclusiveSecondaryVertexFinderTagInfosCHS*
    AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS*
    AK8PFImpactParameterTagInfosSJCHS*
    AK8PFInclusiveSecondaryVertexFinderTagInfosSJCHS*
    AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS*
    ak8PFCHSL1FastL2L3CorrectorChain*
    AK8QGTaggerCHS*
    AK8QGTaggerSubJetsCHS*                
    AK8NjettinessCHS*
    AK8GenJetsCHS*
    AK8FlavorCHS
  )

AK8jetsequenceCHSData = cms.Sequence(
    AK8PFJetsCHS*
    AK8caPFJetsPrunedCHS*
    AK8caPFJetsTrimmedCHS*
    AK8caPFJetsSoftDropCHS*
    AK8PFImpactParameterTagInfosCHS*
    AK8PFInclusiveSecondaryVertexFinderTagInfosCHS*
    AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS*
    AK8PFImpactParameterTagInfosSJCHS*
    AK8PFInclusiveSecondaryVertexFinderTagInfosSJCHS*
    AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS*
    ak8PFCHSL1FastL2L3ResidualCorrectorChain*
    AK8QGTaggerCHS*
    AK8QGTaggerSubJetsCHS*
    AK8NjettinessCHS
  )

def setMiniAODAK8CHS(process) :
    process.AK8PFJetsCHS.src                                                      = cms.InputTag("packedPFCandidates")
    #process.AK8GenJetsCHS.src                                                     = cms.InputTag("packedGenParticles")
    process.AK8PFImpactParameterTagInfosCHS.primaryVertex                         = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.AK8PFImpactParameterTagInfosCHS.candidates                            = cms.InputTag("packedPFCandidates")
    process.AK8PFInclusiveSecondaryVertexFinderTagInfosCHS.extSVCollection        = cms.InputTag('slimmedSecondaryVertices')
    process.AK8PFImpactParameterTagInfosSJCHS.primaryVertex                       = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.AK8PFImpactParameterTagInfosSJCHS.candidates                          = cms.InputTag("packedPFCandidates")
    process.AK8PFInclusiveSecondaryVertexFinderTagInfosSJCHS.extSVCollection      = cms.InputTag('slimmedSecondaryVertices')
    process.QGTaggerAK8                  = process.QGTagger.clone()
    process.QGTaggerAK8.srcJets          = cms.InputTag('slimmedJetsAK8')
    process.QGTaggerAK8.jetsLabel        = cms.string('QGL_AK4PFchs')
    process.QGTaggerAK8.jec              = cms.InputTag('')
    process.QGTaggerAK8.systematicsLabel = cms.string('')

