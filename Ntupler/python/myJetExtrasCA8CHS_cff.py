
import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
CA8GenJetsCHS = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8)
  )

# Jet Flavour
CA8FlavorCHS = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("CA8PFJetsCHS"),
    groomedJets              = cms.InputTag("CA8caPFJetsSoftDropCHS"),
    subjets                  = cms.InputTag("CA8caPFJetsSoftDropCHS", "SubJets"),
    bHadrons                 = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
    cHadrons                 = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
    partons                  = cms.InputTag("selectedHadronsAndPartons","algorithmicPartons"),
    jetAlgorithm             = cms.string("CambridgeAachen"),
    rParam                   = cms.double(0.8),
    ghostRescaling           = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(True)
  )

# take the default AK4 PFJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
CA8PFJetsCHS = ak4PFJets.clone(
    src          = cms.InputTag('pfNoPileUpJME'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8),
    jetPtMin     = cms.double(150)
  )

# Pruned
CA8caPFJetsPrunedCHS = CA8PFJetsCHS.clone(
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
CA8caPFJetsTrimmedCHS = CA8PFJetsCHS.clone(
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
CA8caPFJetsSoftDropCHS = CA8PFJetsCHS.clone(
    useSoftDrop         = cms.bool(True),
    zcut                = cms.double(0.1),
    beta                = cms.double(0.0),
    R0                  = cms.double(0.8),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# q/g discriminator
# Note: need to provide JECs (or corrected jet collection) to QGL calculator
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff import *
#from BaconProd.Ntupler.myCHSCorrections_cff          import *
ca8PFCHSL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFCHSL1FastjetCorrector','ak8PFCHSL2RelativeCorrector','ak8PFCHSL3AbsoluteCorrector')
  )
ca8PFCHSL1FastL2L3CorrectorChain = cms.Sequence(
  ak8PFCHSL1FastjetCorrector * ak8PFCHSL2RelativeCorrector * ak8PFCHSL3AbsoluteCorrector * ca8PFCHSL1FastL2L3Corrector
)

ca8PFCHSL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFCHSL1FastjetCorrector','ak8PFCHSL2RelativeCorrector','ak8PFCHSL3AbsoluteCorrector','ak8PFCHSResidualCorrector')
  )
ca8PFCHSL1FastL2L3ResidualCorrectorChain = cms.Sequence(
  ak8PFCHSL1FastjetCorrector * ak8PFCHSL2RelativeCorrector * ak8PFCHSL3AbsoluteCorrector * ak8PFCHSResidualCorrector * ca8PFCHSL1FastL2L3ResidualCorrector
)

from RecoJets.JetProducers.QGTagger_cfi import *
CA8QGTaggerCHS           = QGTagger.clone()
CA8QGTaggerCHS.srcJets   = cms.InputTag('CA8PFJetsCHS')
CA8QGTaggerCHS.jetsLabel = cms.string('QGL_AK4PFchs')
CA8QGTaggerCHS.jec       = cms.InputTag("ca8PFCHSL1FastL2L3Corrector") # NOTE: use "ca8PFCHSL1FastL2L3Corrector" for MC / "ca8PFCHSL1FastL2L3ResidualCorrector" for Data

CA8QGTaggerSubJetsCHS           = CA8QGTaggerCHS.clone()
CA8QGTaggerSubJetsCHS.srcJets   = cms.InputTag('CA8caPFJetsSoftDropCHS','SubJets')
CA8QGTaggerSubJetsCHS.jetsLabel = cms.string('QGL_AK4PFchs')
CA8QGTaggerSubJetsCHS.jec       = cms.InputTag("ca8PFCHSL1FastL2L3Corrector") # NOTE: use "ca8PFCHSL1FastL2L3Corrector" for MC / "ca8PFCHSL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
CA8NjettinessCHS = Njettiness.clone(
    src   = cms.InputTag('CA8PFJetsCHS'),
    cone  = cms.double(0.8),
    Njets = cms.vuint32(1,2,3,4)
  )

#
# Define sequences
#
CA8genjetsequenceCHS = cms.Sequence(
  CA8GenJetsCHS
  #CA8FlavorCHS
)

CA8jetsequenceCHS = cms.Sequence(
    CA8PFJetsCHS*
    CA8caPFJetsPrunedCHS*
    CA8caPFJetsTrimmedCHS*
    CA8caPFJetsSoftDropCHS*
    ca8PFCHSL1FastL2L3CorrectorChain*
    CA8QGTaggerCHS*
    CA8QGTaggerSubJetsCHS*                
    CA8NjettinessCHS*
    CA8GenJetsCHS*
    CA8FlavorCHS
  )

CA8jetsequenceCHSData = cms.Sequence(
    CA8PFJetsCHS*
    CA8caPFJetsPrunedCHS*
    CA8caPFJetsTrimmedCHS*
    CA8caPFJetsSoftDropCHS*
    ca8PFCHSL1FastL2L3ResidualCorrectorChain*
    CA8QGTaggerCHS*
    CA8QGTaggerSubJetsCHS*
    CA8NjettinessCHS
  )

def setMiniAODCA8CHS(process) :
    process.CA8PFJetsCHS.src                                                        = cms.InputTag("packedPFCandidates")
    #process.CA8GenJetsCHS.src                                                       = cms.InputTag("packedGenParticles")
    process.CA8PFImpactParameterTagInfosCHS.primaryVertex                           = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA8PFImpactParameterTagInfosCHS.candidates                              = cms.InputTag("packedPFCandidates")
    process.CA8PFInclusiveSecondaryVertexFinderTagInfosCHS.extSVCollection          = cms.InputTag('slimmedSecondaryVertices')
    process.CA8PFImpactParameterTagInfosSJCHS.primaryVertex                         = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA8PFImpactParameterTagInfosSJCHS.candidates                            = cms.InputTag("packedPFCandidates")
    process.CA8PFInclusiveSecondaryVertexFinderTagInfosSJCHS.extSVCollection        = cms.InputTag('slimmedSecondaryVertices')
