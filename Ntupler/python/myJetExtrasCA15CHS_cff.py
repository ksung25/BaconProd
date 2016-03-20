import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
CA15GenJetsCHS = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(1.5)
  )

# Jet Flavour
CA15FlavorCHS = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("CA15PFJetsCHS"),
    groomedJets              = cms.InputTag("CA15caPFJetsSoftDropCHS"),
    subjets                  = cms.InputTag("CA15caPFJetsSoftDropCHS", "SubJets"),
    bHadrons                 = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
    cHadrons                 = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
    partons                  = cms.InputTag("selectedHadronsAndPartons","algorithmicPartons"),
    jetAlgorithm             = cms.string("CambridgeAachen"),
    rParam                   = cms.double(1.5),
    ghostRescaling           = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(True)
  )

# take the default AK4 PFJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
CA15PFJetsCHS = ak4PFJets.clone(
    src          = cms.InputTag('pfNoPileUpJME'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(1.5),
    jetPtMin     = cms.double(150)
  )

# Pruned
CA15caPFJetsPrunedCHS = CA15PFJetsCHS.clone(
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
CA15caPFJetsTrimmedCHS = CA15PFJetsCHS.clone(
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
CA15caPFJetsSoftDropCHS = CA15PFJetsCHS.clone(
    useSoftDrop         = cms.bool(True),
    zcut                = cms.double(0.15),
    beta                = cms.double(1.0),
    R0                  = cms.double(1.5),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )
  
# q/g discriminator
# Note: need to provide JECs (or corrected jet collection) to QGL calculator
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff import *
from BaconProd.Ntupler.myCHSCorrections_cff          import *
ca15PFCHSL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFCHSL1FastjetCorrector','ak8PFCHSL2RelativeCorrector','ak8PFCHSL3AbsoluteCorrector')
  )
ca15PFCHSL1FastL2L3CorrectorChain = cms.Sequence(
  ak8PFCHSL1FastjetCorrector * ak8PFCHSL2RelativeCorrector * ak8PFCHSL3AbsoluteCorrector * ca15PFCHSL1FastL2L3Corrector
)

ca15PFCHSL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFCHSL1FastjetCorrector','ak8PFCHSL2RelativeCorrector','ak8PFCHSL3AbsoluteCorrector','ak8PFCHSResidualCorrector')
  )
ca15PFCHSL1FastL2L3ResidualCorrectorChain = cms.Sequence(
  ak8PFCHSL1FastjetCorrector * ak8PFCHSL2RelativeCorrector * ak8PFCHSL3AbsoluteCorrector * ak8PFCHSResidualCorrector * ca15PFCHSL1FastL2L3ResidualCorrector
)

from RecoJets.JetProducers.QGTagger_cfi import *
CA15QGTaggerCHS           = QGTagger.clone()
CA15QGTaggerCHS.srcJets   = cms.InputTag('CA15PFJetsCHS')
CA15QGTaggerCHS.jetsLabel = cms.string('QGL_AK4PFchs')
CA15QGTaggerCHS.jec       = cms.InputTag("ca15PFCHSL1FastL2L3Corrector") # NOTE: use "ca15PFCHSL1FastL2L3Corrector" for MC / "ca15PFCHSL1FastL2L3ResidualCorrector" for Data

CA15QGTaggerSubJetsCHS           = CA15QGTaggerCHS.clone()
CA15QGTaggerSubJetsCHS.srcJets   = cms.InputTag('CA15caPFJetsSoftDropCHS','SubJets')
CA15QGTaggerSubJetsCHS.jetsLabel = cms.string('QGL_AK4PFchs')
CA15QGTaggerSubJetsCHS.jec       = cms.InputTag("ca15PFCHSL1FastL2L3Corrector") # NOTE: use "ca15PFCHSL1FastL2L3Corrector" for MC / "ca15PFCHSL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
CA15NjettinessCHS = Njettiness.clone(
    src   = cms.InputTag('CA15PFJetsCHS'),
    cone  = cms.double(1.5),
    Njets = cms.vuint32(1,2,3,4)
  )

#
# Define sequences
#
CA15genjetsequenceCHS = cms.Sequence(
  CA15GenJetsCHS*
  CA15FlavorCHS
)

CA15jetsequenceCHS = cms.Sequence(
    CA15PFJetsCHS*
    CA15caPFJetsPrunedCHS*
    CA15caPFJetsTrimmedCHS*
    CA15caPFJetsSoftDropCHS*
    ca15PFCHSL1FastL2L3CorrectorChain*
    CA15QGTaggerCHS*
    CA15QGTaggerSubJetsCHS*                
    CA15NjettinessCHS*
    CA15GenJetsCHS*
    CA15FlavorCHS
  )

CA15jetsequenceCHSData = cms.Sequence(
    CA15PFJetsCHS*
    CA15caPFJetsPrunedCHS*
    CA15caPFJetsTrimmedCHS*
    CA15caPFJetsSoftDropCHS*
    #ca15PFCHSL1FastL2L3ResidualCorrectorChain*
    ca15chsL1FastL2L3ResidualChain*
    CA15QGTaggerCHS*
    CA15QGTaggerSubJetsCHS*
    CA15NjettinessCHS
  )

def setMiniAODCA15CHS(process) :
    process.CA15PFJetsCHS.src                                                   = cms.InputTag("packedPFCandidates")
    #process.CA15GenJetsCHS.src                                                  = cms.InputTag("packedGenParticles")
    process.CA15PFImpactParameterTagInfosCHS.primaryVertex                      = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA15PFImpactParameterTagInfosCHS.candidates                         = cms.InputTag("packedPFCandidates")
    process.CA15PFInclusiveSecondaryVertexFinderTagInfosCHS.extSVCollection     = cms.InputTag('slimmedSecondaryVertices')
    process.CA15PFImpactParameterTagInfosSJCHS.primaryVertex                    = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA15PFImpactParameterTagInfosSJCHS.candidates                       = cms.InputTag("packedPFCandidates")
    process.CA15PFInclusiveSecondaryVertexFinderTagInfosSJCHS.extSVCollection   = cms.InputTag('slimmedSecondaryVertices')
