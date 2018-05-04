import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
CA15GenJetsPuppi = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(1.5)
  )

# Jet Flavour
CA15FlavorPuppi = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("CA15PFJetsPuppi"),
    groomedJets              = cms.InputTag("CA15caPFJetsSoftDropPuppi"),
    subjets                  = cms.InputTag("CA15caPFJetsSoftDropPuppi", "SubJets"),
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
CA15PFJetsPuppi = ak4PFJets.clone(
    src          = cms.InputTag('puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(1.5),
    jetPtMin     = cms.double(150)
  )

# Pruned
CA15caPFJetsPrunedPuppi = CA15PFJetsPuppi.clone(
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
CA15caPFJetsTrimmedPuppi = CA15PFJetsPuppi.clone(
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
CA15caPFJetsSoftDropPuppi = CA15PFJetsPuppi.clone(
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
from BaconProd.Ntupler.myPUPPICorrections_cff                  import *

from RecoJets.JetProducers.QGTagger_cfi import *
CA15QGTaggerPuppi           = QGTagger.clone()
CA15QGTaggerPuppi.srcJets   = cms.InputTag('CA15PFJetsPuppi')
CA15QGTaggerPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
CA15QGTaggerPuppi.jec       = cms.InputTag("ca15PuppiL1FastL2L3Corrector") # NOTE: use "ca15PFPuppiL1FastL2L3Corrector" for MC / "ca15PFPuppiL1FastL2L3ResidualCorrector" for Data

CA15QGTaggerSubJetsPuppi           = CA15QGTaggerPuppi.clone()
CA15QGTaggerSubJetsPuppi.srcJets   = cms.InputTag('CA15caPFJetsSoftDropPuppi','SubJets')
CA15QGTaggerSubJetsPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
CA15QGTaggerSubJetsPuppi.jec       = cms.InputTag("ca15PuppiL1FastL2L3Corrector") # NOTE: use "ca15PFPuppiL1FastL2L3Corrector" for MC / "ca15PFPuppiL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
CA15NjettinessPuppi = Njettiness.clone(
    src   = cms.InputTag('CA15PFJetsPuppi'),
    cone  = cms.double(1.5),
    Njets = cms.vuint32(1,2,3,4)
  )

CA15genjetsequencePuppi = cms.Sequence(
  CA15GenJetsPuppi* 
  CA15FlavorPuppi
)

#
# Define sequences
#
CA15jetsequencePuppi = cms.Sequence(
    CA15PFJetsPuppi*
    CA15caPFJetsPrunedPuppi*
    CA15caPFJetsTrimmedPuppi*
    CA15caPFJetsSoftDropPuppi*
    ca15PuppiL1FastL2L3Chain*
    CA15QGTaggerPuppi*
    CA15QGTaggerSubJetsPuppi*                
    CA15NjettinessPuppi*
    CA15FlavorPuppi
  )

CA15jetsequencePuppiData = cms.Sequence(
    CA15PFJetsPuppi*
    CA15caPFJetsPrunedPuppi*
    CA15caPFJetsTrimmedPuppi*
    CA15caPFJetsSoftDropPuppi*
    ca15PuppiL1FastL2L3ResidualChain*
    #ca15PuppiL1FastL2L3Chain*
    CA15QGTaggerPuppi*
    CA15QGTaggerSubJetsPuppi*                
    CA15NjettinessPuppi
  )

def setMiniAODCA15Puppi(process) :
    process.CA15PFImpactParameterTagInfosPuppi.primaryVertex                    = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA15PFImpactParameterTagInfosPuppi.candidates                       = cms.InputTag("packedPFCandidates")
    process.CA15PFInclusiveSecondaryVertexFinderTagInfosPuppi.extSVCollection   = cms.InputTag('slimmedSecondaryVertices')
    process.CA15PFImpactParameterTagInfosSJPuppi.primaryVertex                  = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.CA15PFImpactParameterTagInfosSJPuppi.candidates                     = cms.InputTag("packedPFCandidates")
    process.CA15PFInclusiveSecondaryVertexFinderTagInfosSJPuppi.extSVCollection = cms.InputTag('slimmedSecondaryVertices')
