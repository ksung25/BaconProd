import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
AK4GenJetsPuppi = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4)
  )

# Jet Flavour
AK4FlavorPuppi = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("AK4PFJetsPuppi"),
    bHadrons                 = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
    cHadrons                 = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
    partons                  = cms.InputTag("selectedHadronsAndPartons","algorithmicPartons"),
    jetAlgorithm             = cms.string("AntiKt"),
    rParam                   = cms.double(0.4),
    ghostRescaling           = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(True)
  )

# take the default AK4 PFJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
AK4PFJetsPuppi = ak4PFJets.clone(
    src          = cms.InputTag('puppi'),
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    jetPtMin     = cms.double(1)
  )

# Pruned
AK4caPFJetsPrunedPuppi = AK4PFJetsPuppi.clone(
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
AK4caPFJetsTrimmedPuppi = AK4PFJetsPuppi.clone(
    jetAlgorithm        = cms.string("CambridgeAachen"),
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
AK4caPFJetsSoftDropPuppi = AK4PFJetsPuppi.clone(
    jetAlgorithm        = cms.string("CambridgeAachen"),
    useSoftDrop         = cms.bool(True),
    zcut                = cms.double(0.1),
    beta                = cms.double(0.0),
    R0                  = cms.double(0.4),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# q/g discriminator
# Note: need to provide JECs (or corrected jet collection) to QGL calculator
from BaconProd.Ntupler.myPUPPICorrections_cff                  import *
from RecoJets.JetProducers.QGTagger_cfi import *
AK4QGTaggerPuppi           = QGTagger.clone()
AK4QGTaggerPuppi.srcJets   = cms.InputTag('AK4PFJetsPuppi')
AK4QGTaggerPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
AK4QGTaggerPuppi.jec       = cms.InputTag("ak4PuppiL1FastL2L3Corrector") # NOTE: use "ak4PFPuppiL1FastL2L3Corrector" for MC / "ak4PFPuppiL1FastL2L3ResidualCorrector" for Data

AK4QGTaggerSubJetsPuppi           = AK4QGTaggerPuppi.clone()
AK4QGTaggerSubJetsPuppi.srcJets   = cms.InputTag('AK4caPFJetsSoftDropPuppi','SubJets')
AK4QGTaggerSubJetsPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
AK4QGTaggerSubJetsPuppi.jec       = cms.InputTag("ak4PuppiL1FastL2L3Corrector") # NOTE: use "ak4PFPuppiL1FastL2L3Corrector" for MC / "ak4PFPuppiL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
AK4NjettinessPuppi = Njettiness.clone(
    src   = cms.InputTag('AK4PFJetsPuppi'),
    cone  = cms.double(0.4),
    Njets = cms.vuint32(1,2,3,4)
  )

#
# Define sequences
#
AK4genjetsequencePuppi = cms.Sequence(
  AK4GenJetsPuppi* 
  AK4FlavorPuppi
)

AK4jetsequencePuppi = cms.Sequence(
    AK4PFJetsPuppi*   
    AK4caPFJetsPrunedPuppi*
    AK4caPFJetsTrimmedPuppi*
    AK4caPFJetsSoftDropPuppi*
    ak4PuppiL1FastL2L3Chain* #   => using type 1 Met
    AK4QGTaggerPuppi*
    AK4QGTaggerSubJetsPuppi*                
    AK4NjettinessPuppi*
    AK4FlavorPuppi
  )

AK4jetsequencePuppiData = cms.Sequence(
    AK4PFJetsPuppi*   
    AK4caPFJetsPrunedPuppi*
    AK4caPFJetsTrimmedPuppi*
    AK4caPFJetsSoftDropPuppi*
    ak4PuppiL1FastL2L3ResidualChain* #   => using type 1 Met
    #ak4PuppiL1FastL2L3Chain* #   => using type 1 Met
    AK4QGTaggerPuppi*
    AK4QGTaggerSubJetsPuppi*                
    AK4NjettinessPuppi
  )

def setMiniAODAK4Puppi(process) :
    process.AK4PFImpactParameterTagInfosPuppi.primaryVertex                    = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.AK4PFImpactParameterTagInfosPuppi.candidates                       = cms.InputTag("packedPFCandidates")
    process.AK4PFInclusiveSecondaryVertexFinderTagInfosPuppi.extSVCollection   = cms.InputTag('slimmedSecondaryVertices')
    process.AK4PFImpactParameterTagInfosSJPuppi.primaryVertex                  = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.AK4PFImpactParameterTagInfosSJPuppi.candidates                     = cms.InputTag("packedPFCandidates")
    process.AK4PFInclusiveSecondaryVertexFinderTagInfosSJPuppi.extSVCollection = cms.InputTag('slimmedSecondaryVertices')


    
    
