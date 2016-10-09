import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
AK8GenJetsPuppi = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8)
  )

# Jet Flavour
AK8FlavorPuppi = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("AK8PFJetsPuppi"),
    groomedJets              = cms.InputTag("AK8caPFJetsSoftDropPuppi"),
    subjets                  = cms.InputTag("AK8caPFJetsSoftDropPuppi", "SubJets"),
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
AK8PFJetsPuppi = ak4PFJets.clone(
    src          = cms.InputTag('puppi'),
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.8),
    jetPtMin     = cms.double(150)
  )

# Pruned
AK8caPFJetsPrunedPuppi = AK8PFJetsPuppi.clone(
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
AK8caPFJetsTrimmedPuppi = AK8PFJetsPuppi.clone(
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
AK8caPFJetsSoftDropPuppi = AK8PFJetsPuppi.clone(
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
#from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff import *
from BaconProd.Ntupler.myPUPPICorrections_cff                  import *
from RecoJets.JetProducers.QGTagger_cfi import *
AK8QGTaggerPuppi           = QGTagger.clone()
AK8QGTaggerPuppi.srcJets   = cms.InputTag('AK8PFJetsPuppi')
AK8QGTaggerPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
AK8QGTaggerPuppi.jec       = cms.InputTag("ak8PuppiL1FastL2L3Corrector") # NOTE: use "ca8PFPuppiL1FastL2L3Corrector" for MC / "ca8PFPuppiL1FastL2L3ResidualCorrector" for Data

AK8QGTaggerSubJetsPuppi           = AK8QGTaggerPuppi.clone()
AK8QGTaggerSubJetsPuppi.srcJets   = cms.InputTag('AK8caPFJetsSoftDropPuppi','SubJets')
AK8QGTaggerSubJetsPuppi.jetsLabel = cms.string('QGL_AK4PFchs')
AK8QGTaggerSubJetsPuppi.jec       = cms.InputTag("ak8PuppiL1FastL2L3Corrector") # NOTE: use "ca8PFPuppiL1FastL2L3Corrector" for MC / "ca8PFPuppiL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
AK8NjettinessPuppi = Njettiness.clone(
    src   = cms.InputTag('AK8PFJetsPuppi'),
    cone  = cms.double(0.8),
    Njets = cms.vuint32(1,2,3,4)
  )

AK8genjetsequencePuppi = cms.Sequence(
  AK8GenJetsPuppi* 
  AK8FlavorPuppi
)

#
# Define sequences
#
AK8jetsequencePuppi = cms.Sequence(
    AK8PFJetsPuppi*
    AK8caPFJetsPrunedPuppi*
    AK8caPFJetsTrimmedPuppi*
    AK8caPFJetsSoftDropPuppi*
    ak8PuppiL1FastL2L3Chain*
    AK8QGTaggerPuppi*
    AK8QGTaggerSubJetsPuppi*                
    AK8NjettinessPuppi*
    AK8FlavorPuppi
    )

AK8jetsequencePuppiData = cms.Sequence(
    AK8PFJetsPuppi*
    AK8caPFJetsPrunedPuppi*
    AK8caPFJetsTrimmedPuppi*
    AK8caPFJetsSoftDropPuppi*
    ak8PuppiL1FastL2L3ResidualChain*
    #ak8PuppiL1FastL2L3Chain*
    AK8QGTaggerPuppi*
    AK8QGTaggerSubJetsPuppi*                
    AK8NjettinessPuppi
  )

def setMiniAODAK8Puppi(process) :
    process.AK8PFImpactParameterTagInfosPuppi.primaryVertex                    = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.AK8PFImpactParameterTagInfosPuppi.candidates                       = cms.InputTag("puppi")
    process.AK8PFInclusiveSecondaryVertexFinderTagInfosPuppi.extSVCollection   = cms.InputTag('slimmedSecondaryVertices')
    process.AK8PFImpactParameterTagInfosSJPuppi.primaryVertex                  = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.AK8PFImpactParameterTagInfosSJPuppi.candidates                     = cms.InputTag("puppi")
    process.AK8PFInclusiveSecondaryVertexFinderTagInfosSJPuppi.extSVCollection = cms.InputTag('slimmedSecondaryVertices')
