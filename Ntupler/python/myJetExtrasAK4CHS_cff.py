import FWCore.ParameterSet.Config as cms

# take the default AK4 GenJet producer and modify accordingly for cone size and clustering algorithm
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
AK4GenJetsCHS = ak4GenJetsNoNu.clone(
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4)
  )

# Jet Flavour
AK4FlavorCHS = cms.EDProducer("JetFlavourClustering",
    jets                     = cms.InputTag("ak4PFJetsCHS"),
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
AK4PFJetsCHS = ak4PFJets.clone(
    src          = cms.InputTag('pfNoPileUpJME'),
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    jetPtMin     = cms.double(20)
  )

# Pruned
AK4caPFJetsPrunedCHS = AK4PFJetsCHS.clone(
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
AK4caPFJetsTrimmedCHS = AK4PFJetsCHS.clone(
    jetAlgorithm        = cms.string("CambridgeAachen"),
    useTrimming         = cms.bool(True),
    rFilt               = cms.double(0.2),
    trimPtFracMin       = cms.double(0.06),
    useExplicitGhosts   = cms.bool(True),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets")
  )

# Soft-Drop
AK4caPFJetsSoftDropCHS = AK4PFJetsCHS.clone(
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
from JetMETCorrections.Configuration.JetCorrectors_cff import *
from BaconProd.Ntupler.myCHSCorrections_cff          import *
from RecoJets.JetProducers.QGTagger_cfi import *
AK4QGTaggerCHS           = QGTagger.clone()
AK4QGTaggerCHS.srcJets   = cms.InputTag('ak4PFJetsCHS')
AK4QGTaggerCHS.jetsLabel = cms.string('QGL_AK4PFchs')
AK4QGTaggerCHS.jec       = cms.InputTag("ak4PFCHSL1FastL2L3Corrector") # NOTE: use "ak4PFCHSL1FastL2L3Corrector" for MC / "ak4PFCHSL1FastL2L3ResidualCorrector" for Data

AK4QGTaggerSubJetsCHS           = AK4QGTaggerCHS.clone()
AK4QGTaggerSubJetsCHS.srcJets   = cms.InputTag('AK4caPFJetsSoftDropCHS','SubJets')
AK4QGTaggerSubJetsCHS.jetsLabel = cms.string('QGL_AK4PFchs')
AK4QGTaggerSubJetsCHS.jec       = cms.InputTag("ak4PFCHSL1FastL2L3Corrector") # NOTE: use "ak4PFCHSL1FastL2L3Corrector" for MC / "ak4PFCHSL1FastL2L3ResidualCorrector" for Data

# N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
AK4NjettinessCHS = Njettiness.clone(
    src   = cms.InputTag('ak4PFJetsCHS'),
    cone  = cms.double(0.4),
    Njets = cms.vuint32(1,2,3,4)
  )

#
# Define sequences
#
AK4genjetsequenceCHS = cms.Sequence(
  AK4GenJetsCHS 
  #AK4FlavorCHS
)

AK4jetsequenceCHS = cms.Sequence(
    AK4PFJetsCHS*   ### no need to run, already in AOD
    AK4caPFJetsPrunedCHS*
    AK4caPFJetsTrimmedCHS*
    AK4caPFJetsSoftDropCHS*
#    AK4PFImpactParameterTagInfosCHS*                     ### no need to run, already in AOD
#    AK4PFInclusiveSecondaryVertexFinderTagInfosCHS*      ###
#    AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS*  ###
#    AK4PFImpactParameterTagInfosSJCHS*
#    AK4PFInclusiveSecondaryVertexFinderTagInfosSJCHS*
#    AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS*
    AK4FlavorCHS*
    ak4PFCHSL1FastL2L3CorrectorChain*
    AK4QGTaggerCHS*
    AK4QGTaggerSubJetsCHS*                
    AK4NjettinessCHS
  )

AK4jetsequenceCHSData = cms.Sequence(
#    AK4PFJetsCHS*   ### no need to run, already in AOD
#    AK4caPFJetsPrunedCHS*
#    AK4caPFJetsTrimmedCHS*
#    AK4caPFJetsSoftDropCHS*
#    AK4PFImpactParameterTagInfosCHS*                     ### no need to run, already in AOD
#    AK4PFInclusiveSecondaryVertexFinderTagInfosCHS*      ###
#    AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS*  ###
#    AK4PFImpactParameterTagInfosSJCHS*
#    AK4PFInclusiveSecondaryVertexFinderTagInfosSJCHS*
#    AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS*
#    ak4PFCHSL1FastL2L3ResidualCorrectorChain*
    ak4chsL1FastL2L3ResidualChain* 
    AK4QGTaggerCHS#*
#    AK4QGTaggerSubJetsCHS*                
#    AK4NjettinessCHS
  )

def setMiniAODAK4CHS(process) :
    process.load("RecoJets/JetProducers/QGTagger_cfi")
    process.QGTagger.srcJets          = cms.InputTag('slimmedJets')
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')
    process.QGTagger.jec              = cms.InputTag('')
    process.QGTagger.systematicsLabel = cms.string('')
    process.AK4FlavorCHS.jets         = cms.InputTag('slimmedJets')
    
