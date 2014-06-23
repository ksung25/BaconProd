import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.Configuration.RecoGenJets_cff     import ak5GenJets
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(
         1000022,
         1000012, 1000014, 1000016,
         2000012, 2000014, 2000016,
         1000039, 5100039,
         4000012, 4000014, 4000016,
         9900012, 9900014, 9900016,
         39),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)
genParticlesForJetsNoNu = genParticlesForJets.clone()
genParticlesForJetsNoNu.ignoreParticleIDs += cms.vuint32( 12,14,16)

# Flavour byReference
partons  = cms.EDProducer("PartonSelector",
                          withLeptons = cms.bool(False),
                          src = cms.InputTag("genParticles")
                          )

AK4GenJets = ak5GenJets.clone(
    rParam = cms.double(0.4)
    )

# Flavour byValue PhysDef
AK4byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK4byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK4byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK4PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK4byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK4byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK4jetFlavor    = cms.Sequence(AK4byRef*AK4byValPhys*AK4byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK4PFJets = ak5PFJets.clone(
    rParam = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK4caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK4jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK4jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK4PFJets')
AK4jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK4jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK4jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK4caPFJetsPruned','SubJets')
AK4jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK4jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK4jetImpactParameterTagInfos.jetTracks            = "AK4jetTracksAssociatorAtVertex"
AK4jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK4jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK4jetImpactParameterTagInfos"
AK4jetCombinedSecondaryVertexBJetTags           = combinedSecondaryVertexBJetTags.clone()
AK4jetCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK4jetImpactParameterTagInfos"), cms.InputTag("AK4jetSecondaryVertexTagInfos") )

AK4jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK4jetImpactParameterTagInfosSJ.jetTracks         = "AK4jetTracksAssociatorAtVertexSJ"
AK4jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK4jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK4jetImpactParameterTagInfosSJ"
AK4jetCombinedSecondaryVertexBJetTagsSJ          = combinedSecondaryVertexBJetTags.clone()
AK4jetCombinedSecondaryVertexBJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK4jetImpactParameterTagInfosSJ"), cms.InputTag("AK4jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK4QGTagger                                       = QGTagger.clone()
AK4QGTagger.srcJets                               = cms.InputTag('AK4PFJets')
AK4QGTaggerSubJets                                = AK4QGTagger.clone()
AK4QGTaggerSubJets.srcJets                        = cms.InputTag('AK4caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK4Njettiness                                     = Njettiness.clone()       
AK4Njettiness.src                                 =  cms.InputTag('AK4PFJets')


genjetsequence = cms.Sequence(
    genParticlesForJets            *
    genParticlesForJetsNoNu        *
    partons *
    goodOfflinePrimaryVerticesQG   *
    kt6PFJetsQG                    *
    kt6PFJetsIsoQG                 *
    ak5GenJets)

AK4genjetsequence = cms.Sequence(
    AK4GenJets                     * 
    AK4jetFlavor                   
)
AK4jetsequence = cms.Sequence(
    AK4PFJets                      *
    AK4caPFJetsPruned              *
    AK4jetTracksAssociatorAtVertex    *
    AK4jetImpactParameterTagInfos     *
    AK4jetSecondaryVertexTagInfos     *
    AK4jetTracksAssociatorAtVertexSJ  *
    AK4jetImpactParameterTagInfosSJ   *
    AK4jetSecondaryVertexTagInfosSJ   *
    AK4jetCombinedSecondaryVertexBJetTags * 
    AK4jetCombinedSecondaryVertexBJetTagsSJ  *
    AK4QGTagger                       *
    AK4QGTaggerSubJets                *                
    AK4Njettiness                     
    )
