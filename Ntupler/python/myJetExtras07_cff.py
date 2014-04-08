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

AK7GenJets = ak5GenJets.clone(
    rParam = cms.double(0.7)
    )

# Flavour byValue PhysDef
AK7byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK7byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK7byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK7PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK7byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK7byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK7jetFlavor    = cms.Sequence(AK7byRef*AK7byValPhys*AK7byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK7PFJets = ak5PFJets.clone(
    rParam = cms.double(0.7),
    jetPtMin = cms.double(20)
    )

AK7caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.7),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK7jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK7jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK7PFJets')
AK7jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK7jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK7jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK7caPFJetsPruned','SubJets')
AK7jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK7jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK7jetImpactParameterTagInfos.jetTracks            = "AK7jetTracksAssociatorAtVertex"
AK7jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK7jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK7jetImpactParameterTagInfos"
AK7jetCombinedSecondaryVertexMVABJetTags           = combinedSecondaryVertexMVABJetTags.clone()
AK7jetCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK7jetImpactParameterTagInfos"), cms.InputTag("AK7jetSecondaryVertexTagInfos") )

AK7jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK7jetImpactParameterTagInfosSJ.jetTracks         = "AK7jetTracksAssociatorAtVertexSJ"
AK7jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK7jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK7jetImpactParameterTagInfosSJ"
AK7jetCombinedSecondaryVertexMVABJetTagsSJ          = combinedSecondaryVertexMVABJetTags.clone()
AK7jetCombinedSecondaryVertexMVABJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK7jetImpactParameterTagInfosSJ"), cms.InputTag("AK7jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK7QGTagger                                       = QGTagger.clone()
AK7QGTagger.srcJets                               = cms.InputTag('AK7PFJets')
AK7QGTaggerSubJets                                = AK7QGTagger.clone()
AK7QGTaggerSubJets.srcJets                        = cms.InputTag('AK7caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK7Njettiness                                     = Njettiness.clone()       
AK7Njettiness.src                                 =  cms.InputTag('AK7PFJets')


genjetsequence = cms.Sequence(
    genParticlesForJets            *
    genParticlesForJetsNoNu        *
    partons *
    ak5GenJets)

recojetsequence = cms.Sequence( 
    goodOfflinePrimaryVerticesQG   *
    kt6PFJetsQG                    *
    kt6PFJetsIsoQG                 
)

AK7genjetsequence = cms.Sequence(
    AK7GenJets                     * 
    AK7jetFlavor                   
)
AK7jetsequence = cms.Sequence(
    AK7PFJets                      *
    AK7caPFJetsPruned              *
    AK7jetTracksAssociatorAtVertex    *
    AK7jetImpactParameterTagInfos     *
    AK7jetSecondaryVertexTagInfos     *
    AK7jetTracksAssociatorAtVertexSJ  *
    AK7jetImpactParameterTagInfosSJ   *
    AK7jetSecondaryVertexTagInfosSJ   *
    AK7jetCombinedSecondaryVertexMVABJetTags * 
    AK7jetCombinedSecondaryVertexMVABJetTagsSJ  *
    AK7QGTagger                       *
    AK7QGTaggerSubJets                *                
    AK7Njettiness                     
    )
