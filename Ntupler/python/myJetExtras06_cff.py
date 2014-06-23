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

AK6GenJets = ak5GenJets.clone(
    rParam = cms.double(0.6)
    )

# Flavour byValue PhysDef
AK6byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK6byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK6byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK6PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK6byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK6byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK6jetFlavor    = cms.Sequence(AK6byRef*AK6byValPhys*AK6byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK6PFJets = ak5PFJets.clone(
    rParam = cms.double(0.6),
    jetPtMin = cms.double(20)
    )

AK6caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.6),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK6jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK6jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK6PFJets')
AK6jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK6jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK6jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK6caPFJetsPruned','SubJets')
AK6jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK6jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK6jetImpactParameterTagInfos.jetTracks            = "AK6jetTracksAssociatorAtVertex"
AK6jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK6jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK6jetImpactParameterTagInfos"
AK6jetCombinedSecondaryVertexBJetTags           = combinedSecondaryVertexBJetTags.clone()
AK6jetCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK6jetImpactParameterTagInfos"), cms.InputTag("AK6jetSecondaryVertexTagInfos") )

AK6jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK6jetImpactParameterTagInfosSJ.jetTracks         = "AK6jetTracksAssociatorAtVertexSJ"
AK6jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK6jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK6jetImpactParameterTagInfosSJ"
AK6jetCombinedSecondaryVertexBJetTagsSJ          = combinedSecondaryVertexBJetTags.clone()
AK6jetCombinedSecondaryVertexBJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK6jetImpactParameterTagInfosSJ"), cms.InputTag("AK6jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK6QGTagger                                       = QGTagger.clone()
AK6QGTagger.srcJets                               = cms.InputTag('AK6PFJets')
AK6QGTaggerSubJets                                = AK6QGTagger.clone()
AK6QGTaggerSubJets.srcJets                        = cms.InputTag('AK6caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK6Njettiness                                     = Njettiness.clone()       
AK6Njettiness.src                                 =  cms.InputTag('AK6PFJets')


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

AK6genjetsequence = cms.Sequence(
    AK6GenJets                     * 
    AK6jetFlavor                   
)

AK6jetsequence = cms.Sequence(
    AK6PFJets                      *
    AK6caPFJetsPruned              *
    AK6jetTracksAssociatorAtVertex    *
    AK6jetImpactParameterTagInfos     *
    AK6jetSecondaryVertexTagInfos     *
    AK6jetTracksAssociatorAtVertexSJ  *
    AK6jetImpactParameterTagInfosSJ   *
    AK6jetSecondaryVertexTagInfosSJ   *
    AK6jetCombinedSecondaryVertexBJetTags * 
    AK6jetCombinedSecondaryVertexBJetTagsSJ  *
    AK6QGTagger                       *
    AK6QGTaggerSubJets                *                
    AK6Njettiness                     
    )
