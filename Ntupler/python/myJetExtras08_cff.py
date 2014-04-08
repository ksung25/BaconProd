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

AK8GenJets = ak5GenJets.clone(
    rParam = cms.double(0.8)
    )

# Flavour byValue PhysDef
AK8byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK8byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK8byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK8PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK8byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK8byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK8jetFlavor    = cms.Sequence(AK8byRef*AK8byValPhys*AK8byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK8PFJets = ak5PFJets.clone(
    rParam = cms.double(0.8),
    jetPtMin = cms.double(20)
    )

AK8caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK8jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK8jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK8PFJets')
AK8jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK8jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK8jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK8caPFJetsPruned','SubJets')
AK8jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK8jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK8jetImpactParameterTagInfos.jetTracks            = "AK8jetTracksAssociatorAtVertex"
AK8jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK8jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK8jetImpactParameterTagInfos"
AK8jetCombinedSecondaryVertexMVABJetTags           = combinedSecondaryVertexMVABJetTags.clone()
AK8jetCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK8jetImpactParameterTagInfos"), cms.InputTag("AK8jetSecondaryVertexTagInfos") )

AK8jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK8jetImpactParameterTagInfosSJ.jetTracks         = "AK8jetTracksAssociatorAtVertexSJ"
AK8jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK8jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK8jetImpactParameterTagInfosSJ"
AK8jetCombinedSecondaryVertexMVABJetTagsSJ          = combinedSecondaryVertexMVABJetTags.clone()
AK8jetCombinedSecondaryVertexMVABJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK8jetImpactParameterTagInfosSJ"), cms.InputTag("AK8jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK8QGTagger                                       = QGTagger.clone()
AK8QGTagger.srcJets                               = cms.InputTag('AK8PFJets')
AK8QGTaggerSubJets                                = AK8QGTagger.clone()
AK8QGTaggerSubJets.srcJets                        = cms.InputTag('AK8caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK8Njettiness                                     = Njettiness.clone()       
AK8Njettiness.src                                 =  cms.InputTag('AK8PFJets')


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

AK8genjetsequence = cms.Sequence(
    AK8GenJets                     * 
    AK8jetFlavor                   
)
AK8jetsequence = cms.Sequence(
    AK8PFJets                      *
    AK8caPFJetsPruned              *
    AK8jetTracksAssociatorAtVertex    *
    AK8jetImpactParameterTagInfos     *
    AK8jetSecondaryVertexTagInfos     *
    AK8jetTracksAssociatorAtVertexSJ  *
    AK8jetImpactParameterTagInfosSJ   *
    AK8jetSecondaryVertexTagInfosSJ   *
    AK8jetCombinedSecondaryVertexMVABJetTags * 
    AK8jetCombinedSecondaryVertexMVABJetTagsSJ  *
    AK8QGTagger                       *
    AK8QGTaggerSubJets                *                
    AK8Njettiness                     
    )
