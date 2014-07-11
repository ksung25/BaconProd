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

AK5GenJets = ak5GenJets.clone(
    rParam = cms.double(0.5)
    )

# Flavour byValue PhysDef
AK5byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK5byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK5byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK5PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK5byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK5byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK5jetFlavor    = cms.Sequence(AK5byRef*AK5byValPhys*AK5byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK5PFJets = ak5PFJets.clone(
    rParam = cms.double(0.5),
    jetPtMin = cms.double(20)
    )

AK5caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.5),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK5jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK5jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK5PFJets')
AK5jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK5jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK5jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK5caPFJetsPruned','SubJets')
AK5jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK5jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK5jetImpactParameterTagInfos.jetTracks            = "AK5jetTracksAssociatorAtVertex"
AK5jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK5jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK5jetImpactParameterTagInfos"
AK5jetCombinedSecondaryVertexBJetTags           = combinedSecondaryVertexBJetTags.clone()
AK5jetCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK5jetImpactParameterTagInfos"), cms.InputTag("AK5jetSecondaryVertexTagInfos") )

AK5jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK5jetImpactParameterTagInfosSJ.jetTracks         = "AK5jetTracksAssociatorAtVertexSJ"
AK5jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK5jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK5jetImpactParameterTagInfosSJ"
AK5jetCombinedSecondaryVertexBJetTagsSJ          = combinedSecondaryVertexBJetTags.clone()
AK5jetCombinedSecondaryVertexBJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK5jetImpactParameterTagInfosSJ"), cms.InputTag("AK5jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK5QGTagger                                       = QGTagger.clone()
AK5QGTagger.srcJets                               = cms.InputTag('AK5PFJets')
AK5QGTaggerSubJets                                = AK5QGTagger.clone()
AK5QGTaggerSubJets.srcJets                        = cms.InputTag('AK5caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK5Njettiness                                     = Njettiness.clone()       
AK5Njettiness.src                                 =  cms.InputTag('AK5PFJets')


genjetsequence = cms.Sequence(
    genParticlesForJets            *
    genParticlesForJetsNoNu        *
    partons *
    ak5GenJets)

recojetsequence = cms.Sequence( 
    goodOfflinePrimaryVerticesQG   *
    kt6PFJetsQG                    *
    kt6PFJetsIsoQG                 
    #kt6PFJets                   
)

AK5genjetsequence = cms.Sequence(
    AK5GenJets                     * 
    AK5jetFlavor                   
)

AK5jetsequence = cms.Sequence(
    AK5PFJets                      *
    AK5caPFJetsPruned              *
    AK5jetTracksAssociatorAtVertex    *
    AK5jetImpactParameterTagInfos     *
    AK5jetSecondaryVertexTagInfos     *
    AK5jetTracksAssociatorAtVertexSJ  *
    AK5jetImpactParameterTagInfosSJ   *
    AK5jetSecondaryVertexTagInfosSJ   *
    AK5jetCombinedSecondaryVertexBJetTags * 
    AK5jetCombinedSecondaryVertexBJetTagsSJ  *
    AK5QGTagger                       *
    AK5QGTaggerSubJets                *                
    AK5Njettiness                     
    )
