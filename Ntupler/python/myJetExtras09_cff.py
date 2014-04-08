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

AK9GenJets = ak5GenJets.clone(
    rParam = cms.double(0.9)
    )

# Flavour byValue PhysDef
AK9byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK9byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK9byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK9PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK9byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK9byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK9jetFlavor    = cms.Sequence(AK9byRef*AK9byValPhys*AK9byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK9PFJets = ak5PFJets.clone(
    rParam = cms.double(0.9),
    jetPtMin = cms.double(20)
    )

AK9caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.9),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK9jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK9jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK9PFJets')
AK9jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK9jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK9jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK9caPFJetsPruned','SubJets')
AK9jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK9jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK9jetImpactParameterTagInfos.jetTracks            = "AK9jetTracksAssociatorAtVertex"
AK9jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK9jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK9jetImpactParameterTagInfos"
AK9jetCombinedSecondaryVertexMVABJetTags           = combinedSecondaryVertexMVABJetTags.clone()
AK9jetCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK9jetImpactParameterTagInfos"), cms.InputTag("AK9jetSecondaryVertexTagInfos") )

AK9jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK9jetImpactParameterTagInfosSJ.jetTracks         = "AK9jetTracksAssociatorAtVertexSJ"
AK9jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK9jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK9jetImpactParameterTagInfosSJ"
AK9jetCombinedSecondaryVertexMVABJetTagsSJ          = combinedSecondaryVertexMVABJetTags.clone()
AK9jetCombinedSecondaryVertexMVABJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK9jetImpactParameterTagInfosSJ"), cms.InputTag("AK9jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK9QGTagger                                       = QGTagger.clone()
AK9QGTagger.srcJets                               = cms.InputTag('AK9PFJets')
AK9QGTaggerSubJets                                = AK9QGTagger.clone()
AK9QGTaggerSubJets.srcJets                        = cms.InputTag('AK9caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK9Njettiness                                     = Njettiness.clone()       
AK9Njettiness.src                                 =  cms.InputTag('AK9PFJets')


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

AK9genjetsequence = cms.Sequence(
    AK9GenJets                     * 
    AK9jetFlavor                   
)
AK9jetsequence = cms.Sequence(
    AK9PFJets                      *
    AK9caPFJetsPruned              *
    AK9jetTracksAssociatorAtVertex    *
    AK9jetImpactParameterTagInfos     *
    AK9jetSecondaryVertexTagInfos     *
    AK9jetTracksAssociatorAtVertexSJ  *
    AK9jetImpactParameterTagInfosSJ   *
    AK9jetSecondaryVertexTagInfosSJ   *
    AK9jetCombinedSecondaryVertexMVABJetTags * 
    AK9jetCombinedSecondaryVertexMVABJetTagsSJ  *
    AK9QGTagger                       *
    AK9QGTaggerSubJets                *                
    AK9Njettiness                     
    )
