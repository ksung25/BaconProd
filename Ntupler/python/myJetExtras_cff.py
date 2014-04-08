import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.Configuration.RecoGenJets_cff     import ak5GenJets
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

antiktGenJets = ak5GenJets.clone(
    rParam = cms.double(0.5)
    )

# Flavour byValue PhysDef
AK5byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK5byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              
# Flavour byReference
partons  = cms.EDProducer("PartonSelector",
                          withLeptons = cms.bool(False),
                          src = cms.InputTag("genParticles")
                          )

AK5byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("ak5PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK5byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK5byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

jetFlavor    = cms.Sequence(partons*AK5byRef*AK5byValPhys*AK5byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
ca5PFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.5),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
jetTracksAssociatorAtVertex  .jets   = cms.InputTag('ak5PFJets')
jetTracksAssociatorAtVertex  .tracks = "generalTracks"

jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('ca5PFJetsPruned','SubJets')
jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
jetImpactParameterTagInfos.jetTracks            = "jetTracksAssociatorAtVertex"
jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
jetSecondaryVertexTagInfos.trackIPTagInfos      = "jetImpactParameterTagInfos"
jetCombinedSecondaryVertexMVABJetTags           = combinedSecondaryVertexMVABJetTags.clone()
jetCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("jetImpactParameterTagInfos"), cms.InputTag("jetSecondaryVertexTagInfos") )

jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
jetImpactParameterTagInfosSJ.jetTracks         = "jetTracksAssociatorAtVertexSJ"
jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "jetImpactParameterTagInfosSJ"
jetCombinedSecondaryVertexMVABJetTagsSJ          = combinedSecondaryVertexMVABJetTags.clone()
jetCombinedSecondaryVertexMVABJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("jetImpactParameterTagInfosSJ"), cms.InputTag("jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
QGTagger.srcJets                               = cms.InputTag('ak5PFJets')
QGTaggerSubJets                                = QGTagger.clone()
QGTaggerSubJets.srcJets                        = cms.InputTag('ca5PFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
Njettiness.src                                 =  cms.InputTag('ak5PFJets')


genjetsequence = cms.Sequence(
    genParticlesForJets            *
    genParticlesForJetsNoNu        *
    ak5GenJets                     * 
    jetFlavor                      )

jetsequence = cms.Sequence(
    ca5PFJetsPruned                *
    jetTracksAssociatorAtVertex    *
    jetImpactParameterTagInfos     *
    jetSecondaryVertexTagInfos     *
    jetTracksAssociatorAtVertexSJ  *
    jetImpactParameterTagInfosSJ   *
    jetSecondaryVertexTagInfosSJ   *
    jetCombinedSecondaryVertexMVABJetTags * 
    jetCombinedSecondaryVertexMVABJetTagsSJ  *
    goodOfflinePrimaryVerticesQG   *
    kt6PFJetsQG                    *
    kt6PFJetsIsoQG                 *
    QGTagger                       *
    QGTaggerSubJets                *                
    Njettiness                     
    )
