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

AK12GenJets = ak5GenJets.clone(
    rParam = cms.double(1.2)
    )

# Flavour byValue PhysDef
AK12byValPhys = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK12byRef"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK12byRef = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK12PFJets"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK12byValAlgo = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK12byRef"),
                              physicsDefinition = cms.bool(False),
                              leptonInfo = cms.bool(True))

AK12jetFlavor    = cms.Sequence(AK12byRef*AK12byValPhys*AK12byValAlgo)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK12PFJets = ak5PFJets.clone(
    rParam = cms.double(1.2),
    jetPtMin = cms.double(20)
    )

AK12caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(1.2),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK12jetTracksAssociatorAtVertex          = ic5JetTracksAssociatorAtVertex.clone()
AK12jetTracksAssociatorAtVertex  .jets   = cms.InputTag('AK12PFJets')
AK12jetTracksAssociatorAtVertex  .tracks = "generalTracks"

AK12jetTracksAssociatorAtVertexSJ        = ic5JetTracksAssociatorAtVertex.clone()
AK12jetTracksAssociatorAtVertexSJ.jets   = cms.InputTag('AK12caPFJetsPruned','SubJets')
AK12jetTracksAssociatorAtVertexSJ.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK12jetImpactParameterTagInfos                      = impactParameterTagInfos.clone()
AK12jetImpactParameterTagInfos.jetTracks            = "AK12jetTracksAssociatorAtVertex"
AK12jetSecondaryVertexTagInfos                      = secondaryVertexTagInfos.clone()
AK12jetSecondaryVertexTagInfos.trackIPTagInfos      = "AK12jetImpactParameterTagInfos"
AK12jetCombinedSecondaryVertexBJetTags           = combinedSecondaryVertexBJetTags.clone()
AK12jetCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("AK12jetImpactParameterTagInfos"), cms.InputTag("AK12jetSecondaryVertexTagInfos") )

AK12jetImpactParameterTagInfosSJ                   = impactParameterTagInfos.clone()
AK12jetImpactParameterTagInfosSJ.jetTracks         = "AK12jetTracksAssociatorAtVertexSJ"
AK12jetSecondaryVertexTagInfosSJ                   = secondaryVertexTagInfos.clone()
AK12jetSecondaryVertexTagInfosSJ.trackIPTagInfos     = "AK12jetImpactParameterTagInfosSJ"
AK12jetCombinedSecondaryVertexBJetTagsSJ          = combinedSecondaryVertexBJetTags.clone()
AK12jetCombinedSecondaryVertexBJetTagsSJ.tagInfos = cms.VInputTag( cms.InputTag("AK12jetImpactParameterTagInfosSJ"), cms.InputTag("AK12jetSecondaryVertexTagInfosSJ") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK12QGTagger                                       = QGTagger.clone()
AK12QGTagger.srcJets                               = cms.InputTag('AK12PFJets')
AK12QGTaggerSubJets                                = AK12QGTagger.clone()
AK12QGTaggerSubJets.srcJets                        = cms.InputTag('AK12caPFJetsPruned','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK12Njettiness                                     = Njettiness.clone()       
AK12Njettiness.src                                 =  cms.InputTag('AK12PFJets')


genjetsequence = cms.Sequence(
    genParticlesForJets            *
    genParticlesForJetsNoNu        *
    partons *
    goodOfflinePrimaryVerticesQG   *
    kt6PFJetsQG                    *
    kt6PFJetsIsoQG                 *
    ak5GenJets)

AK12genjetsequence = cms.Sequence(
    AK12GenJets                     * 
    AK12jetFlavor                   
)
AK12jetsequence = cms.Sequence(
    AK12PFJets                      *
    AK12caPFJetsPruned              *
    AK12jetTracksAssociatorAtVertex    *
    AK12jetImpactParameterTagInfos     *
    AK12jetSecondaryVertexTagInfos     *
    AK12jetTracksAssociatorAtVertexSJ  *
    AK12jetImpactParameterTagInfosSJ   *
    AK12jetSecondaryVertexTagInfosSJ   *
    AK12jetCombinedSecondaryVertexBJetTags * 
    AK12jetCombinedSecondaryVertexBJetTagsSJ  *
    AK12QGTagger                       *
    AK12QGTaggerSubJets                *                
    AK12Njettiness                     
    )
