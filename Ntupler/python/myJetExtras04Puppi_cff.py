import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK4byRefPuppi = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK4PFJetsPuppi"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
AK4byValPhysPuppi = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK4byRefPuppi"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK4byValAlgoPuppi = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK4byRefPuppi"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK4jetFlavorPuppi    = cms.Sequence(AK4byRefPuppi*AK4byValPhysPuppi*AK4byValAlgoPuppi)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK4PFJetsPuppi = ak5PFJets.clone(
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK4caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK4jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK4jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK4PFJetsPuppi')
AK4jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK4jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK4jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK4caPFJetsPrunedPuppi','SubJets')
AK4jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK4jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK4jetImpactParameterTagInfosPuppi.jetTracks            = "AK4jetTracksAssociatorAtVertexPuppi"
AK4jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK4jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK4jetImpactParameterTagInfosPuppi"
AK4jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK4jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK4jetImpactParameterTagInfosPuppi"), cms.InputTag("AK4jetSecondaryVertexTagInfosPuppi") )

AK4jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK4jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK4jetTracksAssociatorAtVertexSJPuppi"
AK4jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK4jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK4jetImpactParameterTagInfosSJPuppi"
AK4jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK4jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK4jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK4jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK4QGTaggerPuppi                                       = QGTagger.clone()
AK4QGTaggerPuppi.srcJets                               = cms.InputTag('AK4PFJetsPuppi')
AK4QGTaggerSubJetsPuppi                                = AK4QGTaggerPuppi.clone()
AK4QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK4caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK4NjettinessPuppi                                     = Njettiness.clone()       
AK4NjettinessPuppi.src                                 =  cms.InputTag('AK4PFJetsPuppi')


AK4jetsequencePuppi = cms.Sequence(
    AK4PFJetsPuppi                       *
    AK4caPFJetsPrunedPuppi               *
    AK4jetTracksAssociatorAtVertexPuppi  *
    AK4jetImpactParameterTagInfosPuppi   *
    AK4jetSecondaryVertexTagInfosPuppi   *
    AK4jetTracksAssociatorAtVertexSJPuppi*
    AK4jetImpactParameterTagInfosSJPuppi  *
    AK4jetSecondaryVertexTagInfosSJPuppi  *
    AK4jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK4jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK4QGTaggerPuppi                     *
    AK4QGTaggerSubJetsPuppi              *                
    AK4NjettinessPuppi                   *
    AK4jetFlavorPuppi                       
    )
