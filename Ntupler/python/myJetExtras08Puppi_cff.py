import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK8byRefPuppi = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK8PFJetsPuppi"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
AK8byValPhysPuppi = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK8byRefPuppi"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK8byValAlgoPuppi = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK8byRefPuppi"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK8jetFlavorPuppi    = cms.Sequence(AK8byRefPuppi*AK8byValPhysPuppi*AK8byValAlgoPuppi)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK8PFJetsPuppi = ak5PFJets.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(0.8),
    jetPtMin = cms.double(20)
    )

AK8caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK8jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK8jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK8PFJetsPuppi')
AK8jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK8jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK8jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK8caPFJetsPrunedPuppi','SubJets')
AK8jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK8jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK8jetImpactParameterTagInfosPuppi.jetTracks            = "AK8jetTracksAssociatorAtVertexPuppi"
AK8jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK8jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK8jetImpactParameterTagInfosPuppi"
AK8jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK8jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK8jetImpactParameterTagInfosPuppi"), cms.InputTag("AK8jetSecondaryVertexTagInfosPuppi") )

AK8jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK8jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK8jetTracksAssociatorAtVertexSJPuppi"
AK8jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK8jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK8jetImpactParameterTagInfosSJPuppi"
AK8jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK8jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK8jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK8jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK8QGTaggerPuppi                                       = QGTagger.clone()
AK8QGTaggerPuppi.srcJets                               = cms.InputTag('AK8PFJetsPuppi')
AK8QGTaggerSubJetsPuppi                                = AK8QGTaggerPuppi.clone()
AK8QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK8caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK8NjettinessPuppi                                     = Njettiness.clone()       
AK8NjettinessPuppi.src                                 =  cms.InputTag('AK8PFJetsPuppi')


AK8jetsequencePuppi = cms.Sequence(
    AK8PFJetsPuppi                       *
    AK8caPFJetsPrunedPuppi               *
    AK8jetTracksAssociatorAtVertexPuppi  *
    AK8jetImpactParameterTagInfosPuppi   *
    AK8jetSecondaryVertexTagInfosPuppi   *
    AK8jetTracksAssociatorAtVertexSJPuppi*
    AK8jetImpactParameterTagInfosSJPuppi  *
    AK8jetSecondaryVertexTagInfosSJPuppi  *
    AK8jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK8jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK8QGTaggerPuppi                     *
    AK8QGTaggerSubJetsPuppi              *                
    AK8NjettinessPuppi                   *
    AK8jetFlavorPuppi                       
    )
