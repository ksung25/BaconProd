import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK5byRefPuppi = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK5PFJetsPuppi"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )

AK5byValPhysPuppi = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK5byRefPuppi"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK5byValAlgoPuppi = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK5byRefPuppi"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK5jetFlavorPuppi    = cms.Sequence(AK5byRefPuppi*AK5byValPhysPuppi*AK5byValAlgoPuppi)


#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK5PFJetsPuppi = ak5PFJets.clone(
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK5caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK5jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK5jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK5PFJetsPuppi')
AK5jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK5jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK5jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK5caPFJetsPrunedPuppi','SubJets')
AK5jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK5jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK5jetImpactParameterTagInfosPuppi.jetTracks            = "AK5jetTracksAssociatorAtVertexPuppi"
AK5jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK5jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK5jetImpactParameterTagInfosPuppi"
AK5jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK5jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK5jetImpactParameterTagInfosPuppi"), cms.InputTag("AK5jetSecondaryVertexTagInfosPuppi") )

AK5jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK5jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK5jetTracksAssociatorAtVertexSJPuppi"
AK5jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK5jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK5jetImpactParameterTagInfosSJPuppi"
AK5jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK5jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK5jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK5jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK5QGTaggerPuppi                                       = QGTagger.clone()
AK5QGTaggerPuppi.srcJets                               = cms.InputTag('AK5PFJetsPuppi')
AK5QGTaggerSubJetsPuppi                                = AK5QGTaggerPuppi.clone()
AK5QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK5caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK5NjettinessPuppi                                     = Njettiness.clone()       
AK5NjettinessPuppi.src                                 =  cms.InputTag('AK5PFJetsPuppi')


AK5jetsequencePuppi = cms.Sequence(
    AK5PFJetsPuppi                       *
    AK5caPFJetsPrunedPuppi               *
    AK5jetTracksAssociatorAtVertexPuppi  *
    AK5jetImpactParameterTagInfosPuppi   *
    AK5jetSecondaryVertexTagInfosPuppi   *
    AK5jetTracksAssociatorAtVertexSJPuppi*
    AK5jetImpactParameterTagInfosSJPuppi  *
    AK5jetSecondaryVertexTagInfosSJPuppi  *
    AK5jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK5jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK5QGTaggerPuppi                     *
    AK5QGTaggerSubJetsPuppi              *                
    AK5NjettinessPuppi                   *
    AK5jetFlavorPuppi           
    )
