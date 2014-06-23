import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK12byRefPuppi = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK12PFJetsPuppi"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
AK12byValPhysPuppi = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK12byRefPuppi"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK12byValAlgoPuppi = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK12byRefPuppi"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK12jetFlavorPuppi    = cms.Sequence(AK12byRefPuppi*AK12byValPhysPuppi*AK12byValAlgoPuppi)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK12PFJetsPuppi = ak5PFJets.clone(
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(1.2),
    jetPtMin = cms.double(20)
    )

AK12caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(1.2),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK12jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK12jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK12PFJetsPuppi')
AK12jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK12jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK12jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK12caPFJetsPrunedPuppi','SubJets')
AK12jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK12jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK12jetImpactParameterTagInfosPuppi.jetTracks            = "AK12jetTracksAssociatorAtVertexPuppi"
AK12jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK12jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK12jetImpactParameterTagInfosPuppi"
AK12jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK12jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK12jetImpactParameterTagInfosPuppi"), cms.InputTag("AK12jetSecondaryVertexTagInfosPuppi") )

AK12jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK12jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK12jetTracksAssociatorAtVertexSJPuppi"
AK12jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK12jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK12jetImpactParameterTagInfosSJPuppi"
AK12jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK12jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK12jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK12jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK12QGTaggerPuppi                                       = QGTagger.clone()
AK12QGTaggerPuppi.srcJets                               = cms.InputTag('AK12PFJetsPuppi')
AK12QGTaggerSubJetsPuppi                                = AK12QGTaggerPuppi.clone()
AK12QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK12caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK12NjettinessPuppi                                     = Njettiness.clone()       
AK12NjettinessPuppi.src                                 =  cms.InputTag('AK12PFJetsPuppi')


AK12jetsequencePuppi = cms.Sequence(
    AK12PFJetsPuppi                       *
    AK12caPFJetsPrunedPuppi               *
    AK12jetTracksAssociatorAtVertexPuppi  *
    AK12jetImpactParameterTagInfosPuppi   *
    AK12jetSecondaryVertexTagInfosPuppi   *
    AK12jetTracksAssociatorAtVertexSJPuppi*
    AK12jetImpactParameterTagInfosSJPuppi  *
    AK12jetSecondaryVertexTagInfosSJPuppi  *
    AK12jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK12jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK12QGTaggerPuppi                     *
    AK12QGTaggerSubJetsPuppi              *                
    AK12NjettinessPuppi                   *
    AK12jetFlavorPuppi                       
    )
