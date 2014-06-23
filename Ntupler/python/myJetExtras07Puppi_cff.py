import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK7PFJetsPuppi = ak5PFJets.clone(
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK7caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK7jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK7jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK7PFJetsPuppi')
AK7jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK7jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK7jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK7caPFJetsPrunedPuppi','SubJets')
AK7jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK7jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK7jetImpactParameterTagInfosPuppi.jetTracks            = "AK7jetTracksAssociatorAtVertexPuppi"
AK7jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK7jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK7jetImpactParameterTagInfosPuppi"
AK7jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK7jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK7jetImpactParameterTagInfosPuppi"), cms.InputTag("AK7jetSecondaryVertexTagInfosPuppi") )

AK7jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK7jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK7jetTracksAssociatorAtVertexSJPuppi"
AK7jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK7jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK7jetImpactParameterTagInfosSJPuppi"
AK7jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK7jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK7jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK7jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK7QGTaggerPuppi                                       = QGTagger.clone()
AK7QGTaggerPuppi.srcJets                               = cms.InputTag('AK7PFJetsPuppi')
AK7QGTaggerSubJetsPuppi                                = AK7QGTaggerPuppi.clone()
AK7QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK7caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK7NjettinessPuppi                                     = Njettiness.clone()       
AK7NjettinessPuppi.src                                 =  cms.InputTag('AK7PFJetsPuppi')


AK7jetsequencePuppi = cms.Sequence(
    AK7PFJetsPuppi                       *
    AK7caPFJetsPrunedPuppi               *
    AK7jetTracksAssociatorAtVertexPuppi  *
    AK7jetImpactParameterTagInfosPuppi   *
    AK7jetSecondaryVertexTagInfosPuppi   *
    AK7jetTracksAssociatorAtVertexSJPuppi*
    AK7jetImpactParameterTagInfosSJPuppi  *
    AK7jetSecondaryVertexTagInfosSJPuppi  *
    AK7jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK7jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK7QGTaggerPuppi                     *
    AK7QGTaggerSubJetsPuppi              *                
    AK7NjettinessPuppi                    
    )
