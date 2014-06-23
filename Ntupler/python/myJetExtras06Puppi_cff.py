import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK6PFJetsPuppi = ak5PFJets.clone(
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK6caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK6jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK6jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK6PFJetsPuppi')
AK6jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK6jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK6jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK6caPFJetsPrunedPuppi','SubJets')
AK6jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK6jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK6jetImpactParameterTagInfosPuppi.jetTracks            = "AK6jetTracksAssociatorAtVertexPuppi"
AK6jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK6jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK6jetImpactParameterTagInfosPuppi"
AK6jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK6jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK6jetImpactParameterTagInfosPuppi"), cms.InputTag("AK6jetSecondaryVertexTagInfosPuppi") )

AK6jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK6jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK6jetTracksAssociatorAtVertexSJPuppi"
AK6jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK6jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK6jetImpactParameterTagInfosSJPuppi"
AK6jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK6jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK6jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK6jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK6QGTaggerPuppi                                       = QGTagger.clone()
AK6QGTaggerPuppi.srcJets                               = cms.InputTag('AK6PFJetsPuppi')
AK6QGTaggerSubJetsPuppi                                = AK6QGTaggerPuppi.clone()
AK6QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK6caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK6NjettinessPuppi                                     = Njettiness.clone()       
AK6NjettinessPuppi.src                                 =  cms.InputTag('AK6PFJetsPuppi')


AK6jetsequencePuppi = cms.Sequence(
    AK6PFJetsPuppi                       *
    AK6caPFJetsPrunedPuppi               *
    AK6jetTracksAssociatorAtVertexPuppi  *
    AK6jetImpactParameterTagInfosPuppi   *
    AK6jetSecondaryVertexTagInfosPuppi   *
    AK6jetTracksAssociatorAtVertexSJPuppi*
    AK6jetImpactParameterTagInfosSJPuppi  *
    AK6jetSecondaryVertexTagInfosSJPuppi  *
    AK6jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK6jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK6QGTaggerPuppi                     *
    AK6QGTaggerSubJetsPuppi              *                
    AK6NjettinessPuppi                    
    )
