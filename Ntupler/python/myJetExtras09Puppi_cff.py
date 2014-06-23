import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK9PFJetsPuppi = ak5PFJets.clone(
    src      = cms.InputTag('puppi','Puppi'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK9caPFJetsPrunedPuppi = ak5PFJetsPruned.clone(
    src      = cms.InputTag('puppi','Puppi'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK9jetTracksAssociatorAtVertexPuppi          = ic5JetTracksAssociatorAtVertex.clone()
AK9jetTracksAssociatorAtVertexPuppi  .jets   = cms.InputTag('AK9PFJetsPuppi')
AK9jetTracksAssociatorAtVertexPuppi  .tracks = "generalTracks"

AK9jetTracksAssociatorAtVertexSJPuppi        = ic5JetTracksAssociatorAtVertex.clone()
AK9jetTracksAssociatorAtVertexSJPuppi.jets   = cms.InputTag('AK9caPFJetsPrunedPuppi','SubJets')
AK9jetTracksAssociatorAtVertexSJPuppi.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK9jetImpactParameterTagInfosPuppi                      = impactParameterTagInfos.clone()
AK9jetImpactParameterTagInfosPuppi.jetTracks            = "AK9jetTracksAssociatorAtVertexPuppi"
AK9jetSecondaryVertexTagInfosPuppi                      = secondaryVertexTagInfos.clone()
AK9jetSecondaryVertexTagInfosPuppi.trackIPTagInfos      = "AK9jetImpactParameterTagInfosPuppi"
AK9jetCombinedSecondaryVertexBJetTagsPuppi           = combinedSecondaryVertexBJetTags.clone()
AK9jetCombinedSecondaryVertexBJetTagsPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK9jetImpactParameterTagInfosPuppi"), cms.InputTag("AK9jetSecondaryVertexTagInfosPuppi") )

AK9jetImpactParameterTagInfosSJPuppi                     = impactParameterTagInfos.clone()
AK9jetImpactParameterTagInfosSJPuppi.jetTracks           = "AK9jetTracksAssociatorAtVertexSJPuppi"
AK9jetSecondaryVertexTagInfosSJPuppi                     = secondaryVertexTagInfos.clone()
AK9jetSecondaryVertexTagInfosSJPuppi.trackIPTagInfos     = "AK9jetImpactParameterTagInfosSJPuppi"
AK9jetCombinedSecondaryVertexBJetTagsSJPuppi          = combinedSecondaryVertexBJetTags.clone()
AK9jetCombinedSecondaryVertexBJetTagsSJPuppi.tagInfos = cms.VInputTag( cms.InputTag("AK9jetImpactParameterTagInfosSJPuppi"), cms.InputTag("AK9jetSecondaryVertexTagInfosSJPuppi") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK9QGTaggerPuppi                                       = QGTagger.clone()
AK9QGTaggerPuppi.srcJets                               = cms.InputTag('AK9PFJetsPuppi')
AK9QGTaggerSubJetsPuppi                                = AK9QGTaggerPuppi.clone()
AK9QGTaggerSubJetsPuppi.srcJets                        = cms.InputTag('AK9caPFJetsPrunedPuppi','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK9NjettinessPuppi                                     = Njettiness.clone()       
AK9NjettinessPuppi.src                                 =  cms.InputTag('AK9PFJetsPuppi')


AK9jetsequencePuppi = cms.Sequence(
    AK9PFJetsPuppi                       *
    AK9caPFJetsPrunedPuppi               *
    AK9jetTracksAssociatorAtVertexPuppi  *
    AK9jetImpactParameterTagInfosPuppi   *
    AK9jetSecondaryVertexTagInfosPuppi   *
    AK9jetTracksAssociatorAtVertexSJPuppi*
    AK9jetImpactParameterTagInfosSJPuppi  *
    AK9jetSecondaryVertexTagInfosSJPuppi  *
    AK9jetCombinedSecondaryVertexBJetTagsPuppi* 
    AK9jetCombinedSecondaryVertexBJetTagsSJPuppi*
    AK9QGTaggerPuppi                     *
    AK9QGTaggerSubJetsPuppi              *                
    AK9NjettinessPuppi                    
    )
