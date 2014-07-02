import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK6PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('pfNoElectron'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK6caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src      = cms.InputTag('pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK6jetTracksAssociatorAtVertexCHS          = ic5JetTracksAssociatorAtVertex.clone()
AK6jetTracksAssociatorAtVertexCHS  .jets   = cms.InputTag('AK6PFJetsCHS')
AK6jetTracksAssociatorAtVertexCHS  .tracks = "generalTracks"

AK6jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
AK6jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('AK6caPFJetsPrunedCHS','SubJets')
AK6jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK6jetImpactParameterTagInfosCHS                      = impactParameterTagInfos.clone()
AK6jetImpactParameterTagInfosCHS.jetTracks            = "AK6jetTracksAssociatorAtVertexCHS"
AK6jetSecondaryVertexTagInfosCHS                      = secondaryVertexTagInfos.clone()
AK6jetSecondaryVertexTagInfosCHS.trackIPTagInfos      = "AK6jetImpactParameterTagInfosCHS"
AK6jetCombinedSecondaryVertexBJetTagsCHS           = combinedSecondaryVertexBJetTags.clone()
AK6jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("AK6jetImpactParameterTagInfosCHS"), cms.InputTag("AK6jetSecondaryVertexTagInfosCHS") )

AK6jetImpactParameterTagInfosSJCHS                     = impactParameterTagInfos.clone()
AK6jetImpactParameterTagInfosSJCHS.jetTracks           = "AK6jetTracksAssociatorAtVertexSJCHS"
AK6jetSecondaryVertexTagInfosSJCHS                     = secondaryVertexTagInfos.clone()
AK6jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos     = "AK6jetImpactParameterTagInfosSJCHS"
AK6jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
AK6jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("AK6jetImpactParameterTagInfosSJCHS"), cms.InputTag("AK6jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK6QGTaggerCHS                                       = QGTagger.clone()
AK6QGTaggerCHS.srcJets                               = cms.InputTag('AK6PFJetsCHS')
AK6QGTaggerSubJetsCHS                                = AK6QGTaggerCHS.clone()
AK6QGTaggerSubJetsCHS.srcJets                        = cms.InputTag('AK6caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK6NjettinessCHS                                     = Njettiness.clone()       
AK6NjettinessCHS.src                                 =  cms.InputTag('AK6PFJetsCHS')


AK6jetsequenceCHS = cms.Sequence(
    AK6PFJetsCHS                      *
    AK6caPFJetsPrunedCHS              *
    AK6jetTracksAssociatorAtVertexCHS    *
    AK6jetImpactParameterTagInfosCHS     *
    AK6jetSecondaryVertexTagInfosCHS     *
    AK6jetTracksAssociatorAtVertexSJCHS  *
    AK6jetImpactParameterTagInfosSJCHS   *
    AK6jetSecondaryVertexTagInfosSJCHS   *
    AK6jetCombinedSecondaryVertexBJetTagsCHS * 
    AK6jetCombinedSecondaryVertexBJetTagsSJCHS  *
    AK6QGTaggerCHS                       *
    AK6QGTaggerSubJetsCHS                *                
    AK6NjettinessCHS                     
    )
