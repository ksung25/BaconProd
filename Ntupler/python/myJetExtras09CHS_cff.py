import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK9PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('pfNoElectron'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK9caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src      = cms.InputTag('pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK9jetTracksAssociatorAtVertexCHS          = ic5JetTracksAssociatorAtVertex.clone()
AK9jetTracksAssociatorAtVertexCHS  .jets   = cms.InputTag('AK9PFJetsCHS')
AK9jetTracksAssociatorAtVertexCHS  .tracks = "generalTracks"

AK9jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
AK9jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('AK9caPFJetsPrunedCHS','SubJets')
AK9jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK9jetImpactParameterTagInfosCHS                      = impactParameterTagInfos.clone()
AK9jetImpactParameterTagInfosCHS.jetTracks            = "AK9jetTracksAssociatorAtVertexCHS"
AK9jetSecondaryVertexTagInfosCHS                      = secondaryVertexTagInfos.clone()
AK9jetSecondaryVertexTagInfosCHS.trackIPTagInfos      = "AK9jetImpactParameterTagInfosCHS"
AK9jetCombinedSecondaryVertexBJetTagsCHS           = combinedSecondaryVertexBJetTags.clone()
AK9jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("AK9jetImpactParameterTagInfosCHS"), cms.InputTag("AK9jetSecondaryVertexTagInfosCHS") )

AK9jetImpactParameterTagInfosSJCHS                     = impactParameterTagInfos.clone()
AK9jetImpactParameterTagInfosSJCHS.jetTracks           = "AK9jetTracksAssociatorAtVertexSJCHS"
AK9jetSecondaryVertexTagInfosSJCHS                     = secondaryVertexTagInfos.clone()
AK9jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos     = "AK9jetImpactParameterTagInfosSJCHS"
AK9jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
AK9jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("AK9jetImpactParameterTagInfosSJCHS"), cms.InputTag("AK9jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK9QGTaggerCHS                                       = QGTagger.clone()
AK9QGTaggerCHS.srcJets                               = cms.InputTag('AK9PFJetsCHS')
AK9QGTaggerSubJetsCHS                                = AK9QGTaggerCHS.clone()
AK9QGTaggerSubJetsCHS.srcJets                        = cms.InputTag('AK9caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK9NjettinessCHS                                     = Njettiness.clone()       
AK9NjettinessCHS.src                                 =  cms.InputTag('AK9PFJetsCHS')


AK9jetsequenceCHS = cms.Sequence(
    AK9PFJetsCHS                      *
    AK9caPFJetsPrunedCHS              *
    AK9jetTracksAssociatorAtVertexCHS    *
    AK9jetImpactParameterTagInfosCHS     *
    AK9jetSecondaryVertexTagInfosCHS     *
    AK9jetTracksAssociatorAtVertexSJCHS  *
    AK9jetImpactParameterTagInfosSJCHS   *
    AK9jetSecondaryVertexTagInfosSJCHS   *
    AK9jetCombinedSecondaryVertexBJetTagsCHS * 
    AK9jetCombinedSecondaryVertexBJetTagsSJCHS  *
    AK9QGTaggerCHS                       *
    AK9QGTaggerSubJetsCHS                *                
    AK9NjettinessCHS                     
    )
