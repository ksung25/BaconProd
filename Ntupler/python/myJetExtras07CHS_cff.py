import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK7PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('PFBRECO','pfNoElectron'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK7caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src      = cms.InputTag('PFBRECO','pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK7jetTracksAssociatorAtVertexCHS          = ic5JetTracksAssociatorAtVertex.clone()
AK7jetTracksAssociatorAtVertexCHS  .jets   = cms.InputTag('AK7PFJetsCHS')
AK7jetTracksAssociatorAtVertexCHS  .tracks = "generalTracks"

AK7jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
AK7jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('AK7caPFJetsPrunedCHS','SubJets')
AK7jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK7jetImpactParameterTagInfosCHS                      = impactParameterTagInfos.clone()
AK7jetImpactParameterTagInfosCHS.jetTracks            = "AK7jetTracksAssociatorAtVertexCHS"
AK7jetSecondaryVertexTagInfosCHS                      = secondaryVertexTagInfos.clone()
AK7jetSecondaryVertexTagInfosCHS.trackIPTagInfos      = "AK7jetImpactParameterTagInfosCHS"
AK7jetCombinedSecondaryVertexBJetTagsCHS           = combinedSecondaryVertexBJetTags.clone()
AK7jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("AK7jetImpactParameterTagInfosCHS"), cms.InputTag("AK7jetSecondaryVertexTagInfosCHS") )

AK7jetImpactParameterTagInfosSJCHS                     = impactParameterTagInfos.clone()
AK7jetImpactParameterTagInfosSJCHS.jetTracks           = "AK7jetTracksAssociatorAtVertexSJCHS"
AK7jetSecondaryVertexTagInfosSJCHS                     = secondaryVertexTagInfos.clone()
AK7jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos     = "AK7jetImpactParameterTagInfosSJCHS"
AK7jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
AK7jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("AK7jetImpactParameterTagInfosSJCHS"), cms.InputTag("AK7jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK7QGTaggerCHS                                       = QGTagger.clone()
AK7QGTaggerCHS.srcJets                               = cms.InputTag('AK7PFJetsCHS')
AK7QGTaggerSubJetsCHS                                = AK7QGTaggerCHS.clone()
AK7QGTaggerSubJetsCHS.srcJets                        = cms.InputTag('AK7caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK7NjettinessCHS                                     = Njettiness.clone()       
AK7NjettinessCHS.src                                 =  cms.InputTag('AK7PFJetsCHS')


AK7jetsequenceCHS = cms.Sequence(
    AK7PFJetsCHS                      *
    AK7caPFJetsPrunedCHS              *
    AK7jetTracksAssociatorAtVertexCHS    *
    AK7jetImpactParameterTagInfosCHS     *
    AK7jetSecondaryVertexTagInfosCHS     *
    AK7jetTracksAssociatorAtVertexSJCHS  *
    AK7jetImpactParameterTagInfosSJCHS   *
    AK7jetSecondaryVertexTagInfosSJCHS   *
    AK7jetCombinedSecondaryVertexBJetTagsCHS * 
    AK7jetCombinedSecondaryVertexBJetTagsSJCHS  *
    AK7QGTaggerCHS                       *
    AK7QGTaggerSubJetsCHS                *                
    AK7NjettinessCHS                     
    )
