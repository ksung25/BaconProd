import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK8byRefCHS = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK8PFJetsCHS"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
AK8byValPhysCHS = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK8byRefCHS"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK8byValAlgoCHS = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK8byRefCHS"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK8jetFlavorCHS    = cms.Sequence(AK8byRefCHS*AK8byValPhysCHS*AK8byValAlgoCHS)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK8PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('pfNoElectron'),
    rParam   = cms.double(0.8),
    jetPtMin = cms.double(20)
    )

AK8caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src      = cms.InputTag('pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK8jetTracksAssociatorAtVertexCHS          = ic5JetTracksAssociatorAtVertex.clone()
AK8jetTracksAssociatorAtVertexCHS  .jets   = cms.InputTag('AK8PFJetsCHS')
AK8jetTracksAssociatorAtVertexCHS  .tracks = "generalTracks"

AK8jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
AK8jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('AK8caPFJetsPrunedCHS','SubJets')
AK8jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK8jetImpactParameterTagInfosCHS                      = impactParameterTagInfos.clone()
AK8jetImpactParameterTagInfosCHS.jetTracks            = "AK8jetTracksAssociatorAtVertexCHS"
AK8jetSecondaryVertexTagInfosCHS                      = secondaryVertexTagInfos.clone()
AK8jetSecondaryVertexTagInfosCHS.trackIPTagInfos      = "AK8jetImpactParameterTagInfosCHS"
AK8jetCombinedSecondaryVertexBJetTagsCHS           = combinedSecondaryVertexBJetTags.clone()
AK8jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("AK8jetImpactParameterTagInfosCHS"), cms.InputTag("AK8jetSecondaryVertexTagInfosCHS") )

AK8jetImpactParameterTagInfosSJCHS                     = impactParameterTagInfos.clone()
AK8jetImpactParameterTagInfosSJCHS.jetTracks           = "AK8jetTracksAssociatorAtVertexSJCHS"
AK8jetSecondaryVertexTagInfosSJCHS                     = secondaryVertexTagInfos.clone()
AK8jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos     = "AK8jetImpactParameterTagInfosSJCHS"
AK8jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
AK8jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("AK8jetImpactParameterTagInfosSJCHS"), cms.InputTag("AK8jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK8QGTaggerCHS                                       = QGTagger.clone()
AK8QGTaggerCHS.srcJets                               = cms.InputTag('AK8PFJetsCHS')
AK8QGTaggerSubJetsCHS                                = AK8QGTaggerCHS.clone()
AK8QGTaggerSubJetsCHS.srcJets                        = cms.InputTag('AK8caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK8NjettinessCHS                                     = Njettiness.clone()       
AK8NjettinessCHS.src                                 =  cms.InputTag('AK8PFJetsCHS')

AK8jetsequenceCHS = cms.Sequence(
    AK8PFJetsCHS                      *
    AK8caPFJetsPrunedCHS              *
    AK8jetTracksAssociatorAtVertexCHS    *
    AK8jetImpactParameterTagInfosCHS     *
    AK8jetSecondaryVertexTagInfosCHS     *
    AK8jetTracksAssociatorAtVertexSJCHS  *
    AK8jetImpactParameterTagInfosSJCHS   *
    AK8jetSecondaryVertexTagInfosSJCHS   *
    AK8jetCombinedSecondaryVertexBJetTagsCHS * 
    AK8jetCombinedSecondaryVertexBJetTagsSJCHS  *
    AK8QGTaggerCHS                       *
    AK8QGTaggerSubJetsCHS                *                
    AK8NjettinessCHS                     *
    AK8jetFlavorCHS                   
    )
