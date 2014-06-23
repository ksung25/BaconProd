import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK4byRefCHS = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK4PFJetsCHS"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
AK4byValPhysCHS = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK4byRefCHS"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK4byValAlgoCHS = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK4byRefCHS"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK4jetFlavorCHS    = cms.Sequence(AK4byRefCHS*AK4byValPhysCHS*AK4byValAlgoCHS)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK4PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('pfNoElectron'),
    rParam   = cms.double(0.4),
    jetPtMin = cms.double(20)
    )

AK4caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src      = cms.InputTag('pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.4),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK4jetTracksAssociatorAtVertexCHS          = ic5JetTracksAssociatorAtVertex.clone()
AK4jetTracksAssociatorAtVertexCHS  .jets   = cms.InputTag('AK4PFJetsCHS')
AK4jetTracksAssociatorAtVertexCHS  .tracks = "generalTracks"

AK4jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
AK4jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('AK4caPFJetsPrunedCHS','SubJets')
AK4jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK4jetImpactParameterTagInfosCHS                      = impactParameterTagInfos.clone()
AK4jetImpactParameterTagInfosCHS.jetTracks            = "AK4jetTracksAssociatorAtVertexCHS"
AK4jetSecondaryVertexTagInfosCHS                      = secondaryVertexTagInfos.clone()
AK4jetSecondaryVertexTagInfosCHS.trackIPTagInfos      = "AK4jetImpactParameterTagInfosCHS"
AK4jetCombinedSecondaryVertexBJetTagsCHS           = combinedSecondaryVertexBJetTags.clone()
AK4jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("AK4jetImpactParameterTagInfosCHS"), cms.InputTag("AK4jetSecondaryVertexTagInfosCHS") )

AK4jetImpactParameterTagInfosSJCHS                     = impactParameterTagInfos.clone()
AK4jetImpactParameterTagInfosSJCHS.jetTracks           = "AK4jetTracksAssociatorAtVertexSJCHS"
AK4jetSecondaryVertexTagInfosSJCHS                     = secondaryVertexTagInfos.clone()
AK4jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos     = "AK4jetImpactParameterTagInfosSJCHS"
AK4jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
AK4jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("AK4jetImpactParameterTagInfosSJCHS"), cms.InputTag("AK4jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK4QGTaggerCHS                                       = QGTagger.clone()
AK4QGTaggerCHS.srcJets                               = cms.InputTag('AK4PFJetsCHS')
AK4QGTaggerSubJetsCHS                                = AK4QGTaggerCHS.clone()
AK4QGTaggerSubJetsCHS.srcJets                        = cms.InputTag('AK4caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK4NjettinessCHS                                     = Njettiness.clone()       
AK4NjettinessCHS.src                                 =  cms.InputTag('AK4PFJetsCHS')

AK4jetsequenceCHS = cms.Sequence(
    AK4PFJetsCHS                      *
    AK4caPFJetsPrunedCHS              *
    AK4jetTracksAssociatorAtVertexCHS    *
    AK4jetImpactParameterTagInfosCHS     *
    AK4jetSecondaryVertexTagInfosCHS     *
    AK4jetTracksAssociatorAtVertexSJCHS  *
    AK4jetImpactParameterTagInfosSJCHS   *
    AK4jetSecondaryVertexTagInfosSJCHS   *
    AK4jetCombinedSecondaryVertexBJetTagsCHS * 
    AK4jetCombinedSecondaryVertexBJetTagsSJCHS  *
    AK4QGTaggerCHS                       *
    AK4QGTaggerSubJetsCHS                *                
    AK4NjettinessCHS                     *
    AK4jetFlavorCHS                   
    )
