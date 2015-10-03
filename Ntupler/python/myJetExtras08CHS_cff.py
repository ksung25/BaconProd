import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoGenJets_cff     import ak5GenJets
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
CA8byRefCHS = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("CA8PFJetsCHS"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
CA8byValPhysCHS = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("CA8byRefCHS"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

CA8byValAlgoCHS = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("CA8byRefCHS"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

CA8jetFlavorCHS    = cms.Sequence(CA8byRefCHS*CA8byValPhysCHS*CA8byValAlgoCHS)

CA8GenJetsCHS = ak5GenJets.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8)
    )


#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
CA8PFJetsCHS = ak5PFJets.clone(
    src          = cms.InputTag('pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(0.8),
    jetPtMin     = cms.double(20)
    )

CA8caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src                 = cms.InputTag('pfNoElectron'),
    jetAlgorithm        = cms.string("CambridgeAachen"),
    rParam              = cms.double(0.8),
    doAreaFastjet       = cms.bool(False),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets"),
    jetPtMin            = cms.double(20)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
CA8jetTracksAssociatorAtVertexCHS        = ic5JetTracksAssociatorAtVertex.clone()
CA8jetTracksAssociatorAtVertexCHS.jets   = cms.InputTag('CA8PFJetsCHS')
CA8jetTracksAssociatorAtVertexCHS.tracks = "generalTracks"

CA8jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
CA8jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('CA8caPFJetsPrunedCHS','SubJets')
CA8jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
CA8jetImpactParameterTagInfosCHS                  = impactParameterTagInfos.clone()
CA8jetImpactParameterTagInfosCHS.jetTracks        = "CA8jetTracksAssociatorAtVertexCHS"
CA8jetSecondaryVertexTagInfosCHS                  = secondaryVertexTagInfos.clone()
CA8jetSecondaryVertexTagInfosCHS.trackIPTagInfos  = "CA8jetImpactParameterTagInfosCHS"
CA8jetCombinedSecondaryVertexBJetTagsCHS          = combinedSecondaryVertexBJetTags.clone()
CA8jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("CA8jetImpactParameterTagInfosCHS"), cms.InputTag("CA8jetSecondaryVertexTagInfosCHS") )

CA8jetImpactParameterTagInfosSJCHS                  = impactParameterTagInfos.clone()
CA8jetImpactParameterTagInfosSJCHS.jetTracks        = "CA8jetTracksAssociatorAtVertexSJCHS"
CA8jetSecondaryVertexTagInfosSJCHS                  = secondaryVertexTagInfos.clone()
CA8jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos  = "CA8jetImpactParameterTagInfosSJCHS"
CA8jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
CA8jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("CA8jetImpactParameterTagInfosSJCHS"), cms.InputTag("CA8jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
CA8QGTaggerCHS                = QGTagger.clone()
CA8QGTaggerCHS.srcJets        = cms.InputTag('CA8PFJetsCHS')
CA8QGTaggerSubJetsCHS         = CA8QGTaggerCHS.clone()
CA8QGTaggerSubJetsCHS.srcJets = cms.InputTag('CA8caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
CA8NjettinessCHS     = Njettiness.clone()       
CA8NjettinessCHS.src = cms.InputTag('CA8PFJetsCHS')


CA8genjetsequenceCHS = cms.Sequence(
    CA8GenJetsCHS    * 
    CA8jetFlavorCHS                   
)

CA8jetsequenceCHS = cms.Sequence(
    CA8PFJetsCHS                      *
    CA8caPFJetsPrunedCHS              *
    CA8jetTracksAssociatorAtVertexCHS    *
    CA8jetImpactParameterTagInfosCHS     *
    CA8jetSecondaryVertexTagInfosCHS     *
    CA8jetTracksAssociatorAtVertexSJCHS  *
    CA8jetImpactParameterTagInfosSJCHS   *
    CA8jetSecondaryVertexTagInfosSJCHS   *
    CA8jetCombinedSecondaryVertexBJetTagsCHS * 
    CA8jetCombinedSecondaryVertexBJetTagsSJCHS  *
    CA8QGTaggerCHS                       *
    CA8QGTaggerSubJetsCHS                *                
    CA8NjettinessCHS                   
    )
