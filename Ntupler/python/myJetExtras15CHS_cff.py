import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoGenJets_cff     import ak5GenJets
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
CA15byRefCHS = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("CA15PFJetsCHS"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
CA15byValPhysCHS = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("CA15byRefCHS"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

CA15byValAlgoCHS = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("CA15byRefCHS"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

CA15jetFlavorCHS    = cms.Sequence(CA15byRefCHS*CA15byValPhysCHS*CA15byValAlgoCHS)

CA15GenJetsCHS = ak5GenJets.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam       = cms.double(1.5)
    )

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
CA15PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('pfNoElectron'),
    rParam   = cms.double(1.5),
    jetPtMin = cms.double(100)
    )

CA15caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src                 = cms.InputTag('pfNoElectron'),
    jetAlgorithm        = cms.string("CambridgeAachen"),
    rParam              = cms.double(1.5),
    doAreaFastjet       = cms.bool(False),
    writeCompound       = cms.bool(True),
    jetCollInstanceName = cms.string("SubJets"),
    jetPtMin            = cms.double(100)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
CA15jetTracksAssociatorAtVertexCHS        = ic5JetTracksAssociatorAtVertex.clone()
CA15jetTracksAssociatorAtVertexCHS.jets   = cms.InputTag('CA15PFJetsCHS')
CA15jetTracksAssociatorAtVertexCHS.tracks = "generalTracks"

CA15jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
CA15jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('CA15caPFJetsPrunedCHS','SubJets')
CA15jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
CA15jetImpactParameterTagInfosCHS                  = impactParameterTagInfos.clone()
CA15jetImpactParameterTagInfosCHS.jetTracks        = "CA15jetTracksAssociatorAtVertexCHS"
CA15jetSecondaryVertexTagInfosCHS                  = secondaryVertexTagInfos.clone()
CA15jetSecondaryVertexTagInfosCHS.trackIPTagInfos  = "CA15jetImpactParameterTagInfosCHS"
CA15jetCombinedSecondaryVertexBJetTagsCHS          = combinedSecondaryVertexBJetTags.clone()
CA15jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("CA15jetImpactParameterTagInfosCHS"), cms.InputTag("CA15jetSecondaryVertexTagInfosCHS") )

CA15jetImpactParameterTagInfosSJCHS                  = impactParameterTagInfos.clone()
CA15jetImpactParameterTagInfosSJCHS.jetTracks        = "CA15jetTracksAssociatorAtVertexSJCHS"
CA15jetSecondaryVertexTagInfosSJCHS                  = secondaryVertexTagInfos.clone()
CA15jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos  = "CA15jetImpactParameterTagInfosSJCHS"
CA15jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
CA15jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("CA15jetImpactParameterTagInfosSJCHS"), cms.InputTag("CA15jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
CA15QGTaggerCHS                = QGTagger.clone()
CA15QGTaggerCHS.srcJets        = cms.InputTag('CA15PFJetsCHS')
CA15QGTaggerSubJetsCHS         = CA15QGTaggerCHS.clone()
CA15QGTaggerSubJetsCHS.srcJets = cms.InputTag('CA15caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
CA15NjettinessCHS     = Njettiness.clone()       
CA15NjettinessCHS.src = cms.InputTag('CA15PFJetsCHS')

CA15genjetsequenceCHS = cms.Sequence(
    CA15GenJetsCHS    * 
    CA15jetFlavorCHS                   
)

CA15jetsequenceCHS = cms.Sequence(
    CA15PFJetsCHS                      *
    CA15caPFJetsPrunedCHS              *
    CA15jetTracksAssociatorAtVertexCHS    *
    CA15jetImpactParameterTagInfosCHS     *
    CA15jetSecondaryVertexTagInfosCHS     *
    CA15jetTracksAssociatorAtVertexSJCHS  *
    CA15jetImpactParameterTagInfosSJCHS   *
    CA15jetSecondaryVertexTagInfosSJCHS   *
    CA15jetCombinedSecondaryVertexBJetTagsCHS * 
    CA15jetCombinedSecondaryVertexBJetTagsSJCHS  *
    CA15QGTaggerCHS                       *
    CA15QGTaggerSubJetsCHS                *                
    CA15NjettinessCHS                   
    )
