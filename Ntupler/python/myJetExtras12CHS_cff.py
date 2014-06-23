import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
from RecoJets.JetProducers.ak5PFJetsPruned_cfi  import ak5PFJetsPruned

# Flavour byValue PhysDef
AK12byRefCHS = cms.EDProducer("JetPartonMatcher",
                          jets = cms.InputTag("AK12PFJetsCHS"),
                          coneSizeToAssociate = cms.double(0.3),
                          partons = cms.InputTag("partons")
                          )
AK12byValPhysCHS = cms.EDProducer("JetFlavourIdentifier",
                              srcByReference = cms.InputTag("AK12byRefCHS"),
                              physicsDefinition = cms.bool(True),
                              leptonInfo = cms.bool(True)
                              )                              

AK12byValAlgoCHS = cms.EDProducer("JetFlavourIdentifier",
                                 srcByReference = cms.InputTag("AK12byRefCHS"),
                                 physicsDefinition = cms.bool(False),
                                 leptonInfo = cms.bool(True))

AK12jetFlavorCHS    = cms.Sequence(AK12byRefCHS*AK12byValPhysCHS*AK12byValAlgoCHS)

#for each jet collection run Pruning, subjet b-tagging, quark gluon discrimination,n-subjettiness and subjet quark gluon discrimination
AK12PFJetsCHS = ak5PFJets.clone(
    src      = cms.InputTag('pfNoElectron'),
    rParam   = cms.double(1.2),
    jetPtMin = cms.double(100)
    )

AK12caPFJetsPrunedCHS = ak5PFJetsPruned.clone(
    src      = cms.InputTag('pfNoElectron'),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(1.2),
    doAreaFastjet = cms.bool(False),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(100)
    )

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import ic5JetTracksAssociatorAtVertex
AK12jetTracksAssociatorAtVertexCHS          = ic5JetTracksAssociatorAtVertex.clone()
AK12jetTracksAssociatorAtVertexCHS  .jets   = cms.InputTag('AK12PFJetsCHS')
AK12jetTracksAssociatorAtVertexCHS  .tracks = "generalTracks"

AK12jetTracksAssociatorAtVertexSJCHS        = ic5JetTracksAssociatorAtVertex.clone()
AK12jetTracksAssociatorAtVertexSJCHS.jets   = cms.InputTag('AK12caPFJetsPrunedCHS','SubJets')
AK12jetTracksAssociatorAtVertexSJCHS.tracks = "generalTracks"

from RecoBTag.Configuration.RecoBTag_cff import *
AK12jetImpactParameterTagInfosCHS                      = impactParameterTagInfos.clone()
AK12jetImpactParameterTagInfosCHS.jetTracks            = "AK12jetTracksAssociatorAtVertexCHS"
AK12jetSecondaryVertexTagInfosCHS                      = secondaryVertexTagInfos.clone()
AK12jetSecondaryVertexTagInfosCHS.trackIPTagInfos      = "AK12jetImpactParameterTagInfosCHS"
AK12jetCombinedSecondaryVertexBJetTagsCHS           = combinedSecondaryVertexBJetTags.clone()
AK12jetCombinedSecondaryVertexBJetTagsCHS.tagInfos = cms.VInputTag( cms.InputTag("AK12jetImpactParameterTagInfosCHS"), cms.InputTag("AK12jetSecondaryVertexTagInfosCHS") )

AK12jetImpactParameterTagInfosSJCHS                     = impactParameterTagInfos.clone()
AK12jetImpactParameterTagInfosSJCHS.jetTracks           = "AK12jetTracksAssociatorAtVertexSJCHS"
AK12jetSecondaryVertexTagInfosSJCHS                     = secondaryVertexTagInfos.clone()
AK12jetSecondaryVertexTagInfosSJCHS.trackIPTagInfos     = "AK12jetImpactParameterTagInfosSJCHS"
AK12jetCombinedSecondaryVertexBJetTagsSJCHS          = combinedSecondaryVertexBJetTags.clone()
AK12jetCombinedSecondaryVertexBJetTagsSJCHS.tagInfos = cms.VInputTag( cms.InputTag("AK12jetImpactParameterTagInfosSJCHS"), cms.InputTag("AK12jetSecondaryVertexTagInfosSJCHS") )

from JetTools.AnalyzerToolbox.QGTagger_RecoJets_cff import *
AK12QGTaggerCHS                                       = QGTagger.clone()
AK12QGTaggerCHS.srcJets                               = cms.InputTag('AK12PFJetsCHS')
AK12QGTaggerSubJetsCHS                                = AK12QGTaggerCHS.clone()
AK12QGTaggerSubJetsCHS.srcJets                        = cms.InputTag('AK12caPFJetsPrunedCHS','SubJets')

from JetTools.AnalyzerToolbox.njettinessadder_cfi import *
AK12NjettinessCHS                                     = Njettiness.clone()       
AK12NjettinessCHS.src                                 =  cms.InputTag('AK12PFJetsCHS')

AK12jetsequenceCHS = cms.Sequence(
    AK12PFJetsCHS                      *
    AK12caPFJetsPrunedCHS              *
    AK12jetTracksAssociatorAtVertexCHS    *
    AK12jetImpactParameterTagInfosCHS     *
    AK12jetSecondaryVertexTagInfosCHS     *
    AK12jetTracksAssociatorAtVertexSJCHS  *
    AK12jetImpactParameterTagInfosSJCHS   *
    AK12jetSecondaryVertexTagInfosSJCHS   *
    AK12jetCombinedSecondaryVertexBJetTagsCHS * 
    AK12jetCombinedSecondaryVertexBJetTagsSJCHS  *
    AK12QGTaggerCHS                       *
    AK12QGTaggerSubJetsCHS                *                
    AK12NjettinessCHS                     *
    AK12jetFlavorCHS                   
    )
