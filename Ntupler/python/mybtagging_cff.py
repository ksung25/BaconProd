import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jet and tracks association
myJetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
myJetTracksAssociatorAtVertex.jets   = cms.InputTag("ak5PFJets")
myJetTracksAssociatorAtVertex.tracks = cms.InputTag("generalTracks")


# impact parameter b-tag
myImpactParameterTagInfos = impactParameterTagInfos.clone()
myImpactParameterTagInfos.jetTracks = cms.InputTag("myJetTracksAssociatorAtVertex")
myTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
myTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos") )
myTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
myTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos") )
myJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
myJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos") )
myJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
myJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos") )


# secondary vertex b-tag
mySecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
mySecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag("myImpactParameterTagInfos")

mySimpleSecondaryVertexHighEffBJetTags = simpleSecondaryVertexHighEffBJetTags.clone()
mySimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySecondaryVertexTagInfos") )
mySimpleSecondaryVertexHighPurBJetTags = simpleSecondaryVertexHighPurBJetTags.clone()
mySimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySecondaryVertexTagInfos") )
myCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
myCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos"), cms.InputTag("mySecondaryVertexTagInfos") )
myCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
myCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos"), cms.InputTag("mySecondaryVertexTagInfos") )


# ghost track b-tag
myGhostTrackVertexTagInfos = ghostTrackVertexTagInfos.clone()
myGhostTrackVertexTagInfos.trackIPTagInfos =cms.InputTag( "myImpactParameterTagInfos")
myGhostTrackBJetTags = ghostTrackBJetTags.clone()
myGhostTrackBJetTags.tagInfos = cms.VInputTag( cms.InputTag("myImpactParameterTagInfos"), cms.InputTag("myGhostTrackVertexTagInfos") )


# soft electron b-tag
mySoftElectronTagInfos = softElectronTagInfos.clone()
mySoftElectronTagInfos.jets = cms.InputTag("ak5PFJets")
mySoftElectronByIP3dBJetTags = softElectronByIP3dBJetTags.clone()
mySoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySoftElectronTagInfos") )
mySoftElectronByPtBJetTags = softElectronByPtBJetTags.clone()
mySoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySoftElectronTagInfos") )


# soft muon b-tag
mySoftMuonTagInfos = softMuonTagInfos.clone()
mySoftMuonTagInfos.jets = cms.InputTag("ak5PFJets")
mySoftMuonBJetTags = softMuonBJetTags.clone()
mySoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySoftMuonTagInfos") )
mySoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
mySoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySoftMuonTagInfos") )
mySoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
mySoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("mySoftMuonTagInfos") )


# prepare a path running the new modules
myJetTracksAssociator = cms.Sequence(
  myJetTracksAssociatorAtVertex
)

myJetBtaggingIP = cms.Sequence(
  myImpactParameterTagInfos * (
    myTrackCountingHighEffBJetTags +
    myTrackCountingHighPurBJetTags +
    myJetProbabilityBJetTags +
    myJetBProbabilityBJetTags
  )
)

myJetBtaggingSV = cms.Sequence(
  myImpactParameterTagInfos *
  mySecondaryVertexTagInfos * (
    mySimpleSecondaryVertexHighEffBJetTags +
    mySimpleSecondaryVertexHighPurBJetTags +
    myCombinedSecondaryVertexBJetTags +
    myCombinedSecondaryVertexBJetTags
  )
)

myJetBtaggingGhostTrack = cms.Sequence(
  myGhostTrackVertexTagInfos *
  myGhostTrackBJetTags
)

myJetBtaggingEle = cms.Sequence(
  softElectronCands *
  mySoftElectronTagInfos *
  mySoftElectronByIP3dBJetTags *
  mySoftElectronByPtBJetTags
)

myJetBtaggingMu = cms.Sequence(
  mySoftMuonTagInfos * (
    mySoftMuonBJetTags +
    mySoftMuonByIP3dBJetTags +
    mySoftMuonByPtBJetTags
  )
)

mybtagging = cms.Sequence(
  myJetTracksAssociator * (
    myJetBtaggingIP +
    myJetBtaggingSV +
    myJetBtaggingGhostTrack +
    myJetBtaggingEle +
    myJetBtaggingMu
  )
)
