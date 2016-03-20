from RecoBTag.Configuration.RecoBTag_cff import *

def addBTagging(process,jets='ak4PFJetsCHS',cone=0.4,head='AK4',tail='CHS',useMiniAOD=True,dropSub=False):
    setattr(process, head+'PFImpactParameterTagInfos'+tail, 
            pfImpactParameterTagInfos.clone(
            jets      = cms.InputTag(jets),
            maxDeltaR = cms.double(cone)
            ))
    setattr(process, head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail,
            pfInclusiveSecondaryVertexFinderTagInfos.clone(
            trackIPTagInfos = cms.InputTag(head+"PFImpactParameterTagInfos"+tail)
            ))
    setattr(process, head+'PFSecondaryVertexTagInfos'+tail,
            pfSecondaryVertexTagInfos.clone(
            trackIPTagInfos = cms.InputTag(head+"PFImpactParameterTagInfos"+tail)
            ))
    setattr(process, head+'PFInclusiveSecondaryVertexFinderCvsLTagInfos'+tail, 
            pfInclusiveSecondaryVertexFinderCvsLTagInfos.clone(
            trackIPTagInfos = cms.InputTag(head+"PFImpactParameterTagInfos"+tail)
            ))
    #Soft Lepton Info
    setattr(process, head+'PFSoftPFMuonsTagInfos'+tail,     softPFMuonsTagInfos    .clone(jets=cms.InputTag(jets)))
    setattr(process, head+'PFSoftPFElectronsTagInfos'+tail, softPFElectronsTagInfos.clone(jets=cms.InputTag(jets)))
    setattr(process, head+'PFSoftPFMuonBJetTags'+tail,      softPFMuonBJetTags     .clone(tagInfos=cms.VInputTag(head+'PFSoftPFMuonsTagInfos'+tail)))
    setattr(process, head+'PFSoftPFElectronBJetTags'+tail,  softPFElectronBJetTags .clone(tagInfos=cms.VInputTag(head+'PFSoftPFElectronsTagInfos'+tail)))
    #B-taggers
    setattr(process, head+'PFCombinedInclusiveSecondaryVertexV2BJetTags'+tail,
            pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
            tagInfos = cms.VInputTag( cms.InputTag(head+'PFImpactParameterTagInfos'+tail), 
                                      cms.InputTag(head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail))
            ))
    setattr(process, head+'PFCombinedMVAV2BJetTags'+tail,
            pfCombinedMVAV2BJetTags.clone(
            tagInfos = cms.VInputTag( cms.InputTag(head+"PFImpactParameterTagInfos"+tail),
                                      cms.InputTag(head+"PFInclusiveSecondaryVertexFinderTagInfos"+tail),
                                      cms.InputTag(head+"PFSecondaryVertexTagInfos"+tail),
                                      cms.InputTag(head+"PFSoftPFMuonsTagInfos"+tail),
                                      cms.InputTag(head+"PFSoftPFElectronsTagInfos"+tail)
                                      )))
    #C-taggers
    setattr(process, head+'PFCombinedCvsLJetTags'+tail,
            pfCombinedCvsLJetTags.clone(
            tagInfos = cms.VInputTag( cms.InputTag(head+"PFImpactParameterTagInfos"+tail),
                                      cms.InputTag(head+"PFInclusiveSecondaryVertexFinderCvsLTagInfos"+tail),
                                      cms.InputTag(head+"PFSoftPFMuonsTagInfos"+tail),
                                      cms.InputTag(head+"PFSoftPFElectronsTagInfos"+tail)
                                      )))
    setattr(process, head+'PFCombinedCvsBJetTags'+tail,
            getattr(process,head+'PFCombinedCvsLJetTags'+tail).clone(
            jetTagComputer = cms.string('charmTagsComputerCvsB')))

    # subjet b-tagging
    if not dropSub:
        setattr(process, head+'PFImpactParameterTagInfosSJ'+tail,
                pfImpactParameterTagInfos.clone(
                jets      = cms.InputTag(head+'caPFJetsSoftDrop'+tail,'SubJets'),
                maxDeltaR = cms.double(0.4)
                ))

        setattr(process, head+'PFInclusiveSecondaryVertexFinderTagInfosSJ'+tail,
                pfInclusiveSecondaryVertexFinderTagInfos.clone(
                trackIPTagInfos = cms.InputTag(head+"PFImpactParameterTagInfosSJ"+tail)
                ))
        setattr(process, head+'PFCombinedInclusiveSecondaryVertexV2BJetTagsSJ'+tail,
                pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
                tagInfos = cms.VInputTag( cms.InputTag(head+"PFImpactParameterTagInfosSJ"+tail), 
                                          cms.InputTag(head+"PFInclusiveSecondaryVertexFinderTagInfosSJ"+tail) )
                ))
     #double b-tagging
        setattr(process, head+'PFBoostedDoubleSecondaryVertexBJetTags'+tail,
                pfBoostedDoubleSecondaryVertexAK8BJetTags.clone(
                tagInfos = cms.VInputTag(cms.InputTag(head+"PFImpactParameterTagInfos"+tail), 
                                         cms.InputTag(head+"PFInclusiveSecondaryVertexFinderTagInfos"+tail) )
                ))    

    process.btagging *= getattr(process,head+'PFImpactParameterTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSecondaryVertexTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFInclusiveSecondaryVertexFinderCvsLTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFMuonsTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFElectronsTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFMuonBJetTags'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFElectronBJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedInclusiveSecondaryVertexV2BJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedMVAV2BJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedCvsLJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedCvsBJetTags'+tail)
    if not dropSub:
        process.btagging *= getattr(process,head+'PFImpactParameterTagInfosSJ'+tail)
        process.btagging *= getattr(process,head+'PFInclusiveSecondaryVertexFinderTagInfosSJ'+tail)
        process.btagging *= getattr(process,head+'PFCombinedInclusiveSecondaryVertexV2BJetTagsSJ'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDoubleSecondaryVertexBJetTags'+tail)
    if useMiniAOD:
        getattr(process,head+'PFImpactParameterTagInfos'+tail).primaryVertex                   = cms.InputTag("offlineSlimmedPrimaryVertices")
        getattr(process,head+'PFImpactParameterTagInfos'+tail).candidates                      = cms.InputTag("packedPFCandidates")
        getattr(process,head+'PFSecondaryVertexTagInfos'+tail).extSVCollection                 = cms.InputTag('slimmedSecondaryVertices')
        getattr(process,head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail).extSVCollection  = cms.InputTag('slimmedSecondaryVertices')
        getattr(process,head+'PFInclusiveSecondaryVertexFinderCvsLTagInfos'+tail).extSVCollection = cms.InputTag('slimmedSecondaryVertices')
        getattr(process,head+'PFSoftPFMuonsTagInfos'+tail).muons                         = cms.InputTag("slimmedMuons")
        getattr(process,head+'PFSoftPFMuonsTagInfos'+tail).primaryVertex                 = cms.InputTag("offlineSlimmedPrimaryVertices")
        getattr(process,head+'PFSoftPFElectronsTagInfos'+tail).electrons                 = cms.InputTag("slimmedElectrons")
        getattr(process,head+'PFSoftPFElectronsTagInfos'+tail).primaryVertex             = cms.InputTag("offlineSlimmedPrimaryVertices")
        if not dropSub:
            getattr(process,head+'PFImpactParameterTagInfosSJ'+tail).primaryVertex                 = cms.InputTag("offlineSlimmedPrimaryVertices")
            getattr(process,head+'PFImpactParameterTagInfosSJ'+tail).candidates                    = cms.InputTag("packedPFCandidates")
            getattr(process,head+'PFInclusiveSecondaryVertexFinderTagInfosSJ'+tail).extSVCollection   = cms.InputTag('slimmedSecondaryVertices')

