from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.SoftLeptonByMVAComputers_cff import *
from RecoBTag.TensorFlow.pfDeepDoubleX_cff import *
#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

def addBTaggingAK4CHS(process,jets='ak4PFJetsCHS',cone=0.4,head='AK4',tail='CHS',useMiniAOD=True,dropSub=False):
    setattr(process, 'pfImpactParameterTagInfos', 
            pfImpactParameterTagInfos.clone(
            jets      = cms.InputTag(jets),
            maxDeltaR = cms.double(cone)
            ))
    setattr(process, 'pfInclusiveSecondaryVertexFinderTagInfos',
            pfInclusiveSecondaryVertexFinderTagInfos.clone(
            trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfos")
            ))
    setattr(process, 'pfSecondaryVertexTagInfos',
            pfSecondaryVertexTagInfos.clone(
            trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfos")
            ))
    setattr(process, 'pfInclusiveSecondaryVertexFinderCvsLTagInfos', 
            pfInclusiveSecondaryVertexFinderCvsLTagInfos.clone(
            trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfos")
            ))
    #Soft Lepton Info
    setattr(process, 'softPFMuonsTagInfos',     softPFMuonsTagInfos    .clone(jets=cms.InputTag(jets)))
    setattr(process, 'softPFElectronsTagInfos', softPFElectronsTagInfos.clone(jets=cms.InputTag(jets)))
    setattr(process, 'softPFMuonBJetTags',      softPFMuonBJetTags     .clone(tagInfos=cms.VInputTag('softPFMuonsTagInfos')))
    setattr(process, 'softPFElectronBJetTags',  softPFElectronBJetTags .clone(tagInfos=cms.VInputTag('softPFElectronsTagInfos')))
    #B-taggers
    setattr(process, 'pfDeepCSVTagInfos',
            pfDeepCSVTagInfos.clone(
            svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfos')
            )),
    setattr(process, 'pfDeepCSVJetTags',
            pfDeepCSVJetTags.clone(
            src = cms.InputTag('pfDeepCSVTagInfos'),
            ))
    setattr(process, 'pfDeepCMVATagInfos',
            pfDeepCMVATagInfos.clone(
            deepNNTagInfos = cms.InputTag('pfDeepCSVTagInfos'),
            ipInfoSrc = cms.InputTag("pfImpactParameterTagInfos"),
            muInfoSrc = cms.InputTag("softPFMuonsTagInfos"),
            elInfoSrc = cms.InputTag("softPFElectronsTagInfos"),
            ))
    setattr(process, 'pfDeepCMVAJetTags',  
            pfDeepCMVAJetTags.clone(
            src = cms.InputTag('pfDeepCMVATagInfos')
            ))
    process.btagging *= getattr(process,'pfImpactParameterTagInfos')
    process.btagging *= getattr(process,'pfSecondaryVertexTagInfos')
    process.btagging *= getattr(process,'pfInclusiveSecondaryVertexFinderTagInfos')
    process.btagging *= getattr(process,'pfInclusiveSecondaryVertexFinderCvsLTagInfos')
    process.btagging *= getattr(process,'softPFMuonsTagInfos')
    process.btagging *= getattr(process,'softPFElectronsTagInfos')
    process.btagging *= getattr(process,'softPFMuonBJetTags')
    process.btagging *= getattr(process,'softPFElectronBJetTags')
    process.btagging *= getattr(process,'pfDeepCSVTagInfos')
    process.btagging *= getattr(process,'pfDeepCSVJetTags')
    process.btagging *= getattr(process,'pfDeepCMVATagInfos')
    process.btagging *= getattr(process,'pfDeepCMVAJetTags')
    if useMiniAOD:
        getattr(process,'pfImpactParameterTagInfos').primaryVertex                   = cms.InputTag("offlineSlimmedPrimaryVertices")
        getattr(process,'pfImpactParameterTagInfos').candidates                      = cms.InputTag("packedPFCandidates")
        getattr(process,'pfSecondaryVertexTagInfos').extSVCollection                 = cms.InputTag('slimmedSecondaryVertices')
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfos').extSVCollection  = cms.InputTag('slimmedSecondaryVertices')
        getattr(process,'pfInclusiveSecondaryVertexFinderCvsLTagInfos').extSVCollection = cms.InputTag('slimmedSecondaryVertices')
        getattr(process,'softPFMuonsTagInfos').muons                         = cms.InputTag("slimmedMuons")
        getattr(process,'softPFMuonsTagInfos').primaryVertex                 = cms.InputTag("offlineSlimmedPrimaryVertices")
        getattr(process,'softPFElectronsTagInfos').electrons                 = cms.InputTag("slimmedElectrons")
        getattr(process,'softPFElectronsTagInfos').primaryVertex             = cms.InputTag("offlineSlimmedPrimaryVertices")
 #   updateJetCollection(
 #       process,
 #       labelName = 'UpdatedJEC',
 #       jetSource = cms.InputTag(jets),
 #       jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
 #       btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probc','pfDeepCSVJetTags:probudsg','pfDeepCSVJetTags:probbb']
 #       )


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
    setattr(process, head+'PFDeepCSVTagInfos'+tail,
            pfDeepCSVTagInfos.clone(
            svTagInfos = cms.InputTag(head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail)
            )),
    setattr(process, head+'PFDeepCSVJetTags'+tail,
            pfDeepCSVJetTags.clone(
            src = cms.InputTag(head+'PFDeepCSVTagInfos'+tail),
            ))
    setattr(process, head+'PFDeepCMVATagInfos'+tail,
            pfDeepCMVATagInfos.clone(
            deepNNTagInfos = cms.InputTag(head+'PFDeepCSVTagInfos'+tail),
            ipInfoSrc = cms.InputTag(head+"PFImpactParameterTagInfos"+tail),
            muInfoSrc = cms.InputTag(head+"PFSoftPFMuonsTagInfos"+tail),
            elInfoSrc = cms.InputTag(head+"PFSoftPFElectronsTagInfos"+tail),
            ))
    setattr(process, head+'PFDeepCMVAJetTags'+tail,  
            pfDeepCMVAJetTags.clone(
            src = cms.InputTag(head+'PFDeepCMVATagInfos'+tail)
            ))
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
        if cone < 1.0:
 	    setattr(process, head+'PFBoostedDoubleSVTagInfos'+tail,
                    pfBoostedDoubleSVAK8TagInfos.clone(
                    svTagInfos = cms.InputTag(head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail)
                    )) 
            setattr(process, head+'PFBoostedDoubleSecondaryVertexBJetTags'+tail,
                    pfBoostedDoubleSecondaryVertexAK8BJetTags.clone(
                    tagInfos = cms.VInputTag(cms.InputTag(head+"PFBoostedDoubleSVTagInfos"+tail))
                    #tagInfos = cms.VInputTag(
#cms.InputTag(head+"PFImpactParameterTagInfos"+tail), 
                    ))    
        else:
  	    setattr(process, head+'PFBoostedDoubleSVTagInfos'+tail,
                   pfBoostedDoubleSVCA15TagInfos.clone(
                   svTagInfos = cms.InputTag(head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail)
            	    ))
            setattr(process, head+'PFBoostedDoubleSecondaryVertexBJetTags'+tail,
                    pfBoostedDoubleSecondaryVertexCA15BJetTags.clone(
                    tagInfos = cms.VInputTag(
#cms.InputTag(head+"PFImpactParameterTagInfos"+tail), #pfImpactParameterAK8TagInfos
                                             cms.InputTag(head+"PFBoostedDoubleSVTagInfos"+tail) )
                    ))    
	#Deep Double B
	setattr(process, head+'PFBoostedDeepDoubleBTagInfos'+tail,
            pfDeepDoubleXTagInfos.clone(
	    shallow_tag_infos = cms.InputTag(head+'PFBoostedDoubleSVTagInfos'+tail),
            jets = cms.InputTag(jets)
            ))
        setattr(process, head+'PFBoostedDeepDoubleBJetTags'+tail, 
            pfDeepDoubleBJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDB_mass_independent.pb')
            ))
        setattr(process, head+'PFBoostedDeepDoubleBvLJetTags'+tail,
            pfDeepDoubleBvLJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDB_mass_independent.pb')
            ))
	#Deep Double B - No Mass Sculpting Penalty
        setattr(process, head+'PFBoostedDeepDoubleBvLNoMassSculptPenJetTags'+tail,
            pfDeepDoubleBvLJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDB.pb')
            #graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleB/V01/constant_graph_PtCut.pb')
            ))
        #Deep Double C
        setattr(process, head+'PFBoostedDeepDoubleCvLJetTags'+tail,
            pfDeepDoubleCvLJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDC_mass_independent.pb')
            ))
	#Deep Double C - No Mass Sculpting Penalty
        setattr(process, head+'PFBoostedDeepDoubleCvLNoMassSculptPenJetTags'+tail,
            pfDeepDoubleCvLJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDC.pb')
            ))
        #Deep Double CvB
        setattr(process, head+'PFBoostedDeepDoubleCvBJetTags'+tail,
            pfDeepDoubleCvBJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDCvB_mass_independent.pb')
            ))
	#Deep Double CvB - No Mass Sculpting Penalty
        setattr(process, head+'PFBoostedDeepDoubleCvBNoMassSculptPenJetTags'+tail,
            pfDeepDoubleCvBJetTags.clone(
            src = cms.InputTag(head+"PFBoostedDeepDoubleBTagInfos"+tail),
            graph_path = cms.FileInPath('RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDCvB.pb')
            ))

    process.btagging *= getattr(process,head+'PFImpactParameterTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSecondaryVertexTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFInclusiveSecondaryVertexFinderTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFInclusiveSecondaryVertexFinderCvsLTagInfos'+tail)
    if not dropSub:
        process.btagging *= getattr(process,head+'PFBoostedDoubleSVTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFMuonsTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFElectronsTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFMuonBJetTags'+tail)
    process.btagging *= getattr(process,head+'PFSoftPFElectronBJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedInclusiveSecondaryVertexV2BJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedMVAV2BJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedCvsLJetTags'+tail)
    process.btagging *= getattr(process,head+'PFCombinedCvsBJetTags'+tail)
    process.btagging *= getattr(process,head+'PFDeepCSVTagInfos'+tail)
    process.btagging *= getattr(process,head+'PFDeepCSVJetTags'+tail)
    process.btagging *= getattr(process,head+'PFDeepCMVATagInfos'+tail)
    process.btagging *= getattr(process,head+'PFDeepCMVAJetTags'+tail)
    if not dropSub:
        process.btagging *= getattr(process,head+'PFImpactParameterTagInfosSJ'+tail)
        process.btagging *= getattr(process,head+'PFInclusiveSecondaryVertexFinderTagInfosSJ'+tail)
        process.btagging *= getattr(process,head+'PFCombinedInclusiveSecondaryVertexV2BJetTagsSJ'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDoubleSecondaryVertexBJetTags'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleBTagInfos'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleBvLJetTags'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleBvLNoMassSculptPenJetTags'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleCvLJetTags'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleCvLNoMassSculptPenJetTags'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleCvBJetTags'+tail)
        process.btagging *= getattr(process,head+'PFBoostedDeepDoubleCvBNoMassSculptPenJetTags'+tail)
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
            getattr(process,head+'PFBoostedDeepDoubleBTagInfos'+tail).vertices                     = cms.InputTag("offlineSlimmedPrimaryVertices")
            getattr(process,head+'PFBoostedDeepDoubleBTagInfos'+tail).secondary_vertices           = cms.InputTag("slimmedSecondaryVertices")

