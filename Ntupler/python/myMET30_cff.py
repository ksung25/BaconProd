import FWCore.ParameterSet.Config as cms

def setMet30(process,etaCut,isMiniAOD=True):
    if isMiniAOD:
        process.packedPFCandidates30 = cms.EDFilter("CandPtrSelector",
                                                    src = cms.InputTag("packedPFCandidates"),
                                                    cut = cms.string("abs(eta) < %f"%etaCut))
    else : 
        process.packedPFCandidates30 = cms.EDFilter("CandPtrSelector",
                                                    src = cms.InputTag("particleFlow"),
                                                    cut = cms.string("abs(eta) < %f"%etaCut))
    process.puppi30              = cms.EDFilter("CandPtrSelector",
                                                src = cms.InputTag("puppi"),
                                                cut = cms.string("abs(eta) < %f"%etaCut))
    process.puppinolep30         = cms.EDFilter("CandPtrSelector",
                                                src = cms.InputTag("puppinolep"),
                                                cut = cms.string("abs(eta) < %f"%etaCut))
    process.AK4PFJetsPuppi30             = process.ak4PFJets.clone(src = cms.InputTag("puppi30"))
    process.pfMet30                      = process.pfMetPuppi.clone();
    process.pfMet30.src                  = cms.InputTag('packedPFCandidates30')
    process.pfMetPuppi30                 = process.pfMetPuppi.clone();
    process.pfMetPuppi30.src             = cms.InputTag('puppinolep30')
    process.pfJetMETcorrPuppi30          = process.pfJetMETcorrPuppi.clone(src = 'AK4PFJetsPuppi30')
    process.pfType1PuppiCorrectedMet30   = process.pfType1PuppiCorrectedMet.clone(src='pfMetPuppi30')
    process.pfType1PuppiCorrectedMet30.srcType1Corrections = cms.VInputTag(cms.InputTag('pfJetMETcorrPuppi30', 'type1'))
    process.ak4PFJets30                  = process.ak4PFJets.clone(src = cms.InputTag("packedPFCandidates30"))
    process.pfJetMETcorr30               = process.pfJetMETcorr.clone(src = 'ak4PFJets30')
    process.pfType1CorrectedMet30        = process.pfType1CorrectedMet.clone(src='pfMet30') 
    process.pfType1CorrectedMet30.srcType1Corrections = cms.VInputTag(cms.InputTag('pfJetMETcorr30', 'type1'))
    process.allMET30           = cms.Sequence(process.pfMet30*process.pfMetPuppi30*process.pfJetMETcorr30*process.pfJetMETcorrPuppi30*process.pfType1CorrectedMet30*process.pfType1PuppiCorrectedMet30 )

    #MVA Met 3.0
    process.calibratedAK4PFJetsForPFMVAMEt30 = process.calibratedAK4PFJetsForPFMVAMEt.clone( src = cms.InputTag('ak4PFJets30'))
    process.puJetIdForPFMVAMEt30            = process.puJetIdForPFMVAMEt.clone(jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt30"))     
    process.pfMVAMEt30                      = process.pfMVAMEt.clone()
    process.pfMVAMEt30.srcMVAPileupJetId    = cms.InputTag('puJetIdForPFMVAMEt30','fullDiscriminant')
    process.pfMVAMEt30.srcPFCandidates      = cms.InputTag("packedPFCandidates30")
    process.pfMVAMEt30.srcUncorrJets        = "ak4PFJets30"
    process.pfMVAMEt30.srcCorrJets          = "calibratedAK4PFJetsForPFMVAMEt30" 
    process.pfMVAMEtSequenceNoLep30         = cms.Sequence(process.calibratedAK4PFJetsForPFMVAMEt30*process.puJetIdForPFMVAMEt30*process.pfMVAMEt30)
