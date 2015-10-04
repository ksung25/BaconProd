import FWCore.ParameterSet.Config as cms
from RecoMET.METPUSubtraction.mvaPFMET_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3
from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3Residual

pfMVAMEt.srcLeptons = cms.VInputTag("isomuons","gedGsfElectrons","isotaus","gedPhotons")

MVAMetSeq = cms.Sequence(ak4PFJetsL1FastL2L3 * pfMVAMEtSequence)
MVAMetSeqData = cms.Sequence(ak4PFJetsL1FastL2L3Residual * pfMVAMEtSequence)



from RecoMET.METPUSubtraction.objectSelection_miniAOD_cff import addLeptons

def setMiniAODMVAMet(process):
    addLeptons(process)
    process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
    process.ak4PFJets.doAreaFastjet = cms.bool(True)
    process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
    process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
    process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.pfMVAMEt.srcLeptons  = cms.VInputTag("slimmedMuonsTight","slimmedElectrons","slimmedTausLoose","slimmedPhotons")
    process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
    process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")

