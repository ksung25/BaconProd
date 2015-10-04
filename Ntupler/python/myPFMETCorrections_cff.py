import FWCore.ParameterSet.Config as cms

# Set pfJetMETcorr.jetCorrLabel accordingly in the top level config file
# Use producePFMETCorrections(MC/Data) accordingly in the top level config file

# load jet energy correctors
from JetMETCorrections.Configuration.JetCorrectors_cff  import *
from BaconProd.Ntupler.myPUPPICorrections_cff           import *

#--------------------------------------------------------------------------------
# produce Type 1 MET corrections for PFJets
pfJetMETcorr = cms.EDProducer("PFJetMETcorrInputProducer",
    src = cms.InputTag('ak4PFJets'),
    offsetCorrLabel = cms.InputTag("ak4PFL1FastjetCorrector"),
    jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3Corrector"), # NOTE: use "ak4PFL1FastL2L3Corrector" for MC / "ak4PFL1FastL2L3ResidualCorrector" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)					  
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# use MET corrections to produce Type 1 corrected PFMET objects
pfType1CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(False),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
	cms.InputTag('pfJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define sequence to run all modules
producePFMETCorrectionsMC = cms.Sequence(
  ak4PFL1FastL2L3CorrectorChain +
  pfJetMETcorr +
  pfType1CorrectedMet
)

producePFMETCorrectionsData = cms.Sequence(
  ak4PFL1FastL2L3ResidualCorrectorChain +
  pfJetMETcorr +
  pfType1CorrectedMet
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#Puppi
pfJetMETcorrPuppi = pfJetMETcorr.clone(
    src = 'AK4PFJetsPuppi',
    jetCorrLabel = 'ak4PuppiL1FastL2L3Corrector',
    offsetCorrLabel = 'ak4PuppiL1FastjetCorrector',
    type1JetPtThreshold = cms.double(20)
    )

pfType1PuppiCorrectedMet = pfType1CorrectedMet.clone(
    src = cms.InputTag('pfMetPuppi'),
    srcType1Corrections = [ cms.InputTag("pfJetMETcorrPuppi","type1") ]
    )

#--------------------------------------------------------------------------------
# define sequence to run all modules
producePFMETCorrectionsPuppiMC = cms.Sequence(
  ak4PuppiL1FastL2L3Chain +
  pfJetMETcorrPuppi  +
  pfType1PuppiCorrectedMet
)
producePFMETCorrectionsPuppiData = cms.Sequence(
  ak4PuppiL1FastL2L3ResidualChain +
  pfJetMETcorrPuppi + 
  pfType1PuppiCorrectedMet
)
