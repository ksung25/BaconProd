import FWCore.ParameterSet.Config as cms
from RecoMET.METPUSubtraction.mvaPFMET_leptons_cff import *
from RecoTauTag.Configuration.RecoPFTauTag_cff import *

pfCandsNotInJet = cms.EDProducer("TPPFJetsOnPFCandidates",
  bottomCollection = cms.InputTag("particleFlow"),
  enable	   = cms.bool(True),
  topCollection    = cms.InputTag("ak5PFJets"),
  name  	   = cms.untracked.string('noJet'),
  verbose	   = cms.untracked.bool(False)
)

pfCandsUCDownTmp = cms.EDProducer("ShiftedPFCandidateProducer",
  src         = cms.InputTag("pfCandsNotInJet"),
  uncertainty = cms.double(0.1),
  shiftBy     = cms.double(1.0)
)

pfCandsUCDown = cms.EDProducer("ShiftedPFCandidateProducerForPFMEtMVA",
  dRmatch_PFCandidate = cms.double(0.01),
  srcPFCandidates     = cms.InputTag("particleFlow"),
  srcShiftedObjects   = cms.InputTag("pfCandsUCDownTmp"),
  srcUnshiftedObjects = cms.InputTag("particleFlow")
) 

smearedJets = cms.EDProducer("SmearedPFJetProducer",
  jetResolutions = cms.PSet(
    ptresolthreshold = cms.double(10.0),
    
    resolutionsAlgo = cms.string('AK5PF'),
    resolutionsEra  = cms.string('Spring10'),
    
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    HE_EtResPar = cms.vdouble(0.0, 1.3,  0.05),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    HO_EtResPar = cms.vdouble(0.0, 1.3,  0.005),
     
    EE_PhiResPar = cms.vdouble(0.02511),
    EB_PhiResPar = cms.vdouble(0.00502),
    HB_PhiResPar = cms.vdouble(0.02511),
    HE_PhiResPar = cms.vdouble(0.02511),
    HF_PhiResPar = cms.vdouble(0.05022),
    HO_PhiResPar = cms.vdouble(0.02511),

    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType4 = cms.vdouble(0.042, 0.1, 0.0),
    PF_EtResType5 = cms.vdouble(0.41, 0.52, 0.25),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
        
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    PF_PhiResType4 = cms.vdouble(0.0028, 0.0, 0.0022),
    PF_PhiResType5 = cms.vdouble(0.1, 0.1, 0.13),
    PF_PhiResType6 = cms.vdouble(0.02511),
    PF_PhiResType7 = cms.vdouble(0.02511),    
    
    jdpt0 = cms.vdouble(0.749, 0.829, 1.099, 1.355, 1.584, 1.807, 2.035, 2.217, 2.378, 2.591),
    jdpt1 = cms.vdouble(0.718, 0.813, 1.133, 1.384, 1.588, 1.841, 2.115, 2.379, 2.508, 2.772),
    jdpt2 = cms.vdouble(0.841, 0.937, 1.316, 1.605, 1.919, 2.295, 2.562, 2.722, 2.943, 3.293),
    jdpt3 = cms.vdouble(0.929, 1.04,  1.46,  1.74,  2.042, 2.289, 2.639, 2.837, 2.946, 2.971),
    jdpt4 = cms.vdouble(0.85,  0.961, 1.337, 1.593, 1.854, 2.005, 2.209, 2.533, 2.812, 3.047),
    jdpt5 = cms.vdouble(1.049, 1.149, 1.607, 1.869, 2.012, 2.219, 2.289, 2.412, 2.695, 2.865),
    jdpt6 = cms.vdouble(1.213, 1.298, 1.716, 2.015, 2.191, 2.612, 2.863, 2.879, 2.925, 2.902),
    jdpt7 = cms.vdouble(1.094, 1.139, 1.436, 1.672, 1.831, 2.05,  2.267, 2.549, 2.785, 2.86 ),
    jdpt8 = cms.vdouble(0.889, 0.939, 1.166, 1.365, 1.553, 1.805, 2.06,  2.22,  2.268, 2.247),
    jdpt9 = cms.vdouble(0.843, 0.885, 1.245, 1.665, 1.944, 1.981, 1.972, 2.875, 3.923, 7.51 ),    
    
    jdphi8 = cms.vdouble(0.059, 0.057, 0.051, 0.044, 0.038, 0.035, 0.037, 0.032, 0.028, 0.028),
    jdphi9 = cms.vdouble(0.062, 0.059, 0.053, 0.047, 0.042, 0.045, 0.036, 0.032, 0.034, 0.044),
    jdphi4 = cms.vdouble(0.042, 0.042, 0.043, 0.042, 0.038, 0.036, 0.036, 0.033, 0.031, 0.031),
    jdphi2 = cms.vdouble(0.04,  0.04,  0.04,  0.04,  0.04,  0.038, 0.036, 0.035, 0.034, 0.033),
    jdphi1 = cms.vdouble(0.034, 0.035, 0.035, 0.035, 0.035, 0.034, 0.031, 0.03,  0.029, 0.027),
    jdphi0 = cms.vdouble(0.034, 0.034, 0.034, 0.034, 0.032, 0.031, 0.028, 0.027, 0.027, 0.027),
    jdphi7 = cms.vdouble(0.077, 0.072, 0.059, 0.05,  0.045, 0.042, 0.039, 0.039, 0.037, 0.031),
    jdphi6 = cms.vdouble(0.084, 0.08,  0.072, 0.065, 0.066, 0.06,  0.051, 0.049, 0.045, 0.045),
    jdphi5 = cms.vdouble(0.069, 0.069, 0.064, 0.058, 0.053, 0.049, 0.049, 0.043, 0.039, 0.04),                    
    jdphi3 = cms.vdouble(0.042, 0.043, 0.044, 0.043, 0.041, 0.039, 0.039, 0.036, 0.034, 0.031)
  ),
  src                    = cms.InputTag("ak5PFJets"),
  srcGenJets             = cms.InputTag("ak5GenJets"),
  skipCorrJetPtThreshold = cms.double(0.01),
  skipRawJetPtThreshold  = cms.double(10.0),
  shiftBy                = cms.double(0.0),
  inputFileName          = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
  dRmaxGenJetMatch       = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*(genJetPt - 10.)))'),
  lutName                = cms.string('pfJetResolutionMCtoDataCorrLUT'),
  sigmaMaxGenJetMatch    = cms.double(5.0),
  jetCorrLabel           = cms.string('')
)

smearedCorrJetsInput = cms.EDProducer('PFJetCorrectionProducer',
  src = cms.InputTag('smearedJets'),
  correctors = cms.vstring("ak5PFL1FastL2L3") # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
)

pfCandidatesSmeared = cms.EDProducer("ShiftedPFCandidateProducerForPFMEtMVA",
  srcShiftedObjects   = cms.InputTag("smearedJets"),
  srcPFCandidates     = cms.InputTag("particleFlow"),
  dRmatch_PFCandidate = cms.double(0.5),
  srcUnshiftedObjects = cms.InputTag("ak5PFJets")
)


smearedCorrJets     = smearedCorrJetsInput.clone()
smearedCorrJets.src = 'smearedJets'

pfMEtMVA.srcCorrJets     = 'smearedCorrJets'
pfMEtMVA.srcUncorrJets   = 'smearedJets'
pfMEtMVA.srcPFCandidates = 'pfCandidatesSmeared'
pfMEtMVAUnity            = pfMEtMVA.clone()
pfMEtMVAUnity.inputFileNames = cms.PSet(
  U	= cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v2.root'),
  DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
  CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root'),
  CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root')
)

     
MVAMetSeq = cms.Sequence(pfCandsNotInJet*
                         smearedJets*
			 smearedCorrJetsInput*
			 smearedCorrJets*
			 pfCandidatesSmeared* 
                         PFTau*
			 pfMEtMVAsequence*
			 pfMEtMVAUnity)

