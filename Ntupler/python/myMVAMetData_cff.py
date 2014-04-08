import FWCore.ParameterSet.Config as cms
from RecoMET.METPUSubtraction.mvaPFMET_leptons_cff import *
from RecoTauTag.Configuration.RecoPFTauTag_cff import *

calibratedAK5PFJetsForPFMEtMVA.correctors =  cms.vstring("ak5PFL1FastL2L3Residual") 
pfMEtMVAUnity            = pfMEtMVA.clone()
pfMEtMVAUnity.inputFileNames = cms.PSet(
  U	= cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_53_Sep2013_type1_UnityResponse_v2.root'),
  DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
  CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root'),
  CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root')
)

     
MVAMetSeq = cms.Sequence(
                         PFTau*
			 pfMEtMVAsequence*
			 pfMEtMVAUnity)

