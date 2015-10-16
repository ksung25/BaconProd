import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *

#puppijec =  cms.ESSource("PoolDBESSource",
#                         DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
#                         timetype = cms.string('runnumber'),
#                         toGet = cms.VPSet(
#                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
#                                    tag     = cms.string('JetCorrectorParametersCollection_PY8_RunIISpring15DR74_bx50_MC_AK4PFchs'),
#                                    label   = cms.untracked.string('AK4puppi')
#                                    ),
#                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
#                                    tag     = cms.string('JetCorrectorParametersCollection_PY8_RunIISpring15DR74_bx50_MC_AK8PFchs'),
#                                    label   = cms.untracked.string('AK8puppi')
#                                    )
#                           ),
                         #connect = cms.string('sqlite:PY8_RunIISpring15DR74_bx50_MC.dbX'),
                         #connect = cms.string('sqlite:///BaconProd/Utils/data/PY8_RunIISpring15DR74_bx50_MC.db'),
#                         )                                        

#Puppi Sequence AK4
puppilabel='PFPuppi'
ak4PuppiL1FastjetCorrector  = ak4PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK4'+puppilabel))
ak4PuppiL2RelativeCorrector = ak4PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK4'+puppilabel))
ak4PuppiL3AbsoluteCorrector = ak4PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK4'+puppilabel))
ak4PuppiResidualCorrector   = ak4PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK4'+puppilabel))

ak4PuppiL1FastL2L3Corrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4PuppiL1FastjetCorrector", "ak4PuppiL2RelativeCorrector", "ak4PuppiL3AbsoluteCorrector")
    )
ak4PuppiL1FastL2L3ResidualCorrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4PuppiL1FastjetCorrector", "ak4PuppiL2RelativeCorrector", "ak4PuppiL3AbsoluteCorrector",'ak4PuppiResidualCorrector')
    )
ak4PuppiL1FastL2L3Chain = cms.Sequence(
    ak4PuppiL1FastjetCorrector * ak4PuppiL2RelativeCorrector * ak4PuppiL3AbsoluteCorrector * ak4PuppiL1FastL2L3Corrector
)
ak4PuppiL1FastL2L3ResidualChain = cms.Sequence(
    ak4PuppiL1FastjetCorrector * ak4PuppiL2RelativeCorrector * ak4PuppiL3AbsoluteCorrector * ak4PuppiResidualCorrector * ak4PuppiL1FastL2L3ResidualCorrector
)
#Puppi sequence CA8
#Puppi Sequence               
ak8PuppiL1FastjetCorrector  = ak8PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK8'+puppilabel))
ak8PuppiL2RelativeCorrector = ak8PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK8'+puppilabel))
ak8PuppiL3AbsoluteCorrector = ak8PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK8'+puppilabel))
ak8PuppiResidualCorrector   = ak8PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK8'+puppilabel))
ak8PuppiL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PuppiL1FastjetCorrector','ak8PuppiL2RelativeCorrector','ak8PuppiL3AbsoluteCorrector')
)
ak8PuppiL1FastL2L3Chain = cms.Sequence(
  ak8PuppiL1FastjetCorrector * ak8PuppiL2RelativeCorrector * ak8PuppiL3AbsoluteCorrector * ak8PuppiL1FastL2L3Corrector
)
ak8PuppiL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PuppiL1FastjetCorrector','ak8PuppiL2RelativeCorrector','ak8PuppiL3AbsoluteCorrector','ak8PuppiResidualCorrector')
)
ak8PuppiL1FastL2L3ResidualChain = cms.Sequence(
  ak8PuppiL1FastjetCorrector * ak8PuppiL2RelativeCorrector * ak8PuppiL3AbsoluteCorrector * ak8PuppiResidualCorrector * ak8PuppiL1FastL2L3ResidualCorrector
)
#Puppi Sequnce CA15
ca15PuppiL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PuppiL1FastjetCorrector','ak8PuppiL2RelativeCorrector','ak8PuppiL3AbsoluteCorrector')
)
ca15PuppiL1FastL2L3Chain = cms.Sequence(
  ak8PuppiL1FastjetCorrector * ak8PuppiL2RelativeCorrector * ak8PuppiL3AbsoluteCorrector * ca15PuppiL1FastL2L3Corrector
)
ca15PuppiL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PuppiL1FastjetCorrector','ak8PuppiL2RelativeCorrector','ak8PuppiL3AbsoluteCorrector','ak8PuppiResidualCorrector')
)
ca15PuppiL1FastL2L3ResidualChain = cms.Sequence(
  ak8PuppiL1FastjetCorrector * ak8PuppiL2RelativeCorrector * ak8PuppiL3AbsoluteCorrector * ak8PuppiResidualCorrector * ca15PuppiL1FastL2L3ResidualCorrector
)
