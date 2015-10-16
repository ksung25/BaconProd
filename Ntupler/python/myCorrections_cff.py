import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *

#puppijec =  cms.ESSource("PoolDBESSource",
#                         DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
#                         timetype = cms.string('runnumber'),
#                         toGet = cms.VPSet(
#                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
#                                    tag     = cms.string('JetCorrectorParametersCollection_PY8_RunIISpring15DR74_bx50_MC_AK4PF'),
#                                    label   = cms.untracked.string('AK4puppi')
#                                    ),
#                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
#                                    tag     = cms.string('JetCorrectorParametersCollection_PY8_RunIISpring15DR74_bx50_MC_AK8PF'),
#                                    label   = cms.untracked.string('AK8puppi')
#                                    )
#                           ),
                         #connect = cms.string('sqlite:PY8_RunIISpring15DR74_bx50_MC.dbX'),
                         #connect = cms.string('sqlite:///BaconProd/Utils/data/PY8_RunIISpring15DR74_bx50_MC.db'),
#                         )                                        

# Sequence AK4
label='PFchs'
ak4L1FastjetCorrector  = ak4PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK4'+label))
ak4L2RelativeCorrector = ak4PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK4'+label))
ak4L3AbsoluteCorrector = ak4PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK4'+label))
ak4ResidualCorrector   = ak4PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK4'+label))

ak4L1FastL2L3Corrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4L1FastjetCorrector", "ak4L2RelativeCorrector", "ak4L3AbsoluteCorrector")
    )
ak4L1FastL2L3ResidualCorrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4L1FastjetCorrector", "ak4L2RelativeCorrector", "ak4L3AbsoluteCorrector",'ak4ResidualCorrector')
    )
ak4L1FastL2L3Chain = cms.Sequence(
    ak4L1FastjetCorrector * ak4L2RelativeCorrector * ak4L3AbsoluteCorrector * ak4L1FastL2L3Corrector
)
ak4L1FastL2L3ResidualChain = cms.Sequence(
    ak4L1FastjetCorrector * ak4L2RelativeCorrector * ak4L3AbsoluteCorrector * ak4ResidualCorrector * ak4L1FastL2L3ResidualCorrector
)
# sequence CA8
# Sequence               
ak8L1FastjetCorrector  = ak8PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK8'+label))
ak8L2RelativeCorrector = ak8PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK8'+label))
ak8L3AbsoluteCorrector = ak8PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK8'+label))
ak8ResidualCorrector   = ak8PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK8'+label))
ak8L1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8L1FastjetCorrector','ak8L2RelativeCorrector','ak8L3AbsoluteCorrector')
)
ak8L1FastL2L3Chain = cms.Sequence(
  ak8L1FastjetCorrector * ak8L2RelativeCorrector * ak8L3AbsoluteCorrector * ak8L1FastL2L3Corrector
)
ak8L1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8L1FastjetCorrector','ak8L2RelativeCorrector','ak8L3AbsoluteCorrector','ak8ResidualCorrector')
)
ak8L1FastL2L3ResidualChain = cms.Sequence(
  ak8L1FastjetCorrector * ak8L2RelativeCorrector * ak8L3AbsoluteCorrector * ak8ResidualCorrector * ak8L1FastL2L3ResidualCorrector
)
# Sequnce CA15
ca15L1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8L1FastjetCorrector','ak8L2RelativeCorrector','ak8L3AbsoluteCorrector')
)
ca15L1FastL2L3Chain = cms.Sequence(
  ak8L1FastjetCorrector * ak8L2RelativeCorrector * ak8L3AbsoluteCorrector * ca15L1FastL2L3Corrector
)
ca15L1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8L1FastjetCorrector','ak8L2RelativeCorrector','ak8L3AbsoluteCorrector','ak8ResidualCorrector')
)
ca15L1FastL2L3ResidualChain = cms.Sequence(
  ak8L1FastjetCorrector * ak8L2RelativeCorrector * ak8L3AbsoluteCorrector * ak8ResidualCorrector * ca15L1FastL2L3ResidualCorrector
)
