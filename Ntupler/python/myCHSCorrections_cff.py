import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *

#chs Sequence AK4
chslabel='chs'
ak4chsL1FastjetCorrector  = ak4PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK4'+chslabel))
ak4chsL2RelativeCorrector = ak4PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK4'+chslabel))
ak4chsL3AbsoluteCorrector = ak4PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK4'+chslabel))
ak4chsResidualCorrector   = ak4PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK4'+chslabel))

ak4chsL1FastL2L3Corrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4chsL1FastjetCorrector", "ak4chsL2RelativeCorrector", "ak4chsL3AbsoluteCorrector")
    )
ak4chsL1FastL2L3ResidualCorrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4chsL1FastjetCorrector", "ak4chsL2RelativeCorrector", "ak4chsL3AbsoluteCorrector",'ak4chsResidualCorrector')
    )
ak4chsL1FastL2L3Chain = cms.Sequence(
    ak4chsL1FastjetCorrector * ak4chsL2RelativeCorrector * ak4chsL3AbsoluteCorrector * ak4chsL1FastL2L3Corrector
)
ak4chsL1FastL2L3ResidualChain = cms.Sequence(
    ak4chsL1FastjetCorrector * ak4chsL2RelativeCorrector * ak4chsL3AbsoluteCorrector * ak4chsResidualCorrector * ak4chsL1FastL2L3ResidualCorrector
)
#chs sequence CA8
#chs Sequence               
ak8chsL1FastjetCorrector  = ak8PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK8'+chslabel))
ak8chsL2RelativeCorrector = ak8PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK8'+chslabel))
ak8chsL3AbsoluteCorrector = ak8PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK8'+chslabel))
ak8chsResidualCorrector   = ak8PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK8'+chslabel))
ak8chsL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8chsL1FastjetCorrector','ak8chsL2RelativeCorrector','ak8chsL3AbsoluteCorrector')
)
ak8chsL1FastL2L3Chain = cms.Sequence(
  ak8chsL1FastjetCorrector * ak8chsL2RelativeCorrector * ak8chsL3AbsoluteCorrector * ak8chsL1FastL2L3Corrector
)
ak8chsL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8chsL1FastjetCorrector','ak8chsL2RelativeCorrector','ak8chsL3AbsoluteCorrector','ak8chsResidualCorrector')
)
ak8chsL1FastL2L3ResidualChain = cms.Sequence(
  ak8chsL1FastjetCorrector * ak8chsL2RelativeCorrector * ak8chsL3AbsoluteCorrector * ak8chsResidualCorrector * ak8chsL1FastL2L3ResidualCorrector
)
#chs Sequnce CA15
ca15chsL1FastL2L3Corrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8chsL1FastjetCorrector','ak8chsL2RelativeCorrector','ak8chsL3AbsoluteCorrector')
)
ca15chsL1FastL2L3Chain = cms.Sequence(
  ak8chsL1FastjetCorrector * ak8chsL2RelativeCorrector * ak8chsL3AbsoluteCorrector * ca15chsL1FastL2L3Corrector
)
ca15chsL1FastL2L3ResidualCorrector = cms.EDProducer('ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8chsL1FastjetCorrector','ak8chsL2RelativeCorrector','ak8chsL3AbsoluteCorrector','ak8chsResidualCorrector')
)
ca15chsL1FastL2L3ResidualChain = cms.Sequence(
  ak8chsL1FastjetCorrector * ak8chsL2RelativeCorrector * ak8chsL3AbsoluteCorrector * ak8chsResidualCorrector * ca15chsL1FastL2L3ResidualCorrector
)
