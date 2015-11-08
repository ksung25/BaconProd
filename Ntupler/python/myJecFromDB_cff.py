import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *

jec =  cms.ESSource("PoolDBESSource",
                         DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
                         timetype = cms.string('runnumber'),
                         toGet = cms.VPSet(
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK4PFPuppi'),
                                   label   = cms.untracked.string('AK4Puppi')
                                   ),
                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                    tag     = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK8PFPuppi'),
                                    label   = cms.untracked.string('AK8Puppi')
                                    ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK4PFchs'),
                                   label   = cms.untracked.string('AK4chs')
                                   ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK8PFchs'),
                                   label   = cms.untracked.string('AK8chs')
                                   ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK4PF'),
                                   label   = cms.untracked.string('AK4')
                                   ),
                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                    tag     = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK8PF'),
                                    label   = cms.untracked.string('AK8')
                                    )
                           ),

                    )                                        
