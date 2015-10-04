import FWCore.ParameterSet.Config as cms

process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.jec =  cms.ESSource("PoolDBESSource",
#                            DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
#                            timetype = cms.string('runnumber'),
#                            toGet = cms.VPSet(
#                                      cms.PSet(record  = cms.string('JetCorrectionsRecord'),
#                                               tag     = cms.string('JetCorrectorParametersCollection_PY8_RunIISpring15DR74_bx50_MC_AK4PFchs'),
#                                               label   = cms.untracked.string('AK4puppi')
#                                              ),
#                                      cms.PSet(record  = cms.string('JetCorrectionsRecord'),
#                                               tag     = cms.string('JetCorrectorParametersCollection_PY8_RunIISpring15DR74_bx50_MC_AK8PFchs'),
#                                               label   = cms.untracked.string('AK8puppi')
#                                              )
#                                      ),
#                            connect = cms.string('sqlite:PY8_RunIISpring15DR74_bx50_MC.db'),
#                            )                                        
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')



# define your favorite global tag
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string('74X_mcRun2_asymptotic_v2')

#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag.globaltag = '74X_mcRun2_design_v2::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK4PF    = cms.EDAnalyzer('JetCorrectorDBReader',
                                      # below is the communication to the database
                                      #payloadName    = cms.untracked.string('AK4PFPuppi'),
                                      payloadName    = cms.untracked.string('AK4PFchs'),
                                      globalTag      = cms.untracked.string('T2'),
                                      printScreen    = cms.untracked.bool(False),
                                      createTextFile = cms.untracked.bool(True)
                                      )
process.readAK8PF    = cms.EDAnalyzer('JetCorrectorDBReader',
                                      # below is the communication to the database
                                      payloadName    = cms.untracked.string('AK8PFPuppi'),
                                      #payloadName    = cms.untracked.string('AK8PFchs'),
                                      globalTag      = cms.untracked.string('T2'),
                                      printScreen    = cms.untracked.bool(False),
                                      createTextFile = cms.untracked.bool(True)
                                      )

process.readAK4PFCHS = cms.EDAnalyzer('JetCorrectorDBReader',  
      # below is the communication to the database 
      payloadName    = cms.untracked.string('AK4PFchs'),
      # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
      # but it is recommended to use the GT name that you retrieved the files from.
      globalTag      = cms.untracked.string('MCRUN2_74_V9A'),  
      printScreen    = cms.untracked.bool(False),
      createTextFile = cms.untracked.bool(True)
)
#process.readAK4Calo = process.readAK4PFCHS.clone(payloadName = 'AK4Calo')
#process.readAK5JPT = process.readAK5PF.clone(payloadName = 'AK5JPT')
#process.readAK4CHS  = process.readAK4PF.clone(payloadName = 'AK4PFchs')
#process.p = cms.Path(process.readAK4PF*process.readAK8PF)
process.p = cms.Path(process.readAK4PF*process.readAK8PF)
#*process.readAK4PF)# * process.readAK4CHS)

