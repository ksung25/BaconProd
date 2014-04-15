import FWCore.ParameterSet.Config as cms

process = cms.Process('MakingJetMETExpert')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load('RecoParticleFlow/PFClusterProducer/particleFlowCluster_cff')
process.load('RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi')
process.load('RecoLocalCalo/HcalRecAlgos/hcalRecAlgoESProd_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'START53_V7G::All'


process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
process.load("BaconProd/Ntupler/myRho_cff")

# trigger filter
import os
cmssw_base = os.environ['CMSSW_BASE']
hlt_filename = "BaconAna/DataFormats/data/HLTFile_v0"
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(False)
process.hltHighLevel.HLTPaths = cms.vstring()
hlt_file = open(cmssw_base + "/src/" + hlt_filename, "r")
for line in hlt_file.readlines():
  line = line.strip()              # strip preceding and trailing whitespaces
  if (line[0:3] == 'HLT'):         # assumes typical lines begin with HLT path name (e.g. HLT_Mu15_v1)
    hlt_path = line.split()[0]
    process.hltHighLevel.HLTPaths.extend(cms.untracked.vstring(hlt_path))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
#  firstRun   = cms.untracked.uint32(0),
#  firstLumi  = cms.untracked.uint32(0),
#  firstEvent = cms.untracked.uint32(0),
  fileNames  = cms.untracked.vstring('/store/relval/CMSSW_5_3_6/RelValTTbar/GEN-SIM-RECO/PU_START53_V14-v1/0003/3E3EDF4A-E92C-E211-A1BF-003048D2BD66.root')
#/store/relval/CMSSW_5_3_14/RelValTTbar/GEN-SIM-RECO/START53_LV6_Jan31-v1/00000/348BFC23-468B-E311-AAC5-0025905A6118.root')
#/store/cmst3/group/cmgtools/CMG/ggH125_HTobb_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PFAOD_99.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

is_data_flag = False
do_hlt_filter = False
process.ntupler = cms.EDAnalyzer('ExpertMod',
  skipOnHLTFail = cms.untracked.bool(do_hlt_filter),
  outputName    = cms.untracked.string('expert.root'),
  TriggerFile   = cms.untracked.string(hlt_filename), 

  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    edmPFCandName        = cms.untracked.string('particleFlow'),
    edmPileupInfoName    = cms.untracked.string('addPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot'),
    edmPFMETName         = cms.untracked.string('pfMet'),
    edmPFMETCorrName     = cms.untracked.string('pfType0p1CorrectedMet'),
    edmMVAMETName        = cms.untracked.string('pfMEtMVA'),
    edmMVAMETUnityName   = cms.untracked.string('pfMEtMVAUnity'),
    edmMVAMETNoSmearName = cms.untracked.string('pfMEtMVANoSmear'),
    edmRhoForIsoName     = cms.untracked.string('kt6PFJets'),
    edmRhoForJetEnergy   = cms.untracked.string('kt6PFJets'),
    doFillMET            = cms.untracked.bool(False)
  ),

  GenInfo = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmGenEventInfoName = cms.untracked.string('generator'),
    edmGenParticlesName = cms.untracked.string('genParticles'),
    fillAllGen          = cms.untracked.bool(False)
  ),
  
  PV = cms.untracked.PSet(
    isActive      = cms.untracked.bool(True),   
    edmName       = cms.untracked.string('offlinePrimaryVertices'),
    minNTracksFit = cms.untracked.uint32(0),
    minNdof       = cms.untracked.double(4),
    maxAbsZ       = cms.untracked.double(24),
    maxRho        = cms.untracked.double(2)
  ),
  
  PFCand = cms.untracked.PSet(
    isActive       = cms.untracked.bool(True),
    edmName        = cms.untracked.string('particleFlow'),
    edmPVName      = cms.untracked.string('offlinePrimaryVertices'),
    doAddDepthTime = cms.untracked.bool(False)
  )
)

process.baconSequence = cms.Sequence(process.particleFlowCluster*
                                     process.kt6PFJets*
                                     process.ntupler)
				     
if do_hlt_filter:
  process.p = cms.Path(process.hltHighLevel*process.baconSequence)
else:
  process.p = cms.Path(process.baconSequence)

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                      
                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                        
                                  fileName       = cms.untracked.string ("test.root")                                                                                                                      
)                                                                                                                                                                                                          

# schedule definition                                                                                                                                                                                      
#process.schedule = cms.Schedule(process.bambu_step)
#process.outpath  = cms.EndPath(process.output)                                                                                                                                                             

#
# simple checks to catch some mistakes...
#
if is_data_flag:
  assert process.ntupler.GenInfo.isActive == cms.untracked.bool(False)
