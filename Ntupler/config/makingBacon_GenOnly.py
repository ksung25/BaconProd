import FWCore.ParameterSet.Config as cms

process = cms.Process('MakingBacon')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')

process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(
         1000022,
         1000012, 1000014, 1000016,
         2000012, 2000014, 2000016,
         1000039, 5100039,
         4000012, 4000014, 4000016,
         9900012, 9900014, 9900016,
         39),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)
process.genParticlesForJetsNoNu = process.genParticlesForJets.clone()
process.genParticlesForJetsNoNu.ignoreParticleIDs += cms.vuint32(12,14,16,18)

# Select hadrons and partons for Jet Flavour                                                                                                                                                                  
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.GlobalTag.globaltag = 'START53_V7G::All'

# import custom configurations
#process.load('BaconProd/Ntupler/myGenJets_cff')            # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtrasAK4CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtrasAK8CHS_cff')    # include gen jets and b-tagging

# trigger filter
import os
cmssw_base = os.environ['CMSSW_BASE']

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
  fileNames  = cms.untracked.vstring('file:test.root')
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
hlt_filename = "BaconAna/DataFormats/data/HLTFile_v0"
process.ntupler = cms.EDAnalyzer('NtuplerMod',
  useTrigger    = cms.untracked.bool(False),
  skipOnHLTFail = cms.untracked.bool(do_hlt_filter),
  outputName    = cms.untracked.string('MonoJ.root'),
  TriggerFile   = cms.untracked.string(hlt_filename),
  edmPVName     = cms.untracked.string('offlinePrimaryVertices'),
  edmPFCandName = cms.untracked.string('particleFlow'),
  
  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(False),
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
    doFillMET            = cms.untracked.bool(True)
  ),
  
  GenInfo = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmGenEventInfoName = cms.untracked.string('generator'),
    edmLHEEventInfoName = cms.untracked.string('source'),
    edmGenParticlesName = cms.untracked.string('genParticles'),
    fillLHEWeights      = cms.untracked.bool(True),
    fillAllGen          = cms.untracked.bool(True)
  ),

  GenJet  = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    isActiveFatJet      = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmGenParticlesName = cms.untracked.string('genParticles'),
    genJetName          = cms.untracked.string('AK4GenJetsCHS'),
    genFatJetName       = cms.untracked.string('AK8GenJetsCHS'),
    fillAllGen          = cms.untracked.bool(True)
  ),
)
#process.genParticlesForJets.src     = 'source'
#process.genParticlesForJetsNoNu.src = 'source'
#process.genParticles.src            = 'source'
process.baconSequence = cms.Sequence(
    #                                     process.GeneInfo*
                                     #process.genjetsequence*
                                     process.genParticlesForJetsNoNu*
                                     process.AK4genjetsequenceCHS*
                                     process.AK8GenJetsCHS*
                                     #process.AK8genjetsequenceCHS*                        
                                     process.ntupler)
				     
process.p = cms.Path(process.baconSequence)


