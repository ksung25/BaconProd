import FWCore.ParameterSet.Config as cms

process = cms.Process('MakingBacon')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'START53_V7G::All'

process.load('CommonTools/ParticleFlow/PFBRECO_cff')
process.load('Dummy/Puppi/Puppi_cff')
process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

# import custom configurations
process.load('BaconProd/Ntupler/myJetExtras04_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras05_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras06_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras07_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras08_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras09_cff')    # include gen jets and b-tagging

process.load('BaconProd/Ntupler/myJetExtras04CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras05CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras06CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras07CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras08CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras09CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras12CHS_cff')    # include gen jets and b-tagging

process.load('BaconProd/Ntupler/myJetExtras04Puppi_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras05Puppi_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras06Puppi_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras07Puppi_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras08Puppi_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras09Puppi_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras12Puppi_cff')    # include gen jets and b-tagging

process.load('BaconProd/Ntupler/myMETFilters_cff')        # apply MET filters set to tagging mode
process.load('BaconProd/Ntupler/myMVAMet_cff')            # MVA MET
process.load("BaconProd/Ntupler/myPFMETCorrections_cff")  # PF MET corrections
process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")

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
  fileNames  = cms.untracked.vstring('AAA')
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
process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail = cms.untracked.bool(do_hlt_filter),
  outputName    = cms.untracked.string('BBB'),
  TriggerFile   = cms.untracked.string(hlt_filename),
  edmPVName     = cms.untracked.string('offlinePrimaryVertices'),
  edmPFCandName = cms.untracked.string('particleFlow'),
  
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
    doFillMET            = cms.untracked.bool(True)
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
  
  Electron = cms.untracked.PSet(
    isActive                  = cms.untracked.bool(True),
    minPt                     = cms.untracked.double(5),
    edmName                   = cms.untracked.string('gsfElectrons'),
    edmPFCandName             = cms.untracked.string('particleFlow'),
    edmTrackName              = cms.untracked.string('generalTracks'),
    edmConversionName         = cms.untracked.string('allConversions'),
    edmRhoForEnergyRegression = cms.untracked.string('kt6PFJets'),
    edmEBSuperClusterName     = cms.untracked.string('correctedHybridSuperClusters'),
    edmEESuperCClusterName    = cms.untracked.string('correctedMulti5x5SuperClustersWithPreshower'),
    edmEBRecHitName           = cms.untracked.string('reducedEcalRecHitsEB'),
    edmEERecHitName           = cms.untracked.string('reducedEcalRecHitsEE'),
    
    # ORDERED list of weight files for electron MVA ID (specify paths relative to $CMSSW_BASE/src)
    eleIDFilenames = cms.untracked.vstring('BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml',
                                           'BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml',
					   'BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml',
					   'BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml',
					   'BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml',
					   'BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml'),
    
    # inputs for HZZ4l electron momentum corrections (specify paths relative to $CMSSW_BASE/src)
    doHZZ4lCorr         = cms.untracked.bool(True),    
    isData              = ( cms.untracked.bool(True) if is_data_flag else cms.untracked.bool(False) ),
    doRandForMC         = cms.untracked.bool(True),
    eleEnergyRegWeights = cms.untracked.string('BaconProd/Utils/data/eleEnergyRegWeights_WithSubClusters_VApr15.root'),
    scalesCorr          = cms.untracked.string('BaconProd/Utils/data/scalesCorr.csv'),
    smearsCorrType1     = cms.untracked.string('BaconProd/Utils/data/smearsCorrType1.csv'),
    smearsCorrType2     = cms.untracked.string('BaconProd/Utils/data/smearsCorrType2.csv'),
    smearsCorrType3     = cms.untracked.string('BaconProd/Utils/data/smearsCorrType3.csv'),
    linearityCorr       = cms.untracked.string('BaconProd/Utils/data/linearityNewReg-May2013.csv')
  ),
  
  Muon = cms.untracked.PSet(
    isActive      = cms.untracked.bool(True),
    minPt         = cms.untracked.double(0),
    edmName       = cms.untracked.string('muons'),
    edmPFCandName = cms.untracked.string('particleFlow'),
    
    # inputs for HZZ4l muon momentum corrections (specify paths relative to $CMSSW_BASE/src)
    doHZZ4lCorr  = cms.untracked.bool(True),
    isData       = ( cms.untracked.bool(True) if is_data_flag else cms.untracked.bool(False) ),
    doRandForMC  = cms.untracked.bool(True),
    muScleFitDir = cms.untracked.string('MuScleFit/Calibration/data'),
    
    # save general tracker tracks in our muon collection (used in tag-and-probe for muons)
    doSaveTracks = cms.untracked.bool(False),
    minTrackPt   = cms.untracked.double(20),
    edmTrackName = cms.untracked.string('generalTracks')
    
    ),
  
  Photon = cms.untracked.PSet(
    isActive              = cms.untracked.bool(True),
    minPt                 = cms.untracked.double(0),
    edmName               = cms.untracked.string('photons'),
    edmPFCandName         = cms.untracked.string('particleFlow'),
    edmElectronName       = cms.untracked.string('gsfElectrons'),
    edmConversionName     = cms.untracked.string('allConversions'),
    edmEBSuperClusterName = cms.untracked.string('correctedHybridSuperClusters'),
    edmEESuperClusterName = cms.untracked.string('correctedMulti5x5SuperClustersWithPreshower'),
    edmEBRecHitName       = cms.untracked.string('reducedEcalRecHitsEB'),
    edmEERecHitName       = cms.untracked.string('reducedEcalRecHitsEE'),
    edmRhoForEnergyRegression = cms.untracked.string('kt6PFJets'),
    edmPVName                 = cms.untracked.string('offlinePrimaryVertices')
  ),
  
  Tau = cms.untracked.PSet(
    isActive = cms.untracked.bool(True),
    minPt    = cms.untracked.double(15),
    edmName  = cms.untracked.string('hpsPFTauProducer'),
    ringIsoFile      = cms.untracked.string('BaconProd/Utils/data/gbrfTauIso_apr29a.root'),
    ringIso2File     = cms.untracked.string('BaconProd/Utils/data/gbrfTauIso_v2.root'),
    edmRhoForRingIso = cms.untracked.string('kt6PFJets')
  ),
  
  Jet = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    
    coneSizes = cms.untracked.vdouble(0.5),#,0.8,1.2),
    postFix   = cms.untracked.vstring("CHS"),
    
    edmPVName = cms.untracked.string('offlinePrimaryVertices'),
    
    # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt',
                                       'BaconProd/Utils/data/Summer13_V1_DATA_L2Relative_AK5PF.txt',
				       'BaconProd/Utils/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt',
				       'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt')
                 if is_data_flag else 
                 cms.untracked.vstring('BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt',
                                       'BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt',
				       'BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt')
	       ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/Summer13_V1_MC_Uncertainty_AK5PF.txt')
                  ),
    jecFilesForID = ( cms.untracked.vstring('BaconProd/Utils/data/FT_53_V21_AN3_L1FastJet_AK5PF.txt',
                                            'BaconProd/Utils/data/FT_53_V21_AN3_L2Relative_AK5PF.txt',
				            'BaconProd/Utils/data/FT_53_V21_AN3_L3Absolute_AK5PF.txt',
					    'BaconProd/Utils/data/FT_53_V21_AN3_L2L3Residual_AK5PF.txt')
                      if is_data_flag else
                      cms.untracked.vstring('BaconProd/Utils/data/START53_V15_L1FastJet_AK5PF.txt',
                                            'BaconProd/Utils/data/START53_V15_L2Relative_AK5PF.txt',
				            'BaconProd/Utils/data/START53_V15_L3Absolute_AK5PF.txt')
		    ),
    edmRhoName = cms.untracked.string('kt6PFJets'),
    
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml'),
    
    # names of various jet-related collections WITHOUT prefix (e.g. 'PFJets' instead of 'AK5PFJets')
    # prefix string will be determined by ntupler module based on cone size
    jetName            = cms.untracked.string('PFJets'),
    genJetName         = cms.untracked.string('GenJets'),
    jetFlavorName      = cms.untracked.string('byValAlgo'),
    jetFlavorPhysName  = cms.untracked.string('byValPhys'),
    pruneJetName       = cms.untracked.string('caPFJetsPruned'),
    subJetName         = cms.untracked.string('caPFJetsPruned'),
    csvBTagName        = cms.untracked.string('jetCombinedSecondaryVertexBJetTags'),
    csvBTagSubJetName  = cms.untracked.string('jetCombinedSecondaryVertexBJetTagsSJ'),
    jettiness          = cms.untracked.string('Njettiness'),
    qgLikelihood       = cms.untracked.string('QGTagger'),
    qgLikelihoodSubjet = cms.untracked.string('QGTaggerSubJets')
  ),
  
  PFCand = cms.untracked.PSet(
    isActive       = cms.untracked.bool(False),
    edmName        = cms.untracked.string('particleFlow'),
    edmPVName      = cms.untracked.string('offlinePrimaryVertices'),
    doAddDepthTime = cms.untracked.bool(False)
  )
)

process.baconSequence = cms.Sequence(process.PFBRECO*
                                     #process.puppi*
                                     process.metFilters*
                                     process.producePFMETCorrections*
                                     process.recojetsequence*
                                     process.genjetsequence*
                                     process.AK5jetsequenceCHS*
                                     process.AK5genjetsequence*
                                     process.recoTau*   ### must come after antiktGenJets otherwise conflict on RecoJets/JetProducers/plugins
				     process.MVAMetSeq*
				     process.ntupler)
				     
if do_hlt_filter:
  process.p = cms.Path(process.hltHighLevel*process.baconSequence)
else:
  process.p = cms.Path(process.baconSequence)

#process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
#                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                      
#                                  fileName       = cms.untracked.string ("BBB")                                                                                                                    
#)

# schedule definition                                                                                                       
#process.outpath  = cms.EndPath(process.output)                                                                                                                                                

#
# simple checks to catch some mistakes...
#
if is_data_flag:
  assert process.ntupler.GenInfo.isActive == cms.untracked.bool(False)
  assert process.ntupler.Electron.isData == cms.untracked.bool(True)
  assert process.ntupler.Muon.isData == cms.untracked.bool(True)
  assert process.ntupler.Jet.doGenJet == cms.untracked.bool(False)
else:
  assert process.ntupler.Electron.isData == cms.untracked.bool(False)
  assert process.ntupler.Muon.isData == cms.untracked.bool(False)
