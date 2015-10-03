import FWCore.ParameterSet.Config as cms
import os

process = cms.Process('MakingBacon')

is_data_flag  = False                                    # flag for if process data
do_hlt_filter = False                                    # flag to skip events that fail relevant triggers
hlt_filename  = "BaconAna/DataFormats/data/HLTFile_v2"   # list of relevant triggers

cmssw_base = os.environ['CMSSW_BASE']

#--------------------------------------------------------------------------------
# Import of standard configurations
#================================================================================
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

process.GlobalTag.globaltag = 'START53_V27::All'

process.load('BaconProd/Ntupler/PFBRECO_v2_cff')
process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

#--------------------------------------------------------------------------------
# Import custom configurations
#================================================================================

process.load('BaconProd/Ntupler/myJetExtras05_cff')       # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras08CHS_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras15CHS_cff')    # include gen jets and b-tagging

process.load('BaconProd/Ntupler/myMETFilters_cff')        # apply MET filters set to tagging mode
process.load('BaconProd/Ntupler/myMVAMet_cff')            # MVA MET
process.load("BaconProd/Ntupler/myPFMETCorrections_cff")  # PF MET corrections
if is_data_flag:
  process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
else:
  process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")

#--------------------------------------------------------------------------------
# input settings
#================================================================================
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
  fileNames  = cms.untracked.vstring('file:/afs/cern.ch/work/k/ksung/private/HZZ4lAna/temp/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A_FEEEEFFF-7FFB-E111-8FE2-002618943810.root')
#  fileNames  = cms.untracked.vstring('/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v1/20000/00277FF2-7B84-E211-9475-782BCB27B958.root')
#  fileNames  = cms.untracked.vstring('/store/mc/Summer12_DR53X/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0000/00065798-2704-E211-B308-0025901D4C44.root')
#  fileNames  = cms.untracked.vstring('/store/mc/Summer12_DR53X/DYJetsToLL_PtZ-100_TuneZ2star_8TeV_ext-madgraph-tarball/AODSIM/PU_S10_START53_V7C-v1/00000/001B91CE-7639-E211-B7D1-00261894385A.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

#--------------------------------------------------------------------------------
# Reporting
#================================================================================
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

#--------------------------------------------------------------------------------
# Bacon making settings
#================================================================================
process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail = cms.untracked.bool(do_hlt_filter),
  outputName    = cms.untracked.string('Output.root'),
  TriggerFile   = cms.untracked.string(hlt_filename),
  edmPVName     = cms.untracked.string('offlinePrimaryVertices'),
  edmPFCandName = cms.untracked.string('particleFlow'),
  
  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    edmPFCandName        = cms.untracked.string('particleFlow'),
    edmPileupInfoName    = cms.untracked.string('addPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot'),
    edmPFMETName         = cms.untracked.string('pfMet'),
    edmPFMETCorrName     = cms.untracked.string('pfType1CorrectedMet'),
    edmMVAMETName        = cms.untracked.string('pfMEtMVA'),
    edmMVAMETUnityName   = cms.untracked.string('pfMEtMVAUnity'),
    edmMVAMETNoSmearName = cms.untracked.string('pfMEtMVANoSmear'),
    edmRhoForIsoName     = cms.untracked.string('kt6PFJets'),
    edmRhoForJetEnergy   = cms.untracked.string('kt6PFJets'),
    doFillMET            = cms.untracked.bool(True),
    doFillMETFilters     = cms.untracked.bool(True),
    addSusyGen           = cms.untracked.bool(False)
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
    minPt    = cms.untracked.double(10),
    edmName  = cms.untracked.string('hpsPFTauProducer'),
    ringIsoFile      = cms.untracked.string('BaconProd/Utils/data/gbrfTauIso_apr29a.root'),
    ringIso2File     = cms.untracked.string('BaconProd/Utils/data/gbrfTauIso_v2.root'),
    edmRhoForRingIso = cms.untracked.string('kt6PFJets')
  ),
 
  AK5 = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    coneSize             = cms.untracked.double(0.5),
    doComputeFullJetInfo = cms.untracked.bool(False),
    topTagType           = cms.untracked.string('none'),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

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
                                         'RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml'),

    # names of various jet-related collections
    jetName            = cms.untracked.string('AK5PFJets'),
    genJetName         = cms.untracked.string('AK5GenJets'),
    jetFlavorName      = cms.untracked.string('AK5byValAlgo'),
    jetFlavorPhysName  = cms.untracked.string('AK5byValPhys'),
    pruneJetName       = cms.untracked.string('AK5caPFJetsPruned'),
    subJetName         = cms.untracked.string('AK5caPFJetsPruned'),
    csvBTagName        = cms.untracked.string('AK5jetCombinedSecondaryVertexBJetTags'),
    csvBTagSubJetName  = cms.untracked.string('AK5jetCombinedSecondaryVertexBJetTagsSJ'),
    jettiness          = cms.untracked.string('AK5Njettiness'),
    qgLikelihood       = cms.untracked.string('AK5QGTagger'),
    qgLikelihoodSubjet = cms.untracked.string('AK5QGTaggerSubJets')
  ),

  CA8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(150),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    topTagType           = cms.untracked.string('CMS'),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

    # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L2L3Residual_AK7PFchs.txt')
                 if is_data_flag else
                 cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt')
               ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_Uncertainty_AK7PFchs.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_Uncertainty_AK7PFchs.txt')
                  ),
    jecFilesForID = ( cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L2L3Residual_AK7PFchs.txt')
                      if is_data_flag else
                      cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt')
                    ),
    edmRhoName = cms.untracked.string('kt6PFJets'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),

    # names of various jet-related collections
    jetName            = cms.untracked.string('CA8PFJetsCHS'),
    genJetName         = cms.untracked.string('CA8GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA8byValAlgoCHS'),
    jetFlavorPhysName  = cms.untracked.string('CA8byValPhysCHS'),
    pruneJetName       = cms.untracked.string('CA8caPFJetsPrunedCHS'),
    subJetName         = cms.untracked.string('CA8caPFJetsPrunedCHS'),
    csvBTagName        = cms.untracked.string('CA8jetCombinedSecondaryVertexBJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('CA8jetCombinedSecondaryVertexBJetTagsSJCHS'),
    jettiness          = cms.untracked.string('CA8NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('CA8QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('CA8QGTaggerSubJetsCHS')
  ),

  CA15CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(150),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(True),
    topTagType           = cms.untracked.string('HEP'),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

    # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L2L3Residual_AK7PFchs.txt')
                 if is_data_flag else
                 cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                       'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt')
               ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_Uncertainty_AK7PFchs.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_Uncertainty_AK7PFchs.txt')
                  ),
    jecFilesForID = ( cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L2L3Residual_AK7PFchs.txt')
                      if is_data_flag else
                      cms.untracked.vstring('BaconProd/Utils/data/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L2Relative_AK7PFchs.txt',
                                            'BaconProd/Utils/data/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt')
                    ),
    edmRhoName = cms.untracked.string('kt6PFJets'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),

    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsCHS'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15byValAlgoCHS'),
    jetFlavorPhysName  = cms.untracked.string('CA15byValPhysCHS'),
    pruneJetName       = cms.untracked.string('CA15caPFJetsPrunedCHS'),
    subJetName         = cms.untracked.string('CA15caPFJetsPrunedCHS'),
    csvBTagName        = cms.untracked.string('CA15jetCombinedSecondaryVertexBJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('CA15jetCombinedSecondaryVertexBJetTagsSJCHS'),
    jettiness          = cms.untracked.string('CA15NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsCHS')
  ),

  PFCand = cms.untracked.PSet(
    isActive       = cms.untracked.bool(False),
    edmName        = cms.untracked.string('particleFlow'),
    edmPVName      = cms.untracked.string('offlinePrimaryVertices'),
    doAddDepthTime = cms.untracked.bool(False)
  )
)

process.baconSequence = cms.Sequence(process.PFBRECO*
                                     process.metFilters*
                                     process.producePFMETCorrections*
                                     process.recojetsequence*
                                     process.genjetsequence*
                                     process.AK5jetsequence*
                                     process.AK5genjetsequence*
                                     process.CA8jetsequenceCHS*
                                     process.CA8genjetsequenceCHS*
                                     process.CA15jetsequenceCHS*
                                     process.CA15genjetsequenceCHS*
                                     process.PFTau*   ### must come after antiktGenJets otherwise conflict on RecoJets/JetProducers/plugins
				     process.MVAMetSeq*
				     process.ntupler)

#--------------------------------------------------------------------------------
# apply trigger filter, if necessary
#================================================================================
if do_hlt_filter:
  process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
  process.hltHighLevel.throw = cms.bool(False)
  process.hltHighLevel.HLTPaths = cms.vstring()
  hlt_file = open(cmssw_base + "/src/" + hlt_filename, "r")
  for line in hlt_file.readlines():
    line = line.strip()              # strip preceding and trailing whitespaces
    if (line[0:3] == 'HLT'):         # assumes typical lines begin with HLT path name (e.g. HLT_Mu15_v1)
      hlt_path = line.split()[0]
      process.hltHighLevel.HLTPaths.extend(cms.untracked.vstring(hlt_path))
  process.p = cms.Path(process.hltHighLevel*process.baconSequence)
else:
  process.p = cms.Path(process.baconSequence)

#--------------------------------------------------------------------------------
# simple checks to catch some mistakes...
#================================================================================
if is_data_flag:
  assert process.ntupler.GenInfo.isActive == cms.untracked.bool(False)
  assert process.ntupler.Electron.isData  == cms.untracked.bool(True)
  assert process.ntupler.Muon.isData      == cms.untracked.bool(True)
  assert process.ntupler.AK5.doGenJet     == cms.untracked.bool(False)
  assert process.ntupler.CA8CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA15CHS.doGenJet == cms.untracked.bool(False)
else:
  assert process.ntupler.Electron.isData == cms.untracked.bool(False)
  assert process.ntupler.Muon.isData     == cms.untracked.bool(False)
