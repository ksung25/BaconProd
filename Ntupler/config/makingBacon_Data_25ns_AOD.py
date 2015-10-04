import FWCore.ParameterSet.Config as cms
import os

process = cms.Process('MakingBacon')

is_data_flag  = True                                      # flag for if process data
do_hlt_filter = True                                      # flag to skip events that fail relevant triggers
hlt_filename  = "BaconAna/DataFormats/data/HLTFile_25ns"  # list of relevant triggers
do_alpaca     = True

cmssw_base = os.environ['CMSSW_BASE']

#--------------------------------------------------------------------------------
# Import of standard configurations
#================================================================================
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
process.load("CommonTools/ParticleFlow/pfNoPileUpJME_cff")
process.load('BaconProd/Ntupler/myPUPPICorrections_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
if is_data_flag:
  process.GlobalTag.globaltag = cms.string('74X_dataRun2_v2')
else:
  process.GlobalTag.globaltag = cms.string('74X_mcRun2_asymptotic_v2')
#process.puppijec.connect = cms.string('sqlite:////'+cmssw_base+'/src/BaconProd/Utils/data/PY8_RunIISpring15DR74_bx50_MC.db')

#--------------------------------------------------------------------------------
# Import custom configurations
#================================================================================
# custom jet stuff (incl. GenJets, b-tagging, grooming, njettiness)
process.load('BaconProd/Ntupler/myGenJets_cff')
process.load('BaconProd/Ntupler/myJetExtrasAK4CHS_cff')
process.load('BaconProd/Ntupler/myJetExtrasAK8CHS_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA8CHS_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA15CHS_cff')

process.load('BaconProd/Ntupler/myJetExtrasAK4Puppi_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA8Puppi_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA15Puppi_cff')

# apply MET filters set to tagging mode
process.load('BaconProd/Ntupler/myMETFilters_cff') 

# MVA MET
process.load('BaconProd/Ntupler/myMVAMet_cff')     
from RecoMET.METPUSubtraction.objectSelection_AOD_cff import addLeptons
addLeptons(process)

# PF MET corrections
process.load("BaconProd/Ntupler/myPFMETCorrections_cff")
process.pfJetMETcorr.jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3Corrector")
process.producePFMETCorrections = cms.Sequence(process.producePFMETCorrectionsMC)
if is_data_flag:
  process.pfJetMETcorr.jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3ResidualCorrector")
  process.producePFMETCorrections = cms.Sequence(process.producePFMETCorrectionsData)
  process.AK4QGTaggerCHS.jec = cms.InputTag("ak4PFL1FastL2L3ResidualCorrector")
  process.CA8QGTaggerCHS.jec  = cms.InputTag("ca8PFCHSL1FastL2L3ResidualCorrector")
  process.AK8QGTaggerCHS.jec  = cms.InputTag("ak8PFCHSL1FastL2L3ResidualCorrector")
  process.CA15QGTaggerCHS.jec = cms.InputTag("ca15PFCHSL1FastL2L3ResidualCorrector")
  process.AK4QGTaggerSubJetsCHS.jec  = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector")
  process.CA8QGTaggerSubJetsCHS.jec  = cms.InputTag("ca8PFCHSL1FastL2L3ResidualCorrector")
  process.AK8QGTaggerSubJetsCHS.jec  = cms.InputTag("ak8PFCHSL1FastL2L3ResidualCorrector")
  process.CA15QGTaggerSubJetsCHS.jec = cms.InputTag("ca15PFCHSL1FastL2L3ResidualCorrector")

# produce photon isolation with proper footprint removal
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

# PUPPI
from RecoMET.METProducers.PFMET_cfi import pfMet
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.pfCandNoLep = cms.EDFilter("CandPtrSelector", src = cms.InputTag("particleFlow"), cut = cms.string("abs(pdgId) != 13 && abs(pdgId) != 11 && abs(pdgId) != 15"))
process.pfCandLep   = cms.EDFilter("CandPtrSelector", src = cms.InputTag("particleFlow"), cut = cms.string("abs(pdgId) == 13 || abs(pdgId) == 11 || abs(pdgId) == 15"))
process.puppinolep = process.puppi.clone()
process.puppinolep.candName = 'pfCandNoLep'
process.puppiForMET = cms.EDProducer("CandViewMerger",src = cms.VInputTag( 'puppinolep','pfCandLep'))     
process.pfMetPuppi = pfMet.clone();
process.pfMetPuppi.src = cms.InputTag('puppiForMET')
process.pfMetPuppi.calculateSignificance = False
process.pfJetMETcorrPuppi.jetCorrLabel = cms.InputTag("ak4PuppiL1FastL2L3Corrector")
process.producePFMETCorrectionsPuppi = cms.Sequence(process.producePFMETCorrectionsPuppiMC)
if is_data_flag:
  process.pfJetMETcorrPuppi.jetCorrLabel = cms.InputTag("ak4PuppiL1FastL2L3ResidualCorrector")
  process.producePFMETCorrectionsPuppi   = cms.Sequence(process.producePFMETCorrectionsPuppiData)
  process.AK4QGTaggerPuppi.jec           = cms.InputTag("ak4PuppiL1FastL2L3ResidualCorrector")

# ALPACA
process.load('BaconProd/Ntupler/myAlpacaCorrections_cff')
alpacaMet = ''
alpacaPuppiMet = ''
if do_alpaca: 
  alpacaMet      = ('pfMetAlpacaData'        if is_data_flag else 'pfMetAlpacaMC' )
  alpacaPuppiMet = ('pfMetPuppiAlpacaData'   if is_data_flag else 'pfMetPuppiAlpacaMC' ) 

#MET 3.0
from BaconProd.Ntupler.myMET30_cff import setMet30
setMet30(process,3.0,False)

#JEC
JECTag='Summer15_25nsV5_DATA'
if not is_data_flag: 
  JECTag='Summer15_25nsV2_MC'
ak4CHSJEC = cms.untracked.vstring('BaconProd/Utils/data/'+JECTag+'_L1FastJet_AK4PFchs.txt',
                                  'BaconProd/Utils/data/'+JECTag+'_L2Relative_AK4PFchs.txt',
                                  'BaconProd/Utils/data/'+JECTag+'_L3Absolute_AK4PFchs.txt',
                                  'BaconProd/Utils/data/'+JECTag+'_L2L3Residual_AK4PFchs.txt')

ak8CHSJEC = cms.untracked.vstring('BaconProd/Utils/data/'+JECTag+'_L1FastJet_AK8PFchs.txt',
                                  'BaconProd/Utils/data/'+JECTag+'_L2Relative_AK8PFchs.txt',
                                  'BaconProd/Utils/data/'+JECTag+'_L3Absolute_AK8PFchs.txt',
                                  'BaconProd/Utils/data/'+JECTag+'_L2L3Residual_AK8PFchs.txt')

ca15CHSJEC = ak8CHSJEC

ak4PUPPIJEC = cms.untracked.vstring('BaconProd/Utils/data/'+JECTag+'_L1FastJet_AK4PFPuppi.txt',
                                    'BaconProd/Utils/data/'+JECTag+'_L2Relative_AK4PFPuppi.txt',
                                    'BaconProd/Utils/data/'+JECTag+'_L3Absolute_AK4PFPuppi.txt',
                                    'BaconProd/Utils/data/'+JECTag+'_L2L3Residual_AK4PFPuppi.txt')

ak8PUPPIJEC = cms.untracked.vstring('BaconProd/Utils/data/'+JECTag+'_L1FastJet_AK8PFPuppi.txt',
                                    'BaconProd/Utils/data/'+JECTag+'_L2Relative_AK8PFPuppi.txt',
                                    'BaconProd/Utils/data/'+JECTag+'_L3Absolute_AK8PFPuppi.txt',
                                    'BaconProd/Utils/data/'+JECTag+'_L2L3Residual_AK8PFPuppi.txt')

ca15PUPPIJEC = ak8PUPPIJEC

ak4CHSUnc    = cms.untracked.vstring('BaconProd/Utils/data/'+JECTag+'_Uncertainty_AK4PFchs.txt')
ak8CHSUnc    = ak4CHSUnc
ca15CHSUnc   = ak4CHSUnc

ak4PUPPIUnc  = ak4CHSUnc
ak8PUPPIUnc  = ak4PUPPIUnc
ca15PUPPIUnc = ak4PUPPIUnc

#--------------------------------------------------------------------------------
# input settings
#================================================================================
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/96EF1A5F-8115-E511-AF17-02163E0125CE.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

#--------------------------------------------------------------------------------
# Reporting
#================================================================================
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

#--------------------------------------------------------------------------------
# Bacon making settings
#================================================================================
process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail     = cms.untracked.bool(do_hlt_filter),
  useAOD            = cms.untracked.bool(True),
  outputName        = cms.untracked.string('Output.root'),
  TriggerFile       = cms.untracked.string(hlt_filename),
  edmPVName         = cms.untracked.string('offlinePrimaryVertices'),
  edmGenRunInfoName = cms.untracked.string('generator'),
  
  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    edmPFCandName        = cms.untracked.string('particleFlow'),
    edmPileupInfoName    = cms.untracked.string('addPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot'),
    edmPFMETName         = cms.untracked.string('pfMet'),
    edmPFMETCorrName     = cms.untracked.string('pfType1CorrectedMet'),
    edmMVAMETName        = cms.untracked.string('pfMVAMEt'),
    edmPuppETName        = cms.untracked.string('pfMetPuppi'),
    edmPuppETCorrName    = cms.untracked.string('pfType1PuppiCorrectedMet'),
    edmPFMET30Name       = cms.untracked.string('pfMet30'),
    edmPFMETC30Name      = cms.untracked.string('pfType1CorrectedMet30'),
    edmMVAMET30Name      = cms.untracked.string('pfMVAMEt30'),
    edmPuppET30Name      = cms.untracked.string('pfMetPuppi30'),
    edmPuppET30CorrName  = cms.untracked.string('pfType1PuppiCorrectedMet30'),
    edmAlpacaMETName     = cms.untracked.string(alpacaMet),
    edmPupAlpacaMETName  = cms.untracked.string(alpacaPuppiMet),
    edmRhoForIsoName     = cms.untracked.string('fixedGridRhoFastjetAll'),
    edmRhoForJetEnergy   = cms.untracked.string('fixedGridRhoFastjetAll'),
    doFillMETFilters     = cms.untracked.bool(False),
    doFillMET            = cms.untracked.bool(True)
  ),
  
  GenInfo = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmGenEventInfoName = cms.untracked.string('generator'),
    edmGenParticlesName = cms.untracked.string('genParticles'),
    fillAllGen          = cms.untracked.bool(False),
    fillLHEWeights      = cms.untracked.bool(False)
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
    minPt                     = cms.untracked.double(7),
    edmName                   = cms.untracked.string('gedGsfElectrons'),
    edmBeamspotName           = cms.untracked.string('offlineBeamSpot'),
    edmPFCandName             = cms.untracked.string('particleFlow'),
    edmTrackName              = cms.untracked.string('generalTracks'),
    edmConversionName         = cms.untracked.string('allConversions'),
    edmSuperClusterName       = cms.untracked.string('particleFlowEGamma'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppinolep'),
    usePuppi                  = cms.untracked.bool(True)
  ),
  
  Muon = cms.untracked.PSet(
    isActive                  = cms.untracked.bool(True),
    minPt                     = cms.untracked.double(3),
    edmName                   = cms.untracked.string('muons'),
    edmPFCandName             = cms.untracked.string('particleFlow'),
    # save general tracker tracks in our muon collection (used in tag-and-probe for muons)
    doSaveTracks              = cms.untracked.bool(False),
    minTrackPt                = cms.untracked.double(20),
    edmTrackName              = cms.untracked.string('generalTracks'),
    #puppi
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppinolep'),
    usePuppi                  = cms.untracked.bool(True)    
  ),
  
  Photon = cms.untracked.PSet(
    isActive              = cms.untracked.bool(True),
    minPt                 = cms.untracked.double(10),
    edmName               = cms.untracked.string('gedPhotons'),
    edmPFCandName         = cms.untracked.string('particleFlow'),
    edmElectronName       = cms.untracked.string('gedGsfElectrons'),
    edmConversionName     = cms.untracked.string('allConversions'),
    edmSuperClusterName   = cms.untracked.string('particleFlowEGamma'),
    edmChHadIsoMapTag     = cms.untracked.InputTag("photonIDValueMapProducer:phoChargedIsolation"),        # EGM recommendation not in AOD/MINIAOD
    edmNeuHadIsoMapTag    = cms.untracked.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),  # EGM recommendation not in AOD/MINIAOD
    edmGammaIsoMapTag     = cms.untracked.InputTag("photonIDValueMapProducer:phoPhotonIsolation")          # EGM recommendation not in AOD/MINIAOD
  ),
  
  Tau = cms.untracked.PSet(
    isActive = cms.untracked.bool(True),
    minPt    = cms.untracked.double(10),
    edmName  = cms.untracked.string('hpsPFTauProducer'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppinolep'),
    usePuppi                  = cms.untracked.bool(True),
  ),
  
  AK4CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    coneSize             = cms.untracked.double(0.4),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    
    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ak4CHSJEC, 
    jecUncFiles = ak4CHSUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    # names of various jet-related collections
    jetName            = cms.untracked.string('ak4PFJetsCHS'),
    genJetName         = cms.untracked.string('AK4GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK4FlavorCHS'),
    prunedJetName      = cms.untracked.string('AK4caPFJetsPrunedCHS'),
    trimmedJetName     = cms.untracked.string('AK4caPFJetsTrimmedCHS'),
    softdropJetName    = cms.untracked.string('AK4caPFJetsSoftDropCHS'),
    subJetName         = cms.untracked.string('AK4caPFJetsSoftDropCHS'),
    csvBTagName        = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    csvBTagSubJetName  = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    jettiness          = cms.untracked.string('AK4NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('AK4QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('AK4QGTaggerSubJetsCHS'),
    topTaggerName      = cms.untracked.string('')
  ),

  AK4Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    coneSize             = cms.untracked.double(0.4),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ak4PUPPIJEC,
    jecUncFiles = ak4PUPPIUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    
    # names of various jet-related collections
    jetName            = cms.untracked.string('AK4PFJetsPuppi'),
    genJetName         = cms.untracked.string('AK4GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK4FlavorPuppi'),
    prunedJetName      = cms.untracked.string('AK4caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('AK4caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    jettiness          = cms.untracked.string('AK4NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('AK4QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('AK4QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('')
  ),

  AK8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ak8CHSJEC,
    jecUncFiles = ak8CHSUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),

    # names of various jet-related collections
    jetName            = cms.untracked.string('AK8PFJetsCHS'),
    genJetName         = cms.untracked.string('AK8GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK8FlavorCHS'),
    prunedJetName      = cms.untracked.string('AK8caPFJetsPrunedCHS'),
    trimmedJetName     = cms.untracked.string('AK8caPFJetsTrimmedCHS'),
    softdropJetName    = cms.untracked.string('AK8caPFJetsSoftDropCHS'),
    subJetName         = cms.untracked.string('AK8caPFJetsSoftDropCHS'),
    csvBTagName        = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    jettiness          = cms.untracked.string('AK8NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('AK8QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('AK8QGTaggerSubJetsCHS'),
    topTaggerName      = cms.untracked.string('CMS')
  ),

  CA8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(False),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    
    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ak8CHSJEC,
    jecUncFiles = ak8CHSUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    
    # names of various jet-related collections
    jetName            = cms.untracked.string('CA8PFJetsCHS'),
    genJetName         = cms.untracked.string('CA8GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA8FlavorCHS'),
    prunedJetName      = cms.untracked.string('CA8caPFJetsPrunedCHS'),
    trimmedJetName     = cms.untracked.string('CA8caPFJetsTrimmedCHS'),
    softdropJetName    = cms.untracked.string('CA8caPFJetsSoftDropCHS'),
    subJetName         = cms.untracked.string('CA8caPFJetsSoftDropCHS'),
    csvBTagName        = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    jettiness          = cms.untracked.string('CA8NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('CA8QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('CA8QGTaggerSubJetsCHS'),
    topTaggerName      = cms.untracked.string('CMS')
  ),

  CA8Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    
    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ak8PUPPIJEC,
    jecUncFiles = ak8PUPPIUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    
    # names of various jet-related collections
    jetName            = cms.untracked.string('CA8PFJetsPuppi'),
    genJetName         = cms.untracked.string('CA8GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA8FlavorPuppi'),
    prunedJetName      = cms.untracked.string('CA8caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('CA8caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('CA8caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('CA8caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    jettiness          = cms.untracked.string('CA8NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('CA8QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('CA8QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('CMS')
  ),

  CA15CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ca15CHSJEC,
    jecUncFiles = ca15CHSUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsCHS'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15FlavorCHS'),
    prunedJetName      = cms.untracked.string('CA15caPFJetsPrunedCHS'),
    trimmedJetName     = cms.untracked.string('CA15caPFJetsTrimmedCHS'),
    softdropJetName    = cms.untracked.string('CA15caPFJetsSoftDropCHS'),
    subJetName         = cms.untracked.string('CA15caPFJetsSoftDropCHS'),
    csvBTagName        = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    jettiness          = cms.untracked.string('CA15NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsCHS'),
    topTaggerName      = cms.untracked.string('HEP')
  ),
  CA15Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmPVName   = cms.untracked.string('offlinePrimaryVertices'),
    jecFiles    = ca15PUPPIJEC,
    jecUncFiles = ca15PUPPIUnc,
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsPuppi'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15FlavorPuppi'),
    prunedJetName      = cms.untracked.string('CA15caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('CA15caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('CA15caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('CA15caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    jettiness          = cms.untracked.string('CA15NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('HEP')
  ),
  
  PFCand = cms.untracked.PSet(
    isActive       = cms.untracked.bool(False),
    edmName        = cms.untracked.string('particleFlow'),
    edmPVName      = cms.untracked.string('offlinePrimaryVertices'),
    doAddDepthTime = cms.untracked.bool(False)
  )
)

process.baconSequence = cms.Sequence(#process.puppi*
                                     process.pfNoPileUpJMESequence*
                                     process.packedPFCandidates30*
                                     process.metFilters*
                                     process.producePFMETCorrections*
                                     process.egmGsfElectronIDSequence* 
                                     process.egmPhotonIDSequence*
                                     process.MVAMetSeqData*
                                     process.ak4PFJets30*              #30
                                     process.pfMVAMEtSequenceNoLep30*  #30 Lepton ids done above
                                     process.pfCandNoLep*
                                     process.pfCandLep*
                                     process.puppi*
                                     process.puppinolep*
                                     process.puppiForMET*
                                     process.puppi30*           # 30 
                                     process.puppinolep30*      # 30
                                     process.ak4PuppiL1FastL2L3ResidualChain*
                                     process.AK4PFJetsPuppi30*  # 30
                                     #process.alpacaSequenceData*
                                     process.pfMetPuppi*
                                     process.AK4jetsequenceCHSData*
                                     process.AK4jetsequencePuppiData*
                                     process.producePFMETCorrectionsPuppi*
                                     process.AK8jetsequenceCHSData*
                                     process.CA8jetsequenceCHSData*
                                     process.CA15jetsequenceCHSData*
                                     process.CA8jetsequencePuppiData*
                                     process.CA15jetsequencePuppiData*
                                     process.allMET30*         # 30
                                     process.photonIDValueMapProducer*
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
  assert process.ntupler.AK4CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA8CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA15CHS.doGenJet == cms.untracked.bool(False)


#process.out = cms.OutputModule("PoolOutputModule",                                                                                                                                                   
#                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                      
#                                  fileName       = cms.untracked.string ("test.root")                                                                                                                    
#                                  )   
#
#process.endpath = cms.EndPath(process.out)
