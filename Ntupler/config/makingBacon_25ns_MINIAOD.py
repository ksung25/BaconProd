# combined configuration for all eras for data and MC

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as vp
import os

process = cms.Process('MakingBacon')

options = vp.VarParsing('analysis')

# register command line options
options.register('isData',
                 True, # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )

options.register('doHLTFilter',
                 True, # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )

options.register('doAlpaca',
                 False, # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )

options.register('era',
                 '2016', # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.string, 
                 )

# user input values
options.parseArguments()

hlt_filename  = "BaconAna/DataFormats/data/HLTFile_25ns"   # list of relevant triggers

cmssw_base = os.environ['CMSSW_BASE']
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if options.era == '2016':
    if options.isData:
      process.GlobalTag.globaltag = cms.string('94X_dataRun2_v10')
      JECTag='Summer16_23Sep2016AllV4_DATA'
    else:
      process.GlobalTag.globaltag = cms.string('94X_mcRun2_asymptotic_v3')
      JECTag='Summer16_23Sep2016V4_MC'

elif options.era == '2017':
    if options.isData:
      process.GlobalTag.globaltag = cms.string('101X_dataRun2_Prompt_v11')
      JECTag='Fall17_17Nov2017_V32_102X_DATA'
    else:
      process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v14')
      JECTag='Fall17_17Nov2017_V32_102X_MC'

elif options.era == '2018':
    if options.isData:
      process.GlobalTag.globaltag = cms.string('101X_dataRun2_Prompt_v11')
      JECTag='Fall17_17Nov2017_V32_102X_DATA'
    else:
      process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v15')
      JECTag='Fall17_17Nov2017_V32_102X_MC'

from BaconProd.Ntupler.myJecFromDB_cff    import setupJEC
setupJEC(process,options.isData,JECTag)
if options.isData:
  process.jec.connect = cms.string('sqlite:////'+cmssw_base+'/src/BaconProd/Utils/data/'+JECTag+'.db')
else:
  process.jec.connect = cms.string('sqlite:////'+cmssw_base+'/src/BaconProd/Utils/data/'+JECTag+'.db')

#--------------------------------------------------------------------------------
# Import of standard configurations
#================================================================================
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

if options.era == '2017' or options.era == '2018':
    process.load('RecoBTag.SoftLepton.SoftLeptonByMVAComputers_cff')

process.pfNoPileUpJME = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
process.load('BaconProd/Ntupler/myPUPPICorrections_cff')
process.load('BaconProd/Ntupler/myCHSCorrections_cff')
process.load('BaconProd/Ntupler/myCorrections_cff')

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
process.load('BaconProd/Ntupler/myJetExtrasAK8Puppi_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA15Puppi_cff')

process.load("RecoBTag.ImpactParameter.impactParameter_cff")
process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
process.load("RecoBTag.SoftLepton.softLepton_cff")
process.load("RecoBTag.Combined.combinedMVA_cff")
process.load("RecoBTag.CTagging.cTagging_cff")
process.load("RecoBTag.Combined.deepFlavour_cff")
if options.era == '2016':
    process.pfDeepCSVJetTags.NNConfig = cms.FileInPath('BaconProd/Utils/data/DeepFlavourNoSL.json')

from BaconProd.Ntupler.myBtagging_cff           import addBTagging,addBTaggingAK4CHS
from BaconProd.Ntupler.myGenJets_cff            import setMiniAODGenJets
from BaconProd.Ntupler.myJetExtrasAK4CHS_cff    import setMiniAODAK4CHS
from BaconProd.Ntupler.myJetExtrasAK8CHS_cff    import setMiniAODAK8CHS
from BaconProd.Ntupler.myJetExtrasCA15CHS_cff   import setMiniAODCA15CHS

from BaconProd.Ntupler.myJetExtrasAK4Puppi_cff  import setMiniAODAK4Puppi
from BaconProd.Ntupler.myJetExtrasAK8Puppi_cff  import setMiniAODAK8Puppi
from BaconProd.Ntupler.myJetExtrasCA15Puppi_cff import setMiniAODCA15Puppi

process.btagging = cms.Sequence()
addBTagging(process,'AK4PFJetsPuppi' ,0.4,'AK4' ,'Puppi')
addBTagging(process,'AK8PFJetsPuppi' ,0.8,'AK8' ,'Puppi')
addBTagging(process,'CA15PFJetsPuppi',1.5,'CA15','Puppi')
addBTaggingAK4CHS(process,'updatedPatJets'    ,0.4,'AK4' ,'CHS',True,True)

setMiniAODGenJets(process)
setMiniAODAK4CHS(process)

setMiniAODAK4Puppi (process)
setMiniAODAK8Puppi (process)
setMiniAODCA15Puppi(process)

#METFilters
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#CHS
process.chs = cms.EDFilter("CandPtrSelector",
                           src = cms.InputTag('packedPFCandidates'),
                           cut = cms.string('fromPV')
                           )
if options.isData:
  process.AK4QGTaggerCHS.jec  = cms.InputTag("ak4chsL1FastL2L3ResidualCorrector")
  process.CA8QGTaggerCHS.jec  = cms.InputTag("ca8chsL1FastL2L3ResidualCorrector")
  process.AK8QGTaggerCHS.jec  = cms.InputTag("ak8chsL1FastL2L3ResidualCorrector")
  process.CA15QGTaggerCHS.jec = cms.InputTag("ak8chsL1FastL2L3ResidualCorrector")
  process.AK4QGTaggerSubJetsCHS.jec  = cms.InputTag("ak4chsL1FastL2L3ResidualCorrector")
  process.CA8QGTaggerSubJetsCHS.jec  = cms.InputTag("ca8chsL1FastL2L3ResidualCorrector")
  process.AK8QGTaggerSubJetsCHS.jec  = cms.InputTag("ak8chsL1FastL2L3ResidualCorrector")
  process.CA15QGTaggerSubJetsCHS.jec = cms.InputTag("ak8chsL1FastL2L3ResidualCorrector")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if options.era == '2016':
    setupEgammaPostRecoSeq(process,
                           runVID=True,
                           era='2016-Legacy', #era is new to select between 2016 / 2017,  it defaults to 2017
                           phoIDModules=[]) #bug with default modules for photon VID; off for now
elif options.era == '2017':
    setupEgammaPostRecoSeq(process,
                           runVID=True,
                           era='2017-Nov17ReReco', #era is new to select between 2016 / 2017,  it defaults to 2017
                           phoIDModules=[]) #bug with default modules for photon VID; off for now

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
if options.era == '2016':
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
        DataEra = cms.string("2016BtoH"), 
        UseJetEMPt = cms.bool(False),
        PrefiringRateSystematicUncty = cms.double(0.2),
        SkipWarnings = False)
elif options.era == '2017':
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
        DataEra = cms.string("2017BtoF"), 
        UseJetEMPt = cms.bool(False),
        PrefiringRateSystematicUncty = cms.double(0.2),
        SkipWarnings = False)

# PF MET corrections
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=options.isData,
                           manualJetConfig=True,
                           jetCorLabelL3="ak4chsL1FastL2L3Corrector",
                           jetCorLabelRes="ak4chsL1FastL2L3ResidualCorrector",
                           reclusterJets=True,
                           recoMetFromPFCs=True,
                           postfix="V2"
                           )

# PUPPI Woof Woof
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD (process, True )
process.puppi.useExistingWeights = True
runMetCorAndUncFromMiniAOD(process,
                           isData=options.isData,
                           manualJetConfig=True,
                           metType="Puppi",
                           pfCandColl=cms.InputTag("puppiForMET"),
                           recoMetFromPFCs=True,
                           jetFlavor="AK4PFPuppi",
                           jetCorLabelL3="ak4PuppiL1FastL2L3Corrector",
                           jetCorLabelRes="ak4PuppiL1FastL2L3ResidualCorrector",
                           reclusterJets=True,
                           postfix="Puppi"
                           )

if options.isData:
  process.AK4QGTaggerPuppi.jec           = cms.InputTag("ak4PuppiL1FastL2L3ResidualCorrector")
  process.AK8QGTaggerPuppi.jec           = cms.InputTag("ak8PuppiL1FastL2L3ResidualCorrector")
  process.CA15QGTaggerPuppi.jec          = cms.InputTag("ak8PuppiL1FastL2L3ResidualCorrector")
  process.AK4QGTaggerSubJetsPuppi.jec    = cms.InputTag("ak4PuppiL1FastL2L3ResidualCorrector")
  process.AK8QGTaggerSubJetsPuppi.jec    = cms.InputTag("ak8PuppiL1FastL2L3ResidualCorrector")
  process.CA15QGTaggerSubJetsPuppi.jec   = cms.InputTag("ak8PuppiL1FastL2L3ResidualCorrector")

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
  btagDiscriminators = ['pfDeepCMVAJetTags:probb','pfDeepCMVAJetTags:probc','pfDeepCMVAJetTags:probudsg','pfDeepCMVAJetTags:probbb',
                        'pfDeepCSVJetTags:probb' ,'pfDeepCSVJetTags:probc' ,'pfDeepCSVJetTags:probudsg' ,'pfDeepCSVJetTags:probbb']
  )

# ALPACA
alpacaMet = ''
alpacaPuppiMet = ''
if options.doAlpaca: 
  alpacaMet      = ('pfMetAlpacaData'        if options.isData else 'pfMetAlpacaMC' )
  alpacaPuppiMet = ('pfMetPuppiAlpacaData'   if options.isData else 'pfMetPuppiAlpacaMC' ) 

#--------------------------------------------------------------------------------
# input settings
#================================================================================
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

if options.era == '2016':
    if options.isData:
        test_file = cms.untracked.vstring('/store/data/Run2016C/DoubleEG/MINIAOD/03Feb2017-v1/810000/D8A12591-5DED-E611-BAF8-02163E019C24.root')
    else:
        test_file = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/GluGluHToZG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/702BC677-2AC8-E611-A5B3-02163E011463.root')
elif options.era == '2017':
    if options.isData:
        test_file = cms.untracked.vstring('/store/data/Run2017F/DoubleEG/MINIAOD/31Mar2018-v1/90001/F61EA338-8E37-E811-A203-0025905C4262.root')
    else:
        #test_file = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/46FB5EDE-F708-E811-A50F-0025905C53A4.root')
        test_file = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/005DC030-D3F4-E711-889A-02163E01A62D.root')

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/GluGluHToZG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/702BC677-2AC8-E611-A5B3-02163E011463.root')
                            fileNames = test_file
                            )
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

#--------------------------------------------------------------------------------
# Reporting
#================================================================================
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE'),
)

#--------------------------------------------------------------------------------
# Bacon making settings
#================================================================================
process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail     = cms.untracked.bool(options.doHLTFilter),
  useTrigger        = cms.untracked.bool(True),
  useTriggerObject  = cms.untracked.bool(True),
  TriggerObject     = cms.untracked.string("selectedPatTrigger"),
  TriggerFile       = cms.untracked.string(hlt_filename),
  useAOD            = cms.untracked.bool(False),
  outputName        = cms.untracked.string('Output.root'),
  edmPVName         = cms.untracked.string('offlineSlimmedPrimaryVertices'),
  edmGenRunInfoName = cms.untracked.string('generator'),
  
  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    edmPFCandName        = cms.untracked.string('packedPFCandidates'),
    edmPileupInfoName    = cms.untracked.string('slimmedAddPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot'),
    edmMETName           = cms.untracked.string('slimmedMETs'),
    edmPFMETName         = cms.untracked.InputTag('slimmedMETsV2'),
    edmMVAMETName        = cms.untracked.string(''),
    edmPuppETName        = cms.untracked.InputTag('slimmedMETsPuppi'),
    edmAlpacaMETName     = cms.untracked.string(alpacaMet),
    edmPupAlpacaMETName  = cms.untracked.string(alpacaPuppiMet),
    edmRhoForIsoName     = cms.untracked.string('fixedGridRhoFastjetAll'),
    edmRhoForJetEnergy   = cms.untracked.string('fixedGridRhoFastjetAll'),
    doFillMETFilters     = cms.untracked.bool(True),
    doFillMET            = cms.untracked.bool(True)
  ),
  
  GenInfo = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    edmGenEventInfoName     = cms.untracked.string('generator'),
    edmGenParticlesName     = cms.untracked.string('prunedGenParticles'),
    edmGenPackParticlesName = cms.untracked.string('packedGenParticles'),
    fillAllGen          = cms.untracked.bool(False),
    fillLHEWeights      = cms.untracked.bool(True)
  ),

  GenJet  = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False)),
    isActiveFatJet      = ( cms.untracked.bool(False)),
    #isActive            = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    #isActiveFatJet      = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    edmGenParticlesName = cms.untracked.string('prunedGenParticles'),
    genJetName          = cms.untracked.string('AK4GenJetsCHS'),
    genFatJetName       = cms.untracked.string('AK8GenJetsCHS'),
    fillAllGen          = cms.untracked.bool(False)
  ),
                                   
  PV = cms.untracked.PSet(
    isActive      = cms.untracked.bool(True),   
    edmName       = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    minNTracksFit = cms.untracked.uint32(0),
    minNdof       = cms.untracked.double(4),
    maxAbsZ       = cms.untracked.double(24),
    maxRho        = cms.untracked.double(2)
  ),
  
  Electron = cms.untracked.PSet(
    isActive                  = cms.untracked.bool(True),
    minPt                     = cms.untracked.double(7),
    edmName                   = cms.untracked.string('slimmedElectrons'),
    edmSCName                 = cms.untracked.InputTag('reducedEgamma','reducedSuperClusters'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppiNoLep'),
    usePuppi                  = cms.untracked.bool(True),
    useTriggerObject          = cms.untracked.bool(False),
    edmEcalPFClusterIsoMapTag = cms.untracked.InputTag('electronEcalPFClusterIsolationProducer'),
    edmHcalPFClusterIsoMapTag = cms.untracked.InputTag('electronHcalPFClusterIsolationProducer'),

    edmEleMVASpring16         = cms.untracked.string('ElectronMVAEstimatorRun2Spring16GeneralPurposeV1'),
    edmEleMVAFall17V1Iso        = cms.untracked.string('ElectronMVAEstimatorRun2Fall17IsoV1'),
    edmEleMVAFall17V1NoIso      = cms.untracked.string('ElectronMVAEstimatorRun2Fall17NoIsoV1'),
    edmEleMVAFall17V2Iso        = cms.untracked.string('ElectronMVAEstimatorRun2Fall17IsoV2'),
    edmEleMVAFall17V2NoIso      = cms.untracked.string('ElectronMVAEstimatorRun2Fall17NoIsoV2'),
    edmEleMVASpring16HZZ      = cms.untracked.string('ElectronMVAEstimatorRun2Spring16HZZV1'),

    edmEleMediumMVA           = cms.untracked.string('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
    edmEleTightMVA            = cms.untracked.string('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),

    #edmEleMVA                 = cms.untracked.string('ElectronMVAEstimatorRun2Spring16GeneralPurposeV1'),
    #edmEleMediumMVAIso        = cms.untracked.string(''), 
    #edmEleTightMVAIso         = cms.untracked.string(''),
    #edmEleMVAIso              = cms.untracked.string(''),
    #storeSecondMVA            = cms.untracked.bool(False),
    #storeHZZMVA               = cms.untracked.bool(True),
  ),
  
  Muon = cms.untracked.PSet(
    isActive                  = cms.untracked.bool(True),
    minPt                     = cms.untracked.double(3),
    edmName                   = cms.untracked.string('slimmedMuons'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppiNoLep'),
    usePuppi                  = cms.untracked.bool(True),
    useTriggerObject          = cms.untracked.bool(False),    
  ),
  
  Photon = cms.untracked.PSet(
    isActive              = cms.untracked.bool(True),
    minPt                 = cms.untracked.double(10),
    edmName               = cms.untracked.string('slimmedPhotons'),
    edmSCName             = cms.untracked.InputTag('reducedEgamma','reducedSuperClusters'),
    edmPhoMVASpring16     = cms.untracked.string('PhotonMVAEstimatorRun2Spring16NonTrigV1'),
    edmPhoMVAFall17V1     = cms.untracked.string('PhotonMVAEstimatorRunIIFall17v1'),
    edmPhoMVAFall17V2     = cms.untracked.string('PhotonMVAEstimatorRunIIFall17v2'),
    useTriggerObject      = cms.untracked.bool(False),
  ),
  
  Tau = cms.untracked.PSet(
    isActive = cms.untracked.bool(True),
    minPt    = cms.untracked.double(10),
    edmName  = cms.untracked.string('slimmedTaus'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppiNoLep'),
    usePuppi                  = cms.untracked.bool(True),
    useTriggerObject          = cms.untracked.bool(False),
  ),
  
  AK4CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(False),
    useTriggerObject     = cms.untracked.bool(False),
    minPt                = cms.untracked.double(15),
    coneSize             = cms.untracked.double(0.4),
    addPFCand            = cms.untracked.bool(False),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    showerDecoConf       = cms.untracked.string(''),
    
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak4chsL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak4chsL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK4chs')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),

    # names of various jet-related collections
    jetName              = cms.untracked.string('selectedUpdatedPatJets'),
    genJetName           = cms.untracked.string('slimmedGenJets'),
    csvBTagName          = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    mvaBTagName          = cms.untracked.string('pfCombinedMVAV2BJetTags'),
    cvlcTagName          = cms.untracked.string('pfCombinedCvsLJetTags'),
    cvbcTagName          = cms.untracked.string('pfCombinedCvsBJetTags'),
    qgLikelihood         = cms.untracked.string('QGTagger'),
    deepCSVBTagName      = cms.untracked.string('pfDeepCSVJetTags'),
    deepCMVABTagName     = cms.untracked.string('pfDeepCMVAJetTags'),
    ),

  AK4Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    useTriggerObject     = cms.untracked.bool(False),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    coneSize             = cms.untracked.double(0.4),
    addPFCand            = cms.untracked.bool(True),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    showerDecoConf       = cms.untracked.string(''),
    
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak4PuppiL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak4PuppiL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK4Puppi')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('', 'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    jetBoostedBtaggingFiles = cms.untracked.string('BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),
    # names of various jet-related collections
    jetName            = cms.untracked.string('AK4PFJetsPuppi'),
    genJetName         = cms.untracked.string('AK4GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK4FlavorPuppi'),
    prunedJetName      = cms.untracked.string('AK4caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('AK4caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    mvaBTagName        = cms.untracked.string('AK4PFCombinedMVAV2BJetTagsPuppi'),
    cvlcTagName        = cms.untracked.string('AK4PFCombinedCvsLJetTagsPuppi'),
    cvbcTagName        = cms.untracked.string('AK4PFCombinedCvsBJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    csvDoubleBTagName  = cms.untracked.string('AK4PFBoostedDoubleSecondaryVertexBJetTagsPuppi'),
    deepCSVBTagName    = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi'),
    deepCSVBTagNameb   = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probb'),
    deepCSVBTagNamec   = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probc'),
    deepCSVBTagNamel   = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probudsg'),
    deepCSVBTagNamebb  = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probbb'),
    deepCMVABTagName   = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi'),
    deepCMVABTagNameb  = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probb'),
    deepCMVABTagNamec  = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probc'),
    deepCMVABTagNamel  = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probudsg'),
    deepCMVABTagNamebb = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('AK4PFBoostedDeepDoubleBvLJetTagsPuppi:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('AK4PFBoostedDeepDoubleCvLJetTagsPuppi:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('AK4PFBoostedDeepDoubleCvBJetTagsPuppi:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('AK4PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsPuppi:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('AK4PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsPuppi:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('AK4PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsPuppi:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('AK4PFBoostedDoubleSVTagInfosPuppi'), 
    secVertices        = cms.untracked.string('slimmedSecondaryVertices'),
    edmMuonName        = cms.untracked.string('slimmedMuons'),
    edmElectronName    = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('AK4PFSoftPFMuonsTagInfosPuppi'),
    softPFElectronTagInfoName = cms.untracked.string('AK4PFSoftPFElectronsTagInfosPuppi'),
    jettiness          = cms.untracked.string('AK4NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('AK4QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('AK4QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('')
  ),

  AK8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(False),
    useAOD               = cms.untracked.bool(True),
    useTriggerObject     = cms.untracked.bool(False),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    showerDecoConf       = cms.untracked.string(''),
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    jetBoostedBtaggingFiles = cms.untracked.string('BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),

    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak8chsL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak8chsL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK8chs')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),

    # names of various jet-related collections
    jetName              = cms.untracked.string('AK8PFJetsCHS'),
    genJetName           = cms.untracked.string('AK8GenJetsCHS'),
    jetFlavorName        = cms.untracked.string('AK8FlavorCHS'),
    prunedJetName        = cms.untracked.string('AK8caPFJetsPrunedCHS'),
    trimmedJetName       = cms.untracked.string('AK8caPFJetsTrimmedCHS'),
    softdropJetName      = cms.untracked.string('AK8caPFJetsSoftDropCHS'),
    subJetName           = cms.untracked.string('AK8caPFJetsSoftDropCHS'),
    csvBTagName          = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS'),
    mvaBTagName          = cms.untracked.string('AK8PFCombinedMVAV2BJetTagsCHS'),
    cvlcTagName          = cms.untracked.string('AK8PFCombinedCvsLJetTagsCHS'),
    cvbcTagName          = cms.untracked.string('AK8PFCombinedCvsBJetTagsCHS'),
    csvBTagSubJetName    = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    csvDoubleBTagName    = cms.untracked.string('AK8PFBoostedDoubleSecondaryVertexBJetTagsCHS'),
    deepCSVBTagName      = cms.untracked.string('AK8PFDeepCSVJetTagsCHS'),
    deepCSVBTagNameb     = cms.untracked.string('AK8PFDeepCSVJetTagsCHS:probb'),
    deepCSVBTagNamec     = cms.untracked.string('AK8PFDeepCSVJetTagsCHS:probc'),
    deepCSVBTagNamel     = cms.untracked.string('AK8PFDeepCSVJetTagsCHS:probudsg'),
    deepCSVBTagNamebb    = cms.untracked.string('AK8PFDeepCSVJetTagsCHS:probbb'),
    deepCMVABTagName     = cms.untracked.string('AK8PFDeepCMVAJetTagsCHS'),
    deepCMVABTagNameb    = cms.untracked.string('AK8PFDeepCMVAJetTagsCHS:probb'),
    deepCMVABTagNamec    = cms.untracked.string('AK8PFDeepCMVAJetTagsCHS:probc'),
    deepCMVABTagNamel    = cms.untracked.string('AK8PFDeepCMVAJetTagsCHS:probudsg'),
    deepCMVABTagNamebb   = cms.untracked.string('AK8PFDeepCMVAJetTagsCHS:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('AK8PFBoostedDeepDoubleBvLJetTagsCHS:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('AK8PFBoostedDeepDoubleCvLJetTagsCHS:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('AK8PFBoostedDeepDoubleCvBJetTagsCHS:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('AK8PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsCHS:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('AK8PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsCHS:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('AK8PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsCHS:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('AK8PFBoostedDoubleSVTagInfosCHS'),
    jettiness            = cms.untracked.string('AK8NjettinessCHS'),
    edmMuonName          = cms.untracked.string('slimmedMuons'),
    edmElectronName      = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('AK8PFSoftPFMuonsTagInfosCHS'),
    softPFElectronTagInfoName = cms.untracked.string('AK8PFSoftPFElectronsTagInfosCHS'),
    qgLikelihood         = cms.untracked.string('AK8QGTaggerCHS'),
    qgLikelihoodSubjet   = cms.untracked.string('AK8QGTaggerSubJetsCHS'),
    topTaggerName        = cms.untracked.string('')
  ),

  CA8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(False),
    useAOD               = cms.untracked.bool(False),
    useTriggerObject     = cms.untracked.bool(False),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
        
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak8chsL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak8chsL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK8chs')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    showerDecoConf       = cms.untracked.string(''),
    # names of various jet-related collections
    jetName              = cms.untracked.string('CA8PFJetsCHS'),
    genJetName           = cms.untracked.string('CA8GenJetsCHS'),
    subJetName           = cms.untracked.string('SoftDrop'),
    csvBTagName          = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTags'),
    csvDoubleBTagName    = cms.untracked.string('CA8PFBoostedDoubleSecondaryVertexBJetTagsCHS'),
    deepCSVBTagName      = cms.untracked.string('CA8PFDeepCSVJetTagsCHS'),
    deepCSVBTagNameb     = cms.untracked.string('CA8PFDeepCSVJetTagsCHS:probb'),
    deepCSVBTagNamec     = cms.untracked.string('CA8PFDeepCSVJetTagsCHS:probc'),
    deepCSVBTagNamel     = cms.untracked.string('CA8PFDeepCSVJetTagsCHS:probudsg'),
    deepCSVBTagNamebb    = cms.untracked.string('CA8PFDeepCSVJetTagsCHS:probbb'),
    deepCMVABTagName     = cms.untracked.string('CA8PFDeepCMVAJetTagsCHS'),
    deepCMVABTagNameb    = cms.untracked.string('CA8PFDeepCMVAJetTagsCHS:probb'),
    deepCMVABTagNamec    = cms.untracked.string('CA8PFDeepCMVAJetTagsCHS:probc'),
    deepCMVABTagNamel    = cms.untracked.string('CA8PFDeepCMVAJetTagsCHS:probudsg'),
    deepCMVABTagNamebb   = cms.untracked.string('CA8PFDeepCMVAJetTagsCHS:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('CA8PFBoostedDeepDoubleBvLJetTagsCHS:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('CA8PFBoostedDeepDoubleCvLJetTagsCHS:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('CA8PFBoostedDeepDoubleCvBJetTagsCHS:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('CA8PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsCHS:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('CA8PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsCHS:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('CA8PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsCHS:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('CA8PFBoostedDoubleSVTagInfosCHS'),
    svTagInfoName        = cms.untracked.string('CA8PFSecondaryVertexTagInfosCHS'),
    edmMuonName          = cms.untracked.string('slimmedMuons'),
    edmElectronName      = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('CA8PFSoftPFMuonsTagInfosCHS'),
    softPFElectronTagInfoName = cms.untracked.string('CA8PFSoftPFElectronsTagInfosCHS'),
    qgLikelihood         = cms.untracked.string('CA8QGTaggerCHS'),
    prunedJetName        = cms.untracked.string('CA8PFJetsCHSPruned'),
    trimmedJetName       = cms.untracked.string('CA8PFJetsCHSTrimmed'),
    softdropJetName      = cms.untracked.string('CA8PFJetsCHSSoftDrop'),
    jettiness            = cms.untracked.string('CA8NjettinessCHS'),
    topTaggerName        = cms.untracked.string('CMS')
  ),
                                 
  AK8Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    useTriggerObject     = cms.untracked.bool(False),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doComputeSVInfo      = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak8PuppiL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak8PuppiL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK8Puppi')),
    showerDecoConf       = cms.untracked.string(''),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    jetBoostedBtaggingFiles = cms.untracked.string('BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),

    # names of various jet-related collections
    jetName            = cms.untracked.string('AK8PFJetsPuppi'),
    genJetName         = cms.untracked.string('AK8GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK8FlavorPuppi'),
    prunedJetName      = cms.untracked.string('AK8caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('AK8caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('AK8caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('AK8caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    mvaBTagName        = cms.untracked.string('AK8PFCombinedMVAV2BJetTagsPuppi'),
    cvlcTagName        = cms.untracked.string('AK8PFCombinedCvsLJetTagsPuppi'),
    cvbcTagName        = cms.untracked.string('AK8PFCombinedCvsBJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    csvDoubleBTagName  = cms.untracked.string('AK8PFBoostedDoubleSecondaryVertexBJetTagsPuppi'),
    deepCSVBTagName    = cms.untracked.string('AK8PFDeepCSVJetTagsPuppi'),
    deepCSVBTagNameb   = cms.untracked.string('AK8PFDeepCSVJetTagsPuppi:probb'),
    deepCSVBTagNamec   = cms.untracked.string('AK8PFDeepCSVJetTagsPuppi:probc'),
    deepCSVBTagNamel   = cms.untracked.string('AK8PFDeepCSVJetTagsPuppi:probudsg'),
    deepCSVBTagNamebb  = cms.untracked.string('AK8PFDeepCSVJetTagsPuppi:probbb'),
    deepCMVABTagName   = cms.untracked.string('AK8PFDeepCMVAJetTagsPuppi'),
    deepCMVABTagNameb  = cms.untracked.string('AK8PFDeepCMVAJetTagsPuppi:probb'),
    deepCMVABTagNamec  = cms.untracked.string('AK8PFDeepCMVAJetTagsPuppi:probc'),
    deepCMVABTagNamel  = cms.untracked.string('AK8PFDeepCMVAJetTagsPuppi:probudsg'),
    deepCMVABTagNamebb = cms.untracked.string('AK8PFDeepCMVAJetTagsPuppi:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('AK8PFBoostedDeepDoubleBvLJetTagsPuppi:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('AK8PFBoostedDeepDoubleCvLJetTagsPuppi:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('AK8PFBoostedDeepDoubleCvBJetTagsPuppi:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('AK8PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsPuppi:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('AK8PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsPuppi:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('AK8PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsPuppi:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('AK8PFBoostedDoubleSVTagInfosPuppi'),
    edmMuonName        = cms.untracked.string('slimmedMuons'),
    edmElectronName    = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('AK8PFSoftPFMuonsTagInfosPuppi'),
    softPFElectronTagInfoName = cms.untracked.string('AK8PFSoftPFElectronsTagInfosPuppi'),
    jettiness          = cms.untracked.string('AK8NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('AK8QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('AK8QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('HEP')
  ),

  CA15CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(False),
    useAOD               = cms.untracked.bool(True),
    useTriggerObject     = cms.untracked.bool(False),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    showerDecoConf       = cms.untracked.string(''),
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak8chsL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak8chsL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK8chs')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    jetBoostedBtaggingFiles = cms.untracked.string('BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),

    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsCHS'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15FlavorCHS'),
    prunedJetName      = cms.untracked.string('CA15caPFJetsPrunedCHS'),
    trimmedJetName     = cms.untracked.string('CA15caPFJetsTrimmedCHS'),
    softdropJetName    = cms.untracked.string('CA15caPFJetsSoftDropCHS'),
    subJetName         = cms.untracked.string('CA15caPFJetsSoftDropCHS'),
    csvBTagName        = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS'),
    mvaBTagName        = cms.untracked.string('CA15PFCombinedMVAV2BJetTagsCHS'),
    cvlcTagName        = cms.untracked.string('CA15PFCombinedCvsLJetTagsCHS'),
    cvbcTagName        = cms.untracked.string('CA15PFCombinedCvsBJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    csvDoubleBTagName  = cms.untracked.string('CA15PFBoostedDoubleSecondaryVertexBJetTagsCHS'),
    deepCSVBTagName    = cms.untracked.string('CA15PFDeepCSVJetTagsCHS'),
    deepCSVBTagNameb   = cms.untracked.string('CA15PFDeepCSVJetTagsCHS:probb'),
    deepCSVBTagNamec   = cms.untracked.string('CA15PFDeepCSVJetTagsCHS:probc'),
    deepCSVBTagNamel   = cms.untracked.string('CA15PFDeepCSVJetTagsCHS:probudsg'),
    deepCSVBTagNamebb  = cms.untracked.string('CA15PFDeepCSVJetTagsCHS:probbb'),
    deepCMVABTagName   = cms.untracked.string('CA15PFDeepCMVAJetTagsCHS'),
    deepCMVABTagNameb  = cms.untracked.string('CA15PFDeepCMVAJetTagsCHS:probb'),
    deepCMVABTagNamec  = cms.untracked.string('CA15PFDeepCMVAJetTagsCHS:probc'),
    deepCMVABTagNamel  = cms.untracked.string('CA15PFDeepCMVAJetTagsCHS:probudsg'),
    deepCMVABTagNamebb = cms.untracked.string('CA15PFDeepCMVAJetTagsCHS:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('CA15PFBoostedDeepDoubleBvLJetTagsCHS:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('CA15PFBoostedDeepDoubleCvLJetTagsCHS:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('CA15PFBoostedDeepDoubleCvBJetTagsCHS:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('CA15PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsCHS:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('CA15PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsCHS:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('CA15PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsCHS:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('CA15PFBoostedDoubleSVTagInfosCHS'),
    edmMuonName               = cms.untracked.string('slimmedMuons'),
    edmElectronName           = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('CA15PFSoftPFMuonsTagInfosCHS'),
    softPFElectronTagInfoName = cms.untracked.string('CA15PFSoftPFElectronsTagInfosCHS'),
    jettiness          = cms.untracked.string('CA15NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsCHS'),
    topTaggerName      = cms.untracked.string('HEP')
  ),
  CA15Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useTriggerObject     = cms.untracked.bool(False),
    useAOD               = cms.untracked.bool(True),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doComputeSVInfo      = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak8PuppiL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak8PuppiL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK8Puppi')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    jetBoostedBtaggingFiles = cms.untracked.string('BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),
    showerDecoConf     = cms.untracked.string(''),
    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsPuppi'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15FlavorPuppi'),
    prunedJetName      = cms.untracked.string('CA15caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('CA15caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('CA15caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('CA15caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    mvaBTagName        = cms.untracked.string('CA15PFCombinedMVAV2BJetTagsPuppi'),
    cvlcTagName        = cms.untracked.string('CA15PFCombinedCvsLJetTagsPuppi'),
    cvbcTagName        = cms.untracked.string('CA15PFCombinedCvsBJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    csvDoubleBTagName  = cms.untracked.string('CA15PFBoostedDoubleSecondaryVertexBJetTagsPuppi'),
    deepCSVBTagName    = cms.untracked.string('CA15PFDeepCSVJetTagsPuppi'),
    deepCSVBTagNameb   = cms.untracked.string('CA15PFDeepCSVJetTagsPuppi:probb'),
    deepCSVBTagNamec   = cms.untracked.string('CA15PFDeepCSVJetTagsPuppi:probc'),
    deepCSVBTagNamel   = cms.untracked.string('CA15PFDeepCSVJetTagsPuppi:probudsg'),
    deepCSVBTagNamebb  = cms.untracked.string('CA15PFDeepCSVJetTagsPuppi:probbb'),
    deepCMVABTagName   = cms.untracked.string('CA15PFDeepCMVAJetTagsPuppi'),
    deepCMVABTagNameb  = cms.untracked.string('CA15PFDeepCMVAJetTagsPuppi:probb'),
    deepCMVABTagNamec  = cms.untracked.string('CA15PFDeepCMVAJetTagsPuppi:probc'),
    deepCMVABTagNamel  = cms.untracked.string('CA15PFDeepCMVAJetTagsPuppi:probudsg'),
    deepCMVABTagNamebb = cms.untracked.string('CA15PFDeepCMVAJetTagsPuppi:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('CA15PFBoostedDeepDoubleBvLJetTagsPuppi:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('CA15PFBoostedDeepDoubleCvLJetTagsPuppi:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('CA15PFBoostedDeepDoubleCvBJetTagsPuppi:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('CA15PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsPuppi:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('CA15PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsPuppi:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('CA15PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsPuppi:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('CA15PFBoostedDoubleSVTagInfosPuppi'),
    edmMuonName               = cms.untracked.string('slimmedMuons'),
    edmElectronName           = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('CA15PFSoftPFMuonsTagInfosPuppi'),
    softPFElectronTagInfoName = cms.untracked.string('CA15PFSoftPFElectronsTagInfosPuppi'),
    jettiness          = cms.untracked.string('CA15NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('HEP')
  ),
  
  PFCand = cms.untracked.PSet(
    isActive       = cms.untracked.bool(False),
    edmName        = cms.untracked.string('packedPFCandidates'),
    edmPVName      = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    doAddDepthTime = cms.untracked.bool(False)
  )
)

# overwrite parameters for different eras 
if options.era == '2017' or options.era == '2018':
    process.ntupler.useTriggerObject = cms.untracked.bool(False)
    process.ntupler.TriggerObject = cms.untracked.string("slimmedPatTrigger")

    #electron
    process.ntupler.Electron.edmEleMediumMVA = cms.untracked.string('mvaEleID-Fall17-noIso-V1-wp90')
    process.ntupler.Electron.edmEleTightMVA = cms.untracked.string('mvaEleID-Fall17-noIso-V1-wp80')
    #process.ntupler.Electron.edmEleMVA = cms.untracked.string('ElectronMVAEstimatorRun2Fall17NoIsoV1')
    #process.ntupler.Electron.edmEleMediumMVAIso = cms.untracked.string('mvaEleID-Fall17-iso-V1-wp90')
    #process.ntupler.Electron.edmEleTightMVAIso = cms.untracked.string('mvaEleID-Fall17-iso-V1-wp80')
    #process.ntupler.Electron.edmEleMVAIso = cms.untracked.string('ElectronMVAEstimatorRun2Fall17IsoV1')
    #process.ntupler.Electron.edmEleMVAHZZ = cms.untracked.string('')
    #process.ntupler.Electron.storeSecondMVA = cms.untracked.bool(True)
    #process.ntupler.Electron.storeHZZMVA = cms.untracked.bool(False)

    #AK4CHS
    process.ntupler.AK4CHS.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.AK4CHS.BRegNNMean = cms.untracked.double(1.0610932111740112)
    process.ntupler.AK4CHS.BRegNNStd = cms.untracked.double(0.39077115058898926)
    process.ntupler.AK4CHS.deepCMVABTagNameb = cms.untracked.string('ak4pfDeepCMVAJetTags:probb')
    process.ntupler.AK4CHS.deepCMVABTagNamec = cms.untracked.string('ak4pfDeepCMVAJetTags:probc')
    process.ntupler.AK4CHS.deepCMVABTagNamel = cms.untracked.string('ak4pfDeepCMVAJetTags:probudsg')
    process.ntupler.AK4CHS.deepCMVABTagNamebb = cms.untracked.string('ak4pfDeepCMVAJetTags:probbb')

    #AK4Puppi
    process.ntupler.AK4Puppi.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.AK4Puppi.BRegNNMean = cms.untracked.double(1.0610932111740112)
    process.ntupler.AK4Puppi.BRegNNStd = cms.untracked.double(0.39077115058898926)
    process.ntupler.AK4Puppi.genJetName = cms.untracked.string('slimmedGenJets')

    #AK8CHS
    process.ntupler.AK8CHS.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.AK8CHS.BRegNNMean = cms.untracked.double(1.0610932111740112)
    process.ntupler.AK8CHS.BRegNNStd = cms.untracked.double(0.39077115058898926)
    process.ntupler.AK8CHS.gqlikelihood = cms.untracked.string('AK8QGTaggerCHS') # THIS LOOKS LIKE A BUG!! Difference between files in whether likelihood is capitalized.

    #AK8Puppi 
    process.ntupler.AK8Puppi.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.AK8Puppi.BRegNNMean = cms.untracked.double(1.0610932111740112)
    process.ntupler.AK8Puppi.BRegNNStd = cms.untracked.double(0.39077115058898926)
    process.ntupler.AK8Puppi.genJetName = cms.untracked.string('slimmedGenJetsAK8')

    #CA15CHS
    process.ntupler.CA15CHS.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.CA15CHS.BRegNNMean = cms.untracked.double(1.0610932111740112)
    process.ntupler.BRegNNStd = cms.untracked.double(0.39077115058898926)

    #CA15Puppi
    process.ntupler.CA15Puppi.jetBoostedBtaggingFiles = cms.untracked.string('BaconProd/Utils/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml')
    process.ntupler.CA15Puppi.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.CA15Puppi.BRegNNMean = cms.untracked.double(1.0610932111740112)
    process.ntupler.CA15Puppi.BRegNNStd = cms.untracked.double(0.39077115058898926)
    process.ntupler.CA15Puppi.genJetName = cms.untracked.string('slimmedGenJetsAK8') 

if options.isData:
    process.ntupler.GenInfo.fillLHEWeights = cms.untracked.bool(False)
    process.ntupler.AK4CHS.deepCMVABTagName = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNameb = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNamec = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNamel = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNamebb = cms.untracked.string('')
    

if options.isData:
    process.baconSequence = cms.Sequence(
                                         process.BadPFMuonFilter          *
                                         process.BadChargedCandidateFilter*
                                         process.ak4chsL1FastL2L3ResidualChain*
                                         process.ak4PuppiL1FastL2L3ResidualChain*
                                         process.ak8PuppiL1FastL2L3ResidualChain*
                                         process.ak4chsL1FastL2L3Corrector*
                                         process.ak4PuppiL1FastL2L3Corrector*
                                         process.pfNoPileUpJME            *
                                         process.egammaPostRecoSeq        *
                                         process.prefiringweight          *
                                         process.puppiMETSequence          *
                                         process.AK4jetsequencePuppiData  *
                                         process.AK8jetsequencePuppiData  *
                                         process.CA15jetsequencePuppiData *
                                         process.patJetCorrFactors*
                                         process.updatedPatJets*
                                         process.btagging                 *
                                         process.fullPatMetSequenceV2     *
                                         process.fullPatMetSequencePuppi  *
                                         process.patJetCorrFactorsTransientCorrected*
                                         process.updatedPatJetsTransientCorrected*
                                         process.selectedUpdatedPatJets   *
                                         process.QGTagger                 *
                                         process.ntupler
                                         )
else:
    if options.era == '2016':
        process.baconSequence = cms.Sequence(
                                             process.BadPFMuonFilter          *
                                             process.BadChargedCandidateFilter*
                                             process.ak4chsL1FastL2L3Chain    *
                                             process.ak4PuppiL1FastL2L3Chain  *
                                             process.ak8PuppiL1FastL2L3Chain  *
                                             process.pfNoPileUpJME            *
                                             process.egammaPostRecoSeq        *
                                             process.prefiringweight          *
                                             process.puppiMETSequence         *
                                             process.genjetsequence           *
                                             process.AK4genjetsequenceCHS     *
                                             process.AK4jetsequencePuppi      *
                                             process.AK8jetsequenceCHS        *
                                             process.AK8jetsequencePuppi      *
                                             process.CA15jetsequenceCHS       *
                                             process.CA15jetsequencePuppi     *
                                             process.patJetCorrFactors*
                                             process.updatedPatJets*
                                             process.btagging                 *
                                             process.fullPatMetSequenceV2     *
                                             process.fullPatMetSequencePuppi  *
                                             process.patJetCorrFactorsTransientCorrected*
                                             process.updatedPatJetsTransientCorrected*
                                             process.selectedUpdatedPatJets*
                                             process.QGTagger              *
                                             process.ntupler)
    else:
        process.baconSequence = cms.Sequence(
                                             process.BadPFMuonFilter          *
                                             process.BadChargedCandidateFilter*
                                             process.ak4chsL1FastL2L3Chain    *
                                             process.ak4PuppiL1FastL2L3Chain  *
                                             process.ak8PuppiL1FastL2L3Chain  *
                                             process.pfNoPileUpJME            *
                                             process.egammaPostRecoSeq        *
                                             process.prefiringweight          *
                                             process.puppiMETSequence         *
                                             process.genjetsequence           *
                                             process.AK4genjetsequenceCHS     *
                                             process.AK4jetsequencePuppi      *
                                             process.AK8jetsequenceCHS        *
                                             process.AK8jetsequencePuppi      *
                                             process.CA15jetsequenceCHS       *
                                             process.CA15jetsequencePuppi     *
                                             process.patJetCorrFactors*
                                             process.updatedPatJets*
                                             process.btagging                 *
                                             process.fullPatMetSequenceV2     *
                                             #process.fullPatMetSequencePuppi  *
                                             process.patJetCorrFactorsTransientCorrected*
                                             process.updatedPatJetsTransientCorrected*
                                             process.selectedUpdatedPatJets*
                                             process.QGTagger              *
                                             process.ntupler)

#--------------------------------------------------------------------------------
# apply trigger filter, if necessary
#================================================================================
if options.doHLTFilter:
  process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
  process.hltHighLevel.throw = cms.bool(False)
  process.hltHighLevel.HLTPaths = cms.vstring()
  hlt_file = open(cmssw_base + "/src/" + hlt_filename, "r")
  for line in hlt_file.readlines():
    line = line.strip()              # strip preceding and trailing whitespaces
    if (line[0:3] == 'HLT'):         # assumes typical lines begin with HLT path name (e.g. HLT_Mu15_v1)
      hlt_path = line.split()[0]
      process.hltHighLevel.HLTPaths.extend(cms.untracked.vstring(hlt_path))
  process.p = cms.EndPath(process.hltHighLevel*process.baconSequence)
else:
  process.p = cms.EndPath(process.baconSequence)

#--------------------------------------------------------------------------------
# simple checks to catch some mistakes...
#================================================================================
if options.isData:
  assert process.ntupler.GenInfo.isActive == cms.untracked.bool(False)
  assert process.ntupler.AK4CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA8CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA15CHS.doGenJet == cms.untracked.bool(False)

#process.out = cms.OutputModule("PoolOutputModule",                                                                                                                                                                                                                
#                               outputCommands = cms.untracked.vstring('keep *'),                                                                                                                                                                                
#                               fileName       = cms.untracked.string ("test.root")                                                                                                                                                                              
#                               )                                                                                                                                                                                                                                
#process.endpath = cms.EndPath(process.out)    
