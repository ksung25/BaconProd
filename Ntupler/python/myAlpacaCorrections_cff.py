import FWCore.ParameterSet.Config as cms
from   CommonTools.PileupAlgos.Alpaca_cff import *
from   RecoMET.METProducers.PFMET_cfi     import pfMet

jecFilesPuppiData = ( cms.untracked.vstring('BaconProd/Utils/data/Summer15_25nsV2_DATA_L1FastJet_AK4PFPuppi.txt',
                                            'BaconProd/Utils/data/Summer15_25nsV2_DATA_L2Relative_AK4PFPuppi.txt',
                                            'BaconProd/Utils/data/Summer15_25nsV2_DATA_L3Absolute_AK4PFPuppi.txt',
                                            'BaconProd/Utils/data/Summer15_25nsV2_DATA_L2L3Residual_AK4PFPuppi.txt')
                      )

jecFilesPuppiMC   = ( cms.untracked.vstring('BaconProd/Utils/data/Summer15_25nsV2_MC_L1FastJet_AK4PFPuppi.txt',
                                            'BaconProd/Utils/data/Summer15_25nsV2_MC_L2Relative_AK4PFPuppi.txt',
                                            'BaconProd/Utils/data/Summer15_25nsV2_MC_L3Absolute_AK4PFPuppi.txt'
                                            )
                      )

alpacaMC.candName   = 'pfCandNoLep'
alpacaData.candName = 'pfCandNoLep'

alpacaPuppiMC = alpacaMC.clone()
alpacaPuppiMC.candName   = 'puppi'
alpacaPuppiMC.chJecFiles = jecFilesPuppiMC
alpacaPuppiMC.emJecFiles = jecFilesPuppiMC
alpacaPuppiMC.nhJecFiles = jecFilesPuppiMC

alpacaPuppiNoLepMC = alpacaPuppiMC.clone()
alpacaPuppiNoLepMC.candName   = 'puppinolep'

alpacaPuppiData = alpacaData.clone()
alpacaPuppiData.candName   = 'puppi'
alpacaPuppiData.chJecFiles = jecFilesPuppiData
alpacaPuppiData.emJecFiles = jecFilesPuppiData
alpacaPuppiData.nhJecFiles = jecFilesPuppiData

alpacaPuppiNoLepData = alpacaData.clone()
alpacaPuppiNoLepData.candName   = 'puppinolep'

alpacaFullMC         = cms.EDProducer("CandViewMerger",src = cms.VInputTag( 'alpacaMC','pfCandLep'))     
pfMetAlpacaMC        = pfMet.clone()
pfMetAlpacaMC.src    = cms.InputTag('alpacaFullMC')

alpacaFullData       = cms.EDProducer("CandViewMerger",src = cms.VInputTag( 'alpacaData','pfCandLep'))     
pfMetAlpacaData      = pfMet.clone()
pfMetAlpacaData.src  = cms.InputTag('alpacaFullData')

puppiAlpacaForMETMC     = cms.EDProducer("CandViewMerger",
                                         src = cms.VInputTag( 'alpacaPuppiNoLepMC','pfCandLep')
                                         )     
pfMetPuppiAlpacaMC      = pfMet.clone()
pfMetPuppiAlpacaMC.src  = cms.InputTag('puppiAlpacaForMETMC')

puppiAlpacaForMETData     = cms.EDProducer("CandViewMerger",
                                           src = cms.VInputTag( 'alpacaPuppiNoLepData','pfCandLep')
                                           )     
pfMetPuppiAlpacaData      = pfMet.clone()
pfMetPuppiAlpacaData.src  = cms.InputTag('puppiAlpacaForMETData')

alpacaSequenceMC   = cms.Sequence(
    alpacaMC *
    alpacaFullMC *
    pfMetAlpacaMC *
    alpacaPuppiNoLepMC * 
    #alpacaPuppiMC *
    puppiAlpacaForMETMC*
    pfMetPuppiAlpacaMC
    )

alpacaSequenceData   = cms.Sequence(
    alpacaData *
    alpacaFullData *
    pfMetAlpacaData * 
    alpacaPuppiNoLepData *
    #alpacaPuppiData*
    puppiAlpacaForMETData*
    pfMetPuppiAlpacaData
    )


def setMiniAODAlpaca(process) :
    process.alpacaMC.candName         = cms.InputTag('packedPFCandidates')
    process.alpacaData.candName       = cms.InputTag('packedPFCandidates')

