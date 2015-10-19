#
# A modified version of RecoMET/METFilters/python/metFilters_cff.py to set all filters to "tagging" mode
#

import FWCore.ParameterSet.Config as cms

## The iso-based HBHE noise filter ___________________________________________||
from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

## The CSC beam halo tight filter ____________________________________________||
#from RecoMET.METFilters.CSCTightHaloFilter_cfi import * 
from RecoMET.METProducers.CSCHaloData_cfi import *
from RecoMET.METProducers.EcalHaloData_cfi import *
from RecoMET.METProducers.HcalHaloData_cfi import *
from RecoMET.METProducers.GlobalHaloData_cfi import *
from RecoMET.METProducers.BeamHaloSummary_cfi import *
from RecoMET.METFilters.CSCTightHalo2015Filter_cfi
CSCTightHalo2015Filter.taggingMode = cms.bool(True)

## The HCAL laser filter _____________________________________________________||
from RecoMET.METFilters.hcalLaserEventFilter_cfi import *
hcalLaserEventFilter.taggingMode = cms.bool(True)

## The ECAL dead cell trigger primitive filter _______________________________||
from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

## The EE bad SuperCrystal filter ____________________________________________||
from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = cms.bool(True)

## The ECAL laser correction filter
from RecoMET.METFilters.ecalLaserCorrFilter_cfi import *
ecalLaserCorrFilter.taggingMode = cms.bool(True)

## The Good vertices collection needed by the tracking failure filter ________||
goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)

## The tracking failure filter _______________________________________________||
from RecoMET.METFilters.trackingFailureFilter_cfi import *
trackingFailureFilter.taggingMode = cms.bool(True)

## The tracking POG filters __________________________________________________||
from RecoMET.METFilters.trackingPOGFilters_cff import *
manystripclus53X.taggedMode  = cms.untracked.bool(True)
manystripclus53X.forcedValue = cms.untracked.bool(False)
toomanystripclus53X.taggedMode  = cms.untracked.bool(True)
toomanystripclus53X.forcedValue = cms.untracked.bool(False)
logErrorTooManyClusters.taggedMode  = cms.untracked.bool(True)
logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)
## Also the stored boolean for the three filters is opposite to what we usually
## have for other filters, i.e., true means rejected bad events while false means 
## good events.

metFilters = cms.Sequence(
   HBHENoiseFilterResultProducer *
   hcalLaserEventFilter *
   EcalDeadCellTriggerPrimitiveFilter *
   goodVertices * trackingFailureFilter *
   eeBadScFilter *
   ecalLaserCorrFilter *
   trkPOGFilters *
   CSCHaloData*EcalHaloData*HcalHaloData*GlobalHaloData*BeamHaloSummary*
   CSCTightHalo2015Filter
)
