import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *
from CondCore.DBCommon.CondDBSetup_cfi import *

def setupJEC(process,isData) :
    label='MC'
    if isData:
        label='DATA'
    process.jec =  cms.ESSource("PoolDBESSource",
                                CondDBSetup,
                                toGet = cms.VPSet(
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_'+label+'_AK4PFPuppi'),
                                   label   = cms.untracked.string('AK4Puppi')
                                   ),
                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                    tag     = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_'+label+'_AK8PFPuppi'),
                                    label   = cms.untracked.string('AK8Puppi')
                                    ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_'+label+'_AK4PFchs'),
                                   label   = cms.untracked.string('AK4chs')
                                   ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_'+label+'_AK8PFchs'),
                                   label   = cms.untracked.string('AK8chs')
                                   ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_'+label+'_AK4PF'),
                                   label   = cms.untracked.string('AK4')
                                   ),
                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                    tag     = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_'+label+'_AK8PF'),
                                    label   = cms.untracked.string('AK8')
                                    )
                           ),

                    )                                        

#es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

qgDatabaseVersion = 'v1' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDBSetup,
      toGet = cms.VPSet(),
      connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
)

for type in ['AK4PFchs','AK4PFchs_antib']:
    QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
        record = cms.string('QGLikelihoodRcd'),
        tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
        label  = cms.untracked.string('QGL_'+type)
        )))

es_prefer_qgl = cms.ESPrefer("PoolDBESSource","QGPoolDBESSource")
