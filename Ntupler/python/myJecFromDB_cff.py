import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *
#from CondCore.DBCommon.CondDBSetup_cfi import *
from CondCore.CondDB.CondDB_cfi import *

def setupJEC(process,isData,label) :
    process.jec =  cms.ESSource("PoolDBESSource",
                                CondDB,
                                toGet = cms.VPSet(
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_'+label+'_AK4PFPuppi'),
                                   label   = cms.untracked.string('AK4Puppi')
                                   ),
                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                    tag     = cms.string('JetCorrectorParametersCollection_'+label+'_AK8PFPuppi'),
                                    label   = cms.untracked.string('AK8Puppi')
                                    ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_'+label+'_AK4PFchs'),
                                   label   = cms.untracked.string('AK4chs')
                                   ),
                          cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                   tag     = cms.string('JetCorrectorParametersCollection_'+label+'_AK8PFchs'),
                                   label   = cms.untracked.string('AK8chs')
                                   ),
                           ),

                    )                                        
    ##qgDatabaseVersion = 'v2b' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
    #for type in ['AK4PFchs','AK4PFchs_antib']:
    #    process.jec.toGet.extend(cms.VPSet(cms.PSet(
    #                record = cms.string('QGLikelihoodRcd'),
    #                tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
    #                label  = cms.untracked.string('QGL_'+type)
    #                )))

    #es_prefer_qgl = cms.ESPrefer("PoolDBESSource","QGPoolDBESSource")
