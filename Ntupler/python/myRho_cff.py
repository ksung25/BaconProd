import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
kt6PFJets = kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True )                                                                                                                                   
