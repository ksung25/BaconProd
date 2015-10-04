import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
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
genParticlesForJetsNoNu = genParticlesForJets.clone()
genParticlesForJetsNoNu.ignoreParticleIDs += cms.vuint32(12,14,16,18)

# Select hadrons and partons for Jet Flavour
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons


genjetsequence = cms.Sequence(
  genParticlesForJets*
  genParticlesForJetsNoNu*
  selectedHadronsAndPartons
)

def setMiniAODGenJets(process):
    process.genParticlesForJets.src             = 'packedGenParticles'
    process.genParticlesForJetsNoNu.src         = 'packedGenParticles'
    process.selectedHadronsAndPartons.particles = 'prunedGenParticles'
