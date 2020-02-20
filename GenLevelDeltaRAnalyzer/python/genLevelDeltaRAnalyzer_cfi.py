import FWCore.ParameterSet.Config as cms

genLevelDeltaRAnalyzer = cms.EDAnalyzer('GenLevelDeltaRAnalyzer',
                                 packedGenParticlesSrc = cms.InputTag("packedGenParticles"),
                                 prunedGenParticlesSrc = cms.InputTag("prunedGenParticles"),
                                 outputPath = cms.untracked.string("output.root"),
                                 verbosity = cms.untracked.int32(0)
)
