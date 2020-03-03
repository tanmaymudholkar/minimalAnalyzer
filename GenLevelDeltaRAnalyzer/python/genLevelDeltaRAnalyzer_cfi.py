import FWCore.ParameterSet.Config as cms

genLevelDeltaRAnalyzer = cms.EDAnalyzer('GenLevelDeltaRAnalyzer',
                                        prunedGenParticlesSrc = cms.InputTag("prunedGenParticles"),
                                        genJetsSrc = cms.InputTag("slimmedGenJets"),
                                        outputPath = cms.untracked.string("output.root"),
                                        verbosity = cms.untracked.int32(0)
)
