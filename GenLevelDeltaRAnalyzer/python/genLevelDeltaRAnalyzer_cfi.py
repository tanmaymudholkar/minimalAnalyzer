import FWCore.ParameterSet.Config as cms

genLevelDeltaRAnalyzer = cms.EDAnalyzer('GenLevelDeltaRAnalyzer',
                                        prunedGenParticlesSrc = cms.InputTag("prunedGenParticles"),
                                        genJetsSrc = cms.InputTag("slimmedGenJets"),
                                        eventProgenitor = cms.untracked.string("none"),
                                        verbosity = cms.untracked.int32(0)
)
