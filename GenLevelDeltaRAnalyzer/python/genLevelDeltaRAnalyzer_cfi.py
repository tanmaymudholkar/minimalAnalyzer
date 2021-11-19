import FWCore.ParameterSet.Config as cms

genLevelDeltaRAnalyzer = cms.EDAnalyzer('GenLevelDeltaRAnalyzer',
                                        LHEInfoSrc = cms.InputTag("externalLHEProducer"),
                                        prunedGenParticlesSrc = cms.InputTag("prunedGenParticles"),
                                        genJetsSrc = cms.InputTag("slimmedGenJets"),
                                        # eventProgenitor = cms.untracked.string("none"),
                                        verbosity = cms.untracked.int32(0)
)
