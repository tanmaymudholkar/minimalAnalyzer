import FWCore.ParameterSet.Config as cms

stealthTriggerEfficiency = cms.EDAnalyzer('StealthTriggerEfficiency',
                                          photonSrc = cms.InputTag("slimmedPhotons"),
                                          verbosity = cms.untracked.int32(0)
)
