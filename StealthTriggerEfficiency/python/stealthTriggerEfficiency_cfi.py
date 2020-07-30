import FWCore.ParameterSet.Config as cms

stealthTriggerEfficiency = cms.EDAnalyzer('StealthTriggerEfficiency',
                                          rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                                          triggerCollection = cms.InputTag("TriggerResults","","HLT"),
                                          photonSrc = cms.InputTag("slimmedPhotons"),
                                          verbosity = cms.untracked.int32(0)
)
