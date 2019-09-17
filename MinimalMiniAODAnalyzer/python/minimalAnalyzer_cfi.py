import FWCore.ParameterSet.Config as cms

minimalAnalyzer = cms.EDAnalyzer('MinimalMiniAODAnalyzer',
                                 photonSrc = cms.InputTag("slimmedPhotons"),
                                 genParticleSrc = cms.InputTag("prunedGenParticles"),
                                 jetSrc = cms.InputTag("slimmedJets"),
                                 rhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),
                                 outputPath = cms.untracked.string("output.root"),
                                 verbosity = cms.untracked.int32(0),
                                 filterType = cms.untracked.string("none"),
                                 pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
                                 selectJetsNearPhotons = cms.untracked.bool(False),
                                 selectJetsAwayFromPhotons = cms.untracked.bool(False)
)
