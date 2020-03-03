import FWCore.ParameterSet.Config as cms

MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        GenDeltaRAnalyzer = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        optionalPSet = cms.untracked.bool(True),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        noTimeStamps = cms.untracked.bool(False),
        FwkReport = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1000),
            limit = cms.untracked.int32(10000000)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        FwkJob = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        FwkSummary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1),
            limit = cms.untracked.int32(10000000)
        ),
        threshold = cms.untracked.string('INFO')
    ),
    # cerr = cms.untracked.PSet(
    #     default = cms.untracked.PSet(
    #         limit = cms.untracked.int32(-1)
    #     )
    # ),
    FrameworkJobReport = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        FwkJob = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(10000000)
        )
    ),
    statistics = cms.untracked.vstring('cout_stats'),
    cout_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('WARNING'),
        output = cms.untracked.string('cout')
    ),
    infos = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True),
        optionalPSet = cms.untracked.bool(True),
        Root_NoDictionary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
    ),
    # destinations = cms.untracked.vstring('cout', 'cerr'),
    destinations = cms.untracked.vstring('cout'),
    # debugModules = cms.untracked.vstring(),
    categories = cms.untracked.vstring(
        'GenDeltaRAnalyzer',
        'FwkJob',
        'FwkReport',
        'FwkSummary',
        'Root_NoDictionary'
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport')
)
