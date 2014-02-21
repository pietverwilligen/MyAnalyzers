import FWCore.ParameterSet.Config as cms

process = cms.Process("Prod")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_1_PostLS1v1/src/cmsDriverCommands/SingleMuPt100_cfi_RECO.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_cfi_RECO.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyCmsDriverCommands/SingleMuPt100_cfi_1p6_2p5_RECO.root'
        'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyCmsDriverCommands/SingleMuPt100_cfi_1p6_2p5_DIGI-RAW.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyAnalyzers/MyRecHitAnalyzer/MyRecHits.root')
)

# process.MessageLogger = cms.Service("MessageLogger",
#     debugModules = cms.untracked.vstring('*'),
#     cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
#     # cout = cms.untracked.PSet(threshold = cms.untracked.string('DEBUG')),
#     destinations = cms.untracked.vstring('cout')
# ) 

process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
process.p = cms.Path(process.rpcRecHits)
process.ep = cms.EndPath(process.out)
