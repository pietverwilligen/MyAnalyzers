import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_2p4_192_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p1_64_cfi_DIGI-RAW.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_64_cfi_DIGI-RAW.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_192_cfi_DIGI-RAW.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p1_256_cfi_DIGI-RAW.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p5_192_cfi_DIGI.root'
        'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyCmsDriverCommands/SingleMuPt100_cfi_DIGI-RAW-newCond.root'
    )
)

process.demo = cms.EDAnalyzer('MyREn1DigiAnalyzer',
                              # DATA
                              # DigiLabel= cms.untracked.string("muonRPCDigis")
                              # MONTE-CARLO
                              DigiLabel    = cms.untracked.string("simMuonRPCDigis"),
                              # ROOT Filename
                              RootFileName = cms.untracked.string("MyREn1DigiHistograms_1p6_2p5_192.root"),
                              # DEBUG ?
                              Debug        = cms.untracked.bool(False),
)


process.p = cms.Path(process.demo)
