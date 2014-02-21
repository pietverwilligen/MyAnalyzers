import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_2p1_64_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_2p4_256_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

# process.load('FWCore.MessageService.MessageLogger_cfi')
# process.load('Configuration.Geometry.GeometryExtended2023RPCUpscopeReco_cff')
# process.load('Configuration.Geometry.GeometryExtended2023RPCUpscope_cff')



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
            # =====   CMSSW 6 2 0 SLHC5  =====
            # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p1_64_cfi_GEN-SIM.root'
            # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_64_cfi_GEN-SIM.root'
            # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_256_cfi_GEN-SIM.root'
            'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p5_192_cfi_GEN-SIM.root'
            )
)

process.demo = cms.EDAnalyzer('MyREn1SimHitAnalyzer',
                              RootFileName = cms.untracked.string("MyREn1SimHistograms_SingleMuPt100_eta1p6_eta2p5_192.root"),
                              Debug        = cms.untracked.bool(False),
                              SimPtCut     = cms.untracked.double(5.0),
)


process.p = cms.Path(process.demo)
