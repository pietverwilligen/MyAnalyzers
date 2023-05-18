import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCEtaUpscopeXML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
            # =====   CMSSW 6 2 0 SLHC5   =====
            'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_cfi_GEN-SIM.root'
            # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_cfi_GEN-SIM_withGEM_withoutMuonNumberingChanges.root'
            # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_cfi_GEN-SIM_withGEM.root'
            # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/FourMuPt_1_200_cfi_GEN-SIM.root'
            )
)

process.demo = cms.EDAnalyzer('MySimHitAnalyzer',
                              RootFileName = cms.untracked.string("MySimHistograms.root"),
)


process.p = cms.Path(process.demo)
