import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCEtaUpscopeXML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometry2015XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Upgrade
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_256_cfi_DIGI-RAW.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyCmsDriverCommands/SingleMuPt100_cfi_DIGI-RAW-newCond.root'
        # Run 234029
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00097CE4-A2B2-E411-8FEB-02163E011D8A.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/0049B70A-8EB2-E411-8979-02163E011D1C.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/004DEB4A-9DB2-E411-AB38-02163E011D06.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00D04605-A5B2-E411-A496-02163E0129F4.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00D21ADC-A2B2-E411-8861-02163E011D57.root',
    )
)

process.demo = cms.EDAnalyzer('MyDigiAnalyzer',
                              # DATA
                              DigiLabel= cms.untracked.string("muonRPCDigis"),
                              # MONTE-CARLO
                              # DigiLabel    = cms.untracked.string("simMuonRPCDigis"),
                              # ROOT Filename
                              RootFileName = cms.untracked.string("MyDigiHistograms.root"),
)


process.p = cms.Path(process.demo)
