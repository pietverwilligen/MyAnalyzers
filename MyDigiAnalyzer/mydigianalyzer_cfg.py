import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCEtaUpscopeXML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2015XML_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Upgrade
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_256_cfi_DIGI-RAW.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyCmsDriverCommands/SingleMuPt100_cfi_DIGI-RAW-newCond.root'

        # 2015 --- Cruzet --- Run 234556
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/556/00000/002EA6F7-A0B6-E411-9FBF-02163E011D47.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/556/00000/00996C86-9AB6-E411-B8B2-02163E012720.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/556/00000/0240143B-A4B6-E411-A064-02163E01263E.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/556/00000/025AA905-9EB6-E411-B40B-02163E011CF0.root',
        '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/556/00000/02A5200C-9EB6-E411-9BF2-02163E012120.root',

        # 2015 --- Cruzet --- Run 234029
        # '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00097CE4-A2B2-E411-8FEB-02163E011D8A.root',
        # '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/0049B70A-8EB2-E411-8979-02163E011D1C.root',
        # '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/004DEB4A-9DB2-E411-AB38-02163E011D06.root',
        # '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00D04605-A5B2-E411-A496-02163E0129F4.root',
        # '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00D21ADC-A2B2-E411-8861-02163E011D57.root',

        # 2015 --- MWGR10 --- Run 232956 and 233121
        # '/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/232/956/00001/000ADBF9-50A8-E411-8DE8-02163E012B06.root'

        # 2012 --- Run2012D --- Run 209283
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/0E265F20-0F8B-E211-A022-0025901248FA.root',
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/1A655E1C-0F8B-E211-B57D-001E67397D64.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/3652CC13-0F8B-E211-970C-001E67398110.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/3ACEEB0F-0F8B-E211-88AE-003048D45F34.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/42DE5920-0F8B-E211-981B-0025902008E4.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/46870E0A-0F8B-E211-8164-003048673F3A.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/5070F7FB-0E8B-E211-A97F-003048D476B8.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/581987FC-0E8B-E211-88A4-003048D45F50.root',  
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/5A336D13-0F8B-E211-9ADA-001E67397620.root', 
        # '/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/5C7FBF18-0F8B-E211-8D9F-002481E14E04.root', 

        # in case of problems with opening remote root files :: issue the command "voms-proxy-init -voms cms"
    )
)

process.demo = cms.EDAnalyzer('MyDigiAnalyzer',
                              Debug = cms.untracked.bool(False),
                              # RunNumber = cms.untracked.int32(209283),
                              # RunNumber = cms.untracked.int32(234029),
                              RunNumber = cms.untracked.int32(234556),
                              # DATA / MONTE-CARLO
                              DigiLabel= cms.untracked.string("muonRPCDigis"),
                              # DigiLabel    = cms.untracked.string("simMuonRPCDigis"),
                              # ROOT Filename
                              RootFileName = cms.untracked.string("MyDigiHistograms.root"),

)


process.p = cms.Path(process.demo)
