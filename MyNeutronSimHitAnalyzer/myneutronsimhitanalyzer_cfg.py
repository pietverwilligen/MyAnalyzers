import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Configuration.Geometry.GeometryExtended_cff')
# process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
# process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
# process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

# process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
# process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
# process.load("Geometry.GEMGeometry.gemGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

# Load Events from python file
# ---------------------------------------------------------------------------------------
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# readFiles = cms.untracked.vstring()
# secFiles = cms.untracked.vstring() 
# source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
# ---------------------------------------------------------------------------------------
process.load("NeutronBackgroundFiles_Eta8_1000s_97k")
# process.load("NeutronBackgroundFiles_Eta8_100s_100k")
# process.load("NeutronBackgroundFiles_Eta8_10s_100k")
# process.load("NeutronBackgroundFiles_Eta8_100ms_100k")
# process.load("NeutronBackgroundFiles_Eta7_100ms_100k")
# process.load("NeutronBackgroundFiles_Eta8_10s_100k_test")
# ---------------------------------------------------------------------------------------

process.demo = cms.EDAnalyzer('MyNeutronSimHitAnalyzer',
                              # ---------
                              PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_1000s_97k_76X"),
                              RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_1000s_97k_76X.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_100s_100k_76X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_100s_100k_76X.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_10s_100k_76X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_10s_100k_76X.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_100ms_100k_76X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta8_100ms_100k_76X.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta7_100ms_100k_76X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_xs_eta7_100ms_100k_76X.root"),
                              # ---------
                              BunchSpacing = cms.untracked.double(25.0),
                              COMEnergy    = cms.untracked.double(13.0),
                              MaxSimTime   = cms.untracked.double(1000000000000.0),   # 1000s = 10^12 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(100000000000.0),  # 100s  = 10^11 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(10000000000.0),   # 10s   = 10^10 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(100000000.0),     # 100ms = 10^8  ns [in ns]
                              PhysicsDebug = cms.untracked.bool(False),
                              TechnicDebug = cms.untracked.bool(False),
                              )


process.p = cms.Path(process.demo)
