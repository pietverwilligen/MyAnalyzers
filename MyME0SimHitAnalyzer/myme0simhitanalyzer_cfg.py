import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")



process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Configuration.Geometry.GeometryExtended_cff')
# process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')
# process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015_cff')
# process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
# process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
# process.load("Configuration.Geometry.GeometryExtended2021_cff")
# process.load("Configuration.Geometry.GeometryExtended2021Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2026D99_cff")
process.load("Configuration.Geometry.GeometryExtended2026D99Reco_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
# process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
# process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
# process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
# process.load("Geometry.CSCGeometry.cscGeometry_cfi")
# process.load("Geometry.DTGeometry.dtGeometry_cfi")
# process.load("Geometry.GEMGeometry.gemGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

# Load Events from python file
# ---------------------------------------------------------------------------------------
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# readFiles = cms.untracked.vstring()
# secFiles = cms.untracked.vstring() 
# source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
# ---------------------------------------------------------------------------------------
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
# process.source = cms.Source("PoolSource",
#      fileNames = cms.untracked.vstring('file:/lustre/home/piet/Experiment_CMS/GEMDPG/ME0-fluka-geant-study/MinBias_1100pre13_23240.0_14TeV_2026D49.root')
# )
# ---------------------------------------------------------------------------------------
# process.load("MinBias_2026D49_RelVal_CMSSW_11_0_0_pre3")
# process.load("MinBias_13TeV_GEN_SIM_2021_100k_11X_ProdCutsDefault_v1_allfiles")
# process.load("ZMM_14TeV_GEN_SIM_2021_100k_11X_ProdCutsDefault_v1_allfiles")
# process.load("MinBias_14TeV_GEN_SIM_XS_2026D49_100k_11X_ProdCutsDefault_v1_allfiles")
# process.load("MinBias_14TeV_GEN_SIM_NoBkg_2026D49_100k_11X_ProdCutsDefault_v2_allfiles")
# process.load("MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p1_allfiles")
# ---------------------------------------------------------------------------------------
# CMSSW 13X
# ---------------------------------------------------------------------------------------
# process.load("MinBias_Phase1_14TeV_TuneCP5_100k_Neutron_XS_2023DB_1E4s")
process.load("MinBias_Phase1_14TeV_TuneCP5_1M_Neutron_XS_2023DB_1E4s")
# process.load("MinBias_Phase2_14TeV_TuneCP5_100k_Neutron_XS_2026D99_1E4s")
# ---------------------------------------------------------------------------------------

process.demo = cms.EDAnalyzer('MyME0SimHitAnalyzer',
                              # ---------
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M_v1"),
                              # RootFileName = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M_v1.root"),
                              # PdfFileNameBase = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M__150eV_v2"),
                              # RootFileName = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M_150eV_v2.root"),
                              # PdfFileNameBase = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M__300eV_v2"),
                              # RootFileName = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M_300eV_v2.root"),
                              PdfFileNameBase = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M__450eV_v2"),
                              RootFileName = cms.untracked.string("MyME0Histograms_2023DB_TuneCP5_14TeV_Neutron_XS_1M_450eV_v2.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyME0Histograms_2026D99_TuneCP5_14TeV_Neutron_XS_100k_v1"),
                              # RootFileName = cms.untracked.string("MyME0Histograms_2026D99_TuneCP5_14TeV_Neutron_XS_100k_v1.root"),
                              # ---------
                              BunchSpacing = cms.untracked.double(25.0),
                              # COMEnergy    = cms.untracked.double(13.0),
                              COMEnergy    = cms.untracked.double(14.0),
                              # MaxSimTime   = cms.untracked.double(1000000000000.0), # 1000s = 10^12 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(100000000000.0),  # 100s  = 10^11 ns [in ns]
                              MaxSimTime   = cms.untracked.double(10000000000.0),     # 10s   = 10^10 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(100000000.0),     # 100ms = 10^8  ns [in ns]
                              # GEMOnlyGE11  = cms.untracked.bool(True), # not in use anymore
                              Edep30eV     = cms.untracked.bool(True),    # cut on minimal Energy deposition for ionisation
                              # MinHitEnergy = cms.untracked.double(0.030), # minimum Energy for hit creation (30eV)
                              # MinHitEnergy = cms.untracked.double(0.150), # minimum Energy for hit creation (30eV)
                              # MinHitEnergy = cms.untracked.double(0.300), # minimum Energy for hit creation (30eV)
                              MinHitEnergy = cms.untracked.double(0.450), # minimum Energy for hit creation (30eV)
                              MinHIPEnergy = cms.untracked.double(030.0), # minimum Energy for HIP creation (150keV = 5000 e- | 30keV = 1000 e-)
                              PhysicsDebug = cms.untracked.bool(False),
                              TechnicDebug = cms.untracked.bool(False),
                              )


process.p = cms.Path(process.demo)
