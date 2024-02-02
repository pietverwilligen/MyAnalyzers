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
# process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff') # ... assume not necessary anymore ... 
# process.load('Configuration.Geometry.GeometryExtended2018Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2018_cff')
process.load('Configuration.Geometry.GeometryExtended2026D99Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D99_cff')


# process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
# process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")                            # ... needed? see if I can get rid of it ...
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
# process.load("Geometry.GEMGeometry.gemGeometry_cfi")                          # ... does not exist ...
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

# Load Events from python file
# ---------------------------------------------------------------------------------------
# option A
# ---------------------------------------------------------------------------------------
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# readFiles = cms.untracked.vstring()
# secFiles = cms.untracked.vstring() 
# source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
# ---------------------------------------------------------------------------------------
# option B
# ---------------------------------------------------------------------------------------
# process.load("MinBias_Phase2_14TeV_TuneCP5_100k_Neutron_XS_2026D99_1E4s")
process.load("Test_MinBias_Phase2_14TeV_TuneCP5_100k_Neutron_XS_2026D99_1E4s")
# ---------------------------------------------------------------------------------------
# option C
# ---------------------------------------------------------------------------------------
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
# process.source = cms.Source ("PoolSource",
#     fileNames = cms.untracked.vstring('/store/user/piet/NeutronBackground/MinBias_Phase2_14TeV_GEN_SIM_XS_2026D99mod_100k_1E4s_13X_v1/crab_MinBias_Phase2_14TeV_100k_1E4s_XS_13X_v1/230504_162117/0000/step1_228.root'))
# ---------------------------------------------------------------------------------------

process.demo = cms.EDAnalyzer('MyNeutronSimHitAnalyzer',
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_Run2_Neutron_XS_1E4s"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_Run2_Neutron_XS_1E4s.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_Run2_Neutron_XS_1E4s_SH30eV"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_Run2_Neutron_XS_1E4s_SH30eV.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_Run2_Neutron_XS_1E4s_Test"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_Run2_Neutron_XS_1E4s_Test.root"),
                              # ---------
                              PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_Phase2_2026D99_Neutron_XS_1E4s"),
                              RootFileName = cms.untracked.string("MyNeutronSimHistograms_Phase2_2026D99_Neutron_XS_1E4s.root"),
                              # ---------
                              BunchSpacing = cms.untracked.double(25.0),
                              COMEnergy    = cms.untracked.double(13.0),
                              MaxSimTime   = cms.untracked.double(10000000000000.0),    # 10000s = 10^13 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(1000000000000.0),   # 1000s  = 10^12 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(100000000000.0),    # 100s   = 10^11 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(10000000000.0),     # 10s    = 10^10 ns [in ns]
                              # MaxSimTime   = cms.untracked.double(100000000.0),       # 100ms  = 10^8  ns [in ns]
                              EDepCut30eV  = cms.untracked.bool(True),
                              PhysicsDebug = cms.untracked.bool(True),
                              TechnicDebug = cms.untracked.bool(True),
                              PDFOutput    = cms.untracked.bool(False), # set false for large number of events
                              )


process.p = cms.Path(process.demo)
