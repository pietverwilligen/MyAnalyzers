import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoMuon.Configuration.RecoMuon_cff")
# process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load("Configuration.StandardSequences.Geometry_cff")

process.load('Configuration.Geometry.GeometryExtended2023_cff')
process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

# process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

        ### 2023 ###
        # '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/128D04DF-9AD4-E311-A75F-0025905A60CE.root',
        # '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/40B0C80D-A2D4-E311-AEFE-0025905A48C0.root',
        ### 2019 ###
        # '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES19_62_V8_UPG2019withGEM-v1/00000/52A87B54-A8D4-E311-A7E9-003048D2BF34.root',
        # '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES19_62_V8_UPG2019withGEM-v1/00000/1074FA37-99D4-E311-AA75-001D09F24FEC.root',

        ### GEM Displaced Muon Study ###
        # 'file:/cmshome/piet/SLC6/GEM_DisplacedMuon_2014/CMSSW_6_2_0_SLHC21/src/MyCmsDriverCommands/SingleMuPt10_GEN-SIM-DIGI-RAW-RECO_upgrade2019.root'
        # 'file:/cmshome/piet/SLC6/GEM_DisplacedMuon_2014/CMSSW_6_2_0_SLHC21/src/MyCmsDriverCommands/DarkSUSY_mH_125_mGammaD_0400_ctau_05_8TeV_RECO.root'
        'file:/cmshome/piet/SLC6/GEM_DisplacedMuon_2014/CMSSW_6_2_0_SLHC21/src/MyCmsDriverCommands/DarkSUSY_mH_125_mGammaD_0400_ctau_50p0_8TeV_RECO.root'
    )
)

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.demo = cms.EDAnalyzer('MyDisplacedMuonAnalyzer',
                              PhysicsDebug = cms.untracked.bool(False),
                              TechnicDebug = cms.untracked.bool(False),
                              GenPartDebug = cms.untracked.bool(False),
                              MCTruthMatching = cms.untracked.bool(True),
                              RecoDRMatching  = cms.untracked.bool(True), # instead of fixed Delta R Match (Reco Muon, Sim Track) < 0.15, assign closesth by Muon SimTrack to Reco Muon
                              RootFileName = cms.untracked.string("DisplacedStandAloneMuon_62X_2019.root"),
                              StandAloneTrackCollectionLabel1 = cms.InputTag("standAloneMuons",        ""),
                              StandAloneTrackCollectionLabel2 = cms.InputTag("refittedStandAloneMuons",""),
                              StandAloneTrackCollectionLabel3 = cms.InputTag("standAloneMuons",        "UpdatedAtVtx"),
                              StandAloneTrackCollectionLabel4 = cms.InputTag("refittedStandAloneMuons","UpdatedAtVtx"),
)


process.p = cms.Path(process.demo)
