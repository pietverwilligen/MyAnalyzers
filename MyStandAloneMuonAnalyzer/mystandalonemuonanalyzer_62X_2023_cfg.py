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
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC12/src/MyCmsDriverCommands/SingleMuPt100_cfi_1p6_2p5_RECO.root',
        # '/store/relval/RelValZMM_14TeV/CMSSW_6_2_0_SLHC13-DES23_62_V1_UPG2023Muon-v1/GEN-SIM-RECO/...',
        # eos ls /eos/store/relval
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/328018D6-5696-E311-859E-003048FFD7D4.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/34D74464-5296-E311-BD5E-0025905A60DA.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/58091D33-2896-E311-8B27-0025905A6066.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/68768EE5-3196-E311-9B1F-003048FFD730.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/6C189AB2-2B96-E311-9962-0025905A4964.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/8ACA0D0C-3496-E311-B3A5-0025905A60B8.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/8E1FDDE5-3696-E311-A5BC-003048FFCC0A.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/ACDCEF24-3C96-E311-BF02-002618943958.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/BA23B9A8-2B96-E311-8A20-0026189438DB.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/C6F5DC64-3196-E311-B546-0025905A6126.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/DE6FE81A-2E96-E311-B494-0025905A60F2.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/EE7B8520-2E96-E311-A5FC-0025905A48D0.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/FA554D43-2696-E311-A8E8-0025905A6104.root',
        # '/store/relval/CMSSW_7_0_0_pre13/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START70_V6-v1/00000/FEC66CC9-2996-E311-AE82-0025905A6104.root',

        ### 2023 ###
        '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/128D04DF-9AD4-E311-A75F-0025905A60CE.root',
        '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/40B0C80D-A2D4-E311-AEFE-0025905A48C0.root',
        ### 2019 ###
        # '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES19_62_V8_UPG2019withGEM-v1/00000/52A87B54-A8D4-E311-A7E9-003048D2BF34.root',
        # '/store/relval/CMSSW_6_2_0_SLHC12/RelValSingleMuPt100/GEN-SIM-RECO/DES19_62_V8_UPG2019withGEM-v1/00000/1074FA37-99D4-E311-AA75-001D09F24FEC.root',
    )
)

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.demo = cms.EDAnalyzer('MyStandAloneMuonAnalyzer',
                              PhysicsDebug = cms.untracked.bool(False),
                              TechnicDebug = cms.untracked.bool(False),
                              GenPartDebug = cms.untracked.bool(False),
                              MCTruthMatching = cms.untracked.bool(False),
                              TightID      = cms.untracked.bool(True),
                              LooseID      = cms.untracked.bool(False),
                              RootFileName = cms.untracked.string("STAMuon_62X_2023.root"),
                              # RootFileName = cms.untracked.string("STAMuon_62X_2019.root"),
                              StandAloneTrackCollectionLabel = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
                              GlobalTrackCollectionLabel     = cms.InputTag("globalMuons",""),
                              SegmentsTrackAssociatorParameters = cms.PSet(
                                  segmentsDt = cms.untracked.InputTag("dt4DSegments"),
                                  SelectedSegments = cms.untracked.InputTag("SelectedSegments"),
                                  segmentsCSC = cms.untracked.InputTag("cscSegments")
                              ),
)


process.p = cms.Path(process.demo)
