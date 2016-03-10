import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("RecoMuon.Configuration.RecoMuon_cff")
# process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load("Configuration.StandardSequences.Geometry_cff")

# process.load('Configuration.Geometry.GeometryExtended2023_cff')
# process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_dataRun2_16Dec2015_v0', '')   # 763 data reprocessing

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
        'file:/afs/cern.ch/work/p/piet/Analysis/SLC6/Data2015/CMSSW_7_6_3/src/MyData/StaMuAtEta6_AOD.root',
        # 'file:/afs/cern.ch/work/p/piet/Analysis/SLC6/Data2015/CMSSW_7_6_3/src/MyData/StaMuAtEta6_RECO.root',
    )
)

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.demo = cms.EDAnalyzer('MyStandAloneMuonAnalyzer',
                              MuonServiceProxy,
                              PhysicsDebug = cms.untracked.bool(True),
                              TechnicDebug = cms.untracked.bool(True),
                              GenPartDebug = cms.untracked.bool(False),
                              MCTruthMatching = cms.untracked.bool(False),
                              TightID      = cms.untracked.bool(False),
                              LooseID      = cms.untracked.bool(True),
                              RootFileName = cms.untracked.string("STAMuon_763_2015.root"),
                              # MuonLabel1   = cms.untracked.string("standAloneMuons"),
                              # MuonLabel2   = cms.untracked.string("UpdatedAtVtx"),
                              # Maybe try to pass by input tag ???
                              StandAloneTrackCollectionLabel = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
                              GlobalTrackCollectionLabel     = cms.InputTag("globalMuons",""),
                              SegmentsTrackAssociatorParameters = cms.PSet(
                                  segmentsDt = cms.untracked.InputTag("dt4DSegments"),
                                  SelectedSegments = cms.untracked.InputTag("SelectedSegments"),
                                  segmentsCSC = cms.untracked.InputTag("cscSegments")
                              ),
)


process.p = cms.Path(process.demo)
