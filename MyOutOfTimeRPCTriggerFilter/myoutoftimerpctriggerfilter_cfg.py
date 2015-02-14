import FWCore.ParameterSet.Config as cms
process = cms.Process("Filter")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(
'/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00097CE4-A2B2-E411-8FEB-02163E011D8A.root',
'/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/0049B70A-8EB2-E411-8979-02163E011D1C.root',
'/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/004DEB4A-9DB2-E411-AB38-02163E011D06.root',
'/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00D04605-A5B2-E411-A496-02163E0129F4.root',
'/store/express/Commissioning2015/ExpressCosmics/FEVT/Express-v1/000/234/029/00000/00D21ADC-A2B2-E411-8861-02163E011D57.root',
# RelVal ZMM
# '/store/relval/CMSSW_7_3_0/RelValZMM_13/GEN-SIM-RECO/MCRUN2_73_V9_71XGENSIM_FIXGT-v1/00000/0AC5D8E4-5DA2-E411-8B1D-0025905B8590.root',
# '/store/relval/CMSSW_7_3_0/RelValZMM_13/GEN-SIM-RECO/MCRUN2_73_V9_71XGENSIM_FIXGT-v1/00000/AC3E5968-5CA2-E411-B77E-0025905938AA.root',
# '/store/relval/CMSSW_7_3_0/RelValZMM_13/GEN-SIM-RECO/MCRUN2_73_V9_71XGENSIM_FIXGT-v1/00000/C69684E2-5DA2-E411-ADF0-0025905964C2.root',


)
)
process.Out = cms.OutputModule("PoolOutputModule",
SelectEvents = cms.untracked.PSet(
SelectEvents = cms.vstring('path')
),
# fileName = cms.untracked.string ("MyFilteredEvents_234029_BX0_PointingMuons.root")
fileName = cms.untracked.string ("MyFilteredEvents_234029_BX0_NoDTSegments.root")
)
process.filter = cms.EDFilter('MyOutOfTimeRPCTriggerFilter',
Debug              = cms.untracked.bool(False),
GTReadoutRcd       = cms.InputTag("gtDigis"),
GMTReadoutRcd      = cms.InputTag("gtDigis" ),
STAMuonTrackCollectionLabel = cms.InputTag("standAloneMuons",""),         # UpdatedAtVtx
TrackerTrackCollectionLabel = cms.InputTag("ctfWithMaterialTracksP5",""), # globalMuons or globalCosmicMuons or globalCosmicMuons1Leg
# TrackerTrackCollectionLabel = cms.InputTag("generalTracks",""), # for normal PP / ZMM MC
SelectBX           = cms.untracked.bool(True),
# config 1 :: RPC in bx=1, DT in bx=0
bxRPC             = cms.untracked.int32(0),
bxDT              = cms.untracked.int32(0),
# config 2 :: RPC in bx=0, DT in bx=1
# bxRPC              = cms.untracked.int32(0),
# bxDT               = cms.untracked.int32(1),
AnalyzeTRK         = cms.untracked.bool(False),
SelectTRK          = cms.untracked.bool(False),
SelectAND          = cms.untracked.bool(False),
SelectOR           = cms.untracked.bool(False),
SelectNoDTSegments = cms.untracked.bool(True),
# RootFileName       = cms.untracked.string("MyFilteredHistograms_234029_BX0_PointingMuons.root"),
RootFileName       = cms.untracked.string("MyFilteredHistograms_234029_BX0_NoDTSegments.root"),
)
process.path = cms.Path(process.filter)
process.end = cms.EndPath(process.Out)
