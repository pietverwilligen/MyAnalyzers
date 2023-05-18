import FWCore.ParameterSet.Config as cms

process = cms.Process("RecoMuon")
# Messages
#process.load("RecoMuon.Configuration.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

# Muon Reco
process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")
process.load("RecoMuon.Configuration.RecoMuon_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

# old CMSSW 31X test on 31X RelVal:
# '/store/relval/CMSSW_3_1_0_pre7/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0004/B89FA4AB-CC41-DE11-8348-000423D98B28.root'
# process.load("Configuration.StandardSequences.Geometry_cff")
# process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.GlobalTag.globaltag = 'IDEAL_31X::All'

process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')                                                                                                                                      
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V54::All', '') # 74X data taking
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V50::All', '') # 73X data taking


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:MyFilteredEvents_243484_Collisions.root'
        )
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RecoMuons_243484_Collisions.root')
)

process.p = cms.Path(process.muonrecoComplete)
process.this_is_the_end = cms.EndPath(process.out)




