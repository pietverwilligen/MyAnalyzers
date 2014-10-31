import FWCore.ParameterSet.Config as cms

process = cms.Process("Timing")

# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
# process.load("Configuration.Geometry.GeometryExtended2023_cff")
# process.load("Configuration.Geometry.GeometryExtended2023Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2015Reco_RPC4RE11_cff")
# process.load("Configuration.Geometry.GeometryExtended2023RPCEtaUpscope_cff")
# process.load('Configuration.Geometry.GeometryExtended2023RPCEtaUpscopeReco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2023', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'MC_72_V1', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")
# process.source = cms.Source("PoolSource",
#   # replace 'myfile.root' with the source file you want to us
#   fileNames = cms.untracked.vstring(
#    'file:myfile.root'
#   )
# )
process.timing = cms.EDAnalyzer('RPCTiming'
)


process.p = cms.Path(process.timing)
