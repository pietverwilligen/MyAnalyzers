import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_1_PostLS1v1/src/cmsDriverCommands/SingleMuPt100_cfi_RECO.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC5/src/MyCmsDriverCommands/SingleMuPt100_1p6_2p4_cfi_RECO.root'
        # 'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyCmsDriverCommands/SingleMuPt100_cfi_1p6_2p5_RECO.root'
        'file:/build/piet/Upgrade/Eta_2p4_Releases/CMSSW_6_2_0_SLHC7/src/MyAnalyzers/MyRecHitAnalyzer/MyRecHits.root'
    )
)

process.demo = cms.EDAnalyzer('MyRecHitAnalyzer',
                              RootFileName = cms.untracked.string("MyRecHitHistograms.root"),

)


process.p = cms.Path(process.demo)
