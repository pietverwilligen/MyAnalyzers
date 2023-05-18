import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometryPostLS1XML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_1_PostLS1v1/src/cmsDriverCommands/SingleMuPt100_cfi_RECO.root'
    )
)

process.demo = cms.EDAnalyzer('MyRE4nRecHitAnalyzer',
                              RootFileName = cms.untracked.string("MyRE4nRecHitHistograms.root"),

)


process.p = cms.Path(process.demo)
