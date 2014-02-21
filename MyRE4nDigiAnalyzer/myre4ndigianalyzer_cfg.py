import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

# special stuff from Ianna
# process.load("Configuration.StandardSequences.GeometryDB_cff")
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# from Configuration.AlCa.autoCond import autoCond
# process.GlobalTag.globaltag = autoCond['mc']
# process.XMLFromDBSource.label=''
# process.GlobalTag.toGet = cms.VPSet(cms.PSet(record = cms.string('GeometryFileRcd'),
#                                              tag = cms.string('XMLFILE_Geometry_61YV1_ExtendedPostLS1_mc'),
#                                              connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_GEOMETRY")
#                                              )
#                                     )

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_upscope_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometryPostLS1XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
# process.load("Geometry.CSCGeometry.cscGeometry_cfi")
# process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_1/src/cmsDriverCommands/SingleMuPt100_cfi_DIGI.root'
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_1/src/cmsDriverCommands/SingleMuPt100_cfi_RECO.root'
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_2_0_SLHC1/src/RE11_To_RE41_Configs/Roumyana/SingleMuPt100_cfi_py_GEN_SIM_DIGI.root'
        'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_2_0_SLHC1/src/RE11_To_RE41_Configs/Roumyana/SingleMuPt100_cfi_py_GEN_SIM_DIGI_10evts.root'
    )
)

process.demo = cms.EDAnalyzer('MyRE4nDigiAnalyzer',
                              # DATA
                              # DigiLabel= cms.untracked.string("muonRPCDigis")
                              # MONTE-CARLO
                              DigiLabel    = cms.untracked.string("simMuonRPCDigis"),
                              # DigiLabel    = cms.untracked.string("muonRPCDigis"),
                              RootFileName = cms.untracked.string("MyRE4nDigiHistograms.root"),
)


process.p = cms.Path(process.demo)
