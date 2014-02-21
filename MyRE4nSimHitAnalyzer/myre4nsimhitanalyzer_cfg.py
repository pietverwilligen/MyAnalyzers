import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometryPostLS1XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_upscope_cfi") ### overwrite muon geometry
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
            # =====   CMSSW 6 0 0 SLHC   =====
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_0_SLHCtkpre1/src/cmsDriverCommands/FourMuPt_1_200_cfi_GEN_SIM.root'
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_0_SLHCtkpre1/src/cmsDriverCommands/RE_DIGI.root'
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_0_SLHCtkpre1/src/cmsDriverCommands/dud_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
            # =====   CMSSW 6 1 0 pre3   =====
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_1_0_pre3/src/cmsDriverCommands/FourMuPt_1_200_cfi_GEN_SIM.root'
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_1_0_pre3/src/cmsDriverCommands/SingleMuPt100_cfi_GEN_SIM.root'
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_1_0_pre3/src/cmsDriverCommands/1.0_ProdMinBias+ProdMinBias+DIGIPROD1+RECOPROD1/step1.root'
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_1_0_pre3/src/cmsDriverCommands/SingleMuPt100_cfi_GEN_SIM_25evts.root'
            # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_1_0_pre3/src/cmsDriverCommands/SingleMuPt100_cfi_GEN_SIM.root'
            # =====   CMSSW 6 2 0 SLHC1  =====
            'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_2_0_SLHC1/src/Geometry/RPCGeometry/test/SingleMuPt100_cfi_GEN_SIM.root'
            )
)

process.demo = cms.EDAnalyzer('MyRE4nSimHitAnalyzer',
                              RootFileName = cms.untracked.string("MyRE4nSimHistograms.root"),
)


process.p = cms.Path(process.demo)
