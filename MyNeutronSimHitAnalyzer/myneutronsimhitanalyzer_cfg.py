import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCUpscopeXML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023RPCEtaUpscopeXML_cfi")
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# process.load('Configuration.Geometry.GeometryExtended_cff')
process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')
process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
            # 'file:/build/civanch/CMSSW_7_0_0_pre12/src/mb_8TeV_mu_xs.root'
            # ---------
            # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_test.root'
            # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_wide.root'
            'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_ext.root'
            # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_2015.root'
            )
)

process.demo = cms.EDAnalyzer('MyNeutronSimHitAnalyzer',
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre12"),                              
                              # RootFileName    = cms.untracked.string("MyNeutronSimHistograms_pre12.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13"),                              
                              # RootFileName    = cms.untracked.string("MyNeutronSimHistograms_pre13.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13_wide"),                              
                              # RootFileName    = cms.untracked.string("MyNeutronSimHistograms_pre13_wide.root"),
                              # ---------
                              PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13_ext"),                              
                              RootFileName    = cms.untracked.string("MyNeutronSimHistograms_pre13_ext.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13_2015"),                              
                              # RootFileName    = cms.untracked.string("MyNeutronSimHistograms_pre13_2015.root"),                              
                              Debug           = cms.untracked.bool(False),
)


process.p = cms.Path(process.demo)
