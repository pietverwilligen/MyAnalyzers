import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Configuration.Geometry.GeometryExtended_cff')
# process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')
# process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

# process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
# process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
# process.load("Geometry.GEMGeometry.gemGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2500) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25) )

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
                # 'file:/build/civanch/CMSSW_7_0_0_pre12/src/mb_8TeV_mu_xs.root'
                # ---------
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_test.root'
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_wide.root'
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_ext.root'
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_2015.root'
                # ---------
                # 'file:/afs/cern.ch/user/p/piet/work/NeutronBackgroundFiles/mb_8TeV_mu_hp_2015_eta7.root'
                # 'file:/afs/cern.ch/user/p/piet/work/NeutronBackgroundFiles/mb_8TeV_mu_xs_2015_eta7a.root'
                # 'file:/afs/cern.ch/user/p/piet/work/NeutronBackgroundFiles/mb_8TeV_mu_xs_2012_eta7.root'
                # 'file:/afs/cern.ch/user/p/piet/work/NeutronBackgroundFiles/mb_8TeV_mu_xs_extended_eta7.root'
                # ---------
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_2015_eta7a.root'
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_hp_2015_eta7.root'
                # 'file:/build/civanch/CMSSW_7_0_0_pre13/src/mb_8TeV_mu_xs_2012_eta7.root'
                # 'file:/build/piet/CMSSW_7_0_0_pre12/src/MyCmsDriverCommands/mb_8TeV_mu_xs_2012_eta7.root'
                # 'file:/build/piet/another_pre12/CMSSW_7_0_0_pre12/src/MyCmsDriverCommands/mb_8TeV_mu_xs_2012_eta8.root'    
                # 'file:/build/piet/CMSSW_7_0_0_pre13/src/MyCmsDriverCommands/SingleMuPt100_2015_cfi_GEN-SIM.root'
                # 'file:/build/piet/CMSSW_7_0_0_pre13/src/MyCmsDriverCommands/ZmumuJets_Pt20_300_2015_cfi_GEN-SIM.root'
                # 'file:/build/piet/yet_another_pre13/CMSSW_7_0_0_pre13/src/MyCmsDriverCommands/mb_8TeV_mu_xs_2015_500ns.root'
                # ---------
                # 'file:/build/piet/CMSSW_7_1_0_pre8/src/MyCmsDriverCommands/mb_13TeV_2015_xs_nc_eta7.root'
                # 'file:/build/piet/CMSSW_7_1_0_pre8/src/MyCmsDriverCommands/mb_8TeV_2015_xs_nc_eta7.root'
                # 'file:/build/piet/CMSSW_7_1_0_pre8/src/MyCmsDriverCommands/mb_8TeV_2015_xs_oc_eta7.root'
                # 'file:/build/piet/CMSSW_7_1_0_pre8/src/MyCmsDriverCommands/mb_8TeV_2015_xs_nc_eta8.root'
                # ---------
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_1_1_r62.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_2_1_0kg.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_3_1_RTo.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_4_1_YIU.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_5_1_CKk.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_6_1_TRj.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_7_1_xn2.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_8_2_bGm.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_9_1_39n.root',
                'file:/lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/mb_13TeV_2023Muon_mu_xs_10_1_gka.root',
     )
)

process.demo = cms.EDAnalyzer('MyNeutronSimHitAnalyzer',
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre12"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_pre12.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_pre13.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13_wide"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_pre13_wide.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13_ext"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_pre13_ext.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_pre13_2015"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_pre13_2015.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2015"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2015.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_hp_2015"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_hp_2015.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2012"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2012.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2012_Cavern"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2012_Cavern.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2012_Cavern_Eta8"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2012_Cavern_Eta8.root"),                              
                              # PdfFileNameBase = cms.untracked.string("MySingleMuonSimHistograms_8TeV_xs_2015"),
                              # RootFileName = cms.untracked.string("MySingleMuonSimHistograms_8TeV_xs_2015.root"),
                              # PdfFileNameBase = cms.untracked.string("MyZMuMuJetsSimHistograms_8TeV_xs_2015"),
                              # RootFileName = cms.untracked.string("MyZMuMuJetsSimHistograms_8TeV_xs_2015.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2015_500ns"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_xs_2015_500ns.root"),
                              # ---------
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_2015_xs_nc_eta7_71X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_2015_xs_nc_eta7_71X.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_2015_xs_nc_eta7_71X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_2015_xs_nc_eta7_71X.root"),
                              # PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_8TeV_2015_xs_oc_eta7_71X"),
                              # RootFileName = cms.untracked.string("MyNeutronSimHistograms_8TeV_2015_xs_oc_eta7_71X.root"),
                              # ---------
                              PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_2023Muon_xs_nc_eta7_62X"),
                              RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_2023Muon_xs_nc_eta7_62X.root"),
                              # ---------
                              BunchSpacing = cms.untracked.double(25.0),
                              COMEnergy = cms.untracked.double(13.0),
                              PhysicsDebug = cms.untracked.bool(False),
                              TechnicDebug = cms.untracked.bool(False),
                              )


process.p = cms.Path(process.demo)
