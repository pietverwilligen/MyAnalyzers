import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    # 'file:myfile.root',
    # 'file:/user/piet/CMSSW_Versions/Photon/DATA/PAT_PhotonJets_Pt40to100.root',
    # 'file:/user/piet/CMSSW_Versions/Photon/DATA/PAT_PhotonJets_Pt100to200.root',
    # 'file:/user/piet/CMSSW_Versions/Photon/DATA/PAT_QCDFlat_Pt15to3000.root',
    # 'file:/user/piet/CMSSW_Versions/Photon/CMSSW_3_6_3/src/PATtifier/PAT_Test.root',
    # 'file:/user/piet/CMSSW_Versions/Photon/DATA/QCD_Summer10_23_1_dUA.root',
    # 'file:/tmp/piet/DATA/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'
    'file:/tmp/piet/DATA/PhotonVJets_7TeV-madgraph_Fall10_1_1_YZk.root'
    )
)

process.demo = cms.EDAnalyzer('RootFileScout',
    # Debug = cms.bool(False)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printList = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("genParticles"),
     maxEventsToPrint = cms.untracked.int32(-1)
)
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
     printP4 = cms.untracked.bool(False),
     printPtEtaPhi = cms.untracked.bool(False),
     printVertex = cms.untracked.bool(True),
     printStatus = cms.untracked.bool(False),
     printIndex = cms.untracked.bool(False),
     status = cms.untracked.vint32(1, 2, 3)
)

process.p = cms.Path(
    process.demo #*
    # process.printList *
    # process.printTree
    )
