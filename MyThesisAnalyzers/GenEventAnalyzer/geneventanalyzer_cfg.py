import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/tmp/piet/DATA/PhotonVJets_7TeV-madgraph_Fall10_1_1_YZk.root',
        'file:/tmp/piet/DATA/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'
    )
)

process.demo = cms.EDAnalyzer('GenEventAnalyzer',
                              Debug = cms.bool(False), 
                              Zinv = cms.bool(True),
                              Pat = cms.bool(False),
                              HTCut= cms.double(0.0),
                              RootFileName = cms.string("GenEvent.root"),
                              GenParticles = cms.string("genParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              RecoJets = cms.string("ak5PFJets"),
                              PatJets = cms.string("patJetsAK5PF"),
                              )

process.p = cms.Path(process.demo)
