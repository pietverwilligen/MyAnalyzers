import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/user/piet/CMSSW_Versions/Photon/DATA/PhotonJets_Pt100to200-madgraph_Spring10_1_1_T0L.root'
    )
)

process.demo = cms.EDAnalyzer('Acceptance',
                              RootFileName = cms.string("Acceptance.root")
)


process.p = cms.Path(process.demo)
