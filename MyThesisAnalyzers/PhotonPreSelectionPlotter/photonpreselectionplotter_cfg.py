import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

# process.load("RecoEgamma.Configuration.RecoEgamma_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V14::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")




process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    # 'file:/user/piet/CMSSW_Versions/Photon/CMSSW_3_8_6_patch1/src/DATA/PhotonRun2010B-Nov4ReReco_v1_1_1_b25.root'
    # 'file:/afs/cern.ch/user/p/piet/scratch0/CMSSW_Versions/Photon/DATA/PhotonRun2010B-Nov4ReReco_v1_10Events.root'
    # 'file:/tmp/piet/DATA/PhotonRun2010B-Nov4ReReco_v1_1_1_b25.root' # @ lxplus 437
    'rfio:/castor/cern.ch/user/p/piet/Data/PhotonJets_Pt100to200-madgraph_Spring10_1_1_T0L.root'
    )
)

process.demo = cms.EDAnalyzer('PhotonPreSelectionPlotter'
)


process.p = cms.Path(process.demo)
