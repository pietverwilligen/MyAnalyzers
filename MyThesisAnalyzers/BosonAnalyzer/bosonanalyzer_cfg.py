import FWCore.ParameterSet.Config as cms
import glob
import sys

process = cms.Process("Phanlzr")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V14::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
     # replace 'myfile.root' with the source file you want to use
     fileNames = cms.untracked.vstring(
     # 'file:/user/piet/CMSSW_Versions/Photon/DATA/PhotonJet_Pt300_Spring10.root'
     # 'file:/tmp/piet/DATA/PhotonRun2010B-Nov4ReReco_v1_1_1_b25.root'    # @lxplus 437
     # 'file:/tmp/piet/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'     # @lxplus 445
    'rfio:/castor/cern.ch/user/p/piet/Data/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'
     )
 )

process.phanlzr = cms.EDAnalyzer('BosonAnalyzer',
                              Debug = cms.vint32(0,0,0,0),  # Timing | Matching | Selecting | Luminosity
                              Boson = cms.int32(23),        # PDGId Boson: Y=22 | Z=23 | W=24
                              Data = cms.bool(False),
                              RootFileName = cms.string("Boson.root"),
                              # GenParticles = cms.string("prunedGenParticles"),
                              GenParticles = cms.string("genParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              # RecoPhotons = cms.string("patPhotons"),
                              RecoPhotons = cms.string("photons"),
                              # RecoJets = cms.string("patJetsAK5PF"),
                              RecoJets = cms.string("ak5PFJets"),
                              )

process.lumiJSON = cms.EDAnalyzer("JsonLogger",
                              fileName = cms.untracked.string("mariamaria_json.txt") 
                              )
process.p = cms.Path(process.phanlzr * process.lumiJSON)
