import FWCore.ParameterSet.Config as cms
import glob
import sys

process = cms.Process("Phanlzr")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/user/piet/CMSSW_Versions/Photon/DATA/PhotonJet_Pt300_Spring10.root'
    )
)

process.phanlzr = cms.EDAnalyzer('PhotonAnalyzer',
                              Debug = cms.vint32(0,1,0,0),  # Timing | Matching | Selecting | Luminosity
                              Data = cms.bool(False),
                              SpikeCleaning = cms.bool(False),
                              EGMIsolation = cms.bool(False),
                              RootFileName = cms.string("mSugra1.root"),
                              GenParticles = cms.string("genParticles"),
                              # GenParticles = cms.string("prunedGenParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              RecoPhotons = cms.string("patPhotons"),
                              # RecoPhotons = cms.string("photons"),
                              RecoJets = cms.string("patJetsAK5PF"),
                              )

process.lumiJSON = cms.EDAnalyzer("JsonLogger",
                              fileName = cms.untracked.string("mariamaria_json.txt") 
                              )
process.p = cms.Path(process.ra2PostCleaning * process.phanlzr * process.lumiJSON)
