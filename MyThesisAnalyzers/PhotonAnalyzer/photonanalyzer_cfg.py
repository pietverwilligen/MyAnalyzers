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

process.load('SandBox.Skims.beamHaloFilter_cfi')
process.load('SandBox.Skims.trackingFailureFilter_cfi')
process.ra2PostCleaning = cms.Sequence(
    process.beamHaloFilter+
    process.trackingFailureFilter
)



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
     # replace 'myfile.root' with the source file you want to use
     fileNames = cms.untracked.vstring(
     'file:/user/piet/CMSSW_Versions/Photon/DATA/susypat_353_2_GON.root'
     # 'file:/user/piet/CMSSW_Versions/Photon/DATA/GJets_HT-200_susypat_20events.root'
     )
 )

process.phanlzr = cms.EDAnalyzer('PhotonAnalyzer',
                              Debug = cms.vint32(0,0,0,0),  # Timing | Matching | Selecting | Luminosity
                              Data = cms.bool(True),
                              SpikeCleaning = cms.bool(True),
                              EGMIsolation = cms.bool(False),
                              AmountOfJets = cms.int32(2),
                              RootFileName = cms.string("Photon.root"),
                              GenParticles = cms.string("prunedGenParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              RecoPhotons = cms.string("patPhotons"),
                              # RecoPhotons = cms.string("photons"),
                              RecoJets = cms.string("patJetsAK5PF"),
                              )

process.lumiJSON = cms.EDAnalyzer("JsonLogger",
                              fileName = cms.untracked.string("mariamaria_json.txt") 
                              )
process.p = cms.Path(process.ra2PostCleaning * process.phanlzr * process.lumiJSON)
