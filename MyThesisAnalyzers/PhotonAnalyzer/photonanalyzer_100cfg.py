import FWCore.ParameterSet.Config as cms
import glob
import sys

process = cms.Process("Phanlzr")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START38_V14::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# --- USING * in POOLSOURCE -----------------------------------------------
# -------------------------------------------------------------------------
### mX home dir
###------------
prefix = 'file:'
# dirname = ''
dirname = '/user/piet/MC2010/Fall10MC/GJets_HT-100To200_madgraph_Fall10/'
### mX SE
###------
# prefix = 'dcap://'
# dirname = '/pnfs/iihe/cms/store/user/piet/PAT/EGRun2010A-Sep17ReReco_v3/'
# dirname = '/pnfs/iihe/cms/store/user/piet/Fall10MC/GJets_HT-200_madgraph_Fall10/'
filenames = glob.glob(dirname + 'susypat*.root')
INPUT_FILES = []
#print "List of input files for local submission:  \n"
#print filenames
#print "----------------------------------------------------"
for i in range(0,len(filenames)):
    a = prefix + str(filenames[i])
    INPUT_FILES.append(a)
    #local_outfilename = str(datasetname) + '-madgraphjob-flatntuple.root'

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(INPUT_FILES)
)
                    
# -------------------------------------------------------------------------
# process.source = cms.Source("PoolSource",
#     # replace 'myfile.root' with the source file you want to use
#     fileNames = cms.untracked.vstring(
# 'file:/user/piet/CMSSW_Versions/Photon/DATA/PhotonJet_Pt300_Spring10.root'
#     )
# )

process.load('SandBox.Skims.RA2Leptons_cff')
process.ra2PFElectrons = cms.Sequence(
                  process.patElectronsPFID *
                                process.patElectronsPFIDIso
                                )
process.ra2PFElectronVeto = cms.Sequence(
        ~process.countPFElectronsIDIso
            )


process.phanlzr = cms.EDAnalyzer('PhotonAnalyzer',
                              Debug = cms.vint32(0,1,0,0),  # Timing | Matching | Selecting | Luminosity
                              Data = cms.bool(False),
                              SpikeCleaning = cms.bool(False),
                              EGMIsolation = cms.bool(False),
                              AmountOfJets = cms.int32(2),
                              RootFileName = cms.string("GJets_HT-100To200_madgraph_Fall10.root"),
                              GenParticles = cms.string("prunedGenParticles"),
                              # GenParticles = cms.string("genParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              RecoPhotons = cms.string("patPhotons"),
                              # RecoPhotons = cms.string("photons"),
                              RecoJets = cms.string("patJetsAK5PF"),
                              )

process.lumiJSON = cms.EDAnalyzer("JsonLogger",
                              fileName = cms.untracked.string("mariamaria_json.txt") 
                              )
process.p = cms.Path(process.ra2PFElectrons * process.ra2PFElectronVeto * process.phanlzr * process.lumiJSON)
