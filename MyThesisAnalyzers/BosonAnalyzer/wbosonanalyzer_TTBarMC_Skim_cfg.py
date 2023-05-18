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

process.load('SandBox/Skims/RA2Objects_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# --- USING * in POOLSOURCE -----------------------------------------------
# -------------------------------------------------------------------------
prefix = 'file:'
# dirname = ''
# dirname = '/user/piet/PATTuples/A_Loose_QCD_Pt80_Spring10/'
# prefix = 'dcap://'
dirname = '/data/dalfonso/RA2_38X/moreSTAT/'
filenames = glob.glob(dirname + 'SUSYPAT_TTbar_D6T*.root')
INPUT_FILES = []
# print "List of input files for local submission:  \n"
# print filenames
# print "----------------------------------------------------"
for i in range(0,len(filenames)):
        a = prefix + str(filenames[i])
	INPUT_FILES.append(a)
	#local_outfilename = str(datasetname) + '-madgraphjob-flatntuple.root'

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(INPUT_FILES)
)
                

# process.source = cms.Source("PoolSource",
#     # replace 'myfile.root' with the source file you want to use
#     fileNames = cms.untracked.vstring(
#     # 'file:/user/piet/CMSSW_Versions/Photon/DATA/PhotonJet_Pt300_Spring10.root'
#     # 'file:/tmp/piet/DATA/PhotonRun2010B-Nov4ReReco_v1_1_1_b25.root'    # @lxplus 437
#     # 'file:/tmp/piet/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'     # @lxplus 445
#    'rfio:/castor/cern.ch/user/p/piet/Data/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'
#     )
# )

process.goodVertices = cms.EDFilter(
	  "VertexSelector",
	  filter = cms.bool(False),
	  src = cms.InputTag("offlinePrimaryVertices"),
	  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)


process.phanlzr = cms.EDAnalyzer('BosonAnalyzer',
                              Debug = cms.vint32(0,0,0,0),  # Timing | Matching | Selecting | Luminosity
                              Boson = cms.int32(24),        # PDGId Boson: Y=22 | Z=23 | W=24
                              Data = cms.bool(True),
			      PtCut = cms.double(000.0),
			      AmountOfJets = cms.int32(2),
                              RootFileName = cms.string("SUSYPAT_TTbar_D6T.root"),
                              # GenParticles = cms.string("prunedGenParticles"),
                              GenParticles = cms.string("genParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              RecoPhotons = cms.string("patPhotons"),
                              # RecoPhotons = cms.string("photons"),
                              RecoJets = cms.string("patJetsAK5PF"),
                              # RecoJets = cms.string("ak5PFJets"),
                              )

process.lumiJSON = cms.EDAnalyzer("JsonLogger",
                               fileName = cms.untracked.string("mariamaria_json.txt") 
                              )

process.p = cms.Path(process.goodVertices * process.ra2PFMuons * process.ra2PFElectrons * process.phanlzr * process.lumiJSON)
