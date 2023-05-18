import FWCore.ParameterSet.Config as cms
import glob
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


prefix = 'dcap://'
dirname = '/pnfs/iihe/cms/store/user/piet/MCFactors/GJets_HT-100To200_madgraph_Fall10/'
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
#    # replace 'myfile.root' with the source file you want to use
#     fileNames = cms.untracked.vstring(
#         'file:/tmp/piet/DATA/PhotonVJets_7TeV-madgraph_Fall10_1_1_YZk.root',
#         # 'file:/tmp/piet/DATA/ZinvisibleJets_7TeV-madgraph_Fall10_1_1_gqg.root'
#     )
# )

process.demo = cms.EDAnalyzer('GenEventAnalyzer',
                              Debug = cms.bool(False), 
                              Zinv = cms.bool(False),
                              Pat = cms.bool(False),
                              HTCut= cms.double(0.0),
                              RootFileName = cms.string("PhotonGenEvent_100To200_Jet.root"),
                              GenParticles = cms.string("genParticles"),
                              GenJets = cms.string("ak5GenJets"),
                              RecoJets = cms.string("ak5PFJets"),
                              PatJets = cms.string("patJetsAK5PF"),
                              )

process.p = cms.Path(process.demo)
