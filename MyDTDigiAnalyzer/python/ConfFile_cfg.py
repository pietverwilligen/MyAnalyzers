import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    # fileNames = cms.untracked.vstring(
    #     'file:myfile.root'
    # )
    ################### RAWDIGI1
    # fileNames = cms.untracked.vstring('/store/group/upgrade/muon/RecoFolder/DYToMuMu_2019_2Step_2/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RAW_CMSSW_6_2_0_SLHC23patch1_2019Scenario_2Step_GEMSH/23d4646c3e8a6be200238397ea8208ad/step2_894_1_4Uu.root'),
    ################### RAWDIGI2
    fileNames = cms.untracked.vstring('/store/group/upgrade/muon/RecoFolder/DYToMuMu_2019_2Step/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/calabria_DYToMuMu_GEN-SIM-DIGI-RAW_CMSSW_6_2_0_SLHC23patch1_2019Scenario_2Step/03a26b85950c71e0fe99659a0a214ac3/step2_894_1_57H.root'),
    ################### Filippo                                                                                                                                                                                                                                        
    eventsToProcess = cms.untracked.VEventRange('1:356:29989',),

)

process.demo = cms.EDAnalyzer('MyDTDigiAnalyzer'
)


process.p = cms.Path(process.demo)
