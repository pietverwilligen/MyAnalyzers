import FWCore.ParameterSet.Config as cms

process = cms.Process("Analyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")

# process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDevReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
# CSCGeometry depends on alignment ==> necessary to provide GlobalPositionRecord
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi") 
# process.load("Geometry.CSCGeometry.cscGeometry_cfi")

# process.load("Geometry.CSCGeometry.cscGeometry_cfi")
# process.load("Geometry.DTGeometry.dtGeometry_cfi")
# process.load("Geometry.GEMGeometry.gemGeometry_cfi")

# process.maxEvents and process.source are included in cfgs below
# -------------------------------------------------------------------------------------------------------------------
# process.load("SourceFiles_DYToMuMu_14TeV_HGCALGS_PU140_1ns_500um_1cm_NeutrBkg_7ns_dPhi0p010_dEta0p050_SP4_v2_RECO")
# -------------------------------------------------------------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(
    input           = cms.untracked.int32(-1),
)
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
                            secondaryFileNames = secFiles,
                            # eventsToProcess = cms.untracked.VEventRange('1:1157:115618'), # Run 1, Event 115618, LumiSection 1157
                            # eventsToProcess = cms.untracked.VEventRange('1:2010:200906'), # Run 1, Event 200906, LumiSection 2010
                            # eventsToProcess = cms.untracked.VEventRange('1:22779:540983'), # Run 1 LS 22779 Evt 540983           
                            # eventsToProcess = cms.untracked.VEventRange('1:26043:618522'), # Run 1 LS 26043 Evt 618522           
                            # eventsToProcess = cms.untracked.VEventRange('1:26044:618528'), # Run 1 LS 26044 Evt 618528           
                           )
readFiles.extend(['file:NuGun_RECO_500Events_1ns.root'])
# -------------------------------------------------------------------------------------------------------------------



process.me0timeanalyzer = cms.EDAnalyzer('MyME0InTimePUAnalyzer',
                              # ----------------------------------------------------------------------
                              RootFileName       = cms.untracked.string("ME0InTimeOutOfTimePUtHistograms_NuGun.root"),
                              InvestigateOnlyME0 = cms.untracked.bool(True),    # require at least one signal muon in 2.0 < | eta | < 2.8
                              # ----------------------------------------------------------------------
                              me0DetResX         = cms.untracked.double(1.0),    # [in cm] (single layer resolution)
                              me0DetResY         = cms.untracked.double(5.0),    # [in cm] (single layer resolution)
                              cscDetResX         = cms.untracked.double(0.0150), # [in cm] (chamber resolution :: 75-150um, take here 150um)
                              cscDetResY         = cms.untracked.double(5.0),    # [in cm]
                              dtDetResX          = cms.untracked.double(0.0400), # [in cm] (chamber resolution ::  75-125um in r-phi, take here 125um) ==> seems to fail to often ... take 0.4 mm now
                              dtDetResY          = cms.untracked.double(0.5000), # [in cm] (chamber resolution :: 150-400um in r-z  , take here 400um) ==> seems to fail to often ... take 5.0 mm now
                              nMatchedHitsME0Seg = cms.untracked.int32(3),
                              nMatchedHitsCSCSeg = cms.untracked.int32(3),
                              nMatchedHitsDTSeg  = cms.untracked.int32(6),
                              matchQualityME0    = cms.untracked.double(0.75),   # Percentage of matched hits (matched hits / total hits) >= 75% ==> 3/3 or 3/4, 4/5, 5/6
                              matchQualityReco   = cms.untracked.double(0.75),   # not using this for now ... problem with DTs 2B solved first
                              # ----------------------------------------------------------------------
                              printInfoHepMC           = cms.untracked.bool(False),
                              printInfoSignal          = cms.untracked.bool(False),
                              printInfoPU              = cms.untracked.bool(False),
                              printInfoVtx             = cms.untracked.bool(False),
                              printInfoAll             = cms.untracked.bool(False),
                              printInfoME0Match        = cms.untracked.bool(False),
                              printInfoMuonMatch       = cms.untracked.bool(False),
                              printInfoMuonMatchDetail = cms.untracked.bool(False),
                              printInfoInTimePU        = cms.untracked.bool(False),
                              # ----------------------------------------------------------------------

)

process.p = cms.Path(process.me0timeanalyzer)
