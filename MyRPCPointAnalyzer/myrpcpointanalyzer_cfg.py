import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
# process.load("Geometry.CMSCommonData.cmsExtendedGeometryPostLS1XML_cfi")
process.load("Configuration.Geometry.GeometryExtended2015_cff")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_71_V6::All'

# process.GlobalTag.globaltag = 'MC_31X_V1::All'
# process.GlobalTag.globaltag = 'GR_P_V32::All'





process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        ###  MC  ####
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_6_0_1_PostLS1v1/src/cmsDriverCommands/SingleMuPt100_cfi_RECO.root'
        ### MWGR1 ###
        # MWGR1 :: old reco
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_7_1_1/src/MyFilters/MyME4SegmentFilter/MWGR1_run222608_Skimmed.root',
        # MWGR1 :: new reco
        'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_7_1_1/src/MyFilters/MyME4SegmentFilter/MWGR1_run222608_NewFullReco_Skimmed.root',
        # 'file:/afs/cern.ch/user/p/piet/work/Analysis/CMSSW_7_1_1/src/MyCmsDriverCommands/run222608_newfullreco/run222608_RECO_ls0001.root'    
    )
)

process.rpcPointProducer = cms.EDProducer('RPCPointProducer',

    incldt = cms.untracked.bool(True),
    inclcsc = cms.untracked.bool(True),
    incltrack =  cms.untracked.bool(False),
                                          
    debug = cms.untracked.bool(False),
                                          
    rangestrips = cms.untracked.double(4.),
    rangestripsRB4 = cms.untracked.double(4.),
    MinCosAng = cms.untracked.double(0.85),
    MaxD = cms.untracked.double(120.0),
    MaxDrb4 = cms.untracked.double(150.0),
    ExtrapolatedRegion = cms.untracked.double(0.6), # in stripl/2 in Y and stripw*nstrips/2 in X
    # cscSegments = cms.InputTag('dTandCSCSegmentsinTracks','SelectedCscSegments','OwnParticles'),
    # dt4DSegments = cms.InputTag('dTandCSCSegmentsinTracks','SelectedDtSegments','OwnParticles'),
    cscSegments = cms.InputTag("cscSegments"),
    dt4DSegments = cms.InputTag("dt4DSegments"),                                          
    tracks = cms.InputTag("standAloneMuons"),
    TrackTransformer = cms.PSet(
          DoPredictionsOnly = cms.bool(False),
          Fitter = cms.string('KFFitterForRefitInsideOut'),
          TrackerRecHitBuilder = cms.string('WithTrackAngle'),
          Smoother = cms.string('KFSmootherForRefitInsideOut'),
          MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
          RefitDirection = cms.string('alongMomentum'),
          RefitRPCHits = cms.bool(False),
          Propagator = cms.string('SmartPropagatorAnyRKOpposite')
    )
)


process.demo = cms.EDAnalyzer('MyRPCPointAnalyzer',

      debug        = cms.untracked.bool(False),
      # RootFileName = cms.untracked.string("MyRPCPointHistograms_run222608_ls0001To0013_FirstReco.root"),
      RootFileName = cms.untracked.string("MyRPCPointHistograms_run222608_ls0001To0013_SecndReco.root"),
      cscSegments  = cms.untracked.InputTag('cscSegments'),
      dt4DSegments = cms.untracked.InputTag('dt4DSegments'),
      rpcRecHits   = cms.untracked.InputTag("rpcRecHits"),
      rpcDTPoints  = cms.untracked.InputTag("rpcPointProducer","RPCDTExtrapolatedPoints"),
      rpcCSCPoints = cms.untracked.InputTag("rpcPointProducer","RPCCSCExtrapolatedPoints"),
                              
)

# process.p = cms.Path(process.rpcPointProducer)
process.p = cms.Path(process.rpcPointProducer*process.demo)

