import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

###################################################################################
### VarParsing options ### !!! the order of this configuration file matters !!! ###
###################################################################################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# add a list of strings for events to process
options.register ('eventsToProcess',
                    '',
                    VarParsing.multiplicity.list,
                    VarParsing.varType.string,
                    "Events to process")
options.register ('maxSize',
                    0,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "Maximum (suggested) file size (in Kb)")
options.parseArguments()

process = cms.Process("PickEvent")
process.source = cms.Source ("PoolSource",
      fileNames = cms.untracked.vstring (options.inputFiles),
)

if options.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange (options.eventsToProcess)
    
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32 (options.maxEvents)
)

if options.maxSize:
    process.Out.maxSize = cms.untracked.int32 (options.maxSize)

# process.Out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string (options.outputFile)
# )
    
# process.end = cms.EndPath(process.Out)


###################################################################################
### Do here the rest of your configuration                                      ###
###################################################################################

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Configuration.Geometry.GeometryExtended_cff')
# process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')
# process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.demo = cms.EDAnalyzer('MyNeutronSimHitAnalyzer',
                              # ---------
                              PdfFileNameBase = cms.untracked.string("MyNeutronSimHistograms_13TeV_2023Muon_xs_nc_eta7_62X"),
                              RootFileName = cms.untracked.string("MyNeutronSimHistograms_13TeV_2023Muon_xs_nc_eta7_62X.root"),
                              # ---------
                              BunchSpacing = cms.untracked.double(25.0),
                              COMEnergy = cms.untracked.double(13.0),
                              PhysicsDebug = cms.untracked.bool(False),
                              TechnicDebug = cms.untracked.bool(False),
                              )


process.p = cms.Path(process.demo)
