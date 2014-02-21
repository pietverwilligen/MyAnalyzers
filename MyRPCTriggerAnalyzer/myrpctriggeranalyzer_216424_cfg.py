import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 216424
        'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216424_1435374E-E449-E311-AD72-0025901E3D64.root',
        'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216424_620101F2-E349-E311-BB7B-003048F0E1CC.root',
        # 216425
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216425_D0D06494-C849-E311-B030-003048F23D74.root',
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216425_E84A85EB-C849-E311-8FCB-002481E0DC70.root',
        # 216473
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216473_26E5E9B9-6D49-E311-BACE-0030486780E4.root',
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216473_402B6169-6E49-E311-BAF1-003048D2BF26.root',
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216473_4CFC12AD-6E49-E311-ACB6-003048D2BEE6.root',
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216473_66E057B1-6D49-E311-9CFC-003048C90C08.root',
        # 'file:/afs/cern.ch/user/p/piet/work/public/GRIN_216473_6C7EC9AE-6E49-E311-ADCF-003048D2C020.root',
    )
)

process.demo = cms.EDAnalyzer('MyRPCTriggerAnalyzer',
                              GTReadoutRcd     = cms.InputTag("gtDigis"),
                              GMTReadoutRcd    = cms.InputTag("gtDigis" ),
                              RootFileName = cms.untracked.string("RPCTrigger_GRIN_216424.root"),
                              Debug        = cms.untracked.bool(False),

)


process.p = cms.Path(process.demo)
