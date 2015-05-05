// -*- C++ -*-
//
// Package:    MyAnalyzers/MyOutOfTimeRPCTriggerFilter
// Class:      MyOutOfTimeRPCTriggerFilter
// 
/**\class MyOutOfTimeRPCTriggerFilter MyOutOfTimeRPCTriggerFilter.cc MyAnalyzers/MyOutOfTimeRPCTriggerFilter/plugins/MyOutOfTimeRPCTriggerFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen
//         Created:  Thu, 12 Feb 2015 11:16:03 GMT
//
//


// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// root include files
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// L1 Trigger
#include <DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h>
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include <DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h>

// Magnetic Field and Tracks
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Muon Detector Geometries
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
// #include <Geometry/DTGeometry/interface/DTGeometry.h>
// #include <Geometry/CSCGeometry/interface/CSCGeometry.h>
// #include <Geometry/Records/interface/MuonGeometryRecord.h>
// #include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
// #include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>
// #include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
// #include <Geometry/CommonDetUnit/interface/GeomDet.h>

// RPC Digis and Rechits
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>

// DT RecHits and Segments
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>

// STA and GLB Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Math
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class MyOutOfTimeRPCTriggerFilter : public edm::EDFilter {
   public:
      explicit MyOutOfTimeRPCTriggerFilter(const edm::ParameterSet&);
      ~MyOutOfTimeRPCTriggerFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::InputTag m_gtReadoutLabel;
  edm::InputTag m_gmtReadoutLabel;
  edm::InputTag STAMuLabel, TRACKLabel;
  bool debug;
  bool selectBX, selectTRK;
  bool selectOR, selectAND;
  bool selectNoDTSegments;
  bool selectNoRPCRechits;
  bool selectDTbutNoRPCTrig;
  bool selectRPCbutNoDTTrig;
  bool selectTRKTrack;
  bool selectSTATrack;
  bool selectMUOTrack;
  bool doFilter;
  bool analyzeTRK;
  int select_bx_rpc;
  int select_bx_dt;
  int eventsProcessed;
  int eventsFiltered;

  std::string rootFileName;
  TFile * outputfile;
  TH1F * RPCb_Triggers_Quality;
  TH1F * RPCb_Triggers_ETA_All, * RPCb_Triggers_ETA_Q0, * RPCb_Triggers_ETA_Q1, * RPCb_Triggers_ETA_Q2, * RPCb_Triggers_ETA_Q3;
  TH1F * RPCb_Triggers_PHI_All, * RPCb_Triggers_PHI_Q0, * RPCb_Triggers_PHI_Q1, * RPCb_Triggers_PHI_Q2, * RPCb_Triggers_PHI_Q3;
  TH2F * RPCb_Triggers_ETA_PHI_All;
  TH1F * DTTF_Triggers_Quality;
  TH1F * DTTF_Triggers_ETA_All, * DTTF_Triggers_ETA_Q0, * DTTF_Triggers_ETA_Q1, * DTTF_Triggers_ETA_Q2, * DTTF_Triggers_ETA_Q3;
  TH1F * DTTF_Triggers_PHI_All, * DTTF_Triggers_PHI_Q0, * DTTF_Triggers_PHI_Q1, * DTTF_Triggers_PHI_Q2, * DTTF_Triggers_PHI_Q3;
  TH2F * DTTF_Triggers_ETA_PHI_All;
  TH2F * ETA_Triggers_RPCb_DTTF_All;
  TH2F * PHI_Triggers_RPCb_DTTF_All;
  TH2F * BX_Triggers_RPCb_DTTF_All;

  TH1F * TrackCollections;
  TH1F * D0_TrackExp, * D0_TrackFnd, *D0_TrackEff;
  TH1F * DZ_TrackExp, * DZ_TrackFnd, *DZ_TrackEff;
  TH1F * D0_TrackTrk, * DZ_TrackTrk;
  TH2F * D0DZ_TrackExp, * D0DZ_TrackFnd, * D0DZ_TrackEff; 
  TH2F * D0DZ_TrackTrk;

};

//
// constants, enums and typedefs
//
int n_bx       = 7;
int n_phi      = 144;
int n_eta      = 48;
int n_qua      = 11;
double n_bx_1  = -3.5;
double n_bx_2  = +3.5;
double n_phi_1 = 0.0000;
double n_phi_2 = 6.2832;
double n_eta_1 = -2.40;
double n_eta_2 = +2.40;
double n_eta_exact = 19;
double n_eta_vec[] = {-2.45, -2.4, -2.1, -1.6, -1.2, -1.05, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.05, 1.2, 1.6, 2.1, 2.4, 2.45};
// double n_eta_exact = 64;
// double n_eta_vec[] = {-2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2, -1.95, -1.9, -1.85, -1.8, -1.75, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45};
double n_qua_1 = -0.5;
double n_qua_2 = 10.5;
//
// static data member definitions
//

//
// constructors and destructor
//
MyOutOfTimeRPCTriggerFilter::MyOutOfTimeRPCTriggerFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  m_gtReadoutLabel     = iConfig.getParameter<edm::InputTag>("GTReadoutRcd");
  m_gmtReadoutLabel    = iConfig.getParameter<edm::InputTag>("GMTReadoutRcd");
  debug                = iConfig.getUntrackedParameter<bool>("Debug");
  selectBX             = iConfig.getUntrackedParameter<bool>("SelectBX");
  select_bx_rpc        = iConfig.getUntrackedParameter<int>("bxRPC");
  select_bx_dt         = iConfig.getUntrackedParameter<int>("bxDT");
  analyzeTRK           = iConfig.getUntrackedParameter<bool>("AnalyzeTRK");
  selectTRK            = iConfig.getUntrackedParameter<bool>("SelectTRK");
  selectAND            = iConfig.getUntrackedParameter<bool>("SelectAND");
  selectOR             = iConfig.getUntrackedParameter<bool>("SelectOR");
  selectNoDTSegments   = iConfig.getUntrackedParameter<bool>("SelectNoDTSegments");
  selectNoRPCRechits   = iConfig.getUntrackedParameter<bool>("SelectNoRPCRechits"); 
  selectDTbutNoRPCTrig = iConfig.getUntrackedParameter<bool>("SelectDTbutNoRPCTrig");
  selectRPCbutNoDTTrig = iConfig.getUntrackedParameter<bool>("SelectRPCbutNoDTTrig");
  selectTRKTrack       = iConfig.getUntrackedParameter<bool>("SelectTRKTrack");
  selectSTATrack       = iConfig.getUntrackedParameter<bool>("SelectSTATrack");
  selectMUOTrack       = iConfig.getUntrackedParameter<bool>("SelectMUOTrack");
  doFilter             = iConfig.getUntrackedParameter<bool>("DoFilter");
  rootFileName         = iConfig.getUntrackedParameter<std::string>("RootFileName");
  STAMuLabel           = iConfig.getParameter<edm::InputTag>("STAMuonTrackCollectionLabel");
  TRACKLabel           = iConfig.getParameter<edm::InputTag>("TrackerTrackCollectionLabel");
}


MyOutOfTimeRPCTriggerFilter::~MyOutOfTimeRPCTriggerFilter()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MyOutOfTimeRPCTriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  ++eventsProcessed;

  int evNum = (iEvent.id()).event();
  int rnNum = (iEvent.id()).run();
  int lsNum = (iEvent.id()).luminosityBlock();

  int bx_firstcand_rpcb = -10, bx_firstcand_dttf = -10;
  bool found_first_cand_rpcb = false, found_first_cand_dttf = false;
  bool keepEventBX = false, keepEventTRK = false;
  bool keepEventNoDTSegments = false;
  bool keepEventNoRPCRechits = false;
  bool keepEventRPCNoDTTrigger = false;
  bool keepEventDTNoRPCTrigger = false;
  bool keepEventTRKTrack = false;
  bool keepEventSTATrack = false;
  bool keepEventMUOTrack = false;
  bool keepEvent = false;
  int quality_first_rpcb = -10, quality_first_dttf = -10; 
  double eta_first_rpcb  = -10, eta_first_dttf = -10; 
  double phi_first_rpcb  = -10, phi_first_dttf = -10; 

  bool GlobalTrackExpectedInEvent = false;
  bool GlobalTrackFoundInEvent = true;

  // ==========================
  // === Define Collections ===
  // ==========================
  // Triggers
  edm::Handle<L1MuGMTReadoutCollection> pCollection;
  iEvent.getByLabel(m_gmtReadoutLabel,pCollection);
  const L1MuGMTReadoutCollection * gmtRC = pCollection.product();
  std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr;
  std::vector<L1MuGMTReadoutRecord> gmt_records = gmtRC->getRecords();
    
  // DT Segments
  edm::Handle<DTRecSegment4DCollection> dtSegmentCollection;
  iEvent.getByLabel("dt4DSegments", dtSegmentCollection);
  // edm::Handle<DTRecHit1DPair> dtRecHitCollection;
  // iEvent.getByLabel("dt1DRecHits", dtRecHitCollection);

  // Muons & Tracks
  edm::Handle<reco::TrackCollection> staTracks;
  iEvent.getByLabel(STAMuLabel, staTracks);
  edm::Handle<reco::TrackCollection> trkTracks;
  iEvent.getByLabel(TRACKLabel, trkTracks);
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByLabel("muons", recoMuons);

  // Magnetic Field & Geometries
  edm::ESHandle<MagneticField> theMagneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagneticField);
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);


  // ==================================
  // === get GMT readout collection ===
  // ==================================
  if(debug) { std::cout<<"SELECT BX"<<std::endl; }
  if(selectBX) {
    for( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
      // int BxInEvent = RRItr->getBxInEvent();
      // int BxInEventNew = RRItr->getBxNr();
      // int nrpcB = 0;
      // int ndtB = 0;
      std::vector<L1MuRegionalCand> BrlRpcCands = RRItr->getBrlRPCCands();
      std::vector<L1MuRegionalCand> BrlDtCands = RRItr->getDTBXCands ();
      if(BrlRpcCands.size() != 0 && BrlDtCands.size() == 0) { keepEventRPCNoDTTrigger = true; }
      if(BrlRpcCands.size() == 0 && BrlDtCands.size() != 0) { keepEventDTNoRPCTrigger = true; }
      std::vector<L1MuRegionalCand>::const_iterator RCItr;
      // RPC barrel muon candidates
      for( RCItr = BrlRpcCands.begin(); RCItr !=BrlRpcCands.end(); ++RCItr) {
	if ( !(*RCItr).empty() ) {
	  // Print Information
	  if(debug) {
	    std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | ";
	    std::cout<<"RPCb Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<" bx = "<<RCItr->bx()<<std::endl;
	  }
	  // Select First Candidate and Fill Some Plots --- for now only from first sorted trigger candidate
	  if(!found_first_cand_rpcb) { 
	    // Essential for the filter
	    bx_firstcand_rpcb = RCItr->bx(); found_first_cand_rpcb = true; 
	    // For Filling some histograms later on
	    quality_first_rpcb = RCItr->quality();
	    eta_first_rpcb = RCItr->etaValue();
	    phi_first_rpcb = RCItr->phiValue();

	    // See whether we can match a Stand Alone Muon with those triggers ...
	    reco::MuonCollection::const_iterator  recoMuon;
	    for (recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); ++recoMuon) {
	      if(!recoMuon->isStandAloneMuon()) continue;
	      if(debug) std::cout<<"Trying to match to Stand Alone Muon :: dR = "<<deltaR(eta_first_rpcb,phi_first_rpcb,recoMuon->eta(),recoMuon->phi())<<std::endl;
	      if(deltaR(eta_first_rpcb,phi_first_rpcb,recoMuon->eta(),recoMuon->phi())<0.5) {
		if(debug) {
		  std::cout<<" Matched Stand Alone Muon :: pt = "<<recoMuon->pt()<<" eta = "<<recoMuon->eta()<<" phi = "<<recoMuon->phi();
		  // Piotr Traczyk ::
		  // For the timing calculation OutsideIn and InsideOut really mean along/opposite to the momentum direction. 
		  // So if an outside-in cosmic muon would be reconstructed by the _collision_ algorithm in the upper hemisphere, then the OutsideIn timing calculation would make sense. 
		  // For cosmics reconstructed with the _cosmic_ algorithm, I believe InsideOut is always the one to use.
		  // ---> Cosmics reconstructed with _collision_ algorithm
		  // if(recoMuon->time().direction()<0)  std::cout<<" direction = "<<("OutsideIn")<<" time at IP = "<<recoMuon->time().timeAtIpOutIn<<" +/- "<<recoMuon->time().timeAtIpOutInErr<<" ns";
		  // ---> Cosmics reconstructed with _cosmic_ algorithm:
		  if(recoMuon->time().direction()<0)  std::cout<<" direction = "<<("OutsideIn")<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" +/- "<<recoMuon->time().timeAtIpInOutErr<<" ns";
		  // For cosmics in lower hemisphere nothing changes
		  if(recoMuon->time().direction()>0)  std::cout<<" direction = "<<("InsideOut")<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" +/- "<<recoMuon->time().timeAtIpInOutErr<<" ns";
		  if(recoMuon->time().direction()==0) std::cout<<" direction = "<<("Undefined")<<" time at IP = "<<("Undefined")<<" ns";
		  // std::cout<<" direction = "<<(recoMuon->time().direction()<0?("OutsideIn"):("InsideOut"))<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" ns";
		  // std::cout<<" # segments = "<<"";
		  std::cout<<std::endl;
		  // std::cout<<" d0: "<<recoMuon->outerTrack()->d0()<<" +/- "<<recoMuon->outerTrack()->d0Error()<<" cm dz: "<<recoMuon->outerTrack()->dz()<<" +/- "<<recoMuon->outerTrack()->dzError()<<" cm"<<std::endl; 
		}
		for(trackingRecHit_iterator recHit = recoMuon->outerTrack()->recHitsBegin(); recHit != recoMuon->outerTrack()->recHitsEnd(); ++recHit) {
		  const GeomDet* geomDet = theTrackingGeometry->idToDet((*recHit)->geographicalId());
		  // double r = geomDet->surface().position().perp();
		  double r = geomDet->surface().position().mag();
		  double z = geomDet->toGlobal((*recHit)->localPosition()).z();
		  DetId detid = DetId((*recHit)->geographicalId());
		  int rpcrechits = 0, dtrechits = 0;
		  if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::RPC) {
		    ++rpcrechits;
		    if(debug) {                                                         
		      std::cout<<"RPC Tracking RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm in DetId = "<<detid.rawId();
		      // std::cout<<""<<std::endl; 
		    }
		    // RPC Geometry
		    edm::ESHandle <RPCGeometry> rpcGeom;
		    iSetup.get<MuonGeometryRecord>().get(rpcGeom);
		    // RPC Digis
		    // edm::Handle<RPCDigiCollection> rpcdigis;
		    // iEvent.getByLabel("muonRPCDigis", rpcdigis);
		    // RPCDigiCollection::DigiRangeIterator digiRpc;
		    // Check Digis
		    /*
		    for(digiRpc=rpcdigis->begin(); digiRpc!=rpcdigis->end(); ++digiRpc){
		      if((*digiRpc).first.rawId() != detid.rawId()) continue;
		      RPCDetId detId=(*digiRpc).first;
		      RPCDigiCollection::const_iterator digiRpcItr;
		      std::cout<<"DetId = "<<detId.rawId()<<" digis: "<<std::endl;
		      for (digiRpcItr =(*digiRpc ).second.first; digiRpcItr != (*digiRpc).second.second; ++digiRpcItr){
			int strip= (*digiRpcItr).strip();
			int bx=(*digiRpcItr).bx();
			std::cout<<" RPC Digi: strip = "<<std::noshowpos<<std::setw(2)<<strip<<" bx = "<<std::showpos<<bx<<std::endl;
		      }
		    }
		    */
		    // Print again Tracking RecHits
		    // if(debug) {                                                         
		    //   std::cout<<"RPC Tracking RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm in DetId = "<<detid.rawId();
		    // }
		    int bx_rpc = -10;
		    // RPC Rechits
		    edm::Handle<RPCRecHitCollection> rpcRecHits;
		    iEvent.getByLabel("rpcRecHits","",rpcRecHits);
		    RPCRecHitCollection::const_iterator recHitRpc;
		    // Check Rechits
		    for (recHitRpc = rpcRecHits->begin(); recHitRpc != rpcRecHits->end(); ++recHitRpc) {
		      if((*recHitRpc).rpcId().rawId() != detid.rawId()) continue;
		      RPCDetId rollId = (RPCDetId)(*recHitRpc).rpcId();
		      const RPCRoll* rollasociated = rpcGeom->roll(rollId);
		      const BoundPlane & RPCSurface = rollasociated->surface();
		      double r_rpc = RPCSurface.toGlobal(recHitRpc->localPosition()).mag();
		      double z_rpc = RPCSurface.toGlobal(recHitRpc->localPosition()).z();
		      bx_rpc = -10;
		      // if(debug) std::cout<<"RPC RecHit at "<<"r: "<< r_rpc <<" cm"<<" z: "<<z_rpc<<" cm in DetId = "<<rollId.rawId()<<" with bx = "<<(*recHitRpc).BunchX();
		      if(fabs(r-r_rpc) < 10 && fabs(z-z_rpc) < 10 ){
			if(bx_rpc==-10) bx_rpc = (*recHitRpc).BunchX(); // take the first one ... they are ordered in time ... first comes bx -2,-1,0,1,2,3
			// if(debug) {
			//   std::cout<<" ==> Matched RPC RecHit with "<<"dr: "<< fabs(r-r_rpc) <<" cm"<<" dz: "<<fabs(z-z_rpc)<<" cm in DetId = "<<rollId.rawId();
			//   std::cout<<" with bx = "<<(*recHitRpc).BunchX()<<" first strip = "<<(*recHitRpc).firstClusterStrip()<<" clustersize = "<<(*recHitRpc).clusterSize();
			// }
		      }
		    }
		    if(debug) {
		      if(bx_rpc != -10) std::cout<<" bx = "<<bx_rpc<<std::endl;
		      else std::cout<<""<<std::endl;
		    }
		  }                                                                    
		  if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::DT) {
		    if(debug) {  
		      ++dtrechits;                                                      
		      // std::cout<<"DT Tracking RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm in DetId = "<<detid.rawId(); 
		      // // std::cout<<" and BX = "<<(*recHit)->bx();                     
		      // std::cout<<""<<std::endl;
		    }
		  }
		}
	      }
	    } 
	  }
	}
      }
      // DT barrel muon candidates
      for( RCItr = BrlDtCands.begin(); RCItr !=BrlDtCands.end(); ++RCItr) {
	if ( !(*RCItr).empty() ) {
	  if(debug) {
	    std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | ";
	    std::cout<<"DTTF Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<" bx = "<<RCItr->bx()<<std::endl;
	  }
	  // Select First Candidate and Fill Some Plots --- for now only from first sorted trigger candidate
	  if(!found_first_cand_dttf) { 
	    // Essential for the filter
	    bx_firstcand_dttf = RCItr->bx(); found_first_cand_dttf = true; 
	    // For Filling some histograms later on
	    quality_first_dttf = RCItr->quality();
	    eta_first_dttf = RCItr->etaValue();
	    phi_first_dttf = RCItr->phiValue();

	    // See whether we can match a Stand Alone Muon with those triggers ...
	    reco::MuonCollection::const_iterator  recoMuon;
	    for (recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); ++recoMuon) {
	      if(!recoMuon->isStandAloneMuon()) continue;
	      if(debug) std::cout<<"Trying to match to Stand Alone Muon :: dR = "<<deltaR(eta_first_dttf,phi_first_dttf,recoMuon->eta(),recoMuon->phi())<<std::endl;
	      if(deltaR(eta_first_dttf,phi_first_dttf,recoMuon->eta(),recoMuon->phi())<0.5) {
		if(debug) {
		  std::cout<<" Matched Stand Alone Muon :: pt = "<<recoMuon->pt()<<" eta = "<<recoMuon->eta()<<" phi = "<<recoMuon->phi();
		  std::cout<<" direction = "<<(recoMuon->time().direction()<0?("OutsideIn"):("InsideOut"))<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" ns";
		  std::cout<<std::endl;
		  // std::cout<<" d0: "<<recoMuon->outerTrack()->d0()<<" +/- "<<recoMuon->outerTrack()->d0Error()<<" cm dz: "<<recoMuon->outerTrack()->dz()<<" +/- "<<recoMuon->outerTrack()->dzError()<<" cm"<<std::endl; 
		}
	      }
	    } 
	  }
	}
      }
      // RPC endcap muon candidates
      // CSC endcap muon candidates    
      
      // We are only interessed in the first candidate (trigger candidates are sorted on quality and then on pt)
      // this is not going to work ... vectors are keeping space for triggers at -2, -1, 0, 1, 2 and bx is already assigned ...
      // RCItr = BrlRpcCands.begin();
      // bx_firstcand_rpcb = RCItr->bx(); std::cout<<"bx_firstcand_rpcb = "<<bx_firstcand_rpcb<<std::endl;
      // RCItr = BrlDtCands.begin();
      // bx_firstcand_dttf = RCItr->bx(); std::cout<<"bx_firstcand_dttf = "<<bx_firstcand_dttf<<std::endl;

      /*
      reco::TrackCollection::const_iterator staTrack;
      for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack) {
	if(debug) std::cout<<"Stand Alone Muon :: pt = "<<staTrack->pt()<<" eta = "<<staTrack->eta()<<" phi = "<<staTrack->phi()<<std::endl;
	for(trackingRecHit_iterator recHit = staTrack->recHitsBegin(); recHit != staTrack->recHitsEnd(); ++recHit) {
	  const GeomDet* geomDet = theTrackingGeometry->idToDet((*recHit)->geographicalId());
	  double r = geomDet->surface().position().perp();
	  double z = geomDet->toGlobal((*recHit)->localPosition()).z();
	  DetId detid = DetId((*recHit)->geographicalId());
	  if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::RPC) {
	    if(debug) {
	      std::cout<<"RPC RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm";
	      // std::cout<<" and BX = "<<(*recHit)->bx();
	      std::cout<<""<<std::endl;
	    }
	  }
	  if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::DT) {
	    if(debug) {
	      std::cout<<"DT RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm";
	      // std::cout<<" and BX = "<<(*recHit)->bx();
	      std::cout<<""<<std::endl;
	    }
	  }
	}
      }
      */
    }
    // Fill Some histograms for the BX selected or for the missing triggers
    if((selectBX && bx_firstcand_dttf == select_bx_dt && bx_firstcand_rpcb == select_bx_rpc) || (selectRPCbutNoDTTrig && keepEventRPCNoDTTrigger) || (selectDTbutNoRPCTrig && keepEventDTNoRPCTrigger)) {
      RPCb_Triggers_Quality->Fill(quality_first_rpcb);
      RPCb_Triggers_ETA_PHI_All->Fill(eta_first_rpcb, phi_first_rpcb);
      RPCb_Triggers_ETA_All->Fill(eta_first_rpcb);
      RPCb_Triggers_PHI_All->Fill(phi_first_rpcb);
      switch (quality_first_rpcb) {
      case 0: {
	RPCb_Triggers_ETA_Q0->Fill(eta_first_rpcb);
	RPCb_Triggers_PHI_Q0->Fill(phi_first_rpcb);
      } break;
      case 1: {
	RPCb_Triggers_ETA_Q1->Fill(eta_first_rpcb);
	RPCb_Triggers_PHI_Q1->Fill(phi_first_rpcb);
      } break;
      case 2: {
	RPCb_Triggers_ETA_Q2->Fill(eta_first_rpcb);
	RPCb_Triggers_PHI_Q2->Fill(phi_first_rpcb);
      } break;
      case 3: {
	RPCb_Triggers_ETA_Q3->Fill(eta_first_rpcb);
	RPCb_Triggers_PHI_Q3->Fill(phi_first_rpcb);
      } break;
      }
      DTTF_Triggers_Quality->Fill(quality_first_dttf);
      DTTF_Triggers_ETA_PHI_All->Fill(eta_first_dttf, phi_first_dttf);
      DTTF_Triggers_ETA_All->Fill(eta_first_dttf);
      DTTF_Triggers_PHI_All->Fill(phi_first_dttf);
      // Fix this ... if necessary ... DT quality codes are different w.r.t. RPC
      switch (quality_first_dttf) {
      case 0: {
	DTTF_Triggers_ETA_Q0->Fill(eta_first_dttf);
	DTTF_Triggers_PHI_Q0->Fill(phi_first_dttf);
      } break;
      case 1: {
	DTTF_Triggers_ETA_Q1->Fill(eta_first_dttf);
	DTTF_Triggers_PHI_Q1->Fill(phi_first_dttf);
      } break;
      case 2: {
	DTTF_Triggers_ETA_Q2->Fill(eta_first_dttf);
	DTTF_Triggers_PHI_Q2->Fill(phi_first_dttf);
      } break;
      case 3: {
	DTTF_Triggers_ETA_Q3->Fill(eta_first_dttf);
	DTTF_Triggers_PHI_Q3->Fill(phi_first_dttf);
      } break;
      }
      ETA_Triggers_RPCb_DTTF_All->Fill(eta_first_rpcb,eta_first_dttf);
      PHI_Triggers_RPCb_DTTF_All->Fill(phi_first_rpcb,phi_first_dttf);
      BX_Triggers_RPCb_DTTF_All->Fill(bx_firstcand_rpcb,bx_firstcand_dttf);
    }
  }



  // =====================================
  // === Filter Events with 0 Segments ===
  // =====================================
  if(debug) { std::cout<<"SELECT NO DT SEGMENTS"<<std::endl; }
  if(selectNoDTSegments) {
    for( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
      std::vector<L1MuRegionalCand> BrlDtCands = RRItr->getDTBXCands();
      std::vector<L1MuRegionalCand>::const_iterator RCItr;
      for(RCItr = BrlDtCands.begin(); RCItr !=BrlDtCands.end(); ++RCItr) {
	if ( !(*RCItr).empty() ) {
	  if(debug) {
	    std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | ";
	    std::cout<<"DTTF Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<" bx = "<<RCItr->bx()<<std::endl;
	  }
	  if(dtSegmentCollection->size()==0) keepEventNoDTSegments = true;
	  if(debug) std::cout<<"Amount of DT Segments in Event = "<<dtSegmentCollection->size()<<std::endl;
	  // DTRecSegment4DCollection::const_iterator segment;  
	  // for (segment = dtSegmentCollection->begin(); segment!=dtSegmentCollection->end(); ++segment){} 
	}
      }
    }
  }



  // ========================================
  // === Filter Events with 0 RPC Rechits ===
  // ========================================
  if(debug) { std::cout<<"SELECT NO RPC RECHITS"<<std::endl; }
  if(selectNoRPCRechits) {
    for( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
      std::vector<L1MuRegionalCand> BrlRpcCands = RRItr->getBrlRPCCands();
      std::vector<L1MuRegionalCand>::const_iterator RCItr;
      for( RCItr = BrlRpcCands.begin(); RCItr !=BrlRpcCands.end(); ++RCItr) {
	if ( !(*RCItr).empty() ) {
	  if(debug) {
	    std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | ";
	    std::cout<<"RPCb Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<" bx = "<<RCItr->bx()<<std::endl;
	  }
	  // Strategy :: Keep event as soon as there is a STA-Muon track within |eta| < 1.0 that has no RPC Rechits
	  reco::TrackCollection::const_iterator staTrack;
	  for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack) {
	    if(debug) std::cout<<"Stand Alone Muon :: pt = "<<staTrack->pt()<<" eta = "<<staTrack->eta()<<" phi = "<<staTrack->phi()<<std::endl;
	    for(trackingRecHit_iterator recHit = staTrack->recHitsBegin(); recHit != staTrack->recHitsEnd(); ++recHit) {
	      const GeomDet* geomDet = theTrackingGeometry->idToDet((*recHit)->geographicalId());
	      double r = geomDet->surface().position().perp();
	      double z = geomDet->toGlobal((*recHit)->localPosition()).z();
	      DetId detid = DetId((*recHit)->geographicalId());
	      int rpchitsontrack = 0; 
	      if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::RPC) {
		++rpchitsontrack;
		if(debug) std::cout<<"RPC RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm"<<std::endl;
	      }
	      if(rpchitsontrack==0 && fabs(staTrack->eta())<1.0) {
		keepEventNoRPCRechits = true;
	      }
	    }
	  }
	}
      }
    }
  }



  // =====================================
  // === Analysis of STA and GLB Muons ===
  // =====================================
  if(debug) { std::cout<<"ANALYZE / SELECT TRK"<<std::endl; }
  if(analyzeTRK || selectTRK || selectTRKTrack) {
    // iterators
    reco::TrackCollection::const_iterator staTrack;
    reco::TrackCollection::const_iterator trkTrack;
    reco::MuonCollection::const_iterator  recoMuon;
    if(debug) std::cout<<"Reconstructed STA Muon Tracks: " <<staTracks->size()<< std::endl;
    if(debug) std::cout<<"Reconstructed Tracker  Tracks: " <<trkTracks->size()<<std::endl;
    if(debug) std::cout<<"Reconstructed Muons: "<<recoMuons->size()<<std::endl;
    // Some Statistics of the Tracking:
    if(recoMuons->size()>0) TrackCollections->Fill(0);
    if(staTracks->size()>0) TrackCollections->Fill(1);
    if(trkTracks->size()>0) TrackCollections->Fill(2);

    
    
    // === STA Tracks ======================
    /*
    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack) {
      reco::TransientTrack track(*staTrack,&*theMagneticField,theTrackingGeometry);
      if(debug) {
	std::cout<<" Stand Alone Muon Track :: ";
	std::cout<<" p: "<<track.impactPointTSCP().momentum().mag();
	std::cout<<" pT: "<<track.impactPointTSCP().momentum().perp();
	std::cout<<" eta: "<<track.impactPointTSCP().momentum().eta();
	std::cout<<" phi: "<<track.impactPointTSCP().momentum().phi();
	std::cout<<" chi2: "<<track.chi2();
	// std::cout<<" d0: "<<staTrack->d0();
	// std::cout<<" dz: "<<staTrack->dz();
	std::cout<<" d0: "<<staTrack->d0()<<" +/- "<<staTrack->d0Error()<<" cm dz: "<<staTrack->dz()<<" +/- "<<staTrack->dzError()<<" cm";
	std::cout<<" pos: "<<track.impactPointTSCP().position()<<" cm";
	std::cout<<" with "<<staTrack->recHitsSize()<<" rechits"<<std::endl;
	std::cout<<"--------------------------------------------------------------"<<std::endl;
      }
      keepEventSTATrack = true;
      // double track_eta = track.impactPointTSCP().momentum().eta();
      // double track_phi = track.impactPointTSCP().momentum().phi();
      // double sta_eta = staTrack->eta();
      // double sta_phi = staTrack->phi();
      // double glb_eta;
      // double glb_phi;
      if(fabs(staTrack->d0()) < 80 && fabs(staTrack->dz()) < 100) {
	// in case you want to loop over all inner tracks
	// for (trkTrack = trkTracks->begin(); trkTrack!=trkTracks->end(); ++trkTrack) {
	// }
	// take this as criterion to keep the event
	keepEventTRK = true;
	GlobalTrackExpectedInEvent = true; D0_TrackExp->Fill(fabs(staTrack->d0())); DZ_TrackExp->Fill(fabs(staTrack->dz())); D0DZ_TrackExp->Fill(fabs(staTrack->dz()),fabs(staTrack->d0()));
	if(trkTracks->size()>0) {
	  GlobalTrackFoundInEvent = true;
	  D0_TrackFnd->Fill(fabs(staTrack->d0())); DZ_TrackFnd->Fill(fabs(staTrack->dz())); D0DZ_TrackFnd->Fill(fabs(staTrack->dz()),fabs(staTrack->d0())); 
	  D0_TrackTrk->Fill(fabs(staTrack->d0())); DZ_TrackTrk->Fill(fabs(staTrack->dz())); D0DZ_TrackTrk->Fill(fabs(staTrack->dz()),fabs(staTrack->d0())); 
	}
      }
    }
    // === Analysis ========================
    if(GlobalTrackExpectedInEvent) { 
      TrackCollections->Fill(3);                             
      // std::cout<<"filled bin 4"<<std::endl;
      if(GlobalTrackFoundInEvent) {
	TrackCollections->Fill(4); 
	// std::cout<<"filled bin 5"<<std::endl; 
      }
    }
    */

    // === TRK Tracks ======================
    for (trkTrack = trkTracks->begin(); trkTrack!=trkTracks->end(); ++trkTrack) {
      reco::TransientTrack track(*trkTrack,&*theMagneticField,theTrackingGeometry);
      if(debug) {
	std::cout<<" Tracker Track :: ";
	std::cout<<" p: "<<track.impactPointTSCP().momentum().mag();
	std::cout<<" pT: "<<track.impactPointTSCP().momentum().perp();
	std::cout<<" eta: "<<track.impactPointTSCP().momentum().eta();
	std::cout<<" chi2: "<<track.chi2();
	std::cout<<" d0: "<<trkTrack->d0()<<" +/- "<<trkTrack->d0Error()<<" cm dz: "<<trkTrack->dz()<<" +/- "<<trkTrack->dzError()<<" cm";
	// std::cout<<" pos: "<<track.impactPointTSCP().position()<<" cm";
	std::cout<<" with "<<trkTrack->recHitsSize()<<" rechits"<<std::endl;
	std::cout<<"--------------------------------------------------------------"<<std::endl;
      }
      keepEventTRKTrack = true;
    }
    // === All Muons =======================
    // const reco::Muon * myMuon = &(*recoMuons->begin());
    for (recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); ++recoMuon) {
      if(debug) {
	std::cout<<" Reco Muon :: pt = "<<recoMuon->pt()<<" eta = "<<recoMuon->eta()<<" phi= "<<recoMuon->phi();
	if(recoMuon->isTrackerMuon()) {
	  std::cout<<" d0: "<<recoMuon->innerTrack()->d0()<<" +/- "<<recoMuon->innerTrack()->d0Error()<<" cm dz: "<<recoMuon->innerTrack()->dz()<<" +/- "<<recoMuon->innerTrack()->dzError()<<" cm"; 
	}
	if(recoMuon->isStandAloneMuon()) {
	  std::cout<<" d0: "<<recoMuon->outerTrack()->d0()<<" +/- "<<recoMuon->outerTrack()->d0Error()<<" cm dz: "<<recoMuon->outerTrack()->dz()<<" +/- "<<recoMuon->outerTrack()->dzError()<<" cm"; 
	}
	if(recoMuon->isGlobalMuon()) {
	  std::cout<<" d0: "<<recoMuon->globalTrack()->d0()<<" +/- "<<recoMuon->globalTrack()->d0Error()<<" cm dz: "<<recoMuon->globalTrack()->dz()<<" +/- "<<recoMuon->globalTrack()->dzError()<<" cm"; 
	}
	std::cout<<" is Stand Alone Muon = "<<recoMuon->isStandAloneMuon()<<" is GlobalMuon = "<<recoMuon->isGlobalMuon();
	std::cout<<" is TrackerMuon = "<<recoMuon->isTrackerMuon()<<" is ParticleFlowMuon = "<<recoMuon->isPFMuon()<<std::endl;
	keepEventMUOTrack = true;	
      }
      // recoMuon->innerTrack();
      // recoMuon->outerTrack();
      // recoMuon->globalTrack();
      /*
	reco::TransientTrack muTT(recoMuon->outerTrack(), &*theMagneticField, theTrackingGeometry);
	TrajectoryStateOnSurface innerMuTSOS = muTT.innermostMeasurementState();
	TrajectoryStateOnSurface outerMuTSOS = muTT.outermostMeasurementState();
	GlobalPoint pointMuonIn  = innerMuTSOS.globalPosition();
	GlobalPoint pointMuonOut = outerMuTSOS.globalPosition();
      */
      // Take Stand Alone Muons from here ...
      if(recoMuon->isStandAloneMuon() && (fabs(recoMuon->outerTrack()->d0()) < 120 && fabs(recoMuon->outerTrack()->dz()) < 300)) {
	keepEventTRK = true;
	GlobalTrackExpectedInEvent = true; D0_TrackExp->Fill(fabs(recoMuon->outerTrack()->d0())); DZ_TrackExp->Fill(fabs(recoMuon->outerTrack()->dz())); 
	D0DZ_TrackExp->Fill(fabs(recoMuon->outerTrack()->dz()),fabs(recoMuon->outerTrack()->d0()));
	if(trkTracks->size()>0) {
	  GlobalTrackFoundInEvent = true;
	  D0_TrackFnd->Fill(fabs(recoMuon->outerTrack()->d0())); DZ_TrackFnd->Fill(fabs(recoMuon->outerTrack()->dz())); D0DZ_TrackFnd->Fill(fabs(recoMuon->outerTrack()->dz()),fabs(recoMuon->outerTrack()->d0())); 
	  D0_TrackTrk->Fill(fabs(recoMuon->outerTrack()->d0())); DZ_TrackTrk->Fill(fabs(recoMuon->outerTrack()->dz())); D0DZ_TrackTrk->Fill(fabs(recoMuon->outerTrack()->dz()),fabs(recoMuon->outerTrack()->d0())); 
	}
      }
    }
    // === Analysis ========================
    if(GlobalTrackExpectedInEvent) { 
      TrackCollections->Fill(3);                             
      // std::cout<<"filled bin 4"<<std::endl;
      if(GlobalTrackFoundInEvent) {
	TrackCollections->Fill(4); 
	// std::cout<<"filled bin 5"<<std::endl; 
      }
    }
  }



  // ======================================================
  // === Take decision here to keep or throw away Event ===
  // ======================================================
  if(debug) { std::cout<<"TAKE FINAL FILTER DECISION"<<std::endl; }
  // === Select based on RPC / DT Trigger BX ==============
  if(selectBX && bx_firstcand_dttf == select_bx_dt && bx_firstcand_rpcb == select_bx_rpc) { keepEventBX = true; }
  else keepEventBX = false;
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" first bx RPC Cand = "<<std::showpos<<bx_firstcand_rpcb<<" first bx DT Cand = "<<std::showpos<<bx_firstcand_dttf<<" ==> Event Kept (BX) :: "<<std::noshowpos<<keepEventBX<<std::endl;
  }
  // === Select based on STA Muon Tracks ==================
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" Event Kept (STA) :: "<<std::noshowpos<<keepEventTRK<<std::endl;
  }
  // === Select based on DT Segments ======================
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" Event Kept (BX) :: "<<std::noshowpos<<keepEventBX<<" Event Kept (DT)  :: "<<std::noshowpos<<keepEventNoDTSegments<<" ==> Event Kept :: "<<std::noshowpos<<keepEventBX*keepEventNoDTSegments<<std::endl;
  }
  // === Select based on RPC Rechits ======================
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" Event Kept (BX) :: "<<std::noshowpos<<keepEventBX<<" Event Kept (RPC) :: "<<std::noshowpos<<keepEventNoRPCRechits<<" ==> Event Kept :: "<<std::noshowpos<<keepEventBX*keepEventNoRPCRechits<<std::endl;
  }
  // === Select based on TRK Tracks ======================
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" Event Kept (TRK) :: "<<std::noshowpos<<keepEventTRKTrack<<std::endl;
  }
  // === Select based on STA Tracks ======================
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" Event Kept (STA) :: "<<std::noshowpos<<keepEventSTATrack<<std::endl;
  }
  // === Select based on MUO Tracks ======================
  if(debug) {
    std::cout<<"Run "<<std::setw(9)<<rnNum<<" Luminosity Block "<<std::setw(4)<<lsNum<<" Event "<<std::setw(12)<<evNum<<" || ";
    std::cout<<" Event Kept (MUO) :: "<<std::noshowpos<<keepEventMUOTrack<<std::endl;
  }
  // === Global Decision ==================================
  if(selectAND && ((selectBX && keepEventBX) && (selectTRK && keepEventTRK))) { keepEvent = true; ++eventsFiltered; } 
  if(selectOR  && ((selectBX && keepEventBX) || (selectTRK && keepEventTRK))) { keepEvent = true; ++eventsFiltered; }
  if(selectNoDTSegments && keepEventNoDTSegments) { keepEvent = true; ++eventsFiltered; }
  if(selectNoRPCRechits && keepEventNoRPCRechits)  { keepEvent = true; ++eventsFiltered; }
  if(selectRPCbutNoDTTrig && keepEventRPCNoDTTrigger) { keepEvent = true; ++eventsFiltered; }
  if(selectDTbutNoRPCTrig && keepEventDTNoRPCTrigger) { keepEvent = true; ++eventsFiltered; }
  if(selectSTATrack && keepEventSTATrack) { keepEvent = true; ++eventsFiltered; }
  if(selectTRKTrack && keepEventTRKTrack) { keepEvent = true; ++eventsFiltered; }
  if(selectMUOTrack && keepEventMUOTrack) { keepEvent = true; ++eventsFiltered; }

  if(debug) { 
    std::cout<<"||||||||||||||||||||||||||||||||||| ==> Event Kept :: "<<keepEvent<<std::endl;
  }
  if(doFilter) return keepEvent;
  else return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MyOutOfTimeRPCTriggerFilter::beginJob()
{
  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );
  RPCb_Triggers_Quality = new TH1F("RPCb_Triggers_Quality", "RPCb_Triggers_Quality", n_qua, n_qua_1, n_qua_2);
  RPCb_Triggers_PHI_All = new TH1F("RPCb_Triggers_PHI_All", "RPCb_Triggers_PHI_All", n_phi, n_phi_1, n_phi_2);
  RPCb_Triggers_PHI_Q0  = new TH1F("RPCb_Triggers_PHI_Q0",  "RPCb_Triggers_PHI_Q0",  n_phi, n_phi_1, n_phi_2);
  RPCb_Triggers_PHI_Q1  = new TH1F("RPCb_Triggers_PHI_Q1",  "RPCb_Triggers_PHI_Q1",  n_phi, n_phi_1, n_phi_2);
  RPCb_Triggers_PHI_Q2  = new TH1F("RPCb_Triggers_PHI_Q2",  "RPCb_Triggers_PHI_Q2",  n_phi, n_phi_1, n_phi_2);
  RPCb_Triggers_PHI_Q3  = new TH1F("RPCb_Triggers_PHI_Q3",  "RPCb_Triggers_PHI_Q3",  n_phi, n_phi_1, n_phi_2);
  RPCb_Triggers_ETA_All = new TH1F("RPCb_Triggers_ETA_All", "RPCb_Triggers_ETA_All", n_eta_exact, n_eta_vec);
  RPCb_Triggers_ETA_Q0  = new TH1F("RPCb_Triggers_ETA_Q0",  "RPCb_Triggers_ETA_Q0",  n_eta_exact, n_eta_vec);
  RPCb_Triggers_ETA_Q1  = new TH1F("RPCb_Triggers_ETA_Q1",  "RPCb_Triggers_ETA_Q1",  n_eta_exact, n_eta_vec);
  RPCb_Triggers_ETA_Q2  = new TH1F("RPCb_Triggers_ETA_Q2",  "RPCb_Triggers_ETA_Q2",  n_eta_exact, n_eta_vec);
  RPCb_Triggers_ETA_Q3  = new TH1F("RPCb_Triggers_ETA_Q3",  "RPCb_Triggers_ETA_Q3",  n_eta_exact, n_eta_vec);
  RPCb_Triggers_ETA_PHI_All = new TH2F("RPCb_Triggers_ETA_PHI_All,", "RPCb_Triggers_ETA_PHI_All,", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);

  DTTF_Triggers_Quality = new TH1F("DTTF_Triggers_Quality", "DTTF_Triggers_Quality", n_qua, n_qua_1, n_qua_2);
  DTTF_Triggers_PHI_All = new TH1F("DTTF_Triggers_PHI_All", "DTTF_Triggers_PHI_All", n_phi, n_phi_1, n_phi_2);
  DTTF_Triggers_PHI_Q0  = new TH1F("DTTF_Triggers_PHI_Q0",  "DTTF_Triggers_PHI_Q0",  n_phi, n_phi_1, n_phi_2);
  DTTF_Triggers_PHI_Q1  = new TH1F("DTTF_Triggers_PHI_Q1",  "DTTF_Triggers_PHI_Q1",  n_phi, n_phi_1, n_phi_2);
  DTTF_Triggers_PHI_Q2  = new TH1F("DTTF_Triggers_PHI_Q2",  "DTTF_Triggers_PHI_Q2",  n_phi, n_phi_1, n_phi_2);
  DTTF_Triggers_PHI_Q3  = new TH1F("DTTF_Triggers_PHI_Q3",  "DTTF_Triggers_PHI_Q3",  n_phi, n_phi_1, n_phi_2);
  DTTF_Triggers_ETA_All = new TH1F("DTTF_Triggers_ETA_All", "DTTF_Triggers_ETA_All", n_eta_exact, n_eta_vec);
  DTTF_Triggers_ETA_Q0  = new TH1F("DTTF_Triggers_ETA_Q0",  "DTTF_Triggers_ETA_Q0",  n_eta_exact, n_eta_vec);
  DTTF_Triggers_ETA_Q1  = new TH1F("DTTF_Triggers_ETA_Q1",  "DTTF_Triggers_ETA_Q1",  n_eta_exact, n_eta_vec);
  DTTF_Triggers_ETA_Q2  = new TH1F("DTTF_Triggers_ETA_Q2",  "DTTF_Triggers_ETA_Q2",  n_eta_exact, n_eta_vec);
  DTTF_Triggers_ETA_Q3  = new TH1F("DTTF_Triggers_ETA_Q3",  "DTTF_Triggers_ETA_Q3",  n_eta_exact, n_eta_vec);
  DTTF_Triggers_ETA_PHI_All = new TH2F("DTTF_Triggers_ETA_PHI_All,", "DTTF_Triggers_ETA_PHI_All,", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);

  ETA_Triggers_RPCb_DTTF_All = new TH2F("ETA_Triggers_RPCb_DTTF_All,", "ETA_Triggers_RPCb_DTTF_All,", n_eta_exact, n_eta_vec, n_eta_exact, n_eta_vec);
  PHI_Triggers_RPCb_DTTF_All = new TH2F("PHI_Triggers_RPCb_DTTF_All,", "PHI_Triggers_RPCb_DTTF_All,", n_phi, n_phi_1, n_phi_2, n_phi, n_phi_1, n_phi_2);
  BX_Triggers_RPCb_DTTF_All  = new TH2F("BX_Triggers_RPCb_DTTF_All,",  "BX_Triggers_RPCb_DTTF_All,",  n_bx, n_bx_1, n_bx_2, n_bx, n_bx_1, n_bx_2);

  TrackCollections = new TH1F("TrackCollections", "TrackCollections", 21, -0.5, 20.5);
  D0_TrackExp = new TH1F("D0_TrackExp", "D0_TrackExp",  60, 0, 120);
  DZ_TrackExp = new TH1F("DZ_TrackExp", "DZ_TrackExp", 150, 0, 300);
  D0_TrackFnd = new TH1F("D0_TrackFnd", "D0_TrackFnd",  60, 0, 120);
  DZ_TrackFnd = new TH1F("DZ_TrackFnd", "DZ_TrackFnd", 150, 0, 300);
  D0_TrackEff = new TH1F("D0_TrackEff", "D0_TrackEff",  60, 0, 120);
  DZ_TrackEff = new TH1F("DZ_TrackEff", "DZ_TrackEff", 150, 0, 300);

  D0DZ_TrackExp = new TH2F("D0DZ_TrackExp", "D0DZ_TrackExp", 15, 0, 300, 12, 0, 120);
  D0DZ_TrackFnd = new TH2F("D0DZ_TrackFnd", "D0DZ_TrackFnd", 15, 0, 300, 12, 0, 120);
  D0DZ_TrackEff = new TH2F("D0DZ_TrackEff", "D0DZ_TrackEff", 15, 0, 300, 12, 0, 120);
  D0_TrackTrk = new TH1F("D0_TrackTrk", "D0_TrackTrk", 60, 0,  120);
  DZ_TrackTrk = new TH1F("DZ_TrackTrk", "DZ_TrackTrk", 150, 0, 300);
  D0DZ_TrackTrk = new TH2F("D0DZ_TrackTrk", "D0DZ_TrackTrk", 15, 0, 300, 12, 0, 120);

  eventsProcessed = 0;
  eventsFiltered = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyOutOfTimeRPCTriggerFilter::endJob() {

  outputfile->cd();
  RPCb_Triggers_Quality->Write();
  RPCb_Triggers_ETA_All->Write(); RPCb_Triggers_ETA_Q0->Write(); RPCb_Triggers_ETA_Q1->Write(); RPCb_Triggers_ETA_Q2->Write(); RPCb_Triggers_ETA_Q3->Write();
  RPCb_Triggers_PHI_All->Write(); RPCb_Triggers_PHI_Q0->Write(); RPCb_Triggers_PHI_Q1->Write(); RPCb_Triggers_PHI_Q2->Write(); RPCb_Triggers_PHI_Q3->Write();
  RPCb_Triggers_ETA_PHI_All->Write();
  DTTF_Triggers_Quality->Write();
  DTTF_Triggers_ETA_All->Write(); DTTF_Triggers_ETA_Q0->Write(); DTTF_Triggers_ETA_Q1->Write(); DTTF_Triggers_ETA_Q2->Write(); DTTF_Triggers_ETA_Q3->Write();
  DTTF_Triggers_PHI_All->Write(); DTTF_Triggers_PHI_Q0->Write(); DTTF_Triggers_PHI_Q1->Write(); DTTF_Triggers_PHI_Q2->Write(); DTTF_Triggers_PHI_Q3->Write();
  DTTF_Triggers_ETA_PHI_All->Write();
  ETA_Triggers_RPCb_DTTF_All->Write();
  PHI_Triggers_RPCb_DTTF_All->Write();
  BX_Triggers_RPCb_DTTF_All->Write();


  TrackCollections->GetXaxis()->SetBinLabel(1, "recoMuon");
  TrackCollections->GetXaxis()->SetBinLabel(2, "standAloneMuon");
  TrackCollections->GetXaxis()->SetBinLabel(3, "Inner Tracks");
  TrackCollections->GetXaxis()->SetBinLabel(4, "Track Expected");
  TrackCollections->GetXaxis()->SetBinLabel(5, "Track Found");
  TrackCollections->Write();

  D0_TrackExp->Write();
  DZ_TrackExp->Write();
  D0_TrackFnd->Write();
  DZ_TrackFnd->Write();
  D0_TrackEff->Divide(D0_TrackFnd,D0_TrackExp,1,1,"");
  DZ_TrackEff->Divide(DZ_TrackFnd,DZ_TrackExp,1,1,"");
  D0_TrackEff->Write();
  DZ_TrackEff->Write();

  D0DZ_TrackExp->Write();
  D0DZ_TrackFnd->Write();
  D0DZ_TrackEff->Divide(D0DZ_TrackFnd,D0DZ_TrackExp,1,1,"");
  D0DZ_TrackEff->Write();
  D0_TrackTrk->Write();
  DZ_TrackTrk->Write();
  D0DZ_TrackTrk->Write();


  outputfile->Close();

  /*
  TrackCollections-GetXaxis()->

    vector<reco::Track>                   "beamhaloTracks"            ""                "RECO"    
    vector<reco::Track>                   "ckfInOutTracksFromConversions"   ""                "RECO"    
    vector<reco::Track>                   "ckfOutInTracksFromConversions"   ""                "RECO"    
    vector<reco::Track>                   "cosmicMuons"               ""                "RECO"    
    vector<reco::Track>                   "cosmicMuons1Leg"           ""                "RECO"    
    vector<reco::Track>                   "cosmicMuonsEndCapsOnly"    ""                "RECO"    
    vector<reco::Track>                   "cosmicMuonsNoRPC"          ""                "RECO"    
    vector<reco::Track>                   "cosmicMuonsWitht0Correction"   ""                "RECO"    
    vector<reco::Track>                   "cosmictrackfinderP5"       ""                "RECO"    
    vector<reco::Track>                   "ctfWithMaterialTracksP5"   ""                "RECO"    
    vector<reco::Track>                   "ctfWithMaterialTracksP5LHCNavigation"   ""                "RECO"    
    vector<reco::Track>                   "globalBeamHaloMuonEndCapslOnly"   ""                "RECO"    
    vector<reco::Track>                   "globalCosmicMuons"         ""                "RECO"    
    vector<reco::Track>                   "globalCosmicMuons1Leg"     ""                "RECO"    
    vector<reco::Track>                   "globalCosmicMuonsNoRPC"    ""                "RECO"    
    vector<reco::Track>                   "globalCosmicMuonsWitht0Correction"   ""                "RECO"    
    vector<reco::Track>                   "globalCosmicSplitMuons"    ""                "RECO"    
    vector<reco::Track>                   "regionalCosmicTracks"      ""                "RECO"    
    vector<reco::Track>                   "splittedTracksP5"          ""                "RECO"    
    vector<reco::Track>                   "standAloneMuons"           ""                "RECO"    
    vector<reco::Track>                   "standAloneMuons"           "UpdatedAtVtx"    "RECO"    
    vector<reco::Track>                   "tevMuons"                  "default"         "RECO"    
    vector<reco::Track>                   "tevMuons"                  "dyt"             "RECO"    
    vector<reco::Track>                   "tevMuons"                  "firstHit"        "RECO"    
    vector<reco::Track>                   "tevMuons"                  "picky"           "RECO"    
  */

  std::cout<<"======================================="<<std::endl;
  std::cout<<"=====  Summary   ======================"<<std::endl;
  std::cout<<"======================================="<<std::endl;
  std::cout<<"=== Events Processed :: "<<std::setw(12)<<eventsProcessed<<" ="<<std::endl;
  std::cout<<"=== Events Filtered  :: "<<std::setw(12)<<eventsFiltered<<" ="<<std::endl;
  std::cout<<"======================================="<<std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void
MyOutOfTimeRPCTriggerFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MyOutOfTimeRPCTriggerFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MyOutOfTimeRPCTriggerFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MyOutOfTimeRPCTriggerFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyOutOfTimeRPCTriggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MyOutOfTimeRPCTriggerFilter);
