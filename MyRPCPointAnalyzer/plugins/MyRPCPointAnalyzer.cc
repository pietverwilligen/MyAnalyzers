// -*- C++ -*-
//
// Package:    MyRPCPointAnalyzer
// Class:      MyRPCPointAnalyzer
// 
/**\class MyRPCPointAnalyzer MyRPCPointAnalyzer.cc MyAnalyzers/MyRPCPointAnalyzer/plugins/MyRPCPointAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen,161 R-006,+41227676292,
//         Created:  Wed Oct 24 17:28:30 CEST 2012
// $Id$
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
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
 
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "DataFormats/Provenance/interface/Timestamp.h"


#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>


#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>



//
// class declaration
//

class MyRPCPointAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRPCPointAnalyzer(const edm::ParameterSet&);
      ~MyRPCPointAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::ESHandle <RPCGeometry> rpcGeom;
  edm::ESHandle <DTGeometry>  dtGeom;
  edm::ESHandle <CSCGeometry> cscGeom;

  bool debug;
  int nEvent, nCSCext, nDText;
  int nCSCext_ME1, nCSCext_ME2, nCSCext_ME3, nCSCext_ME4;

  std::string rootFileName;
  TFile * outputfile;
  TH1F * LocalFrame_DeltaX,  * LocalFrame_DeltaY,  * LocalFrame_DeltaZ;
  TH1F * GlobalFrame_DeltaX, * GlobalFrame_DeltaY, * GlobalFrame_DeltaZ;
  TH1F * LocalFrame_ME1_DeltaX, * LocalFrame_ME2_DeltaX, * LocalFrame_ME3_DeltaX, * LocalFrame_ME4_DeltaX;

  edm::InputTag cscSegments;
  edm::InputTag dt4DSegments;
  edm::InputTag rpcRecHitsLabel;
  edm::InputTag rpcDTPointsLabel;
  edm::InputTag rpcCSCPointsLabel;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MyRPCPointAnalyzer::MyRPCPointAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  debug             = iConfig.getUntrackedParameter<bool>("debug",false);

  cscSegments       = iConfig.getUntrackedParameter<edm::InputTag>("cscSegments");
  dt4DSegments      = iConfig.getUntrackedParameter<edm::InputTag>("dt4DSegments");
  rpcRecHitsLabel   = iConfig.getUntrackedParameter<edm::InputTag>("rpcRecHits");
  rpcDTPointsLabel  = iConfig.getUntrackedParameter<edm::InputTag>("rpcDTPoints");
  rpcCSCPointsLabel = iConfig.getUntrackedParameter<edm::InputTag>("rpcCSCPoints");

  rootFileName      = iConfig.getUntrackedParameter<std::string>("RootFileName");

  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

  LocalFrame_DeltaX  = new TH1F("LocalFrame_DeltaX",  "LocalFrame_DeltaX",  200, -10.00, 10.00);
  LocalFrame_DeltaY  = new TH1F("LocalFrame_DeltaY",  "LocalFrame_DeltaY",  200, -10.00, 10.00);
  LocalFrame_DeltaZ  = new TH1F("LocalFrame_DeltaZ",  "LocalFrame_DeltaZ",  200, -10.00, 10.00);
  GlobalFrame_DeltaX = new TH1F("GlobalFrame_DeltaX", "GlobalFrame_DeltaX", 200, -10.00, 10.00);
  GlobalFrame_DeltaY = new TH1F("GlobalFrame_DeltaY", "GlobalFrame_DeltaY", 200, -10.00, 10.00);
  GlobalFrame_DeltaZ = new TH1F("GlobalFrame_DeltaZ", "GlobalFrame_DeltaZ", 200, -10.00, 10.00);

  LocalFrame_ME1_DeltaX  = new TH1F("LocalFrame_ME1_DeltaX",  "LocalFrame_ME1_DeltaX",  200, -10.00, 10.00);
  LocalFrame_ME2_DeltaX  = new TH1F("LocalFrame_ME2_DeltaX",  "LocalFrame_ME2_DeltaX",  200, -10.00, 10.00);
  LocalFrame_ME3_DeltaX  = new TH1F("LocalFrame_ME3_DeltaX",  "LocalFrame_ME3_DeltaX",  200, -10.00, 10.00);
  LocalFrame_ME4_DeltaX  = new TH1F("LocalFrame_ME4_DeltaX",  "LocalFrame_ME4_DeltaX",  200, -10.00, 10.00);

  nEvent  = 0;
  nCSCext = 0; 
  nDText  = 0;
  nCSCext_ME1 = 0; nCSCext_ME2 = 0; nCSCext_ME3 = 0; nCSCext_ME4 = 0;
}


MyRPCPointAnalyzer::~MyRPCPointAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();
  LocalFrame_DeltaX->Write();
  LocalFrame_DeltaY->Write();
  LocalFrame_DeltaZ->Write();
  GlobalFrame_DeltaX->Write();
  GlobalFrame_DeltaY->Write();
  GlobalFrame_DeltaZ->Write();

  LocalFrame_ME1_DeltaX->Write();
  LocalFrame_ME2_DeltaX->Write();
  LocalFrame_ME3_DeltaX->Write();
  LocalFrame_ME4_DeltaX->Write();

  std::cout<<"========================================="<<std::endl;
  std::cout<<"= Amount of Events Processed:    "<<std::setw(6)<<nEvent <<" ="<<std::endl;
  std::cout<<"= Successful  DT Extrapolations: "<<std::setw(6)<<nDText <<" ="<<std::endl;
  std::cout<<"= Successful CSC Extrapolations: "<<std::setw(6)<<nCSCext<<" ="<<std::endl;
  std::cout<<"-----------------------------------------"<<std::endl;
  std::cout<<"=                ME1 Extrap:     "<<std::setw(6)<<nCSCext_ME1<<" ="<<std::endl;
  std::cout<<"=                ME2 Extrap:     "<<std::setw(6)<<nCSCext_ME2<<" ="<<std::endl;
  std::cout<<"=                ME3 Extrap:     "<<std::setw(6)<<nCSCext_ME3<<" ="<<std::endl;
  std::cout<<"=                ME4 Extrap:     "<<std::setw(6)<<nCSCext_ME4<<" ="<<std::endl;
  std::cout<<"========================================="<<std::endl;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyRPCPointAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++nEvent;
  int evtNum = (iEvent.id()).event();
  int lumNum = (iEvent.id()).luminosityBlock();
  int runNum = (iEvent.id()).run();

  if(debug) std::cout<<"----- Run "<<std::setw(7)<<runNum<<" Lumi "<<std::setw(3)<<lumNum<<" Event "<<std::setw(12)<<evtNum<<" ---------------------------------------------------------------------------"<<std::endl;

  edm::Handle<DTRecSegment4DCollection> allDTSegments;
  iEvent.getByLabel(dt4DSegments, allDTSegments);

  edm::Handle<CSCSegmentCollection> allCSCSegments;
  iEvent.getByLabel(cscSegments, allCSCSegments);

  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByLabel(rpcRecHitsLabel,rpcRecHits);

  edm::Handle<RPCRecHitCollection> rpcDTPoints;
  iEvent.getByLabel(rpcDTPointsLabel,rpcDTPoints);

  edm::Handle<RPCRecHitCollection> rpcCSCPoints;
  iEvent.getByLabel(rpcCSCPointsLabel,rpcCSCPoints);

  // ================
  // Collection Sizes
  // ================
  if(debug) {
    std::cout<<"Amount of DT | CSC | RPC Segments/Hits and DT | CSC Extrap in this Event: [ ";
    std::cout<<allDTSegments->size()<<" | "<<allCSCSegments->size()<<" | "<<rpcRecHits->size()<<" | "<<rpcDTPoints->size()<<" | "<<rpcCSCPoints->size()<<" ]"<<std::endl;
  }


  // =========================================
  // Quick Overview Segments in ME4
  // =========================================
  if(debug) {
    CSCSegmentCollection::const_iterator segment;
    for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment){
      // CSCDetId CSCId = segment->cscDetId();
      if((segment->cscDetId()).station()==4) {
	std::cout<<"CSC Segment in "<<segment->cscDetId()<<std::endl;
      }
    }
  }


  // =========================================
  // Quick Overview Extrapolated Points in RE4
  // =========================================
  if(debug) {
    RPCRecHitCollection::const_iterator point;
    for(point = rpcDTPoints->begin();  point != rpcDTPoints->end();  point++) { 
      // RPCGeomServ srv(point->rpcId()); std::cout<<"DT Extapolated Point in "  <<(point->rpcId()).rawId()<<" = "<<srv.name()<<std::endl;
    }
    for(point = rpcCSCPoints->begin(); point != rpcCSCPoints->end(); point++) { 
      if((point->rpcId()).station()==4) {
	RPCGeomServ srv(point->rpcId()); std::cout<<"CSC Extapolated Point in " <<(point->rpcId()).rawId()<<" = "<<srv.name()<<std::endl;}
    }
    std::cout<<"RPC Rechits in { ";
    for(point = rpcRecHits->begin();   point != rpcRecHits->end();   point++) { 
      if((point->rpcId()).station()==4 && (point->rpcId()).region()!=0){
	std::cout<<(point->rpcId()).rawId()<<" ";   
      }
    }
    std::cout<<" }"<<std::endl;
    std::cout<<"RPC Rechits in { ";
    for(point = rpcRecHits->begin();   point != rpcRecHits->end();   point++) { 
      if((point->rpcId()).station()==4 && (point->rpcId()).region()!=0){
	RPCGeomServ srv(point->rpcId()); std::cout<<srv.name()<<" "; 
      }
    }
    std::cout<<" }"<<std::endl;
  }
  // ================
  // RPC recHits
  // ================
  // edm::Handle<RPCRecHitCollection> rpcRecHits; 
  // iEvent.getByLabel("rpcRecHits","",rpcRecHits);
  RPCRecHitCollection::const_iterator recHit;
  for (recHit = rpcRecHits->begin(); recHit != rpcRecHits->end(); recHit++) {
    RPCDetId rollId = (RPCDetId)(*recHit).rpcId();
    RPCGeomServ rollsrv(rollId);
    std::string nameRoll = rollsrv.name();

    LocalPoint RPCLocalPoint=recHit->localPosition();
    // LocalError recHitErr=recHit->localPositionError();
    const RPCRoll* Roll = rpcGeom->roll(rollId);
    const BoundPlane & RPCSurface = Roll->surface(); 
    GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal(RPCLocalPoint);

    // ==========================
    // DT extrapolated points
    // ========================== 
    if(rpcDTPoints.isValid() && rpcDTPoints->begin()!=rpcDTPoints->end()) {

      RPCRecHitCollection::const_iterator rpcPoint;
      for(rpcPoint = rpcDTPoints->begin(); rpcPoint != rpcDTPoints->end(); rpcPoint++){

	// Only details if there is an extrapolated point in this DetId
	RPCDetId  rpcId = rpcPoint->rpcId();
	if(rpcId.rawId() != rollId.rawId()) continue;

	// RPC DetId and Extrapolated DetId match
	RPCGeomServ rpcsrv(rpcId);
	std::string nameRPC = rpcsrv.name();

	const RPCRoll* DTExtrap_Roll           = rpcGeom->roll(rpcId);
	LocalPoint DTExtrap_RPCLocalPoint      = rpcPoint->localPosition();
	const BoundPlane & DTExtrap_RPCSurface = DTExtrap_Roll->surface();
	GlobalPoint DTExtrap_RPCGlobalPoint    = DTExtrap_RPCSurface.toGlobal(DTExtrap_RPCLocalPoint);
	// const float stripPredicted =DTExtrap_Roll->strip(LocalPoint(DTExtrap_RPCLocalPoint.x(),DTExtrap_RPCLocalPoint.y(),0.));
	const float stripPredicted =DTExtrap_Roll->strip(DTExtrap_RPCLocalPoint);

	if(debug) {
	  // Information about original RPC Rechit
	  std::cout<<"RPC Rec Hit in "    <<rollId.rawId()<<" = "<<nameRoll;
	  std::cout<<" Local Position = " <<RPCLocalPoint;//<<" Local Uncrt = "<<
	  std::cout<<" Global Position = "<<RPCGlobalPoint; // <<" Global Uncrt = "<<;
	  std::cout<<" Clustersize = "    <<recHit->clusterSize()<<" First strip of cluster = "<<recHit->firstClusterStrip()<<" BX = "<<recHit->BunchX()<<std::endl;
	  // Information about extrapolated RPC Point
	  std::cout<<"RPC Ext Pnt in "    <<rpcId.rawId()<<" = "<<nameRPC;
	  std::cout<<" Local Position  = "<<DTExtrap_RPCLocalPoint;
	  std::cout<<" Global Position = "<<DTExtrap_RPCGlobalPoint;
	  std::cout<<" Strip Predicted = "<<stripPredicted<<std::endl;
	}

	++nDText;

	LocalFrame_DeltaX->Fill(DTExtrap_RPCLocalPoint.x()-RPCLocalPoint.x());
	LocalFrame_DeltaY->Fill(DTExtrap_RPCLocalPoint.y()-RPCLocalPoint.y());
	LocalFrame_DeltaZ->Fill(DTExtrap_RPCLocalPoint.z()-RPCLocalPoint.z());
	GlobalFrame_DeltaX->Fill(DTExtrap_RPCGlobalPoint.x()-RPCGlobalPoint.x());
	GlobalFrame_DeltaY->Fill(DTExtrap_RPCGlobalPoint.y()-RPCGlobalPoint.y());
	GlobalFrame_DeltaZ->Fill(DTExtrap_RPCGlobalPoint.z()-RPCGlobalPoint.z());

      } // End Loop over DT Extrap Points
    } // End Check whether DT Extrap Points are empty

    // ==========================
    // CSC extrapolated points
    // ==========================
    if(rpcCSCPoints.isValid()) if(rpcCSCPoints->begin()!=rpcCSCPoints->end()) {

      RPCRecHitCollection::const_iterator rpcPoint;
      for(rpcPoint = rpcCSCPoints->begin(); rpcPoint != rpcCSCPoints->end(); rpcPoint++){

	// Only details if there is an extrapolated point in this DetId
	RPCDetId  rpcId = rpcPoint->rpcId();
	if(rpcId.rawId() != rollId.rawId()) continue;

	// RPC DetId and Extrapolated DetId match
	RPCGeomServ rpcsrv(rpcId);
	std::string nameRPC = rpcsrv.name();

	const RPCRoll* CSCExtrap_Roll           = rpcGeom->roll(rpcId);
	LocalPoint CSCExtrap_RPCLocalPoint      = rpcPoint->localPosition();
	const BoundPlane & CSCExtrap_RPCSurface = CSCExtrap_Roll->surface();
	GlobalPoint CSCExtrap_RPCGlobalPoint    = CSCExtrap_RPCSurface.toGlobal(CSCExtrap_RPCLocalPoint);
	// const float stripPredicted =CSCExtrap_Roll->strip(LocalPoint(CSCExtrap_RPCLocalPoint.x(),CSCExtrap_RPCLocalPoint.y(),0.));
	const float stripPredicted =CSCExtrap_Roll->strip(CSCExtrap_RPCLocalPoint);

	if(debug) {
	  // Information about original RPC Rechit
	  std::cout<<"RPC Rec Hit in "    <<rollId.rawId()<<" = "<<nameRoll;
	  std::cout<<" Local Position = " <<RPCLocalPoint;//<<" Local Uncrt = "<<
	  std::cout<<" Global Position = "<<RPCGlobalPoint; // <<" Global Uncrt = "<<;
	  std::cout<<" Clustersize = "    <<recHit->clusterSize()<<" First strip of cluster = "<<recHit->firstClusterStrip()<<" BX = "<<recHit->BunchX()<<std::endl;
	  // Information about extrapolated RPC Point
	  std::cout<<"RPC Ext Pnt in "    <<rpcId.rawId()<<" = "<<nameRPC;
	  std::cout<<" Local Position  = "<<CSCExtrap_RPCLocalPoint;
	  std::cout<<" Global Position = "<<CSCExtrap_RPCGlobalPoint;
	  std::cout<<" Strip Predicted = "<<stripPredicted<<std::endl;
	}

	++nCSCext;

	LocalFrame_DeltaX->Fill(CSCExtrap_RPCLocalPoint.x()-RPCLocalPoint.x());
	LocalFrame_DeltaY->Fill(CSCExtrap_RPCLocalPoint.y()-RPCLocalPoint.y());
	LocalFrame_DeltaZ->Fill(CSCExtrap_RPCLocalPoint.z()-RPCLocalPoint.z());
	GlobalFrame_DeltaX->Fill(CSCExtrap_RPCGlobalPoint.x()-RPCGlobalPoint.x());
	GlobalFrame_DeltaY->Fill(CSCExtrap_RPCGlobalPoint.y()-RPCGlobalPoint.y());
	GlobalFrame_DeltaZ->Fill(CSCExtrap_RPCGlobalPoint.z()-RPCGlobalPoint.z());

	if(rpcId.station()==1) { ++nCSCext_ME1; LocalFrame_ME1_DeltaX->Fill(CSCExtrap_RPCLocalPoint.x()-RPCLocalPoint.x());}
	if(rpcId.station()==2) { ++nCSCext_ME2; LocalFrame_ME2_DeltaX->Fill(CSCExtrap_RPCLocalPoint.x()-RPCLocalPoint.x());}
	if(rpcId.station()==3) { ++nCSCext_ME3; LocalFrame_ME3_DeltaX->Fill(CSCExtrap_RPCLocalPoint.x()-RPCLocalPoint.x());}
	if(rpcId.station()==4) { ++nCSCext_ME4; LocalFrame_ME4_DeltaX->Fill(CSCExtrap_RPCLocalPoint.x()-RPCLocalPoint.x());}

      } // End Loop over CSC Extrap Points
    } // End Check whether CSC Extrap Points are empty

  } // End Loop over RPC Rechits

  if(debug) std::cout<<"------------------------------------------------------------------------------------------------------------------------"<<std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRPCPointAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRPCPointAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyRPCPointAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
}

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyRPCPointAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyRPCPointAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyRPCPointAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRPCPointAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRPCPointAnalyzer);
