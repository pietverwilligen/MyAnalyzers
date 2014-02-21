// -*- C++ -*-
//
// Package:    MyRecHitAnalyzer
// Class:      MyRecHitAnalyzer
// 
/**\class MyRecHitAnalyzer MyRecHitAnalyzer.cc MyAnalyzers/MyRecHitAnalyzer/src/MyRecHitAnalyzer.cc

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

class MyRecHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRecHitAnalyzer(const edm::ParameterSet&);
      ~MyRecHitAnalyzer();

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
  std::string rootFileName;
  TFile * outputfile;

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
MyRecHitAnalyzer::MyRecHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

}


MyRecHitAnalyzer::~MyRecHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyRecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // ================
  // RPC recHits
  // ================
  edm::Handle<RPCRecHitCollection> rpcRecHits; 
  iEvent.getByLabel("rpcRecHits","",rpcRecHits);

  // count the number of RPC rechits
  int nRPC = 0;
  RPCRecHitCollection::const_iterator recHit;
  for (recHit = rpcRecHits->begin(); recHit != rpcRecHits->end(); recHit++) {
    nRPC++;
  }
 
  std::cout<<"The Number of RecHits is "<<nRPC<<std::endl;       
  for (recHit = rpcRecHits->begin(); recHit != rpcRecHits->end(); recHit++) {
    RPCDetId rollId = (RPCDetId)(*recHit).rpcId();
    RPCGeomServ rpcsrv(rollId);
    LocalPoint recHitPos=recHit->localPosition();
    // LocalError recHitErr=recHit->localPositionError();
    const RPCRoll* rollasociated = rpcGeom->roll(rollId);
    const BoundPlane & RPCSurface = rollasociated->surface(); 
    GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal(recHitPos);
    std::cout<<"RPC Rec Hit in "<<rpcsrv.name();
    // std::cout<<" Local Position = "<<recHitPos<<" Local Uncrt = "<<
    std::cout<<" Global Position = "<<RPCGlobalPoint; // <<" Global Uncrt = "<<;
    std::cout<<" Clustersize = "<<recHit->clusterSize()<<" First strip of cluster = "<<recHit->firstClusterStrip()<<" BX = "<<recHit->BunchX()<<std::endl;

  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRecHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRecHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyRecHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
}

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyRecHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyRecHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyRecHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRecHitAnalyzer);
