// -*- C++ -*-
//
// Package:    RPCTiming
// Class:      RPCTiming
// 
/**\class RPCTiming RPCTiming.cc MyAnalyzers/RPCTiming/src/RPCTiming.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piet Verwilligen
//         Created:  Thu Jan  8 10:09:12 CET 2009
// $Id: RPCTiming.cc,v 1.1 2009/06/12 16:08:28 piet Exp $
//
//


// system include files
#include <memory>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <set>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>

//
// class decleration
//

class RPCTiming : public edm::EDAnalyzer {
   public:
      explicit RPCTiming(const edm::ParameterSet&);
      ~RPCTiming();


   private:
  // virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
  // virtual void endJob() ;

      // ----------member data ---------------------------
  std::ofstream fout;

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
RPCTiming::RPCTiming(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  fout.open("RPCTiming_upscope.dat");
  std::cout <<"======================== Opening output file"<< std::endl;

}


RPCTiming::~RPCTiming()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  fout.close();
  std::cout <<"======================== Closing output file"<< std::endl;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
RPCTiming::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // using namespace edm;
  /*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
  */

  edm::ESHandle<RPCGeometry> pDD;
  iSetup.get<MuonGeometryRecord>().get( pDD );

  for(std::vector<RPCRoll*>::const_iterator it = pDD->rolls().begin(); it != pDD->rolls().end(); it++){
    if( dynamic_cast<RPCRoll*>( *it ) != 0 ){ // check if dynamic cast is ok: cast ok => 1
      RPCRoll* roll = dynamic_cast<RPCRoll*>( *it );
      RPCDetId detId = roll->id();
      int rawDetId = detId.rawId();
      RPCGeomServ rpcsrv(detId);
      std::string name = rpcsrv.name();


      // my way
      int strips = roll->nstrips();
      LocalPoint lp = roll->centreOfStrip(strips/2);
      // GlobalPoint gp = roll->toGlobal(lp);
      // float x = gp.x();
      // float y = gp.y();
      // float z = gp.z();
      // GlobalVector* gv = new GlobalVector(x,y,z);
      // all distances are in cm and all time is in ns, therefore c is in cm/ns                                                                                                                                        
      // float c = 29.9792458;
      // float timing = gv->mag()/c;


      // raffas way
      const BoundPlane & RPCSurface = roll->surface();
      GlobalPoint gp = RPCSurface.toGlobal(LocalPoint(0,0,0));
      // all distances are in cm and all time is in ns, therefore c is in cm/ns
      float c = 29.9792458;
      float timing = gp.mag()/c;

      float r = sqrt(pow(gp.x(),2) + pow(gp.y(),2));
      float striplength, stripwidth; 

      if (detId.region()==0){
	const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(roll->topology()));
	striplength = top_->stripLength();
	stripwidth  = top_->pitch();
      }
      else {
	const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&(roll->topology()));
	striplength = top_->stripLength();
	stripwidth  = top_->pitch();
      }

      // RPCTiming
      if(1) std::cout<<"RPCDetId: "<<rawDetId<<" "<<name<<" R: "<<r<<" cm | striplength: "<<striplength<<" cm | stripwidth: "<<stripwidth<<" cm"<<std::endl; // Debug
      if(0) std::cout<<"RPCDetId: "<<rawDetId<<" "<<name<<" LocalPoint: "<<lp<<" GlobalPoint: "<<gp<<" Time: "<<timing<<std::endl; // Debug
      fout<<rawDetId<<"  "<<timing<<std::endl; // RPCTiming.dat

      // RPCDetId_Noise.dat
      // fout<<name<<"  "<<rawDetId<<"  "<<"0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05"<<std::endl; //RPCDetId_Noise.dat

      // RPCDetId_Eff.dat
      // fout<<name<<"  "<<rawDetId<<"  "<<"0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95  0.95"<<std::endl; // RPCDetId_Eff.dat

      // RPC Names
      // fout<<name<<"\t\t"<<rawDetId<<std::endl; // RPCDetId_Name.dat


    }

  }


}


// ------------ method called once each job just before starting event loop  ------------
/*
void 
RPCTiming::beginJob(const edm::EventSetup&)
{
}
*/

// ------------ method called once each job just after ending the event loop  ------------
/*
void 
RPCTiming::endJob() {
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(RPCTiming);
