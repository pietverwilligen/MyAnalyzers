// -*- C++ -*-
//
// Package:    MyDigiAnalyzer
// Class:      MyDigiAnalyzer
// 
/**\class MyDigiAnalyzer MyDigiAnalyzer.cc MyAnalyzers/MyDigiAnalyzer/src/MyDigiAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen,161 R-006,+41227676292,
//         Created:  Sat Oct  6 17:08:32 CEST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


// root include files
#include <TRandom.h>
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
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
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///Data Format
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
///Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
///Log messages
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

//
// class declaration
//

class MyDigiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyDigiAnalyzer(const edm::ParameterSet&);
      ~MyDigiAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::ESHandle<RPCGeometry> rpcGeo;
  std::string digiLabel, rootFileName;
  TFile * outputfile;

  edm::EDGetTokenT<edm::PSimHitContainer> RPCSimHit_Token;
  edm::EDGetTokenT<RPCDigiCollection>     RPCDigi_Token;


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
MyDigiAnalyzer::MyDigiAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   digiLabel     = iConfig.getUntrackedParameter<std::string>("DigiLabel","muonRPCDigis");
   rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");

   RPCSimHit_Token = consumes<edm::PSimHitContainer>(edm::InputTag(std::string("g4SimHits"), std::string("MuonRPCHits")));
   RPCDigi_Token   = consumes<RPCDigiCollection>(edm::InputTag(digiLabel));

   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

}


MyDigiAnalyzer::~MyDigiAnalyzer()
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
MyDigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<RPCDigiCollection> rpcdigis;
  iEvent.getByToken(RPCDigi_Token, rpcdigis);

  //Loop on digi collection
  for( RPCDigiCollection::DigiRangeIterator collectionItr=rpcdigis->begin(); collectionItr!=rpcdigis->end(); ++collectionItr){

    RPCDetId detId=(*collectionItr).first; 
    uint32_t id=detId(); 

    // const GeomDet* gdet=rpcGeo->idToDet(detId);
    // const BoundPlane & surface = gdet->surface();

    // get roll name
    RPCGeomServ RPCname(detId);
    std::string nameRoll = RPCname.name();
    std::stringstream os;
    // get info
    /*
    int region  = detId.region();
    int station = detId.station();
    int ring    = detId.ring();
    std::string wheelOrDiskType;
    int wheelOrDiskNumber;
    std::string stationOrRingType;
    int stationOrRingNumber;
    if(region == 0) {
      wheelOrDiskType     = "Wheel";  
      wheelOrDiskNumber   = (int)detId.ring();
      stationOrRingType   = "Station";
      stationOrRingNumber = (int)detId.station();
    }
    else{
      wheelOrDiskType     =  "Disk";
      wheelOrDiskNumber   = region*(int)detId.station();
      stationOrRingType   = "Ring";
      stationOrRingNumber = (int)detId.ring();
    }
    int sector = detId.sector();
    */
    std::cout<<" "<<std::endl;
    // std::cout<<" RPC DetId: "<<std::setw(12)<<id<<" a.k.a. "<<std::setw(18)<<nameRoll<<" which is in "<<std::setw(5)<<wheelOrDiskType<<" "<<std::setw(2)<<wheelOrDiskNumber<<" ";
    // std::cout<<std::setw(7)<<stationOrRingType<<" "<<std::setw(2)<<stationOrRingNumber<<" sector "<<std::setw(2)<<sector<<std::endl;
    std::cout<<" RPC DetId: "<<std::setw(12)<<id<<" a.k.a. "<<std::setw(18)<<nameRoll<<" details: "<<detId<<std::endl;
    std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;
    RPCDigiCollection::const_iterator digiItr; 
    //loop on digis of given roll
    for (digiItr =(*collectionItr ).second.first;digiItr != (*collectionItr ).second.second; ++digiItr){
      int strip= (*digiItr).strip();
      int bx=(*digiItr).bx();
      std::cout<<"     Digi: strip = "<<std::setw(2)<<strip<<" bx = "<<std::setw(2)<<bx<<std::endl;
      // Fill here your histograms
    }
    std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;

  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyDigiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyDigiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyDigiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{

  iSetup.get<MuonGeometryRecord>().get(rpcGeo);

  // Loop on geometry for fun
  // ... can later be used to book histograms
  for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
    if(dynamic_cast< RPCChamber* >( *it ) != 0 ){
      RPCChamber* ch = dynamic_cast< RPCChamber* >( *it ); 
      std::vector< const RPCRoll*> roles = (ch->rolls());
      for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
	RPCDetId rpcId = (*r)->id();
	int region=rpcId.region();
	RPCGeomServ rpcsrv(rpcId);
	std::string nameRoll = rpcsrv.name();
	int ring;
	if(rpcId.region() == 0) 
	  ring = rpcId.ring();
	else 
	  ring = rpcId.region()*rpcId.station();
	std::pair<int,int> regionRing(region,ring);
      }
    }
  }// end loop on geometry
}

// ------------ method called when ending the processing of a run  ------------
void 
MyDigiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyDigiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyDigiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyDigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyDigiAnalyzer);
