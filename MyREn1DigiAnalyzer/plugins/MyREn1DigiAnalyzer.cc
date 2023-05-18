// -*- C++ -*-
//
// Package:    MyREn1DigiAnalyzer
// Class:      MyREn1DigiAnalyzer
// 
/**\class MyREn1DigiAnalyzer MyREn1DigiAnalyzer.cc MyAnalyzers/MyREn1DigiAnalyzer/src/MyREn1DigiAnalyzer.cc

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
#include "TAxis.h"
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

class MyREn1DigiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyREn1DigiAnalyzer(const edm::ParameterSet&);
      ~MyREn1DigiAnalyzer();

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
  bool debug;
  TFile * outputfile;

  edm::EDGetTokenT<edm::PSimHitContainer> RPCSimHit_Token;
  edm::EDGetTokenT<RPCDigiCollection>     RPCDigi_Token;

  TH2F * ST_Digis_RE_P, * ST_Digis_RE_N;
  TH2F * BX_Digis_RE_P, * BX_Digis_RE_N;
};

//
// constants, enums and typedefs
//
int n_strips = 256; int n_s1 = 1; int n_s2 = 256;
int n_bx = 5; double n_b1 = -2.5; double n_b2 = +2.5;
int n_rolls = 10; double n_r1 = 1; double n_r2 = 11;

std::string rolls_P[10] = {"RE+3/1A","RE+3/1B","RE+3/1C","RE+3/1D","RE+3/1E","RE+4/1A","RE+4/1B","RE+4/1C","RE+4/1D","RE+4/1E"};
std::string rolls_N[10] = {"RE-3/1A","RE-3/1B","RE-3/1C","RE-3/1D","RE-3/1E","RE-4/1A","RE-4/1B","RE-4/1C","RE-4/1D","RE-4/1E"};

//
// static data member definitions
//

//
// constructors and destructor
//
MyREn1DigiAnalyzer::MyREn1DigiAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   digiLabel     = iConfig.getUntrackedParameter<std::string>("DigiLabel","muonRPCDigis");
   rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
   debug         = iConfig.getUntrackedParameter<bool>("Debug");

   RPCSimHit_Token = consumes<edm::PSimHitContainer>(edm::InputTag(std::string("g4SimHits"), std::string("MuonRPCHits")));
   RPCDigi_Token   = consumes<RPCDigiCollection>(edm::InputTag(digiLabel));

   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

   ST_Digis_RE_P = new TH2F("ST_Digis_RE_P", "Strip Occupancy of RE+N/1", n_strips, n_s1, n_s2, n_rolls, n_r1, n_r2);
   BX_Digis_RE_P = new TH2F("BX_Digis_RE_P", "Bunch Crossings of RE+N/1", n_bx, n_b1, n_b2, n_rolls, n_r1, n_r2);
   ST_Digis_RE_N = new TH2F("ST_Digis_RE_N", "Strip Occupancy of RE-N/1", n_strips, n_s1, n_s2, n_rolls, n_r1, n_r2);
   BX_Digis_RE_N = new TH2F("BX_Digis_RE_N", "Bunch Crossings of RE-N/1", n_bx, n_b1, n_b2, n_rolls, n_r1, n_r2);
}


MyREn1DigiAnalyzer::~MyREn1DigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  for(int i=0; i<10; ++i) {
    ST_Digis_RE_P->GetYaxis()->SetBinLabel(i+1,rolls_P[i].c_str());
    BX_Digis_RE_P->GetYaxis()->SetBinLabel(i+1,rolls_P[i].c_str());
    ST_Digis_RE_N->GetYaxis()->SetBinLabel(i+1,rolls_N[i].c_str());
    BX_Digis_RE_N->GetYaxis()->SetBinLabel(i+1,rolls_N[i].c_str());
  }

  ST_Digis_RE_P->Write();
  BX_Digis_RE_P->Write();
  ST_Digis_RE_N->Write();
  BX_Digis_RE_N->Write();

  // outputfile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyREn1DigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    if(debug) {
      std::cout<<" "<<std::endl;
      std::cout<<" RPC DetId: "<<std::setw(12)<<id<<" a.k.a. "<<std::setw(18)<<nameRoll<<" details: "<<detId<<std::endl;
      std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;
    }
    RPCDigiCollection::const_iterator digiItr; 
    //loop on digis of given roll
    for (digiItr =(*collectionItr ).second.first;digiItr != (*collectionItr ).second.second; ++digiItr){
      int strip= (*digiItr).strip();
      int bx=(*digiItr).bx();
      if(debug) {
	std::cout<<"     Digi: strip = "<<std::setw(2)<<strip<<" bx = "<<std::setw(2)<<bx<<std::endl;
      }
      // Fill here your histograms
      if(detId.region()==1 && detId.ring()==1) {
	// Filling only station 3 and 4
	if(detId.station()==3 || detId.station()==4) {
	  // std::cout<<"Fill Positive :: (strip,bx,row) = ("<<strip<<", "<<bx<<", "<<(detId.station()-3)*5+detId.roll()<<")"<<std::endl;
	  ST_Digis_RE_P->Fill(strip,(detId.station()-3)*5+detId.roll());
	  BX_Digis_RE_P->Fill(bx,(detId.station()-3)*5+detId.roll());
	}
      }
      if(detId.region()==-1 && detId.ring()==1) {
	// Filling only station 3 and 4
	if(detId.station()==3 || detId.station()==4) {
	  // std::cout<<"Fill Negative :: (strip,bx,row) = ("<<strip<<", "<<bx<<", "<<(detId.station()-3)*5+detId.roll()<<")"<<std::endl;
	  ST_Digis_RE_N->Fill(strip,(detId.station()-3)*5+detId.roll());
	  BX_Digis_RE_N->Fill(bx,(detId.station()-3)*5+detId.roll());
	}
      }

    }
    if(debug) std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;

  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyREn1DigiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyREn1DigiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyREn1DigiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
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
MyREn1DigiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyREn1DigiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyREn1DigiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyREn1DigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyREn1DigiAnalyzer);
