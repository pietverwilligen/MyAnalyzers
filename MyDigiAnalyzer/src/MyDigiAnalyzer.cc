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
#include "FWCore/Framework/interface/Run.h"
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
  virtual void beginJob();
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

  // 5XY
  edm::InputTag RPCSimHit_Tag, RPCDigi_Tag;

  // 6XY and beyond
  // edm::EDGetTokenT<edm::PSimHitContainer> RPCSimHit_Token;
  // edm::EDGetTokenT<RPCDigiCollection>     RPCDigi_Token;

  TH1F * BXGap_CountPerDigi,   * BXGap_CountPerEvent;
  TH1F * Results_CountPerDigi, * Results_CountPerEvent;
  bool debug;
  int runnumber;
  std::vector< std::pair<uint32_t,int> > DetIdBookKeeping;

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
   debug         = iConfig.getUntrackedParameter<bool>("Debug");
   runnumber     = iConfig.getUntrackedParameter<int>("RunNumber");
   // 5XY  
   digiLabel     = iConfig.getUntrackedParameter<std::string>("DigiLabel","muonRPCDigis");
   // 6XY and beyond
   // RPCDigi_Token   = consumes<RPCDigiCollection>(edm::InputTag(digiLabel));
   rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

  std::stringstream h1ss; h1ss<<"BXGap_CountPerDigi_Run_"<<runnumber;    std::string h1s = h1ss.str();
  std::stringstream h2ss; h2ss<<"BXGap_CountPerEvent_Run_"<<runnumber;   std::string h2s = h2ss.str();
  std::stringstream h3ss; h3ss<<"Results_CountPerDigi_Run_"<<runnumber;  std::string h3s = h3ss.str();
  std::stringstream h4ss; h4ss<<"Results_CountPerEvent_Run_"<<runnumber; std::string h4s = h4ss.str();

  // Book some Histograms
  BXGap_CountPerDigi    = new TH1F(h1s.c_str(), h1s.c_str(),  8, -2.5, 5.5);
  BXGap_CountPerEvent   = new TH1F(h2s.c_str(), h2s.c_str(),  8, -2.5, 5.5);
  Results_CountPerDigi  = new TH1F(h3s.c_str(), h3s.c_str(),  7, -1.5, 5.5);
  Results_CountPerEvent = new TH1F(h4s.c_str(), h4s.c_str(),  7, -1.5, 5.5);

  /*
  if(debug) {
    std::cout<<"[beginJob] Booked: "<<h1s<<" = "<<BXGap_CountPerDigi<<std::endl;
    std::cout<<"[beginJob] Booked: "<<h2s<<" = "<<BXGap_CountPerEvent<<std::endl;
    std::cout<<"[beginJob] Booked: "<<h3s<<" = "<<Results_CountPerDigi<<std::endl;
    std::cout<<"[beginJob] Booked: "<<h4s<<" = "<<Results_CountPerEvent<<std::endl;
  }         
  */
}


MyDigiAnalyzer::~MyDigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  std::cout<<"======================================="<<std::endl;
  std::cout<<"=== Events ::         "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(1)<<" ==="<<std::endl; // Bin 1 = -2 is for all Digis 
  std::cout<<"=== |--> bxgap = 0 :: "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(3)<<" ==="<<std::endl; // Bin 2 = -1 is for the Unique Digis, Bin 3 = 0 is for Copied Digis (BXGap 0)
  std::cout<<"=== |--> bxgap = 1 :: "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(4)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 2 :: "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(5)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 3 :: "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(6)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 4 :: "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(7)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 5 :: "<<std::setw(12)<<BXGap_CountPerEvent->GetBinContent(8)<<" ==="<<std::endl;
  std::cout<<"======================================="<<std::endl;
  std::cout<<"=== Digis ::          "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(1)<<" ==="<<std::endl;  // Bin 1 = -2 is for all Digis 
  std::cout<<"=== |--> bxgap = 0 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(3)<<" ==="<<std::endl;  // Bin 2 = -1 is for the Unique Digis, Bin 3 = 0 is for Copied Digis (BXGap 0)
  std::cout<<"=== |--> bxgap = 1 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(4)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 2 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(5)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 3 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(6)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 4 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(7)<<" ==="<<std::endl;
  std::cout<<"=== |--> bxgap = 5 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(8)<<" ==="<<std::endl;
  std::cout<<"======================================="<<std::endl;

  if(debug) std::cout<<"[endRun] dealing with BXGap_CountPerEvent = "<<BXGap_CountPerEvent<<" and BXGap_CountPerDigi = "<<BXGap_CountPerDigi<<std::endl;
  BXGap_CountPerEvent->SetBinContent(2,BXGap_CountPerEvent->GetBinContent(1));
  BXGap_CountPerDigi->SetBinContent(2,BXGap_CountPerDigi->GetBinContent(1)-BXGap_CountPerDigi->GetBinContent(3));

  BXGap_CountPerEvent->GetYaxis()->SetTitle("Events (-)");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(1,"All Events");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(2,"Unique Digis");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(3,"Copied Digis");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(4,"Digi BX+1");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(5,"Digi BX+2");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(6,"Digi BX+3");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(7,"Digi BX+4");
  BXGap_CountPerEvent->GetXaxis()->SetBinLabel(8,"Digi BX+5");
  BXGap_CountPerEvent->LabelsOption("u");

  BXGap_CountPerDigi->GetYaxis()->SetTitle("Digis (-)");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(1,"All Digis");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(2,"Unique Digis");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(3,"Copied Digis");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(4,"Digi BX+1");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(5,"Digi BX+2");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(6,"Digi BX+3");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(7,"Digi BX+4");
  BXGap_CountPerDigi->GetXaxis()->SetBinLabel(8,"Digi BX+5");
  BXGap_CountPerDigi->LabelsOption("u");

  if(debug) std::cout<<"[endRun] dealing with Results_CountPerEvent = "<<Results_CountPerEvent<<" and Results_CountPerDigi = "<<Results_CountPerDigi<<std::endl;
  // Results_CountPerDigi = new TH1F("Results_CountPerDigi",  "Results_CountPerDigi",  8, -1.5, 6.5);
  int uniquedigis = BXGap_CountPerDigi->GetBinContent(2);
  Results_CountPerDigi->GetYaxis()->SetTitle("% of Unique Digis (-)");
  Results_CountPerDigi->SetBinContent(1,BXGap_CountPerDigi->GetBinContent(2)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(1,"Unique Digis");
  Results_CountPerDigi->SetBinContent(2,BXGap_CountPerDigi->GetBinContent(3)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(2,"Copied Digis");
  Results_CountPerDigi->SetBinContent(3,BXGap_CountPerDigi->GetBinContent(4)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(3,"Digis in BX+1");
  Results_CountPerDigi->SetBinContent(4,BXGap_CountPerDigi->GetBinContent(5)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(4,"Digis in BX+2");
  Results_CountPerDigi->SetBinContent(5,BXGap_CountPerDigi->GetBinContent(6)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(5,"Digis in BX+3");
  Results_CountPerDigi->SetBinContent(6,BXGap_CountPerDigi->GetBinContent(7)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(6,"Digis in BX+4");
  Results_CountPerDigi->SetBinContent(7,BXGap_CountPerDigi->GetBinContent(8)*100.0/uniquedigis); Results_CountPerDigi->GetXaxis()->SetBinLabel(7,"Digis in BX+5");
  Results_CountPerDigi->LabelsOption("u");

  std::cout<<"ok"<<std::endl;

  // Results_CountPerEvent = new TH1F("Results_CountPerEvent",  "Results_CountPerEvent",  8, -1.5, 6.5);
  int events = BXGap_CountPerEvent->GetBinContent(1);
  Results_CountPerEvent->GetYaxis()->SetTitle("% of Events (-)");
  Results_CountPerEvent->SetBinContent(1,BXGap_CountPerEvent->GetBinContent(2)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(1,"Unique Digis");
  Results_CountPerEvent->SetBinContent(2,BXGap_CountPerEvent->GetBinContent(3)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(2,"Copied Digis");
  Results_CountPerEvent->SetBinContent(3,BXGap_CountPerEvent->GetBinContent(4)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(3,"Digis in BX+1");
  Results_CountPerEvent->SetBinContent(4,BXGap_CountPerEvent->GetBinContent(5)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(4,"Digis in BX+2");
  Results_CountPerEvent->SetBinContent(5,BXGap_CountPerEvent->GetBinContent(6)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(5,"Digis in BX+3");
  Results_CountPerEvent->SetBinContent(6,BXGap_CountPerEvent->GetBinContent(7)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(6,"Digis in BX+4");
  Results_CountPerEvent->SetBinContent(7,BXGap_CountPerEvent->GetBinContent(8)*100.0/events); Results_CountPerEvent->GetXaxis()->SetBinLabel(7,"Digis in BX+5");
  Results_CountPerEvent->LabelsOption("u");

  std::cout<<"ok"<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"======================================="<<std::endl;
  std::cout<<"=== DetIds with BXGap = 1,2 or 3    ==="<<std::endl;
  std::cout<<"====== rolls affected ======== "<<std::setw(4)<<DetIdBookKeeping.size()<<" ===";
  std::cout<<"====================================================================================================================="<<std::endl;
  // Sort the rolls according to occurence ...
  std::sort(DetIdBookKeeping.begin(), DetIdBookKeeping.end(),
	    [](const std::pair<uint32_t,int>& a, const std::pair<uint32_t,int>& b) {
	      return a.second> b.second;
	    });
  // Print the sorted rolls ...
  std::vector< std::pair<uint32_t, int> >::const_iterator MyDetIdVectorItr;
  int accumulation = 0;
  int total_BX123 = BXGap_CountPerDigi->GetBinContent(4)+BXGap_CountPerDigi->GetBinContent(5)+BXGap_CountPerDigi->GetBinContent(6);
  for(MyDetIdVectorItr=DetIdBookKeeping.begin(); MyDetIdVectorItr!=DetIdBookKeeping.end(); ++MyDetIdVectorItr) {

    RPCDetId detId((*MyDetIdVectorItr).first);
    RPCGeomServ RPCname(detId);
    std::string nameRoll = RPCname.name();
    accumulation += (*MyDetIdVectorItr).second;
    std::cout<<"=== "<<std::setw(12)<<(*MyDetIdVectorItr).first<<" with "<<std::setw(6)<<(*MyDetIdVectorItr).second<<" counts  = ";
    std::cout<<std::setw(9)<<(*MyDetIdVectorItr).second*100.0/total_BX123<<" % [of total] and "<<std::setw(9)<<accumulation*100.0/total_BX123<<" % [accumulated] ===";
    std::cout<<"   --->   "<<nameRoll<<std::endl;
  }
  std::cout<<"============================================================================================================================================================"<<std::endl;

  outputfile->cd();
  BXGap_CountPerDigi->Write();
  BXGap_CountPerEvent->Write();
  Results_CountPerDigi->Write();
  Results_CountPerEvent->Write();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyDigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // if(runnumber == (int)(iEvent.id()).event() || runnumber == -1) {
  // if(debug) std::cout<<"Processing Run "<<runnumber<<std::endl;
  // int evNum = (iEvent.id()).event();
  // int rnNum = (iEvent.id()).run();
  // int lsNum = (iEvent.id()).luminosityBlock();
  
    
  edm::Handle<RPCDigiCollection> rpcdigis;
  // 5XY
  iEvent.getByLabel(digiLabel, rpcdigis);
  // 6XY
  // iEvent.getByToken(RPCDigi_Token, rpcdigis);
  
  // count per event whether there was 
  // - a digi duplicated (bx-gap 0, i.e. same bx)
  // - a digi in the next bx (bx-gap 1)
  // - a digi in the next-to-next bx (bx-gap 2)
  // - ...
  // up to five
  bool bxgap_counted[] = {0,0,0,0,0,0};
  
  //Loop on digi collection
  for( RPCDigiCollection::DigiRangeIterator collectionItr=rpcdigis->begin(); collectionItr!=rpcdigis->end(); ++collectionItr){
    
    RPCDetId detId=(*collectionItr).first; 
    uint32_t id=detId(); 
    
    // const GeomDet* gdet=rpcGeo->idToDet(detId);
    // const BoundPlane & surface = gdet->surface();
    
    // get roll name
    /*
      RPCGeomServ RPCname(detId);
      std::string nameRoll = RPCname.name();
      std::stringstream os;
    */
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
    // std::cout<<" "<<std::endl;
    // std::cout<<" RPC DetId: "<<std::setw(12)<<id<<" a.k.a. "<<std::setw(18)<<nameRoll<<" which is in "<<std::setw(5)<<wheelOrDiskType<<" "<<std::setw(2)<<wheelOrDiskNumber<<" ";
    // std::cout<<std::setw(7)<<stationOrRingType<<" "<<std::setw(2)<<stationOrRingNumber<<" sector "<<std::setw(2)<<sector<<std::endl;
    if(debug) {
      RPCGeomServ RPCname(detId);
      std::string nameRoll = RPCname.name();
      std::stringstream os;
      std::cout<<" RPC DetId: "<<std::setw(12)<<id<<" a.k.a. "<<std::setw(18)<<nameRoll<<" details: "<<detId<<std::endl;
      std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;
    }
    
    
    RPCDigiCollection::const_iterator digiItr;
    // save digis of one roll in a vector (but save each strip only once, if same strip gives a second digi, fill histogram, but don't save it in the vector)
    std::vector<RPCDigi> MyDigiVector;
    std::vector<RPCDigi>::const_iterator MyDigiVectorItr;
    // loop on digis of given roll
    
    for (digiItr =(*collectionItr ).second.first;digiItr != (*collectionItr ).second.second; ++digiItr){
      if(debug) {
	int strip= (*digiItr).strip();
	int bx=(*digiItr).bx();
	if(debug) {std::cout<<"     Digi: strip = "<<std::setw(2)<<strip<<" bx = "<<std::setw(2)<<bx<<std::endl;}
      }
      // Fill here your histograms
      // -------------------------
      // Analysis of bx gap
      // -------------------------
      // Count All Digis ... commented out for now ...
      BXGap_CountPerDigi->Fill(-2); // std::cout<<"Count Digi"<<std::endl; 
      // If first digi of a Roll then save it in the vector and exit the loop
      bool save_this_digi = true; // start with true, as soon as you find a match, you put it to false ...
      if(MyDigiVector.size()==0) {
	MyDigiVector.push_back((*digiItr));
	if(debug) {std::cout<<"     --> First Digi of Roll :: pushed back and going to the next Digi"<<std::endl;}
      }
      else {
	if(debug) {std::cout<<"     --> Not First Digi of Roll :: loop over all saved Digis to compare"<<std::endl;}
	// If second or more digi of a Roll compare it to the previous digis saved in the vector
	// If there is a digi in the same strip, measure the bx gap
	// If there is no digi with the same strip, save it in the vector
	for(MyDigiVectorItr=MyDigiVector.begin(); MyDigiVectorItr!=MyDigiVector.end(); ++MyDigiVectorItr) {
	  if((*digiItr).strip() == (*MyDigiVectorItr).strip()) {
	    save_this_digi = false;
	    int bxgap = fabs((*MyDigiVectorItr).bx()-(*digiItr).bx());
	    BXGap_CountPerDigi->Fill(bxgap); // std::cout<<"filled "<<bxgap<<std::endl;
	    bxgap_counted[bxgap] = 1;
	    if(debug) {
	      std::cout<<"BX Gap found :: RPC DetId = "<<std::noshowpos<<id<<" ";
	      std::cout<<"first digi = ["<<std::noshowpos<<(*MyDigiVectorItr).strip()<<","<<std::showpos<<(*MyDigiVectorItr).bx()<<"] ";
	      std::cout<<"second digi = ["<<std::noshowpos<<(*digiItr).strip()<<","<<std::showpos<<(*digiItr).bx()<<"] ==> BX Gap = "<<bxgap<<std::endl;
	    }
	    // Store DetId Information for BX Gaps = 1,2,3
	    if(bxgap < 4 && bxgap > 0) { // 1,2,3 should not happen ... store their detids
	      std::vector< std::pair<uint32_t, int> >::iterator MyDetIdVectorItr;
	      // check whether roll was already a source for this type of digis
	      bool detidfound = false;
	      for(MyDetIdVectorItr=DetIdBookKeeping.begin(); MyDetIdVectorItr!=DetIdBookKeeping.end(); ++MyDetIdVectorItr) {
		if((*MyDetIdVectorItr).first == id) { ++(*MyDetIdVectorItr).second; detidfound = true;}
	      }
	      // if roll was not found, then save it now
	      if(detidfound==false) {
		std::pair<uint32_t, int> mypair(id,1);
		DetIdBookKeeping.push_back(mypair); 
	      }
	    }
	  }
	}
	if(save_this_digi) { MyDigiVector.push_back((*digiItr)); /*std::cout<<"Digi added to vector"<<std::endl;*/ }
      }
    }
    if(debug) {
      std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;
    }
    /*
      std::cout<<" "<<std::endl;
      std::cout<<"======================================="<<std::endl;
      std::cout<<"=== Digis ::          "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(1)<<" ==="<<std::endl;
      std::cout<<"=== |--> bxgap = 0 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(2)<<" ==="<<std::endl;
      std::cout<<"=== |--> bxgap = 1 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(3)<<" ==="<<std::endl;
      std::cout<<"=== |--> bxgap = 2 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(4)<<" ==="<<std::endl;
      std::cout<<"=== |--> bxgap = 3 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(5)<<" ==="<<std::endl;
      std::cout<<"=== |--> bxgap = 4 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(6)<<" ==="<<std::endl;
      std::cout<<"=== |--> bxgap = 5 :: "<<std::setw(12)<<BXGap_CountPerDigi->GetBinContent(7)<<" ==="<<std::endl;
      std::cout<<"======================================="<<std::endl;
      std::cout<<" "<<std::endl;
    */
  }
  
  // Now count Event-based
  BXGap_CountPerEvent->Fill(-2); // Count Total Amount of Events
  if(bxgap_counted[0]) BXGap_CountPerEvent->Fill(0);
  if(bxgap_counted[1]) BXGap_CountPerEvent->Fill(1);
  if(bxgap_counted[2]) BXGap_CountPerEvent->Fill(2);
  if(bxgap_counted[3]) BXGap_CountPerEvent->Fill(3);
  if(bxgap_counted[4]) BXGap_CountPerEvent->Fill(4);
  if(bxgap_counted[5]) BXGap_CountPerEvent->Fill(5);
  // }
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
MyDigiAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  /*
  int runnr = iRun.run();

  std::stringstream h1ss; h1ss<<"BXGap_CountPerDigi_Run_"<<runnr;    std::string h1s = h1ss.str();
  std::stringstream h2ss; h2ss<<"BXGap_CountPerEvent_Run_"<<runnr;   std::string h2s = h2ss.str();
  std::stringstream h3ss; h3ss<<"Results_CountPerDigi_Run_"<<runnr;  std::string h3s = h3ss.str();
  std::stringstream h4ss; h4ss<<"Results_CountPerEvent_Run_"<<runnr; std::string h4s = h4ss.str();

  // Book some Histograms
  BXGap_CountPerDigi   = new TH1F(h1s.c_str(), h1s.c_str(),  8, -2.5, 5.5);
  BXGap_CountPerEvent  = new TH1F(h2s.c_str(), h2s.c_str(),  8, -2.5, 5.5);
  Results_CountPerDigi = new TH1F(h3s.c_str(), h3s.c_str(),  7, -1.5, 5.5);
  Results_CountPerEvent = new TH1F(h4s.c_str(), h4s.c_str(),  7, -1.5, 5.5);

  if(debug) {
    std::cout<<"[beginRun] Booked: "<<h1s<<" = "<<BXGap_CountPerDigi<<std::endl;
    std::cout<<"[beginRun] Booked: "<<h2s<<" = "<<BXGap_CountPerEvent<<std::endl;
    std::cout<<"[beginRun] Booked: "<<h3s<<" = "<<Results_CountPerDigi<<std::endl;
    std::cout<<"[beginRun] Booked: "<<h4s<<" = "<<Results_CountPerEvent<<std::endl;
  }
  */

  // GEOMETRY
  // iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  /*
  // Loop on geometry for fun
  // ... can later be used to book histograms
  for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
    // Pre 7XY 
    // if(dynamic_cast< RPCChamber* >( *it ) != 0 ){
    //   RPCChamber* ch = dynamic_cast< RPCChamber* >( *it ); 
    // Post &XY
    if(dynamic_cast< const RPCChamber* >( *it ) != 0 ){
      const RPCChamber* ch = dynamic_cast< const RPCChamber* >( *it ); 
      std::vector< const RPCRoll*> rolls = (ch->rolls());
      for(std::vector<const RPCRoll*>::const_iterator r = rolls.begin();r != rolls.end(); ++r){
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
  */

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
