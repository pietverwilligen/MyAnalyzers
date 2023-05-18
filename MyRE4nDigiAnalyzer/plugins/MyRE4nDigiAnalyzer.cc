// -*- C++ -*-
//
// Package:    MyRE4nDigiAnalyzer
// Class:      MyRE4nDigiAnalyzer
// 
/**\class MyRE4nDigiAnalyzer MyRE4nDigiAnalyzer.cc MyAnalyzers/MyRE4nDigiAnalyzer/src/MyRE4nDigiAnalyzer.cc

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

//
// class declaration
//

class MyRE4nDigiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRE4nDigiAnalyzer(const edm::ParameterSet&);
      ~MyRE4nDigiAnalyzer();

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
  TH1F * BX_Digis_RE4_Plus,  * BX_Digis_RE4_Ring1_Plus,  * BX_Digis_RE4_Ring2_Plus,  * BX_Digis_RE4_Ring3_Plus,  * ST_Digis_RE4_Plus,  * ST_Digis_RE4_Ring1_Plus,  * ST_Digis_RE4_Ring2_Plus,  * ST_Digis_RE4_Ring3_Plus;
  TH1F * BX_Digis_RE4_Minus, * BX_Digis_RE4_Ring1_Minus, * BX_Digis_RE4_Ring2_Minus, * BX_Digis_RE4_Ring3_Minus, * ST_Digis_RE4_Minus, * ST_Digis_RE4_Ring1_Minus, * ST_Digis_RE4_Ring2_Minus, * ST_Digis_RE4_Ring3_Minus; 

  TH1F * BX_Digis_RE3_Plus,  * BX_Digis_RE3_Ring1_Plus,  * BX_Digis_RE3_Ring2_Plus,  * BX_Digis_RE3_Ring3_Plus,  * ST_Digis_RE3_Plus,  * ST_Digis_RE3_Ring1_Plus,  * ST_Digis_RE3_Ring2_Plus,  * ST_Digis_RE3_Ring3_Plus;
  TH1F * BX_Digis_RE3_Minus, * BX_Digis_RE3_Ring1_Minus, * BX_Digis_RE3_Ring2_Minus, * BX_Digis_RE3_Ring3_Minus, * ST_Digis_RE3_Minus, * ST_Digis_RE3_Ring1_Minus, * ST_Digis_RE3_Ring2_Minus, * ST_Digis_RE3_Ring3_Minus; 

  TCanvas * BX_Digis_RE4, * BX_Digis_RE3, * ST_Digis_RE4, * ST_Digis_RE3;

};

//
// constants, enums and typedefs
//
int n_bx  = 11;  double n1_bx  = -5.5,  n2_bx  = 5.5;
int n_st  = 34;  double n1_st  = 0,     n2_st  = 33;
//
// static data member definitions
//

//
// constructors and destructor
//
MyRE4nDigiAnalyzer::MyRE4nDigiAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   digiLabel     = iConfig.getUntrackedParameter<std::string>("DigiLabel","muonRPCDigis");
   rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");

   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

   BX_Digis_RE4_Plus       = new TH1F("BX_Digis_RE4_Plus", "BX_Digis_RE4_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Ring1_Plus = new TH1F("BX_Digis_RE4_Ring1_Plus", "BX_Digis_RE4_Ring1_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Ring2_Plus = new TH1F("BX_Digis_RE4_Ring2_Plus", "BX_Digis_RE4_Ring2_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Ring3_Plus = new TH1F("BX_Digis_RE4_Ring3_Plus", "BX_Digis_RE4_Ring3_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Minus       = new TH1F("BX_Digis_RE4_Minus", "BX_Digis_RE4_Minus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Ring1_Minus = new TH1F("BX_Digis_RE4_Ring1_Minus", "BX_Digis_RE4_Ring1_Minus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Ring2_Minus = new TH1F("BX_Digis_RE4_Ring2_Minus", "BX_Digis_RE4_Ring2_Minus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE4_Ring3_Minus = new TH1F("BX_Digis_RE4_Ring3_Minus", "BX_Digis_RE4_Ring3_Minus", n_bx, n1_bx, n2_bx);

   ST_Digis_RE4_Plus       = new TH1F("ST_Digis_RE4_Plus", "ST_Digis_RE4_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Ring1_Plus = new TH1F("ST_Digis_RE4_Ring1_Plus", "ST_Digis_RE4_Ring1_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Ring2_Plus = new TH1F("ST_Digis_RE4_Ring2_Plus", "ST_Digis_RE4_Ring2_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Ring3_Plus = new TH1F("ST_Digis_RE4_Ring3_Plus", "ST_Digis_RE4_Ring3_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Minus       = new TH1F("ST_Digis_RE4_Minus", "ST_Digis_RE4_Minus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Ring1_Minus = new TH1F("ST_Digis_RE4_Ring1_Minus", "ST_Digis_RE4_Ring1_Minus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Ring2_Minus = new TH1F("ST_Digis_RE4_Ring2_Minus", "ST_Digis_RE4_Ring2_Minus", n_st, n1_st, n2_st);
   ST_Digis_RE4_Ring3_Minus = new TH1F("ST_Digis_RE4_Ring3_Minus", "ST_Digis_RE4_Ring3_Minus", n_st, n1_st, n2_st);

   BX_Digis_RE3_Plus       = new TH1F("BX_Digis_RE3_Plus", "BX_Digis_RE3_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Ring1_Plus = new TH1F("BX_Digis_RE3_Ring1_Plus", "BX_Digis_RE3_Ring1_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Ring2_Plus = new TH1F("BX_Digis_RE3_Ring2_Plus", "BX_Digis_RE3_Ring2_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Ring3_Plus = new TH1F("BX_Digis_RE3_Ring3_Plus", "BX_Digis_RE3_Ring3_Plus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Minus       = new TH1F("BX_Digis_RE3_Minus", "BX_Digis_RE3_Minus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Ring1_Minus = new TH1F("BX_Digis_RE3_Ring1_Minus", "BX_Digis_RE3_Ring1_Minus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Ring2_Minus = new TH1F("BX_Digis_RE3_Ring2_Minus", "BX_Digis_RE3_Ring2_Minus", n_bx, n1_bx, n2_bx);
   BX_Digis_RE3_Ring3_Minus = new TH1F("BX_Digis_RE3_Ring3_Minus", "BX_Digis_RE3_Ring3_Minus", n_bx, n1_bx, n2_bx);

   ST_Digis_RE3_Plus       = new TH1F("ST_Digis_RE3_Plus", "ST_Digis_RE3_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Ring1_Plus = new TH1F("ST_Digis_RE3_Ring1_Plus", "ST_Digis_RE3_Ring1_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Ring2_Plus = new TH1F("ST_Digis_RE3_Ring2_Plus", "ST_Digis_RE3_Ring2_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Ring3_Plus = new TH1F("ST_Digis_RE3_Ring3_Plus", "ST_Digis_RE3_Ring3_Plus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Minus       = new TH1F("ST_Digis_RE3_Minus", "ST_Digis_RE3_Minus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Ring1_Minus = new TH1F("ST_Digis_RE3_Ring1_Minus", "ST_Digis_RE3_Ring1_Minus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Ring2_Minus = new TH1F("ST_Digis_RE3_Ring2_Minus", "ST_Digis_RE3_Ring2_Minus", n_st, n1_st, n2_st);
   ST_Digis_RE3_Ring3_Minus = new TH1F("ST_Digis_RE3_Ring3_Minus", "ST_Digis_RE3_Ring3_Minus", n_st, n1_st, n2_st);



}


MyRE4nDigiAnalyzer::~MyRE4nDigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();
  BX_Digis_RE4_Plus->Write();
  BX_Digis_RE4_Ring1_Plus->Write();
  BX_Digis_RE4_Ring2_Plus->Write();
  BX_Digis_RE4_Ring3_Plus->Write();
  BX_Digis_RE4_Minus->Write();
  BX_Digis_RE4_Ring1_Minus->Write();
  BX_Digis_RE4_Ring2_Minus->Write();
  BX_Digis_RE4_Ring3_Minus->Write();

  ST_Digis_RE4_Plus->Write();
  ST_Digis_RE4_Ring1_Plus->Write();
  ST_Digis_RE4_Ring2_Plus->Write();
  ST_Digis_RE4_Ring3_Plus->Write();
  ST_Digis_RE4_Minus->Write();
  ST_Digis_RE4_Ring1_Minus->Write();
  ST_Digis_RE4_Ring2_Minus->Write();
  ST_Digis_RE4_Ring3_Minus->Write();

  BX_Digis_RE3_Plus->Write();
  BX_Digis_RE3_Ring1_Plus->Write();
  BX_Digis_RE3_Ring2_Plus->Write();
  BX_Digis_RE3_Ring3_Plus->Write();
  BX_Digis_RE3_Minus->Write();
  BX_Digis_RE3_Ring1_Minus->Write();
  BX_Digis_RE3_Ring2_Minus->Write();
  BX_Digis_RE3_Ring3_Minus->Write();

  ST_Digis_RE3_Plus->Write();
  ST_Digis_RE3_Ring1_Plus->Write();
  ST_Digis_RE3_Ring2_Plus->Write();
  ST_Digis_RE3_Ring3_Plus->Write();
  ST_Digis_RE3_Minus->Write();
  ST_Digis_RE3_Ring1_Minus->Write();
  ST_Digis_RE3_Ring2_Minus->Write();
  ST_Digis_RE3_Ring3_Minus->Write();

  BX_Digis_RE4 = new TCanvas("BX_Digis_RE4", "BX_Digis_RE4", 800, 600);
  BX_Digis_RE3 = new TCanvas("BX_Digis_RE3", "BX_Digis_RE3", 800, 600);
  ST_Digis_RE4 = new TCanvas("ST_Digis_RE4", "ST_Digis_RE4", 800, 600); 
  ST_Digis_RE3 = new TCanvas("ST_Digis_RE3", "ST_Digis_RE3", 800, 600);

  BX_Digis_RE4->cd();  BX_Digis_RE4->Divide(4,2);
  BX_Digis_RE4->cd(1); BX_Digis_RE4_Plus->Draw();        BX_Digis_RE4_Plus->GetXaxis()->SetTitle("BX [-]");        BX_Digis_RE4_Plus->GetYaxis()->SetTitle("entries [-]");        BX_Digis_RE4_Plus->SetTitle("RE+4 Digis");
  BX_Digis_RE4->cd(2); BX_Digis_RE4_Ring1_Plus->Draw();  BX_Digis_RE4_Ring1_Plus->GetXaxis()->SetTitle("BX [-]");  BX_Digis_RE4_Ring1_Plus->GetYaxis()->SetTitle("entries [-]");  BX_Digis_RE4_Ring1_Plus->SetTitle("RE+4/1 Digis");
  BX_Digis_RE4->cd(3); BX_Digis_RE4_Ring2_Plus->Draw();  BX_Digis_RE4_Ring2_Plus->GetXaxis()->SetTitle("BX [-]");  BX_Digis_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]");  BX_Digis_RE4_Ring2_Plus->SetTitle("RE+4/2 Digis");
  BX_Digis_RE4->cd(4); BX_Digis_RE4_Ring3_Plus->Draw();  BX_Digis_RE4_Ring3_Plus->GetXaxis()->SetTitle("BX [-]");  BX_Digis_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]");  BX_Digis_RE4_Ring3_Plus->SetTitle("RE+4/3 Digis");
  BX_Digis_RE4->cd(5); BX_Digis_RE4_Minus->Draw();       BX_Digis_RE4_Minus->GetXaxis()->SetTitle("BX [-]");       BX_Digis_RE4_Minus->GetYaxis()->SetTitle("entries [-]");       BX_Digis_RE4_Minus->SetTitle("RE-4 Digis");
  BX_Digis_RE4->cd(6); BX_Digis_RE4_Ring1_Minus->Draw(); BX_Digis_RE4_Ring1_Minus->GetXaxis()->SetTitle("BX [-]"); BX_Digis_RE4_Ring1_Minus->GetYaxis()->SetTitle("entries [-]"); BX_Digis_RE4_Ring1_Minus->SetTitle("RE-4/1 Digis");
  BX_Digis_RE4->cd(7); BX_Digis_RE4_Ring2_Minus->Draw(); BX_Digis_RE4_Ring2_Minus->GetXaxis()->SetTitle("BX [-]"); BX_Digis_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); BX_Digis_RE4_Ring2_Minus->SetTitle("RE-4/2 Digis");
  BX_Digis_RE4->cd(8); BX_Digis_RE4_Ring3_Minus->Draw(); BX_Digis_RE4_Ring3_Minus->GetXaxis()->SetTitle("BX [-]"); BX_Digis_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); BX_Digis_RE4_Ring3_Minus->SetTitle("RE-4/3 Digis");

  BX_Digis_RE3->cd();  BX_Digis_RE3->Divide(4,2);
  BX_Digis_RE3->cd(1); BX_Digis_RE3_Plus->Draw();        BX_Digis_RE3_Plus->GetXaxis()->SetTitle("BX [-]");        BX_Digis_RE3_Plus->GetYaxis()->SetTitle("entries [-]");        BX_Digis_RE3_Plus->SetTitle("RE+3 Digis");
  BX_Digis_RE3->cd(2); BX_Digis_RE3_Ring1_Plus->Draw();  BX_Digis_RE3_Ring1_Plus->GetXaxis()->SetTitle("BX [-]");  BX_Digis_RE3_Ring1_Plus->GetYaxis()->SetTitle("entries [-]");  BX_Digis_RE3_Ring1_Plus->SetTitle("RE+3/2 Digis");
  BX_Digis_RE3->cd(3); BX_Digis_RE3_Ring2_Plus->Draw();  BX_Digis_RE3_Ring2_Plus->GetXaxis()->SetTitle("BX [-]");  BX_Digis_RE3_Ring2_Plus->GetYaxis()->SetTitle("entries [-]");  BX_Digis_RE3_Ring2_Plus->SetTitle("RE+3/2 Digis");
  BX_Digis_RE3->cd(4); BX_Digis_RE3_Ring3_Plus->Draw();  BX_Digis_RE3_Ring3_Plus->GetXaxis()->SetTitle("BX [-]");  BX_Digis_RE3_Ring3_Plus->GetYaxis()->SetTitle("entries [-]");  BX_Digis_RE3_Ring3_Plus->SetTitle("RE+3/3 Digis");
  BX_Digis_RE3->cd(5); BX_Digis_RE3_Minus->Draw();       BX_Digis_RE3_Minus->GetXaxis()->SetTitle("BX [-]");       BX_Digis_RE3_Minus->GetYaxis()->SetTitle("entries [-]");       BX_Digis_RE3_Minus->SetTitle("RE-3 Digis");
  BX_Digis_RE3->cd(6); BX_Digis_RE3_Ring1_Minus->Draw(); BX_Digis_RE3_Ring1_Minus->GetXaxis()->SetTitle("BX [-]"); BX_Digis_RE3_Ring1_Minus->GetYaxis()->SetTitle("entries [-]"); BX_Digis_RE3_Ring1_Minus->SetTitle("RE-3/2 Digis");
  BX_Digis_RE3->cd(7); BX_Digis_RE3_Ring2_Minus->Draw(); BX_Digis_RE3_Ring2_Minus->GetXaxis()->SetTitle("BX [-]"); BX_Digis_RE3_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); BX_Digis_RE3_Ring2_Minus->SetTitle("RE-3/2 Digis");
  BX_Digis_RE3->cd(8); BX_Digis_RE3_Ring3_Minus->Draw(); BX_Digis_RE3_Ring3_Minus->GetXaxis()->SetTitle("BX [-]"); BX_Digis_RE3_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); BX_Digis_RE3_Ring3_Minus->SetTitle("RE-3/3 Digis");

  ST_Digis_RE4->cd();  ST_Digis_RE4->Divide(4,2);
  ST_Digis_RE4->cd(1); ST_Digis_RE4_Plus->Draw();        ST_Digis_RE4_Plus->GetXaxis()->SetTitle("Strip [-]");        ST_Digis_RE4_Plus->GetYaxis()->SetTitle("entries [-]");        ST_Digis_RE4_Plus->SetTitle("RE+4 Digis");
  ST_Digis_RE4->cd(2); ST_Digis_RE4_Ring1_Plus->Draw();  ST_Digis_RE4_Ring2_Plus->GetXaxis()->SetTitle("Strip [-]");  ST_Digis_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]");  ST_Digis_RE4_Ring2_Plus->SetTitle("RE+4/2 Digis");
  ST_Digis_RE4->cd(3); ST_Digis_RE4_Ring2_Plus->Draw();  ST_Digis_RE4_Ring2_Plus->GetXaxis()->SetTitle("Strip [-]");  ST_Digis_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]");  ST_Digis_RE4_Ring2_Plus->SetTitle("RE+4/2 Digis");
  ST_Digis_RE4->cd(4); ST_Digis_RE4_Ring3_Plus->Draw();  ST_Digis_RE4_Ring3_Plus->GetXaxis()->SetTitle("Strip  [-]"); ST_Digis_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]");  ST_Digis_RE4_Ring3_Plus->SetTitle("RE+4/3 Digis");
  ST_Digis_RE4->cd(5); ST_Digis_RE4_Minus->Draw();       ST_Digis_RE4_Minus->GetXaxis()->SetTitle("Strip [-]");       ST_Digis_RE4_Minus->GetYaxis()->SetTitle("entries [-]");       ST_Digis_RE4_Minus->SetTitle("RE-4 Digis");
  ST_Digis_RE4->cd(6); ST_Digis_RE4_Ring1_Minus->Draw(); ST_Digis_RE4_Ring2_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_Digis_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); ST_Digis_RE4_Ring2_Minus->SetTitle("RE-4/2 Digis");
  ST_Digis_RE4->cd(7); ST_Digis_RE4_Ring2_Minus->Draw(); ST_Digis_RE4_Ring2_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_Digis_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); ST_Digis_RE4_Ring2_Minus->SetTitle("RE-4/2 Digis");
  ST_Digis_RE4->cd(8); ST_Digis_RE4_Ring3_Minus->Draw(); ST_Digis_RE4_Ring3_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_Digis_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); ST_Digis_RE4_Ring3_Minus->SetTitle("RE-4/3 Digis");

  ST_Digis_RE3->cd();  ST_Digis_RE3->Divide(4,2);
  ST_Digis_RE3->cd(1); ST_Digis_RE3_Plus->Draw();        ST_Digis_RE3_Plus->GetXaxis()->SetTitle("Strip [-]");        ST_Digis_RE3_Plus->GetYaxis()->SetTitle("entries [-]");        ST_Digis_RE3_Plus->SetTitle("RE+3 Digis");
  ST_Digis_RE3->cd(2); ST_Digis_RE3_Ring1_Plus->Draw();  ST_Digis_RE3_Ring1_Plus->GetXaxis()->SetTitle("Strip [-]");  ST_Digis_RE3_Ring1_Plus->GetYaxis()->SetTitle("entries [-]");  ST_Digis_RE3_Ring1_Plus->SetTitle("RE+3/2 Digis");
  ST_Digis_RE3->cd(3); ST_Digis_RE3_Ring2_Plus->Draw();  ST_Digis_RE3_Ring2_Plus->GetXaxis()->SetTitle("Strip [-]");  ST_Digis_RE3_Ring2_Plus->GetYaxis()->SetTitle("entries [-]");  ST_Digis_RE3_Ring2_Plus->SetTitle("RE+3/2 Digis");
  ST_Digis_RE3->cd(4); ST_Digis_RE3_Ring3_Plus->Draw();  ST_Digis_RE3_Ring3_Plus->GetXaxis()->SetTitle("Strip [-]");  ST_Digis_RE3_Ring3_Plus->GetYaxis()->SetTitle("entries [-]");  ST_Digis_RE3_Ring3_Plus->SetTitle("RE+3/3 Digis");
  ST_Digis_RE3->cd(5); ST_Digis_RE3_Minus->Draw();       ST_Digis_RE3_Minus->GetXaxis()->SetTitle("Strip [-]");       ST_Digis_RE3_Minus->GetYaxis()->SetTitle("entries [-]");       ST_Digis_RE3_Minus->SetTitle("RE-3 Digis");
  ST_Digis_RE3->cd(6); ST_Digis_RE3_Ring1_Minus->Draw(); ST_Digis_RE3_Ring1_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_Digis_RE3_Ring1_Minus->GetYaxis()->SetTitle("entries [-]"); ST_Digis_RE3_Ring1_Minus->SetTitle("RE-3/2 Digis");
  ST_Digis_RE3->cd(7); ST_Digis_RE3_Ring2_Minus->Draw(); ST_Digis_RE3_Ring2_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_Digis_RE3_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); ST_Digis_RE3_Ring2_Minus->SetTitle("RE-3/2 Digis");
  ST_Digis_RE3->cd(8); ST_Digis_RE3_Ring3_Minus->Draw(); ST_Digis_RE3_Ring3_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_Digis_RE3_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); ST_Digis_RE3_Ring3_Minus->SetTitle("RE-3/3 Digis");

  
 BX_Digis_RE4->Write();
 BX_Digis_RE3->Write();
 ST_Digis_RE4->Write();
 ST_Digis_RE3->Write();


}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyRE4nDigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<RPCDigiCollection> rpcdigis;
  iEvent.getByLabel(digiLabel, "", rpcdigis);


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
    std::cout<<" "<<std::endl;
    std::cout<<" RPC DetId: "<<std::setw(12)<<id<<" a.k.a. "<<std::setw(18)<<nameRoll<<" which is in "<<std::setw(5)<<wheelOrDiskType<<" "<<std::setw(2)<<wheelOrDiskNumber<<" ";
    std::cout<<std::setw(7)<<stationOrRingType<<" "<<std::setw(2)<<stationOrRingNumber<<" sector "<<std::setw(2)<<sector<<std::endl;
    std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;
    RPCDigiCollection::const_iterator digiItr; 
    //loop on digis of given roll
    for (digiItr =(*collectionItr ).second.first;digiItr != (*collectionItr ).second.second; ++digiItr){
      int strip= (*digiItr).strip();
      int bx=(*digiItr).bx();
      std::cout<<"     Digi: strip = "<<std::setw(2)<<strip<<" bx = "<<std::setw(2)<<bx<<std::endl;
      // Fill here your histograms

      // Positive Endcap
      if(region == 1 && station == 4) {
	BX_Digis_RE4_Plus->Fill(bx); ST_Digis_RE4_Plus->Fill(strip);
	// Rings
	if(ring==1) { BX_Digis_RE4_Ring1_Plus->Fill(bx);  ST_Digis_RE4_Ring1_Plus->Fill(strip);} 
	if(ring==2) { BX_Digis_RE4_Ring2_Plus->Fill(bx);  ST_Digis_RE4_Ring2_Plus->Fill(strip);} 
	if(ring==3) { BX_Digis_RE4_Ring3_Plus->Fill(bx);  ST_Digis_RE4_Ring3_Plus->Fill(strip);} 
      }
      if(region == 1 && station == 3) {
        BX_Digis_RE3_Plus->Fill(bx); ST_Digis_RE3_Plus->Fill(strip);
        // Rings
        if(ring==1) { BX_Digis_RE3_Ring1_Plus->Fill(bx);  ST_Digis_RE3_Ring1_Plus->Fill(strip);}
        if(ring==2) { BX_Digis_RE3_Ring2_Plus->Fill(bx);  ST_Digis_RE3_Ring2_Plus->Fill(strip);}
        if(ring==3) { BX_Digis_RE3_Ring3_Plus->Fill(bx);  ST_Digis_RE3_Ring3_Plus->Fill(strip);}
      }
      // Negative Endcap
      if(region == -1 && station == 4) {
	BX_Digis_RE4_Minus->Fill(bx); ST_Digis_RE4_Minus->Fill(strip);
	// Rings
	if(ring==1) { BX_Digis_RE4_Ring1_Minus->Fill(bx);  ST_Digis_RE4_Ring1_Minus->Fill(strip);} 
	if(ring==2) { BX_Digis_RE4_Ring2_Minus->Fill(bx);  ST_Digis_RE4_Ring2_Minus->Fill(strip);} 
	if(ring==3) { BX_Digis_RE4_Ring3_Minus->Fill(bx);  ST_Digis_RE4_Ring3_Minus->Fill(strip);} 
      }
      if(region == -1 && station == 3) {
	BX_Digis_RE3_Minus->Fill(bx); ST_Digis_RE3_Minus->Fill(strip);
	// Rings
	if(ring==1) { BX_Digis_RE3_Ring1_Minus->Fill(bx);  ST_Digis_RE3_Ring1_Minus->Fill(strip);} 
	if(ring==2) { BX_Digis_RE3_Ring2_Minus->Fill(bx);  ST_Digis_RE3_Ring2_Minus->Fill(strip);} 
	if(ring==3) { BX_Digis_RE3_Ring3_Minus->Fill(bx);  ST_Digis_RE3_Ring3_Minus->Fill(strip);} 
      }

    }
    std::cout<<" ---------------------------------------------------------------------------------------------"<<std::endl;

  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRE4nDigiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRE4nDigiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyRE4nDigiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
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
MyRE4nDigiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyRE4nDigiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyRE4nDigiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRE4nDigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRE4nDigiAnalyzer);
