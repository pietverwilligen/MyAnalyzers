// -*- C++ -*-
//
// Package:    MyRE4nRecHitAnalyzer
// Class:      MyRE4nRecHitAnalyzer
// 
/**\class MyRE4nRecHitAnalyzer MyRE4nRecHitAnalyzer.cc MyAnalyzers/MyRE4nRecHitAnalyzer/src/MyRE4nRecHitAnalyzer.cc

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

class MyRE4nRecHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRE4nRecHitAnalyzer(const edm::ParameterSet&);
      ~MyRE4nRecHitAnalyzer();

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
  TH1F * BX_RecHits_RE4_Plus,  * BX_RecHits_RE4_Ring2_Plus,  * BX_RecHits_RE4_Ring3_Plus,  * ST_RecHits_RE4_Plus,  * ST_RecHits_RE4_Ring2_Plus,  * ST_RecHits_RE4_Ring3_Plus;
  TH1F * BX_RecHits_RE4_Minus, * BX_RecHits_RE4_Ring2_Minus, * BX_RecHits_RE4_Ring3_Minus, * ST_RecHits_RE4_Minus, * ST_RecHits_RE4_Ring2_Minus, * ST_RecHits_RE4_Ring3_Minus;
  std::vector<double> x_p, y_p, z_p, x_n, y_n, z_n, r_pr, z_pr, r_nr, z_nr;
  TGraph  * RE4_Plus_XY_All, * RE4_Minus_XY_All, * RE_Minus_YZ_All, * RE_Plus_YZ_All;
  TCanvas * Canvas_RE4_Plus_XY, * Canvas_RE4_Minus_XY, * Canvas_RE_Plus_YZ, * Canvas_RE_Minus_YZ;
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
MyRE4nRecHitAnalyzer::MyRE4nRecHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");

  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

  BX_RecHits_RE4_Plus       = new TH1F("BX_RecHits_RE4_Plus", "BX_RecHits_RE4_Plus", n_bx, n1_bx, n2_bx);
  BX_RecHits_RE4_Ring2_Plus = new TH1F("BX_RecHits_RE4_Ring2_Plus", "BX_RecHits_RE4_Ring2_Plus", n_bx, n1_bx, n2_bx);
  BX_RecHits_RE4_Ring3_Plus = new TH1F("BX_RecHits_RE4_Ring3_Plus", "BX_RecHits_RE4_Ring3_Plus", n_bx, n1_bx, n2_bx);
  BX_RecHits_RE4_Minus       = new TH1F("BX_RecHits_RE4_Minus", "BX_RecHits_RE4_Minus", n_bx, n1_bx, n2_bx);
  BX_RecHits_RE4_Ring2_Minus = new TH1F("BX_RecHits_RE4_Ring2_Minus", "BX_RecHits_RE4_Ring2_Minus", n_bx, n1_bx, n2_bx);
  BX_RecHits_RE4_Ring3_Minus = new TH1F("BX_RecHits_RE4_Ring3_Minus", "BX_RecHits_RE4_Ring3_Minus", n_bx, n1_bx, n2_bx);

  ST_RecHits_RE4_Plus       = new TH1F("ST_RecHits_RE4_Plus", "ST_RecHits_RE4_Plus", n_st, n1_st, n2_st);
  ST_RecHits_RE4_Ring2_Plus = new TH1F("ST_RecHits_RE4_Ring2_Plus", "ST_RecHits_RE4_Ring2_Plus", n_st, n1_st, n2_st);   
  ST_RecHits_RE4_Ring3_Plus = new TH1F("ST_RecHits_RE4_Ring3_Plus", "ST_RecHits_RE4_Ring3_Plus", n_st, n1_st, n2_st);
  ST_RecHits_RE4_Minus       = new TH1F("ST_RecHits_RE4_Minus", "ST_RecHits_RE4_Minus", n_st, n1_st, n2_st);
  ST_RecHits_RE4_Ring2_Minus = new TH1F("ST_RecHits_RE4_Ring2_Minus", "ST_RecHits_RE4_Ring2_Minus", n_st, n1_st, n2_st);
  ST_RecHits_RE4_Ring3_Minus = new TH1F("ST_RecHits_RE4_Ring3_Minus", "ST_RecHits_RE4_Ring3_Minus", n_st, n1_st, n2_st);

}


MyRE4nRecHitAnalyzer::~MyRE4nRecHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();
  BX_RecHits_RE4_Plus->Write();
  BX_RecHits_RE4_Ring2_Plus->Write();
  BX_RecHits_RE4_Ring3_Plus->Write();
  BX_RecHits_RE4_Minus->Write();
  BX_RecHits_RE4_Ring2_Minus->Write();
  BX_RecHits_RE4_Ring3_Minus->Write();

  ST_RecHits_RE4_Plus->Write();
  ST_RecHits_RE4_Ring2_Plus->Write();
  ST_RecHits_RE4_Ring3_Plus->Write();
  ST_RecHits_RE4_Minus->Write();
  ST_RecHits_RE4_Ring2_Minus->Write();
  ST_RecHits_RE4_Ring3_Minus->Write();


  TCanvas * BX_RecHits_RE4 = new TCanvas("BX_RecHits_RE4", "BX_RecHits_RE4", 800, 600);
  TCanvas * ST_RecHits_RE4 = new TCanvas("ST_RecHits_RE4", "ST_RecHits_RE4", 800, 600);
  BX_RecHits_RE4->cd(); BX_RecHits_RE4->Divide(3,2);
  BX_RecHits_RE4->cd(1); BX_RecHits_RE4_Plus->Draw(); BX_RecHits_RE4_Plus->GetXaxis()->SetTitle("BX [-]"); BX_RecHits_RE4_Plus->GetYaxis()->SetTitle("entries [-]"); BX_RecHits_RE4_Plus->SetTitle("RE+4 Digis");
  BX_RecHits_RE4->cd(2); BX_RecHits_RE4_Ring2_Plus->Draw(); BX_RecHits_RE4_Ring2_Plus->GetXaxis()->SetTitle("BX [-]"); BX_RecHits_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]"); BX_RecHits_RE4_Ring2_Plus->SetTitle("RE+4/2 Digis");
  BX_RecHits_RE4->cd(3); BX_RecHits_RE4_Ring3_Plus->Draw(); BX_RecHits_RE4_Ring3_Plus->GetXaxis()->SetTitle("BX [-]"); BX_RecHits_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]"); BX_RecHits_RE4_Ring3_Plus->SetTitle("RE+4/3 Digis");
  BX_RecHits_RE4->cd(4); BX_RecHits_RE4_Minus->Draw();      BX_RecHits_RE4_Minus->GetXaxis()->SetTitle("BX [-]"); BX_RecHits_RE4_Minus->GetYaxis()->SetTitle("entries [-]"); BX_RecHits_RE4_Minus->SetTitle("RE-4 Digis");
  BX_RecHits_RE4->cd(5);BX_RecHits_RE4_Ring2_Minus->Draw();BX_RecHits_RE4_Ring2_Minus->GetXaxis()->SetTitle("BX [-]");BX_RecHits_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); BX_RecHits_RE4_Ring2_Minus->SetTitle("RE-4/2 Digis");
  BX_RecHits_RE4->cd(6);BX_RecHits_RE4_Ring3_Minus->Draw();BX_RecHits_RE4_Ring3_Minus->GetXaxis()->SetTitle("BX [-]");BX_RecHits_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); BX_RecHits_RE4_Ring3_Minus->SetTitle("RE-4/3 Digis");

  ST_RecHits_RE4->cd(); ST_RecHits_RE4->Divide(3,2);
  ST_RecHits_RE4->cd(1); ST_RecHits_RE4_Plus->Draw(); ST_RecHits_RE4_Plus->GetXaxis()->SetTitle("Strip [-]"); ST_RecHits_RE4_Plus->GetYaxis()->SetTitle("entries [-]"); ST_RecHits_RE4_Plus->SetTitle("RE+4 Digis");
  ST_RecHits_RE4->cd(2); ST_RecHits_RE4_Ring2_Plus->Draw();  ST_RecHits_RE4_Ring2_Plus->GetXaxis()->SetTitle("Strip [-]"); ST_RecHits_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]"); ST_RecHits_RE4_Ring2_Plus->SetTitle("RE+4/2 Digis");
  ST_RecHits_RE4->cd(3); ST_RecHits_RE4_Ring3_Plus->Draw();  ST_RecHits_RE4_Ring3_Plus->GetXaxis()->SetTitle("Strip  [-]"); ST_RecHits_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]"); ST_RecHits_RE4_Ring3_Plus->SetTitle("RE+4/3 Digis");
  ST_RecHits_RE4->cd(4); ST_RecHits_RE4_Minus->Draw();       ST_RecHits_RE4_Minus->GetXaxis()->SetTitle("Strip [-]"); ST_RecHits_RE4_Minus->GetYaxis()->SetTitle("entries [-]"); ST_RecHits_RE4_Minus->SetTitle("RE-4 Digis");
  ST_RecHits_RE4->cd(5); ST_RecHits_RE4_Ring2_Minus->Draw(); ST_RecHits_RE4_Ring2_Minus->GetXaxis()->SetTitle("Strip [-]");ST_RecHits_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); ST_RecHits_RE4_Ring2_Minus->SetTitle("RE-4/2 Digis");
  ST_RecHits_RE4->cd(6); ST_RecHits_RE4_Ring3_Minus->Draw(); ST_RecHits_RE4_Ring3_Minus->GetXaxis()->SetTitle("Strip [-]");ST_RecHits_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); ST_RecHits_RE4_Ring3_Minus->SetTitle("RE-4/3 Digis");

  const int n_p = x_p.size();  double x_ap[n_p]; double y_ap[n_p];  // double z_ap[n_p];
  const int n_n = x_n.size();  double x_an[n_n]; double y_an[n_n];  // double z_an[n_n];
  const int n_pr = r_pr.size();  double r_apr[n_pr]; double z_apr[n_pr];
  const int n_nr = r_nr.size();  double r_anr[n_nr]; double z_anr[n_nr];

  for(int i=0; i< n_p; ++i) { x_ap[i] = x_p[i]; y_ap[i] = y_p[i]; /*z_ap[i] = z_p[i];*/ }
  for(int i=0; i< n_n; ++i) { x_an[i] = x_n[i]; y_an[i] = y_n[i]; /*z_an[i] = z_n[i];*/ }
  for(int i=0; i< n_pr; ++i) { r_apr[i] = r_pr[i]; z_apr[i] = z_pr[i]; }
  for(int i=0; i< n_nr; ++i) { r_anr[i] = r_nr[i]; z_anr[i] = z_nr[i]; }
  RE4_Plus_XY_All  = new TGraph(n_p, x_ap, y_ap); std::cout<<"RE+4 All SimHits: "<<n_p<<std::endl;
  RE4_Minus_XY_All = new TGraph(n_n, x_an, y_an); std::cout<<"RE-4 All SimHits: "<<n_n<<std::endl;
  RE_Plus_YZ_All   = new TGraph(n_pr, z_apr, r_apr);
  RE_Minus_YZ_All  = new TGraph(n_nr, z_anr, r_anr);

  Canvas_RE4_Plus_XY  = new TCanvas("Canvas_RE4_Plus_XY",  "Canvas_RE4_Plus_XY", 800, 600);
  Canvas_RE4_Minus_XY = new TCanvas("Canvas_RE4_Minus_XY", "Canvas_RE4_Minus_XY", 800, 600);
  Canvas_RE_Plus_YZ   = new TCanvas("Canvas_RE_Plus_YZ",   "Canvas_RE_Plus_YZ", 800, 600);
  Canvas_RE_Minus_YZ  = new TCanvas("Canvas_RE_Minus_YZ",  "Canvas_RE_Minus_YZ", 800, 600);

  // XY and RZ Graphs
  Canvas_RE4_Plus_XY->cd();
  RE4_Plus_XY_All->SetMarkerStyle(5);  RE4_Plus_XY_All->Draw("AP");  RE4_Plus_XY_All->GetXaxis()->SetTitle("X [cm]");  RE4_Plus_XY_All->GetYaxis()->SetTitle("Y [cm]");   RE4_Plus_XY_All->SetTitle("RE+4 RecHits");
  RE4_Plus_XY_All->Draw("AP");
  Canvas_RE4_Minus_XY->cd();
  RE4_Minus_XY_All->SetMarkerStyle(5); RE4_Minus_XY_All->Draw("AP"); RE4_Minus_XY_All->GetXaxis()->SetTitle("X [cm]"); RE4_Minus_XY_All->GetYaxis()->SetTitle("Y [cm]");  RE4_Minus_XY_All->SetTitle("RE-4 RecHits");
  RE4_Minus_XY_All->Draw("AP");  
  Canvas_RE_Plus_YZ->cd();  RE_Plus_YZ_All->SetMarkerStyle(5);    RE_Plus_YZ_All->Draw("AP");    RE_Plus_YZ_All->GetXaxis()->SetTitle("Z [cm]");    RE_Plus_YZ_All->GetYaxis()->SetTitle("R [cm]");     RE_Plus_YZ_All->SetTitle("RE RecHits");
  Canvas_RE_Minus_YZ->cd(); RE_Minus_YZ_All->SetMarkerStyle(5);   RE_Minus_YZ_All->Draw("AP");   RE_Minus_YZ_All->GetXaxis()->SetTitle("Z [cm]");   RE_Minus_YZ_All->GetYaxis()->SetTitle("R [cm]");    RE_Minus_YZ_All->SetTitle("RE RecHits");

  BX_RecHits_RE4->Write();
  ST_RecHits_RE4->Write();

  Canvas_RE4_Plus_XY->Write();
  Canvas_RE4_Minus_XY->Write();
  Canvas_RE_Plus_YZ->Write();
  Canvas_RE_Minus_YZ->Write();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyRE4nRecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   
  // int digis_Barrel=0;
  // int digis_Endcap=0;
 
  // std::cout<<"The Number of RecHits is "<<nRPC<<std::endl;       
  for (recHit = rpcRecHits->begin(); recHit != rpcRecHits->end(); recHit++) {
    RPCDetId rollId = (RPCDetId)(*recHit).rpcId();
    RPCGeomServ rpcsrv(rollId);
    LocalPoint recHitPos=recHit->localPosition();
    const RPCRoll* rollasociated = rpcGeom->roll(rollId);
    const BoundPlane & RPCSurface = rollasociated->surface(); 
    GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal(recHitPos);
    std::cout<<"RPC Rec Hit in "<<rpcsrv.name()<<" bx = "<<recHit->BunchX()<<" Global Position = "<<RPCGlobalPoint<<std::endl;
    int region = rollId.region();
    int station = rollId.station();
    int ring    = rollId.ring();


    int bx = recHit->BunchX();
    // int cl = recHit->clusterSize();
    int st = recHit->firstClusterStrip();


    // Positive Endcap
    if(region == 1) { r_pr.push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2))); z_pr.push_back(RPCGlobalPoint.z()); }
    if(region == 1 && station == 4) {
      BX_RecHits_RE4_Plus->Fill(bx); ST_RecHits_RE4_Plus->Fill(st);
      x_p.push_back(RPCGlobalPoint.x()); y_p.push_back(RPCGlobalPoint.y()); z_p.push_back(RPCGlobalPoint.z());
      // Rings
      if(ring==2) { BX_RecHits_RE4_Ring2_Plus->Fill(bx);  ST_RecHits_RE4_Ring2_Plus->Fill(st);}
      if(ring==3) { BX_RecHits_RE4_Ring3_Plus->Fill(bx);  ST_RecHits_RE4_Ring3_Plus->Fill(st);}
    }
    // Negative Endcap
    if(region == -1) { r_nr.push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2))); z_nr.push_back(RPCGlobalPoint.z()); } 
    if(region == -1 && station == 4) {
      BX_RecHits_RE4_Minus->Fill(bx); ST_RecHits_RE4_Minus->Fill(st);
      x_n.push_back(RPCGlobalPoint.x()); y_n.push_back(RPCGlobalPoint.y()); z_n.push_back(RPCGlobalPoint.z());
      // Rings
      if(ring==2) { BX_RecHits_RE4_Ring2_Minus->Fill(bx);  ST_RecHits_RE4_Ring2_Minus->Fill(st);}
      if(ring==3) { BX_RecHits_RE4_Ring3_Minus->Fill(bx);  ST_RecHits_RE4_Ring3_Minus->Fill(st);}
    }


  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRE4nRecHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRE4nRecHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyRE4nRecHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
}

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyRE4nRecHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyRE4nRecHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyRE4nRecHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRE4nRecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRE4nRecHitAnalyzer);
