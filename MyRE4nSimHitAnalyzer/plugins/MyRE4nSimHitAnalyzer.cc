// -*- C++ -*-
//
// Package:    MyRE4nSimHitAnalyzer
// Class:      MyRE4nSimHitAnalyzer
// 
/**\class MyRE4nSimHitAnalyzer MyRE4nSimHitAnalyzer.cc MyAnalyzers/MyRE4nSimHitAnalyzer/src/MyRE4nSimHitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen,161 R-006,+41227676292,
//         Created:  Wed Oct 10 17:36:38 CEST 2012
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
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// why is this one not included by default
#include "FWCore/Framework/interface/ESHandle.h"


// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/GeometrySurface/interface/Surface.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>


#include "DataFormats/Common/interface/Handle.h"
// Simulation
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "FastSimulation/Tracking/test/FastTrackAnalyzer.h"
// DetIds
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
// Digis
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
// RecHits
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
// Candidates
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
// Tracking
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
// Track
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Math/interface/LorentzVectorFwd.h"

//
// class declaration
//

class MyRE4nSimHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRE4nSimHitAnalyzer(const edm::ParameterSet&);
      ~MyRE4nSimHitAnalyzer();
  edm::ESHandle <RPCGeometry> rpcGeom;
  edm::ESHandle <CSCGeometry> cscGeom;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  std::string rootFileName;
  TFile * outputfile;
  TH1F * TOF_SimHits_RE4_Plus,  * TOF_SimHits_RE4_Ring2_Plus,  * TOF_SimHits_RE4_Ring3_Plus;  
  TH1F * TOF_SimHits_RE4_Minus, * TOF_SimHits_RE4_Ring2_Minus, * TOF_SimHits_RE4_Ring3_Minus;

  TH1F * TOF_XCheck1_SimHits_RE4_Plus,  * TOF_XCheck1_SimHits_RE4_Ring2_Plus,  * TOF_XCheck1_SimHits_RE4_Ring3_Plus;  
  TH1F * TOF_XCheck1_SimHits_RE4_Minus, * TOF_XCheck1_SimHits_RE4_Ring2_Minus, * TOF_XCheck1_SimHits_RE4_Ring3_Minus;
  TH1F * TOF_XCheck2_SimHits_RE4_Plus,  * TOF_XCheck2_SimHits_RE4_Ring2_Plus,  * TOF_XCheck2_SimHits_RE4_Ring3_Plus;  
  TH1F * TOF_XCheck2_SimHits_RE4_Minus, * TOF_XCheck2_SimHits_RE4_Ring2_Minus, * TOF_XCheck2_SimHits_RE4_Ring3_Minus;
  TH1F * TOF_Difference1_SimHits_RE4_Plus,  * TOF_Difference1_SimHits_RE4_Ring2_Plus,  * TOF_Difference1_SimHits_RE4_Ring3_Plus;  
  TH1F * TOF_Difference1_SimHits_RE4_Minus, * TOF_Difference1_SimHits_RE4_Ring2_Minus, * TOF_Difference1_SimHits_RE4_Ring3_Minus;
  TH1F * TOF_Difference2_SimHits_RE4_Plus,  * TOF_Difference2_SimHits_RE4_Ring2_Plus,  * TOF_Difference2_SimHits_RE4_Ring3_Plus;  
  TH1F * TOF_Difference2_SimHits_RE4_Minus, * TOF_Difference2_SimHits_RE4_Ring2_Minus, * TOF_Difference2_SimHits_RE4_Ring3_Minus;

  TH1F * TOF_SimHits_RE3_Plus,  * TOF_SimHits_RE3_Ring2_Plus,  * TOF_SimHits_RE3_Ring3_Plus;
  TH1F * TOF_SimHits_RE3_Minus, * TOF_SimHits_RE3_Ring2_Minus, * TOF_SimHits_RE3_Ring3_Minus;

  TH1F * TOF_SimHits_RE4_Ring2_A_Plus, * TOF_SimHits_RE4_Ring2_B_Plus, * TOF_SimHits_RE4_Ring2_C_Plus, * TOF_SimHits_RE4_Ring3_A_Plus, * TOF_SimHits_RE4_Ring3_B_Plus, * TOF_SimHits_RE4_Ring3_C_Plus;
  TH1F * TOF_SimHits_RE4_Ring2_A_Minus, * TOF_SimHits_RE4_Ring2_B_Minus, * TOF_SimHits_RE4_Ring2_C_Minus, * TOF_SimHits_RE4_Ring3_A_Minus, * TOF_SimHits_RE4_Ring3_B_Minus, * TOF_SimHits_RE4_Ring3_C_Minus;

  TH1F * VTX_PX, * VTX_PY, *VTX_PZ, * VTX_AX, * VTX_AY, *VTX_AZ; 

  std::vector<double> x_p, y_p, z_p, x_n, y_n, z_n, r_r, z_r;
  TGraph  * RE4_Plus_XY, * RE4_Minus_XY, * RE_YZ; 
  TCanvas * Canvas_RE4_Plus_XY, * Canvas_RE4_Minus_XY, * Canvas_RE_YZ, * Canvas_RE4_TOF,  * Canvas_RE4_TOF_XCheck, * Canvas_RE4_TOF_Difference, * Canvas_RE3_TOF, * Canvas_VTX; 

};

//
// constants, enums and typedefs
//
int n_tof  = 75;  double n1_tof  = 30,    n2_tof = 45;
int n_dif  = 75;  double n1_dif  = -7.5,  n2_dif = 7.5; 
int n_v_x  = 40;  double n1_v_x  = -2.0,  n2_v_x = 2.0;
int n_v_z  = 50;  double n1_v_z  = -25.0, n2_v_z = 25.0;
// pu  seems strange
int p_v_x  = 50;  double p1_v_x  = -500.0,  p2_v_x = 500.0;
int p_v_z  = 125;  double p1_v_z  = -1250.0, p2_v_z = 1250.0;


//
// static data member definitions
//

//
// constructors and destructor
//
MyRE4nSimHitAnalyzer::MyRE4nSimHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");

  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

  TOF_SimHits_RE4_Plus        = new TH1F("TOF_SimHits_RE4_Plus",        "TOF_SimHits_RE4_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_Plus  = new TH1F("TOF_SimHits_RE4_Ring2_Plus",  "TOF_SimHits_RE4_Ring2_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_Plus  = new TH1F("TOF_SimHits_RE4_Ring3_Plus",  "TOF_SimHits_RE4_Ring3_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Minus       = new TH1F("TOF_SimHits_RE4_Minus",       "TOF_SimHits_RE4_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_Minus = new TH1F("TOF_SimHits_RE4_Ring2_Minus", "TOF_SimHits_RE4_Ring2_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_Minus = new TH1F("TOF_SimHits_RE4_Ring3_Minus", "TOF_SimHits_RE4_Ring3_Minus", n_tof, n1_tof, n2_tof);

  TOF_XCheck1_SimHits_RE4_Plus        = new TH1F("TOF_XCheck1_SimHits_RE4_Plus",        "TOF_XCheck1_SimHits_RE4_Plus", n_tof, n1_tof, n2_tof);
  TOF_XCheck1_SimHits_RE4_Ring2_Plus  = new TH1F("TOF_XCheck1_SimHits_RE4_Ring2_Plus",  "TOF_XCheck1_SimHits_RE4_Ring2_Plus", n_tof, n1_tof, n2_tof);
  TOF_XCheck1_SimHits_RE4_Ring3_Plus  = new TH1F("TOF_XCheck1_SimHits_RE4_Ring3_Plus",  "TOF_XCheck1_SimHits_RE4_Ring3_Plus", n_tof, n1_tof, n2_tof);
  TOF_XCheck1_SimHits_RE4_Minus       = new TH1F("TOF_XCheck1_SimHits_RE4_Minus",       "TOF_XCheck1_SimHits_RE4_Minus", n_tof, n1_tof, n2_tof);
  TOF_XCheck1_SimHits_RE4_Ring2_Minus = new TH1F("TOF_XCheck1_SimHits_RE4_Ring2_Minus", "TOF_XCheck1_SimHits_RE4_Ring2_Minus", n_tof, n1_tof, n2_tof);
  TOF_XCheck1_SimHits_RE4_Ring3_Minus = new TH1F("TOF_XCheck1_SimHits_RE4_Ring3_Minus", "TOF_XCheck1_SimHits_RE4_Ring3_Minus", n_tof, n1_tof, n2_tof);
  TOF_XCheck2_SimHits_RE4_Plus        = new TH1F("TOF_XCheck2_SimHits_RE4_Plus",        "TOF_XCheck2_SimHits_RE4_Plus", n_tof, n1_tof, n2_tof);
  TOF_XCheck2_SimHits_RE4_Ring2_Plus  = new TH1F("TOF_XCheck2_SimHits_RE4_Ring2_Plus",  "TOF_XCheck2_SimHits_RE4_Ring2_Plus", n_tof, n1_tof, n2_tof);
  TOF_XCheck2_SimHits_RE4_Ring3_Plus  = new TH1F("TOF_XCheck2_SimHits_RE4_Ring3_Plus",  "TOF_XCheck2_SimHits_RE4_Ring3_Plus", n_tof, n1_tof, n2_tof);
  TOF_XCheck2_SimHits_RE4_Minus       = new TH1F("TOF_XCheck2_SimHits_RE4_Minus",       "TOF_XCheck2_SimHits_RE4_Minus", n_tof, n1_tof, n2_tof);
  TOF_XCheck2_SimHits_RE4_Ring2_Minus = new TH1F("TOF_XCheck2_SimHits_RE4_Ring2_Minus", "TOF_XCheck2_SimHits_RE4_Ring2_Minus", n_tof, n1_tof, n2_tof);
  TOF_XCheck2_SimHits_RE4_Ring3_Minus = new TH1F("TOF_XCheck2_SimHits_RE4_Ring3_Minus", "TOF_XCheck2_SimHits_RE4_Ring3_Minus", n_tof, n1_tof, n2_tof);

  TOF_Difference1_SimHits_RE4_Plus        = new TH1F("TOF_Difference1_SimHits_RE4_Plus",        "TOF_Difference1_SimHits_RE4_Plus", n_dif, n1_dif, n2_dif);
  TOF_Difference1_SimHits_RE4_Ring2_Plus  = new TH1F("TOF_Difference1_SimHits_RE4_Ring2_Plus",  "TOF_Difference1_SimHits_RE4_Ring2_Plus", n_dif, n1_dif, n2_dif);
  TOF_Difference1_SimHits_RE4_Ring3_Plus  = new TH1F("TOF_Difference1_SimHits_RE4_Ring3_Plus",  "TOF_Difference1_SimHits_RE4_Ring3_Plus", n_dif, n1_dif, n2_dif);
  TOF_Difference1_SimHits_RE4_Minus       = new TH1F("TOF_Difference1_SimHits_RE4_Minus",       "TOF_Difference1_SimHits_RE4_Minus", n_dif, n1_dif, n2_dif);
  TOF_Difference1_SimHits_RE4_Ring2_Minus = new TH1F("TOF_Difference1_SimHits_RE4_Ring2_Minus", "TOF_Difference1_SimHits_RE4_Ring2_Minus", n_dif, n1_dif, n2_dif);
  TOF_Difference1_SimHits_RE4_Ring3_Minus = new TH1F("TOF_Difference1_SimHits_RE4_Ring3_Minus", "TOF_Difference1_SimHits_RE4_Ring3_Minus", n_dif, n1_dif, n2_dif);
  TOF_Difference2_SimHits_RE4_Plus        = new TH1F("TOF_Difference2_SimHits_RE4_Plus",        "TOF_Difference2_SimHits_RE4_Plus", n_dif, n1_dif, n2_dif);
  TOF_Difference2_SimHits_RE4_Ring2_Plus  = new TH1F("TOF_Difference2_SimHits_RE4_Ring2_Plus",  "TOF_Difference2_SimHits_RE4_Ring2_Plus", n_dif, n1_dif, n2_dif);
  TOF_Difference2_SimHits_RE4_Ring3_Plus  = new TH1F("TOF_Difference2_SimHits_RE4_Ring3_Plus",  "TOF_Difference2_SimHits_RE4_Ring3_Plus", n_dif, n1_dif, n2_dif);
  TOF_Difference2_SimHits_RE4_Minus       = new TH1F("TOF_Difference2_SimHits_RE4_Minus",       "TOF_Difference2_SimHits_RE4_Minus", n_dif, n1_dif, n2_dif);
  TOF_Difference2_SimHits_RE4_Ring2_Minus = new TH1F("TOF_Difference2_SimHits_RE4_Ring2_Minus", "TOF_Difference2_SimHits_RE4_Ring2_Minus", n_dif, n1_dif, n2_dif);
  TOF_Difference2_SimHits_RE4_Ring3_Minus = new TH1F("TOF_Difference2_SimHits_RE4_Ring3_Minus", "TOF_Difference2_SimHits_RE4_Ring3_Minus", n_dif, n1_dif, n2_dif);


  TOF_SimHits_RE3_Plus        = new TH1F("TOF_SimHits_RE3_Plus",        "TOF_SimHits_RE3_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE3_Ring2_Plus  = new TH1F("TOF_SimHits_RE3_Ring2_Plus",  "TOF_SimHits_RE3_Ring2_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE3_Ring3_Plus  = new TH1F("TOF_SimHits_RE3_Ring3_Plus",  "TOF_SimHits_RE3_Ring3_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE3_Minus       = new TH1F("TOF_SimHits_RE3_Minus",       "TOF_SimHits_RE3_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE3_Ring2_Minus = new TH1F("TOF_SimHits_RE3_Ring2_Minus", "TOF_SimHits_RE3_Ring2_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE3_Ring3_Minus = new TH1F("TOF_SimHits_RE3_Ring3_Minus", "TOF_SimHits_RE3_Ring3_Minus", n_tof, n1_tof, n2_tof);

  TOF_SimHits_RE4_Ring2_A_Plus  = new TH1F("TOF_SimHits_RE4_Ring2_A_Plus",  "TOF_SimHits_RE4_Ring2_A_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_A_Plus  = new TH1F("TOF_SimHits_RE4_Ring3_A_Plus",  "TOF_SimHits_RE4_Ring3_A_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_B_Plus  = new TH1F("TOF_SimHits_RE4_Ring2_B_Plus",  "TOF_SimHits_RE4_Ring2_B_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_B_Plus  = new TH1F("TOF_SimHits_RE4_Ring3_B_Plus",  "TOF_SimHits_RE4_Ring3_B_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_C_Plus  = new TH1F("TOF_SimHits_RE4_Ring2_C_Plus",  "TOF_SimHits_RE4_Ring2_C_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_C_Plus  = new TH1F("TOF_SimHits_RE4_Ring3_C_Plus",  "TOF_SimHits_RE4_Ring3_C_Plus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_A_Minus  = new TH1F("TOF_SimHits_RE4_Ring2_A_Minus",  "TOF_SimHits_RE4_Ring2_A_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_A_Minus  = new TH1F("TOF_SimHits_RE4_Ring3_A_Minus",  "TOF_SimHits_RE4_Ring3_A_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_B_Minus  = new TH1F("TOF_SimHits_RE4_Ring2_B_Minus",  "TOF_SimHits_RE4_Ring2_B_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_B_Minus  = new TH1F("TOF_SimHits_RE4_Ring3_B_Minus",  "TOF_SimHits_RE4_Ring3_B_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring2_C_Minus  = new TH1F("TOF_SimHits_RE4_Ring2_C_Minus",  "TOF_SimHits_RE4_Ring2_C_Minus", n_tof, n1_tof, n2_tof);
  TOF_SimHits_RE4_Ring3_C_Minus  = new TH1F("TOF_SimHits_RE4_Ring3_C_Minus",  "TOF_SimHits_RE4_Ring3_C_Minus", n_tof, n1_tof, n2_tof);

  VTX_PX = new TH1F("VTX_PX", "x position of Primary Vertex", n_v_x, n1_v_x, n2_v_x);
  VTX_PY = new TH1F("VTX_PY", "y position of Primary Vertex", n_v_x, n1_v_x, n2_v_x);
  VTX_PZ = new TH1F("VTX_PZ", "z position of Primary Vertex", n_v_z, n1_v_z, n2_v_z);

  VTX_AX = new TH1F("VTX_AX", "x position of All Vertices", p_v_x, p1_v_x, p2_v_x);
  VTX_AY = new TH1F("VTX_AY", "y position of All vertices", p_v_x, p1_v_x, p2_v_x);
  VTX_AZ = new TH1F("VTX_AZ", "z position of All vertices", p_v_z, p1_v_z, p2_v_z);

}


MyRE4nSimHitAnalyzer::~MyRE4nSimHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();
  /*
  TOF_SimHits_RE4_Plus->Write();
  TOF_SimHits_RE4_Ring2_Plus->Write();
  TOF_SimHits_RE4_Ring3_Plus->Write();
  TOF_SimHits_RE4_Minus->Write();
  TOF_SimHits_RE4_Ring2_Minus->Write();
  TOF_SimHits_RE4_Ring3_Minus->Write();

  TOF_SimHits_RE3_Plus->Write();
  TOF_SimHits_RE3_Ring2_Plus->Write();
  TOF_SimHits_RE3_Ring3_Plus->Write();
  TOF_SimHits_RE3_Minus->Write();
  TOF_SimHits_RE3_Ring2_Minus->Write();
  TOF_SimHits_RE3_Ring3_Minus->Write();

  TOF_SimHits_RE4_Ring2_A_Plus->Write();
  TOF_SimHits_RE4_Ring3_A_Plus->Write();
  TOF_SimHits_RE4_Ring2_B_Plus->Write();
  TOF_SimHits_RE4_Ring3_B_Plus->Write();
  TOF_SimHits_RE4_Ring2_C_Plus->Write();
  TOF_SimHits_RE4_Ring3_C_Plus->Write();
  TOF_SimHits_RE4_Ring2_A_Minus->Write();
  TOF_SimHits_RE4_Ring3_A_Minus->Write();
  TOF_SimHits_RE4_Ring2_B_Minus->Write();
  TOF_SimHits_RE4_Ring3_B_Minus->Write();
  TOF_SimHits_RE4_Ring2_C_Minus->Write();
  TOF_SimHits_RE4_Ring3_C_Minus->Write();
  */
  TOF_SimHits_RE4_Ring2_A_Plus->SetFillColor(kRed);
  TOF_SimHits_RE4_Ring3_A_Plus->SetFillColor(kRed);
  TOF_SimHits_RE4_Ring2_B_Plus->SetFillColor(kMagenta);
  TOF_SimHits_RE4_Ring3_B_Plus->SetFillColor(kMagenta);
  TOF_SimHits_RE4_Ring2_C_Plus->SetFillColor(kBlue);
  TOF_SimHits_RE4_Ring3_C_Plus->SetFillColor(kBlue);
  TOF_SimHits_RE4_Ring2_A_Minus->SetFillColor(kRed);
  TOF_SimHits_RE4_Ring3_A_Minus->SetFillColor(kRed);
  TOF_SimHits_RE4_Ring2_B_Minus->SetFillColor(kMagenta);
  TOF_SimHits_RE4_Ring3_B_Minus->SetFillColor(kMagenta);
  TOF_SimHits_RE4_Ring2_C_Minus->SetFillColor(kBlue);
  TOF_SimHits_RE4_Ring3_C_Minus->SetFillColor(kBlue);



  const int n_p = x_p.size();  double x_ap[n_p]; double y_ap[n_p];  // double z_ap[n_p];
  const int n_n = x_n.size();  double x_an[n_n]; double y_an[n_n];  // double z_an[n_n];
  const int n_r = r_r.size();  double r_ar[n_r]; double z_ar[n_r];
  for(int i=0; i< n_p; ++i) { x_ap[i] = x_p[i]; y_ap[i] = y_p[i]; /*z_ap[i] = z_p[i];*/ }
  for(int i=0; i< n_n; ++i) { x_an[i] = x_n[i]; y_an[i] = y_n[i]; /*z_an[i] = z_n[i];*/ }
  for(int i=0; i< n_r; ++i) { r_ar[i] = r_r[i]; z_ar[i] = z_r[i]; }

  RE4_Plus_XY  = new TGraph(n_p, x_ap, y_ap);
  RE4_Minus_XY = new TGraph(n_n, x_an, y_an);
  RE_YZ        = new TGraph(n_r, z_ar, r_ar);

  Canvas_RE4_Plus_XY  = new TCanvas("Canvas_RE4_Plus_XY",  "Canvas_RE4_Plus_XY", 800, 600);
  Canvas_RE4_Minus_XY = new TCanvas("Canvas_RE4_Minus_XY", "Canvas_RE4_Minus_XY", 800, 600);
  Canvas_RE_YZ        = new TCanvas("Canvas_RE_YZ",        "Canvas_RE_YZ", 800, 600);

  Canvas_RE4_TOF            = new TCanvas("Canvas_RE4_TOF",      "Canvas_RE4_TOF", 800, 600);
  Canvas_RE4_TOF_XCheck     = new TCanvas("Canvas_RE4_TOF_XCheck", "Canvas_RE4_TOF_XCheck", 800, 600);
  Canvas_RE4_TOF_Difference = new TCanvas("Canvas_RE4_TOF_Difference", "Canvas_RE4_TOF_Difference", 800, 600);
  Canvas_RE3_TOF            = new TCanvas("Canvas_RE3_TOF",      "Canvas_RE3_TOF", 800, 600);
  Canvas_VTX                = new TCanvas("Canvas_VTX",      "Canvas_VTX", 800, 600);



  THStack * RE4_Ring2_Plus = new THStack("RE4_Ring2_Plus",   "RE+4 Ring2");
  THStack * RE4_Ring3_Plus = new THStack("RE4_Ring3_Plus",   "RE+4 Ring3");
  THStack * RE4_Ring2_Minus = new THStack("RE4_Ring2_Minus", "RE-4 Ring2");
  THStack * RE4_Ring3_Minus = new THStack("RE4_Ring3_Minus", "RE-4 Ring3");
  TLegend *l1 = new TLegend(0.175,0.625,0.375,0.875,NULL,"brNDC");
  l1->SetLineColor(1); l1->SetLineStyle(1); l1->SetLineWidth(2); l1->SetFillColor(4000); l1->SetBorderSize(1);
  l1->AddEntry(TOF_SimHits_RE4_Ring2_A_Plus, "Ch","");
  l1->AddEntry(TOF_SimHits_RE4_Ring2_A_Plus, "A","f");
  l1->AddEntry(TOF_SimHits_RE4_Ring2_B_Plus, "B","f");
  l1->AddEntry(TOF_SimHits_RE4_Ring2_C_Plus, "C","f");

  // Calculate the difference
  /*
  TOF_Difference_SimHits_RE4_Plus->Add(TOF_SimHits_RE4_Plus,1);  TOF_Difference_SimHits_RE4_Plus->Add(TOF_XCheck_SimHits_RE4_Plus,-1);
  TOF_Difference_SimHits_RE4_Ring2_Plus->Add(TOF_SimHits_RE4_Ring2_Plus,1);  TOF_Difference_SimHits_RE4_Ring2_Plus->Add(TOF_XCheck_SimHits_RE4_Ring2_Plus,-1);
  TOF_Difference_SimHits_RE4_Ring3_Plus->Add(TOF_SimHits_RE4_Ring3_Plus,1);  TOF_Difference_SimHits_RE4_Ring3_Plus->Add(TOF_XCheck_SimHits_RE4_Ring3_Plus,-1);
  TOF_Difference_SimHits_RE4_Minus->Add(TOF_SimHits_RE4_Minus,1);  TOF_Difference_SimHits_RE4_Minus->Add(TOF_XCheck_SimHits_RE4_Minus,-1);
  TOF_Difference_SimHits_RE4_Ring2_Minus->Add(TOF_SimHits_RE4_Ring2_Minus,1);  TOF_Difference_SimHits_RE4_Ring2_Minus->Add(TOF_XCheck_SimHits_RE4_Ring2_Minus,-1);
  TOF_Difference_SimHits_RE4_Ring3_Minus->Add(TOF_SimHits_RE4_Ring3_Minus,1);  TOF_Difference_SimHits_RE4_Ring3_Minus->Add(TOF_XCheck_SimHits_RE4_Ring3_Minus,-1);
  */
  // !!! Don't take the difference of the distribution, but take the distribution of the differences !!!

  // XY and RZ Graphs
 Canvas_RE4_Plus_XY->cd();  RE4_Plus_XY->SetMarkerStyle(5);  RE4_Plus_XY->Draw("AP");  RE4_Plus_XY->GetXaxis()->SetTitle("X [cm]");  RE4_Plus_XY->GetYaxis()->SetTitle("Y [cm]");   RE4_Plus_XY->SetTitle("RE+4 SimHits");
 Canvas_RE4_Minus_XY->cd(); RE4_Minus_XY->SetMarkerStyle(5); RE4_Minus_XY->Draw("AP"); RE4_Minus_XY->GetXaxis()->SetTitle("X [cm]"); RE4_Minus_XY->GetYaxis()->SetTitle("Y [cm]");  RE4_Minus_XY->SetTitle("RE-4 SimHits");
 Canvas_RE_YZ->cd();        RE_YZ->SetMarkerStyle(5);        RE_YZ->Draw("AP");        RE_YZ->GetXaxis()->SetTitle("Z [cm]");        RE_YZ->GetYaxis()->SetTitle("R [cm]");         RE_YZ->SetTitle("RE SimHits");

 // VTX Graphs
 Canvas_VTX->cd();  Canvas_VTX->Divide(3,2);
 Canvas_VTX->cd(1); VTX_PX->Draw(); VTX_PX->GetXaxis()->SetTitle("x [cm]"); VTX_PX->GetYaxis()->SetTitle("entries [-]"); VTX_PX->SetTitle("Primary Vertex X Position");
 Canvas_VTX->cd(2); VTX_PY->Draw(); VTX_PY->GetXaxis()->SetTitle("y [cm]"); VTX_PY->GetYaxis()->SetTitle("entries [-]"); VTX_PY->SetTitle("Primary Vertex Y Position");
 Canvas_VTX->cd(3); VTX_PZ->Draw(); VTX_PZ->GetXaxis()->SetTitle("z [cm]"); VTX_PZ->GetYaxis()->SetTitle("entries [-]"); VTX_PZ->SetTitle("Primary Vertex Z Position");
 Canvas_VTX->cd(4); VTX_AX->Draw(); VTX_AX->GetXaxis()->SetTitle("x [cm]"); VTX_AX->GetYaxis()->SetTitle("entries [-]"); VTX_AX->SetTitle("SimVertices in RPC :: X");
 Canvas_VTX->cd(5); VTX_AY->Draw(); VTX_AY->GetXaxis()->SetTitle("y [cm]"); VTX_AY->GetYaxis()->SetTitle("entries [-]"); VTX_AY->SetTitle("SimVertices in RPC :: Y");
 Canvas_VTX->cd(6); VTX_AZ->Draw(); VTX_AZ->GetXaxis()->SetTitle("z [cm]"); VTX_AZ->GetYaxis()->SetTitle("entries [-]"); VTX_AZ->SetTitle("SimVertices in RPC :: Z");

 // RE4 TOF
 Canvas_RE4_TOF->cd();  Canvas_RE4_TOF->Divide(3,2);
 Canvas_RE4_TOF->cd(1); TOF_SimHits_RE4_Plus->Draw(); TOF_SimHits_RE4_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE4_Plus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE4_Plus->SetTitle("RE+4 SimHits");
 Canvas_RE4_TOF->cd(2); TOF_SimHits_RE4_Ring2_Plus->Draw(); TOF_SimHits_RE4_Ring2_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE4_Ring2_Plus->SetTitle("RE+4/2 SimHits");
 RE4_Ring2_Plus->Add(TOF_SimHits_RE4_Ring2_A_Plus); RE4_Ring2_Plus->Add(TOF_SimHits_RE4_Ring2_B_Plus); RE4_Ring2_Plus->Add(TOF_SimHits_RE4_Ring2_C_Plus); RE4_Ring2_Plus->Draw("Same"); l1->Draw(); gPad->RedrawAxis();
 Canvas_RE4_TOF->cd(3); TOF_SimHits_RE4_Ring3_Plus->Draw(); TOF_SimHits_RE4_Ring3_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE4_Ring3_Plus->SetTitle("RE+4/3 SimHits"); 
 RE4_Ring3_Plus->Add(TOF_SimHits_RE4_Ring3_A_Plus); RE4_Ring3_Plus->Add(TOF_SimHits_RE4_Ring3_B_Plus); RE4_Ring3_Plus->Add(TOF_SimHits_RE4_Ring3_C_Plus); RE4_Ring3_Plus->Draw("Same"); l1->Draw(); gPad->RedrawAxis();
 Canvas_RE4_TOF->cd(4); TOF_SimHits_RE4_Minus->Draw(); TOF_SimHits_RE4_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE4_Minus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE4_Minus->SetTitle("RE-4 SimHits"); 
 Canvas_RE4_TOF->cd(5); TOF_SimHits_RE4_Ring2_Minus->Draw(); TOF_SimHits_RE4_Ring2_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE4_Ring2_Minus->SetTitle("RE-4/2 SimHits"); 
 RE4_Ring2_Minus->Add(TOF_SimHits_RE4_Ring2_A_Minus); RE4_Ring2_Minus->Add(TOF_SimHits_RE4_Ring2_B_Minus); RE4_Ring2_Minus->Add(TOF_SimHits_RE4_Ring2_C_Minus); RE4_Ring2_Minus->Draw("Same"); l1->Draw(); gPad->RedrawAxis();
 Canvas_RE4_TOF->cd(6); TOF_SimHits_RE4_Ring3_Minus->Draw(); TOF_SimHits_RE4_Ring3_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE4_Ring3_Minus->SetTitle("RE-4/3 SimHits"); 
 RE4_Ring3_Minus->Add(TOF_SimHits_RE4_Ring3_A_Minus); RE4_Ring3_Minus->Add(TOF_SimHits_RE4_Ring3_B_Minus); RE4_Ring3_Minus->Add(TOF_SimHits_RE4_Ring3_C_Minus); RE4_Ring3_Minus->Draw("Same"); l1->Draw(); gPad->RedrawAxis();

 // RE4 XCheck
 Canvas_RE4_TOF_XCheck->cd();  Canvas_RE4_TOF_XCheck->Divide(3,2);
 Canvas_RE4_TOF_XCheck->cd(1); TOF_XCheck1_SimHits_RE4_Plus->Draw(); TOF_XCheck1_SimHits_RE4_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_XCheck1_SimHits_RE4_Plus->GetYaxis()->SetTitle("entries [-]"); TOF_XCheck1_SimHits_RE4_Plus->SetTitle("RE+4 XCheck");
 TOF_XCheck2_SimHits_RE4_Plus->SetLineColor(kRed); TOF_XCheck2_SimHits_RE4_Plus->Draw("Same");
 Canvas_RE4_TOF_XCheck->cd(2); TOF_XCheck1_SimHits_RE4_Ring2_Plus->Draw(); TOF_XCheck1_SimHits_RE4_Ring2_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_XCheck1_SimHits_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_XCheck1_SimHits_RE4_Ring2_Plus->SetTitle("RE+4/2 Xcheck"); TOF_XCheck2_SimHits_RE4_Ring2_Plus->SetLineColor(kRed); TOF_XCheck2_SimHits_RE4_Ring2_Plus->Draw("Same");
 Canvas_RE4_TOF_XCheck->cd(3); TOF_XCheck1_SimHits_RE4_Ring3_Plus->Draw(); TOF_XCheck1_SimHits_RE4_Ring3_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_XCheck1_SimHits_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_XCheck1_SimHits_RE4_Ring3_Plus->SetTitle("RE+4/3 XCheck"); TOF_XCheck2_SimHits_RE4_Ring3_Plus->SetLineColor(kRed); TOF_XCheck2_SimHits_RE4_Ring3_Plus->Draw("Same");
 Canvas_RE4_TOF_XCheck->cd(4); TOF_XCheck1_SimHits_RE4_Minus->Draw(); TOF_XCheck1_SimHits_RE4_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_XCheck1_SimHits_RE4_Minus->GetYaxis()->SetTitle("entries [-]"); TOF_XCheck1_SimHits_RE4_Minus->SetTitle("RE-4 XCheck"); 
 TOF_XCheck2_SimHits_RE4_Minus->SetLineColor(kRed); TOF_XCheck2_SimHits_RE4_Minus->Draw("Same");
 Canvas_RE4_TOF_XCheck->cd(5); TOF_XCheck1_SimHits_RE4_Ring2_Minus->Draw(); TOF_XCheck1_SimHits_RE4_Ring2_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_XCheck1_SimHits_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_XCheck1_SimHits_RE4_Ring2_Minus->SetTitle("RE-4/2 XCheck"); TOF_XCheck2_SimHits_RE4_Ring2_Minus->SetLineColor(kRed); TOF_XCheck2_SimHits_RE4_Ring2_Minus->Draw("Same");
 Canvas_RE4_TOF_XCheck->cd(6); TOF_XCheck1_SimHits_RE4_Ring3_Minus->Draw(); TOF_XCheck1_SimHits_RE4_Ring3_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_XCheck1_SimHits_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_XCheck1_SimHits_RE4_Ring3_Minus->SetTitle("RE-4/3 XCheck");  TOF_XCheck2_SimHits_RE4_Ring3_Minus->SetLineColor(kRed); TOF_XCheck2_SimHits_RE4_Ring3_Minus->Draw("Same");

 // RE4 Difference
 Canvas_RE4_TOF_Difference->cd();  Canvas_RE4_TOF_Difference->Divide(3,2);
 Canvas_RE4_TOF_Difference->cd(1); TOF_Difference1_SimHits_RE4_Plus->Draw(); TOF_Difference1_SimHits_RE4_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_Difference1_SimHits_RE4_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_Difference1_SimHits_RE4_Plus->SetTitle("RE+4 Difference"); TOF_Difference2_SimHits_RE4_Plus->SetLineColor(kRed); TOF_Difference2_SimHits_RE4_Plus->Draw("Same");
 Canvas_RE4_TOF_Difference->cd(2); TOF_Difference1_SimHits_RE4_Ring2_Plus->Draw(); TOF_Difference1_SimHits_RE4_Ring2_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_Difference1_SimHits_RE4_Ring2_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_Difference1_SimHits_RE4_Ring2_Plus->SetTitle("RE+4/2 Difference"); TOF_Difference2_SimHits_RE4_Ring2_Plus->SetLineColor(kRed); TOF_Difference2_SimHits_RE4_Ring2_Plus->Draw("Same");
 Canvas_RE4_TOF_Difference->cd(3); TOF_Difference1_SimHits_RE4_Ring3_Plus->Draw(); TOF_Difference1_SimHits_RE4_Ring3_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_Difference1_SimHits_RE4_Ring3_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_Difference1_SimHits_RE4_Ring3_Plus->SetTitle("RE+4/3 Difference"); TOF_Difference2_SimHits_RE4_Ring3_Plus->SetLineColor(kRed); TOF_Difference2_SimHits_RE4_Ring3_Plus->Draw("Same");
 Canvas_RE4_TOF_Difference->cd(4); TOF_Difference1_SimHits_RE4_Minus->Draw(); TOF_Difference1_SimHits_RE4_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_Difference1_SimHits_RE4_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_Difference1_SimHits_RE4_Minus->SetTitle("RE-4 Difference"); TOF_Difference2_SimHits_RE4_Minus->SetLineColor(kRed); TOF_Difference2_SimHits_RE4_Minus->Draw("Same");
 Canvas_RE4_TOF_Difference->cd(5); TOF_Difference1_SimHits_RE4_Ring2_Minus->Draw(); TOF_Difference1_SimHits_RE4_Ring2_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_Difference1_SimHits_RE4_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_Difference1_SimHits_RE4_Ring2_Minus->SetTitle("RE-4/2 Difference"); TOF_Difference2_SimHits_RE4_Ring2_Minus->SetLineColor(kRed); TOF_Difference2_SimHits_RE4_Ring2_Minus->Draw("Same");
 Canvas_RE4_TOF_Difference->cd(6); TOF_Difference1_SimHits_RE4_Ring3_Minus->Draw(); TOF_Difference1_SimHits_RE4_Ring3_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_Difference1_SimHits_RE4_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_Difference1_SimHits_RE4_Ring3_Minus->SetTitle("RE-4/3 Difference");  TOF_Difference2_SimHits_RE4_Ring3_Minus->SetLineColor(kRed); TOF_Difference2_SimHits_RE4_Ring3_Minus->Draw("Same");

 // RE3 TOF
 Canvas_RE3_TOF->cd();  Canvas_RE3_TOF->Divide(3,2);
 Canvas_RE3_TOF->cd(1); TOF_SimHits_RE3_Plus->Draw(); TOF_SimHits_RE3_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE3_Plus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE3_Plus->SetTitle("RE+3 SimHits");
 Canvas_RE3_TOF->cd(2); TOF_SimHits_RE3_Ring2_Plus->Draw(); TOF_SimHits_RE3_Ring2_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE3_Ring2_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_SimHits_RE3_Ring2_Plus->SetTitle("RE+3/2 SimHits"); 
 Canvas_RE3_TOF->cd(3); TOF_SimHits_RE3_Ring3_Plus->Draw(); TOF_SimHits_RE3_Ring3_Plus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE3_Ring3_Plus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_SimHits_RE3_Ring3_Plus->SetTitle("RE+3/3 SimHits");
 Canvas_RE3_TOF->cd(4); TOF_SimHits_RE3_Minus->Draw(); TOF_SimHits_RE3_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE3_Minus->GetYaxis()->SetTitle("entries [-]"); TOF_SimHits_RE3_Minus->SetTitle("RE-4 SimHits");
 Canvas_RE3_TOF->cd(5); TOF_SimHits_RE3_Ring2_Minus->Draw(); TOF_SimHits_RE3_Ring2_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE3_Ring2_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_SimHits_RE3_Ring2_Minus->SetTitle("RE-3/2 SimHits");
 Canvas_RE3_TOF->cd(6); TOF_SimHits_RE3_Ring3_Minus->Draw(); TOF_SimHits_RE3_Ring3_Minus->GetXaxis()->SetTitle("tof [ns]"); TOF_SimHits_RE3_Ring3_Minus->GetYaxis()->SetTitle("entries [-]"); 
 TOF_SimHits_RE3_Ring3_Minus->SetTitle("RE-3/3 SimHits");




  Canvas_RE4_Plus_XY->Write();
  Canvas_RE4_Minus_XY->Write();
  Canvas_RE_YZ->Write();
  Canvas_VTX->Write();
  Canvas_RE4_TOF->Write();
  Canvas_RE4_TOF_XCheck->Write();
  Canvas_RE4_TOF_Difference->Write();
  Canvas_RE3_TOF->Write();
  /*
  RE4_Plus_XY->Write();
  RE4_Minus_XY->Write();
  RE_YZ->Write();
  */
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyRE4nSimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // SimHits
  std::cout << " Getting the SimHits " <<std::endl;
  std::vector<edm::Handle<edm::PSimHitContainer> > theSimHitContainers;
  iEvent.getManyByType(theSimHitContainers);
  std::cout << " The number of SimHit Containers is  " << theSimHitContainers.size() <<std::endl;
  std::vector<PSimHit> theSimHits;
  for (int i = 0; i < int(theSimHitContainers.size()); ++i) {
    theSimHits.insert(theSimHits.end(),theSimHitContainers.at(i)->begin(),theSimHitContainers.at(i)->end());
  }
  // SimTracks
  std::vector<SimTrack> theSimTracks;
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel("g4SimHits",SimTk);
  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  std::cout << "This Event has " <<  theSimTracks.size() << " sim tracks " << std::endl;
  // SimVertices
  std::vector<SimVertex> theSimVertices; 
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel("g4SimHits",SimVtx);
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  std::cout << "This Event has " <<  theSimVertices.size() << " sim vertices " << std::endl;


  // ===========================================================
  //      Search for the right SimTrack and right SimVertex
  // ===========================================================
  int vertex_index = -1;
  for (std::vector<SimTrack>::const_iterator iTrack = theSimTracks.begin(); iTrack != theSimTracks.end(); ++iTrack) {
    SimTrack simtrack = (*iTrack);
    // if(fabs(simtrack.type()) != 13 || simtrack.momentum().vect().perp() < 95.0) continue;
    if(fabs(simtrack.type()) != 13) continue;
    else {
      double mu_pt = simtrack.momentum().pt();
      std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<simtrack.momentum().pt()<<" GeV/c from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      // double Px = simtrack.momentum().x(); double Py = simtrack.momentum().y(); double mu_pt = sqrt(pow(Px,2)+pow(Py,2));
      // std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<mu_pt<<" GeV/c from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      // std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      if(mu_pt > 95) {
	vertex_index = simtrack.vertIndex();
      }
    }
  }
  double vtx_x = -99.9, vtx_y = -99.9, vtx_z = -99.9; 
  int counter = 0;
  for (std::vector<SimVertex>::const_iterator iVertex = theSimVertices.begin(); iVertex != theSimVertices.end(); ++iVertex) {
    SimVertex simvertex = (*iVertex);
    if(vertex_index==counter) {
      vtx_x = simvertex.position().x();
      vtx_y = simvertex.position().y();
      vtx_z = simvertex.position().z();
      std::cout<<"Muon SimVertex Found: x = "<<vtx_x<<" y = "<<vtx_y<<" z = "<<vtx_z<<std::endl;
    }
    else {
      VTX_AX->Fill(simvertex.position().x()); VTX_AY->Fill(simvertex.position().y()); VTX_AZ->Fill(simvertex.position().z());
      std::cout<<"PileUp SimVertex Found: x = "<<simvertex.position().x()<<" y = "<<simvertex.position().y()<<" z = "<<simvertex.position().z()<<std::endl;
    }
    ++counter;
  }
  GlobalPoint VtxGlobalPoint(vtx_x,vtx_y,vtx_z);
  VTX_PX->Fill(VtxGlobalPoint.x()); VTX_PY->Fill(VtxGlobalPoint.y()); VTX_PZ->Fill(VtxGlobalPoint.z());



  // ===============================
  //      Loop over the SimHits
  // ===============================
  for (std::vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); ++iHit) {
    DetId theDetUnitId((*iHit).detUnitId());
    DetId simdetid= DetId((*iHit).detUnitId());

    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){ // Only RPCs
      RPCDetId rollId(theDetUnitId);
      RPCGeomServ rpcsrv(rollId);
      //std::cout << " Reading the Roll"<<std::endl;
      const RPCRoll* rollasociated = rpcGeom->roll(rollId);
      //std::cout << " Getting the Surface"<<std::endl;
      const BoundPlane & RPCSurface = rollasociated->surface();
      GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal((*iHit).localPosition());
      // GlobalPoint VtxGlobalPoint; if(theSimVertices.size() == 0) {;}  else {theSimVertices[0].position();}
      // VTX_X->Fill(VtxGlobalPoint.x()); VTX_Y->Fill(VtxGlobalPoint.y()); VTX_Z->Fill(VtxGlobalPoint.z());
      // Only RE4
      if ((rollId.station()==4) && rollId.region() != 0) {
	std::cout<<"RPC SimHit in "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
	std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta();
	std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
	// Timing as calculated from the global point
	// all distances are in cm and all time is in ns, therefore c is in cm/ns
	// GlobalPoint gp = RPCSurface.toGlobal(LocalPoint(0,0,0)); // global point of the middle of the rol
	float c = 29.9792458;
	// float time0 = gp.mag()/c; // time w.r.t. the middle of the roll
	float time1 = RPCGlobalPoint.mag()/c;
	float time2 = sqrt(pow(RPCGlobalPoint.x()-VtxGlobalPoint.x(),2)+pow(RPCGlobalPoint.y()-VtxGlobalPoint.y(),2)+pow(RPCGlobalPoint.z()-VtxGlobalPoint.z(),2))/c;
	// Fill Histograms
	if(rollId.region() == 1) {
	  // XY Graph RE+4
          x_p.push_back(RPCGlobalPoint.x()); y_p.push_back(RPCGlobalPoint.y()); z_p.push_back(RPCGlobalPoint.z());
	  TOF_SimHits_RE4_Plus->Fill((*iHit).timeOfFlight());
	  TOF_XCheck1_SimHits_RE4_Plus->Fill(time1); TOF_XCheck2_SimHits_RE4_Plus->Fill(time2);
	  TOF_Difference1_SimHits_RE4_Plus->Fill((*iHit).timeOfFlight()-time1); TOF_Difference2_SimHits_RE4_Plus->Fill((*iHit).timeOfFlight()-time2);
          if(rollId.ring() == 2) { 
	    TOF_SimHits_RE4_Ring2_Plus->Fill((*iHit).timeOfFlight()); 
	    if(rollId.roll() == 1) { TOF_SimHits_RE4_Ring2_A_Plus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 2) { TOF_SimHits_RE4_Ring2_B_Plus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 3) { TOF_SimHits_RE4_Ring2_C_Plus->Fill((*iHit).timeOfFlight());}
	    TOF_XCheck1_SimHits_RE4_Ring2_Plus->Fill(time1); TOF_XCheck2_SimHits_RE4_Ring2_Plus->Fill(time2);
	    TOF_Difference1_SimHits_RE4_Ring2_Plus->Fill((*iHit).timeOfFlight()-time1); TOF_Difference2_SimHits_RE4_Ring2_Plus->Fill((*iHit).timeOfFlight()-time2);
	  }
          if(rollId.ring() == 3) { 
	    TOF_SimHits_RE4_Ring3_Plus->Fill((*iHit).timeOfFlight()); 
	    if(rollId.roll() == 1) { TOF_SimHits_RE4_Ring3_A_Plus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 2) { TOF_SimHits_RE4_Ring3_B_Plus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 3) { TOF_SimHits_RE4_Ring3_C_Plus->Fill((*iHit).timeOfFlight());}
	    TOF_XCheck1_SimHits_RE4_Ring3_Plus->Fill(time1); TOF_XCheck2_SimHits_RE4_Ring3_Plus->Fill(time2);
	    TOF_Difference1_SimHits_RE4_Ring3_Plus->Fill((*iHit).timeOfFlight()-time1); TOF_Difference2_SimHits_RE4_Ring3_Plus->Fill((*iHit).timeOfFlight()-time2);
	  }
	}
	if(rollId.region() == -1) {
	  // XY Graph RE-4
	  x_n.push_back(RPCGlobalPoint.x()); y_n.push_back(RPCGlobalPoint.y()); z_n.push_back(RPCGlobalPoint.z());
	  TOF_SimHits_RE4_Minus->Fill((*iHit).timeOfFlight());
	  TOF_XCheck1_SimHits_RE4_Minus->Fill(time1); TOF_XCheck2_SimHits_RE4_Minus->Fill(time2);
	  TOF_Difference1_SimHits_RE4_Minus->Fill((*iHit).timeOfFlight()-time1); TOF_Difference2_SimHits_RE4_Minus->Fill((*iHit).timeOfFlight()-time2);
          if(rollId.ring() == 2) { 
	    TOF_SimHits_RE4_Ring2_Minus->Fill((*iHit).timeOfFlight());
	    if(rollId.roll() == 1) { TOF_SimHits_RE4_Ring2_A_Minus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 2) { TOF_SimHits_RE4_Ring2_B_Minus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 3) { TOF_SimHits_RE4_Ring2_C_Minus->Fill((*iHit).timeOfFlight());}
	    TOF_XCheck1_SimHits_RE4_Ring2_Minus->Fill(time1); TOF_XCheck2_SimHits_RE4_Ring2_Minus->Fill(time2);
	    TOF_Difference1_SimHits_RE4_Ring2_Minus->Fill((*iHit).timeOfFlight()-time1); TOF_Difference2_SimHits_RE4_Ring2_Minus->Fill((*iHit).timeOfFlight()-time2);
	  }
          if(rollId.ring() == 3) { 
	    TOF_SimHits_RE4_Ring3_Minus->Fill((*iHit).timeOfFlight());
	    if(rollId.roll() == 1) { TOF_SimHits_RE4_Ring3_A_Minus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 2) { TOF_SimHits_RE4_Ring3_B_Minus->Fill((*iHit).timeOfFlight());}
	    if(rollId.roll() == 3) { TOF_SimHits_RE4_Ring3_C_Minus->Fill((*iHit).timeOfFlight());}
	    TOF_XCheck1_SimHits_RE4_Ring3_Minus->Fill(time1); TOF_XCheck2_SimHits_RE4_Ring3_Minus->Fill(time2);
	    TOF_Difference1_SimHits_RE4_Ring3_Minus->Fill((*iHit).timeOfFlight()-time1);  TOF_Difference2_SimHits_RE4_Ring3_Minus->Fill((*iHit).timeOfFlight()-time2);
	  }
	}
      }
      // Only RE3                                                                                                                                                                                                                  
      if ((rollId.station()==3) && rollId.region() != 0) {
	std::cout<<"RPC SimHit in "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
	std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta();
	std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
        // Fill Histograms                                                                                                                                                                                                         
        if(rollId.region() == 1) {
          TOF_SimHits_RE3_Plus->Fill((*iHit).timeOfFlight());
          if(rollId.ring() == 2) { TOF_SimHits_RE3_Ring2_Plus->Fill((*iHit).timeOfFlight());}
	  if(rollId.ring() == 3) { TOF_SimHits_RE3_Ring3_Plus->Fill((*iHit).timeOfFlight());}
	}
	if(rollId.region() == -1) {
          TOF_SimHits_RE3_Minus->Fill((*iHit).timeOfFlight());
          if(rollId.ring() == 2) {TOF_SimHits_RE3_Ring2_Minus->Fill((*iHit).timeOfFlight());}
	  if(rollId.ring() == 3) {TOF_SimHits_RE3_Ring3_Minus->Fill((*iHit).timeOfFlight());}
	}
      }
      if(rollId.region() != 0) {
	// RZ graph endcap
	r_r.push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2))); z_r.push_back(RPCGlobalPoint.z());
      }
    }
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::CSC){ // Only CSCs
      CSCDetId rollId(theDetUnitId);
      // CSCGeomServ rpcsrv(rollId);
      GlobalPoint CSCGlobalPoint = cscGeom->idToDet(rollId)->toGlobal((*iHit).localPosition());
      // Only ME4
      if (rollId.station()==4) {
	// std::cout<<"CSC SimHit in "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
	std::cout<<"CSC SimHit in "<<std::setw(24)<<rollId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
	std::cout<<" | r = "<<std::setw(12)<<CSCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<CSCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<CSCGlobalPoint.eta();
	std::cout<<" | global position = "<<CSCGlobalPoint<<std::endl;
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRE4nSimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRE4nSimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyRE4nSimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
}


// ------------ method called when ending the processing of a run  ------------
void 
MyRE4nSimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyRE4nSimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyRE4nSimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRE4nSimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRE4nSimHitAnalyzer);
