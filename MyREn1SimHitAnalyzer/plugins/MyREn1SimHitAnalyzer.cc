// -*- C++ -*-
//
// Package:    MyREn1SimHitAnalyzer
// Class:      MyREn1SimHitAnalyzer
// 
/**\class MyREn1SimHitAnalyzer MyREn1SimHitAnalyzer.cc MyAnalyzers/MyREn1SimHitAnalyzer/src/MyREn1SimHitAnalyzer.cc

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
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TEllipse.h"


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

class MyREn1SimHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyREn1SimHitAnalyzer(const edm::ParameterSet&);
      ~MyREn1SimHitAnalyzer();
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
  bool debug;
  double simPtCut;
  TFile * outputfile;

  TH1F * ParticleIDs, * MuonHitEta, * MuonHitPhi, * MuonHitELoss, * EleHitELoss, * OtherHitELoss;

  TH1F * TOF_SimHits_RE1_Ring1_All_Plus,  * TOF_SimHits_RE2_Ring1_All_Plus,  * TOF_SimHits_RE3_Ring1_All_Plus,  * TOF_SimHits_RE4_Ring1_All_Plus;  
  TH1F * TOF_SimHits_RE1_Ring1_All_Minus, * TOF_SimHits_RE2_Ring1_All_Minus, * TOF_SimHits_RE3_Ring1_All_Minus, * TOF_SimHits_RE4_Ring1_All_Minus;
  std::vector<TH1F*> TOF_SimHits_REn_Ring1_All_Plus, TOF_SimHits_REn_Ring1_All_Minus;  

  TF1 * Fit_SimHits_RE1_Ring1_All_Plus,  * Fit_SimHits_RE2_Ring1_All_Plus,  * Fit_SimHits_RE3_Ring1_All_Plus,  * Fit_SimHits_RE4_Ring1_All_Plus;  
  TF1 * Fit_SimHits_RE1_Ring1_All_Minus, * Fit_SimHits_RE2_Ring1_All_Minus, * Fit_SimHits_RE3_Ring1_All_Minus, * Fit_SimHits_RE4_Ring1_All_Minus;
  std::vector<TF1*> Fit_SimHits_REn_Ring1_All_Plus, Fit_SimHits_REn_Ring1_All_Minus;  


  TH1F * TOF_SimHits_RE1_Ring1_A_Plus,  * TOF_SimHits_RE1_Ring1_B_Plus,  * TOF_SimHits_RE1_Ring1_C_Plus,  * TOF_SimHits_RE1_Ring1_D_Plus,  * TOF_SimHits_RE1_Ring1_E_Plus;
  TH1F * TOF_SimHits_RE2_Ring1_A_Plus,  * TOF_SimHits_RE2_Ring1_B_Plus,  * TOF_SimHits_RE2_Ring1_C_Plus,  * TOF_SimHits_RE2_Ring1_D_Plus,  * TOF_SimHits_RE2_Ring1_E_Plus;
  TH1F * TOF_SimHits_RE3_Ring1_A_Plus,  * TOF_SimHits_RE3_Ring1_B_Plus,  * TOF_SimHits_RE3_Ring1_C_Plus,  * TOF_SimHits_RE3_Ring1_D_Plus,  * TOF_SimHits_RE3_Ring1_E_Plus; 
  TH1F * TOF_SimHits_RE4_Ring1_A_Plus,  * TOF_SimHits_RE4_Ring1_B_Plus,  * TOF_SimHits_RE4_Ring1_C_Plus,  * TOF_SimHits_RE4_Ring1_D_Plus,  * TOF_SimHits_RE4_Ring1_E_Plus;
  TH1F * TOF_SimHits_RE1_Ring1_A_Minus, * TOF_SimHits_RE1_Ring1_B_Minus, * TOF_SimHits_RE1_Ring1_C_Minus, * TOF_SimHits_RE1_Ring1_D_Minus, * TOF_SimHits_RE1_Ring1_E_Minus;
  TH1F * TOF_SimHits_RE2_Ring1_A_Minus, * TOF_SimHits_RE2_Ring1_B_Minus, * TOF_SimHits_RE2_Ring1_C_Minus, * TOF_SimHits_RE2_Ring1_D_Minus, * TOF_SimHits_RE2_Ring1_E_Minus;
  TH1F * TOF_SimHits_RE3_Ring1_A_Minus, * TOF_SimHits_RE3_Ring1_B_Minus, * TOF_SimHits_RE3_Ring1_C_Minus, * TOF_SimHits_RE3_Ring1_D_Minus, * TOF_SimHits_RE3_Ring1_E_Minus; 
  TH1F * TOF_SimHits_RE4_Ring1_A_Minus, * TOF_SimHits_RE4_Ring1_B_Minus, * TOF_SimHits_RE4_Ring1_C_Minus, * TOF_SimHits_RE4_Ring1_D_Minus, * TOF_SimHits_RE4_Ring1_E_Minus;
  std::vector<TH1F*> TOF_SimHits_RE1_Ring1_Plus,  TOF_SimHits_RE2_Ring1_Plus,  TOF_SimHits_RE3_Ring1_Plus,  TOF_SimHits_RE4_Ring1_Plus;
  std::vector<TH1F*> TOF_SimHits_RE1_Ring1_Minus, TOF_SimHits_RE2_Ring1_Minus, TOF_SimHits_RE3_Ring1_Minus, TOF_SimHits_RE4_Ring1_Minus;
  std::vector< std::vector<TH1F*> > TOF_SimHits_RE_Ring1_Plus, TOF_SimHits_RE_Ring1_Minus;

  TF1 * Fit_SimHits_RE1_Ring1_A_Plus,  * Fit_SimHits_RE1_Ring1_B_Plus,  * Fit_SimHits_RE1_Ring1_C_Plus,  * Fit_SimHits_RE1_Ring1_D_Plus,  * Fit_SimHits_RE1_Ring1_E_Plus;
  TF1 * Fit_SimHits_RE2_Ring1_A_Plus,  * Fit_SimHits_RE2_Ring1_B_Plus,  * Fit_SimHits_RE2_Ring1_C_Plus,  * Fit_SimHits_RE2_Ring1_D_Plus,  * Fit_SimHits_RE2_Ring1_E_Plus;
  TF1 * Fit_SimHits_RE3_Ring1_A_Plus,  * Fit_SimHits_RE3_Ring1_B_Plus,  * Fit_SimHits_RE3_Ring1_C_Plus,  * Fit_SimHits_RE3_Ring1_D_Plus,  * Fit_SimHits_RE3_Ring1_E_Plus; 
  TF1 * Fit_SimHits_RE4_Ring1_A_Plus,  * Fit_SimHits_RE4_Ring1_B_Plus,  * Fit_SimHits_RE4_Ring1_C_Plus,  * Fit_SimHits_RE4_Ring1_D_Plus,  * Fit_SimHits_RE4_Ring1_E_Plus;
  TF1 * Fit_SimHits_RE1_Ring1_A_Minus, * Fit_SimHits_RE1_Ring1_B_Minus, * Fit_SimHits_RE1_Ring1_C_Minus, * Fit_SimHits_RE1_Ring1_D_Minus, * Fit_SimHits_RE1_Ring1_E_Minus;
  TF1 * Fit_SimHits_RE2_Ring1_A_Minus, * Fit_SimHits_RE2_Ring1_B_Minus, * Fit_SimHits_RE2_Ring1_C_Minus, * Fit_SimHits_RE2_Ring1_D_Minus, * Fit_SimHits_RE2_Ring1_E_Minus;
  TF1 * Fit_SimHits_RE3_Ring1_A_Minus, * Fit_SimHits_RE3_Ring1_B_Minus, * Fit_SimHits_RE3_Ring1_C_Minus, * Fit_SimHits_RE3_Ring1_D_Minus, * Fit_SimHits_RE3_Ring1_E_Minus; 
  TF1 * Fit_SimHits_RE4_Ring1_A_Minus, * Fit_SimHits_RE4_Ring1_B_Minus, * Fit_SimHits_RE4_Ring1_C_Minus, * Fit_SimHits_RE4_Ring1_D_Minus, * Fit_SimHits_RE4_Ring1_E_Minus;
  std::vector<TF1*> Fit_SimHits_RE1_Ring1_Plus,  Fit_SimHits_RE2_Ring1_Plus,  Fit_SimHits_RE3_Ring1_Plus,  Fit_SimHits_RE4_Ring1_Plus;
  std::vector<TF1*> Fit_SimHits_RE1_Ring1_Minus, Fit_SimHits_RE2_Ring1_Minus, Fit_SimHits_RE3_Ring1_Minus, Fit_SimHits_RE4_Ring1_Minus;
  std::vector< std::vector<TF1*> > Fit_SimHits_RE_Ring1_Plus, Fit_SimHits_RE_Ring1_Minus;

  TPaveStats * Box_SimHits_RE1_Ring1_A_Plus,  * Box_SimHits_RE1_Ring1_B_Plus,  * Box_SimHits_RE1_Ring1_C_Plus,  * Box_SimHits_RE1_Ring1_D_Plus,  * Box_SimHits_RE1_Ring1_E_Plus;
  TPaveStats * Box_SimHits_RE2_Ring1_A_Plus,  * Box_SimHits_RE2_Ring1_B_Plus,  * Box_SimHits_RE2_Ring1_C_Plus,  * Box_SimHits_RE2_Ring1_D_Plus,  * Box_SimHits_RE2_Ring1_E_Plus;
  TPaveStats * Box_SimHits_RE3_Ring1_A_Plus,  * Box_SimHits_RE3_Ring1_B_Plus,  * Box_SimHits_RE3_Ring1_C_Plus,  * Box_SimHits_RE3_Ring1_D_Plus,  * Box_SimHits_RE3_Ring1_E_Plus; 
  TPaveStats * Box_SimHits_RE4_Ring1_A_Plus,  * Box_SimHits_RE4_Ring1_B_Plus,  * Box_SimHits_RE4_Ring1_C_Plus,  * Box_SimHits_RE4_Ring1_D_Plus,  * Box_SimHits_RE4_Ring1_E_Plus;
  TPaveStats * Box_SimHits_RE1_Ring1_A_Minus, * Box_SimHits_RE1_Ring1_B_Minus, * Box_SimHits_RE1_Ring1_C_Minus, * Box_SimHits_RE1_Ring1_D_Minus, * Box_SimHits_RE1_Ring1_E_Minus;
  TPaveStats * Box_SimHits_RE2_Ring1_A_Minus, * Box_SimHits_RE2_Ring1_B_Minus, * Box_SimHits_RE2_Ring1_C_Minus, * Box_SimHits_RE2_Ring1_D_Minus, * Box_SimHits_RE2_Ring1_E_Minus;
  TPaveStats * Box_SimHits_RE3_Ring1_A_Minus, * Box_SimHits_RE3_Ring1_B_Minus, * Box_SimHits_RE3_Ring1_C_Minus, * Box_SimHits_RE3_Ring1_D_Minus, * Box_SimHits_RE3_Ring1_E_Minus; 
  TPaveStats * Box_SimHits_RE4_Ring1_A_Minus, * Box_SimHits_RE4_Ring1_B_Minus, * Box_SimHits_RE4_Ring1_C_Minus, * Box_SimHits_RE4_Ring1_D_Minus, * Box_SimHits_RE4_Ring1_E_Minus;
  std::vector<TPaveStats*> Box_SimHits_RE1_Ring1_Plus,  Box_SimHits_RE2_Ring1_Plus,  Box_SimHits_RE3_Ring1_Plus,  Box_SimHits_RE4_Ring1_Plus;
  std::vector<TPaveStats*> Box_SimHits_RE1_Ring1_Minus, Box_SimHits_RE2_Ring1_Minus, Box_SimHits_RE3_Ring1_Minus, Box_SimHits_RE4_Ring1_Minus;
  std::vector< std::vector<TPaveStats*> > Box_SimHits_RE_Ring1_Plus, Box_SimHits_RE_Ring1_Minus;


  TH1F * VTX_PX, * VTX_PY, *VTX_PZ, * VTX_AX, * VTX_AY, *VTX_AZ; 


  std::vector< std::vector< std::vector< double > > > x_p, y_p, z_p, r_p, x_n, y_n, z_n, r_n;
  std::vector< std::vector< double > > r_r, z_r;

  std::vector< std::vector<TGraph*> > RE_Plus_XY, RE_Minus_XY;
  std::vector< std::vector<TEllipse*> > DetBound;
  std::vector<TGraph*> RE_YZ;

  TCanvas * Canvas_RE1_Plus_XY, * Canvas_RE1_Minus_XY, * Canvas_RE1_Plus_TOF,  * Canvas_RE1_Minus_TOF; 
  TCanvas * Canvas_RE2_Plus_XY, * Canvas_RE2_Minus_XY, * Canvas_RE2_Plus_TOF,  * Canvas_RE2_Minus_TOF; 
  TCanvas * Canvas_RE3_Plus_XY, * Canvas_RE3_Minus_XY, * Canvas_RE3_Plus_TOF,  * Canvas_RE3_Minus_TOF; 
  TCanvas * Canvas_RE4_Plus_XY, * Canvas_RE4_Minus_XY, * Canvas_RE4_Plus_TOF,  * Canvas_RE4_Minus_TOF;
  std::vector<TCanvas*>  Canvas_RE_Plus_XY, Canvas_RE_Minus_XY, Canvas_RE_Plus_TOF,  Canvas_RE_Minus_TOF;  
  TCanvas * Canvas_RE_YZ, * Canvas_VTX; 

  TH2F * Histo_XY_RE11_P, * Histo_XY_RE21_P, * Histo_XY_RE31_P, * Histo_XY_RE41_P;
  TH2F * Histo_XY_RE11_N, * Histo_XY_RE21_N, * Histo_XY_RE31_N, * Histo_XY_RE41_N;
  TH2F * Histo_RZ_RE;
  std::vector< TH2F* > Histo_XY_RE_P, Histo_XY_RE_N;

  edm::EDGetTokenT<edm::PSimHitContainer> RPCSimHit_Token;



};

//
// constants, enums and typedefs
//
int n_tof  = 300; double n1_tof  = 15,    n2_tof = 45;                  // 100 ps resolution

int n_tof_re11 = 160; double n1_tof_re11 = 18.00, n2_tof_re11 = 22.00;  //  25 ps resolution
int n_tof_re21 = 160; double n1_tof_re21 = 25.75, n2_tof_re21 = 29.75;  //  25 ps resolution
int n_tof_re31 = 160; double n1_tof_re31 = 31.75, n2_tof_re31 = 35.75;  //  25 ps resolution
int n_tof_re41 = 160; double n1_tof_re41 = 34.75, n2_tof_re41 = 38.75;  //  25 ps resolution

int n_dif  = 75;  double n1_dif  = -7.5,  n2_dif = 7.5; 
// int n_v_x  = 40;   double n1_v_x  = -2.0,  n2_v_x = 2.0;
int n_v_x  = 500;   double n1_v_x  = -0.5,  n2_v_x = 0.5;
int n_v_z  = 500;   double n1_v_z  = -25.0, n2_v_z = 25.0;
// pu
int p_v_x  = 500;   double p1_v_x  = -500.0,  p2_v_x = 500.0;
int p_v_z  = 1250;  double p1_v_z  = -1250.0, p2_v_z = 1250.0;

int n_xy_x = 700; double n_xy_x1 = -350;  double n_xy_x2 = +350; 
int n_xy_y = 700; double n_xy_y1 = -350;  double n_xy_y2 = +350; 
int n_zr_z = 300; double n_zr_z1 = -1500; double n_zr_z2 = +1500;
int n_zr_r = 700; double n_zr_r1 = -350;  double n_zr_r2 = +350;  

const int n_stations = 4;
const int n_rolls = 5;
// colours
// From ROOT RTypes.h: 
// enum EColor { kWhite =0,   kBlack =1,   kGray=920,
//               kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
//               kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
int colours[n_rolls] = {800, 900, 880, 860, 840};
std::string pos_stations[n_stations] = {"RE+1/1 simhits", "RE+2/1 simhits", "RE+3/1 simhits", "RE+4/1 simhits"};
std::string neg_stations[n_stations] = {"RE-1/1 simhits", "RE-2/1 simhits", "RE-3/1 simhits", "RE-4/1 simhits"};


//
// static data member definitions
//

// template <typename T, size_t N> const T* mybegin(const T (&a)[N]) { return a; }    
// template <typename T, size_t N> const T* myend  (const T (&a)[N]) { return a+N; }


//
// constructors and destructor
//
MyREn1SimHitAnalyzer::MyREn1SimHitAnalyzer(const edm::ParameterSet& iConfig)

{
 
  //now do what ever initialization is needed
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
  debug         = iConfig.getUntrackedParameter<bool>("Debug");
  simPtCut      = iConfig.getUntrackedParameter<double>("SimPtCut");
  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

  RPCSimHit_Token = consumes<edm::PSimHitContainer>(edm::InputTag(std::string("g4SimHits"), std::string("MuonRPCHits")));

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::MyREn1SimHitAnalyzer CONSTRUCTOR]"<<std::endl;

  ParticleIDs = new TH1F("ParticleIDs","ParticleIDs", 1001, -500.5, 500.5);
  MuonHitEta = new TH1F("MuonHitEta", "MuonHitEta", 240, -2.4, 2.4);
  MuonHitPhi = new TH1F("MuonHitPhi", "MuonHitPhi", 144, -3.14, 3.14);
  EleHitELoss   = new TH1F("EleHitELoss",    "Electron Hit Energy Loss",    1000, 1e-7, 1e-4);
  MuonHitELoss  = new TH1F("MuonHitELoss",   "Muon Hit Energy Loss",        1000, 1e-7, 1e-4);
  OtherHitELoss = new TH1F("OtherHitELoss",  "Other Part. Hit Energy Loss", 1000, 1e-7, 1e-4);

  // TOF All Hits inside one station :: Hists ============================================================================================
  TOF_SimHits_RE1_Ring1_All_Plus  = new TH1F("TOF_SimHits_RE1_Ring1_All_Plus",  "TOF_SimHits_RE1_Ring1_All_Plus",  n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE2_Ring1_All_Plus  = new TH1F("TOF_SimHits_RE2_Ring1_All_Plus",  "TOF_SimHits_RE2_Ring2_All_Plus",  n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE3_Ring1_All_Plus  = new TH1F("TOF_SimHits_RE3_Ring1_All_Plus",  "TOF_SimHits_RE3_Ring3_All_Plus",  n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE4_Ring1_All_Plus  = new TH1F("TOF_SimHits_RE4_Ring1_All_Plus",  "TOF_SimHits_RE4_Ring4_All_Plus",  n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE1_Ring1_All_Minus = new TH1F("TOF_SimHits_RE1_Ring1_All_Minus", "TOF_SimHits_RE1_Ring1_All_Minus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE2_Ring1_All_Minus = new TH1F("TOF_SimHits_RE2_Ring2_All_Minus", "TOF_SimHits_RE2_Ring2_All_Minus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE3_Ring1_All_Minus = new TH1F("TOF_SimHits_RE3_Ring3_All_Minus", "TOF_SimHits_RE3_Ring3_All_Minus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE4_Ring1_All_Minus = new TH1F("TOF_SimHits_RE4_Ring4_All_Minus", "TOF_SimHits_RE4_Ring4_All_Minus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_REn_Ring1_All_Plus.push_back(TOF_SimHits_RE1_Ring1_All_Plus);   
  TOF_SimHits_REn_Ring1_All_Plus.push_back(TOF_SimHits_RE2_Ring1_All_Plus);   
  TOF_SimHits_REn_Ring1_All_Plus.push_back(TOF_SimHits_RE3_Ring1_All_Plus); 
  TOF_SimHits_REn_Ring1_All_Plus.push_back(TOF_SimHits_RE4_Ring1_All_Plus); 
  TOF_SimHits_REn_Ring1_All_Minus.push_back(TOF_SimHits_RE1_Ring1_All_Minus);
  TOF_SimHits_REn_Ring1_All_Minus.push_back(TOF_SimHits_RE2_Ring1_All_Minus);   
  TOF_SimHits_REn_Ring1_All_Minus.push_back(TOF_SimHits_RE3_Ring1_All_Minus); 
  TOF_SimHits_REn_Ring1_All_Minus.push_back(TOF_SimHits_RE4_Ring1_All_Minus);
  // =====================================================================================================================================  

  // TOF All Hits inside one station :: Fits =============================================================================================
  Fit_SimHits_RE1_Ring1_All_Plus  = new TF1("Fit_SimHits_RE1_Ring1_All_Plus","gaus", n1_tof_re11,n2_tof_re11);
  Fit_SimHits_RE2_Ring1_All_Plus  = new TF1("Fit_SimHits_RE2_Ring1_All_Plus","gaus", n1_tof_re21,n2_tof_re21);
  Fit_SimHits_RE3_Ring1_All_Plus  = new TF1("Fit_SimHits_RE3_Ring1_All_Plus","gaus", n1_tof_re31,n2_tof_re31);
  Fit_SimHits_RE4_Ring1_All_Plus  = new TF1("Fit_SimHits_RE4_Ring1_All_Plus","gaus", n1_tof_re41,n2_tof_re41);
  Fit_SimHits_RE1_Ring1_All_Minus = new TF1("Fit_SimHits_RE1_Ring1_All_Minus","gaus",n1_tof_re11,n2_tof_re11);
  Fit_SimHits_RE2_Ring1_All_Minus = new TF1("Fit_SimHits_RE2_Ring1_All_Minus","gaus",n1_tof_re21,n2_tof_re21);
  Fit_SimHits_RE3_Ring1_All_Minus = new TF1("Fit_SimHits_RE3_Ring1_All_Minus","gaus",n1_tof_re31,n2_tof_re31);
  Fit_SimHits_RE4_Ring1_All_Minus = new TF1("Fit_SimHits_RE4_Ring1_All_Minus","gaus",n1_tof_re41,n2_tof_re41);
  Fit_SimHits_REn_Ring1_All_Plus.push_back(Fit_SimHits_RE1_Ring1_All_Plus);   
  Fit_SimHits_REn_Ring1_All_Plus.push_back(Fit_SimHits_RE2_Ring1_All_Plus);   
  Fit_SimHits_REn_Ring1_All_Plus.push_back(Fit_SimHits_RE3_Ring1_All_Plus); 
  Fit_SimHits_REn_Ring1_All_Plus.push_back(Fit_SimHits_RE4_Ring1_All_Plus); 
  Fit_SimHits_REn_Ring1_All_Minus.push_back(Fit_SimHits_RE1_Ring1_All_Minus);
  Fit_SimHits_REn_Ring1_All_Minus.push_back(Fit_SimHits_RE2_Ring1_All_Minus);   
  Fit_SimHits_REn_Ring1_All_Minus.push_back(Fit_SimHits_RE3_Ring1_All_Minus); 
  Fit_SimHits_REn_Ring1_All_Minus.push_back(Fit_SimHits_RE4_Ring1_All_Minus);  
  // =====================================================================================================================================

  // TOF All Hits inside the different rolls  :: Hists ===================================================================================
  TOF_SimHits_RE1_Ring1_A_Plus  = new TH1F("TOF_SimHits_RE1_Ring1_A_Plus",  "TOF_SimHits_RE1_Ring1_A_Plus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_B_Plus  = new TH1F("TOF_SimHits_RE1_Ring1_B_Plus",  "TOF_SimHits_RE1_Ring1_B_Plus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_C_Plus  = new TH1F("TOF_SimHits_RE1_Ring1_C_Plus",  "TOF_SimHits_RE1_Ring1_C_Plus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_D_Plus  = new TH1F("TOF_SimHits_RE1_Ring1_D_Plus",  "TOF_SimHits_RE1_Ring1_D_Plus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_E_Plus  = new TH1F("TOF_SimHits_RE1_Ring1_E_Plus",  "TOF_SimHits_RE1_Ring1_E_Plus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE2_Ring1_A_Plus  = new TH1F("TOF_SimHits_RE2_Ring1_A_Plus",  "TOF_SimHits_RE2_Ring1_A_Plus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_B_Plus  = new TH1F("TOF_SimHits_RE2_Ring1_B_Plus",  "TOF_SimHits_RE2_Ring1_B_Plus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_C_Plus  = new TH1F("TOF_SimHits_RE2_Ring1_C_Plus",  "TOF_SimHits_RE2_Ring1_C_Plus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_D_Plus  = new TH1F("TOF_SimHits_RE2_Ring1_D_Plus",  "TOF_SimHits_RE2_Ring1_D_Plus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_E_Plus  = new TH1F("TOF_SimHits_RE2_Ring1_E_Plus",  "TOF_SimHits_RE2_Ring1_E_Plus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE3_Ring1_A_Plus  = new TH1F("TOF_SimHits_RE3_Ring1_A_Plus",  "TOF_SimHits_RE3_Ring1_A_Plus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_B_Plus  = new TH1F("TOF_SimHits_RE3_Ring1_B_Plus",  "TOF_SimHits_RE3_Ring1_B_Plus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_C_Plus  = new TH1F("TOF_SimHits_RE3_Ring1_C_Plus",  "TOF_SimHits_RE3_Ring1_C_Plus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_D_Plus  = new TH1F("TOF_SimHits_RE3_Ring1_D_Plus",  "TOF_SimHits_RE3_Ring1_D_Plus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_E_Plus  = new TH1F("TOF_SimHits_RE3_Ring1_E_Plus",  "TOF_SimHits_RE3_Ring1_E_Plus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE4_Ring1_A_Plus  = new TH1F("TOF_SimHits_RE4_Ring1_A_Plus",  "TOF_SimHits_RE4_Ring1_A_Plus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_B_Plus  = new TH1F("TOF_SimHits_RE4_Ring1_B_Plus",  "TOF_SimHits_RE4_Ring1_B_Plus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_C_Plus  = new TH1F("TOF_SimHits_RE4_Ring1_C_Plus",  "TOF_SimHits_RE4_Ring1_C_Plus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_D_Plus  = new TH1F("TOF_SimHits_RE4_Ring1_D_Plus",  "TOF_SimHits_RE4_Ring1_D_Plus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_E_Plus  = new TH1F("TOF_SimHits_RE4_Ring1_E_Plus",  "TOF_SimHits_RE4_Ring1_E_Plus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE1_Ring1_Plus.push_back(TOF_SimHits_RE1_Ring1_A_Plus);   TOF_SimHits_RE1_Ring1_Plus.push_back(TOF_SimHits_RE1_Ring1_B_Plus); 
  TOF_SimHits_RE1_Ring1_Plus.push_back(TOF_SimHits_RE1_Ring1_C_Plus);   TOF_SimHits_RE1_Ring1_Plus.push_back(TOF_SimHits_RE1_Ring1_D_Plus);   TOF_SimHits_RE1_Ring1_Plus.push_back(TOF_SimHits_RE1_Ring1_E_Plus);
  TOF_SimHits_RE2_Ring1_Plus.push_back(TOF_SimHits_RE2_Ring1_A_Plus);   TOF_SimHits_RE2_Ring1_Plus.push_back(TOF_SimHits_RE2_Ring1_B_Plus); 
  TOF_SimHits_RE2_Ring1_Plus.push_back(TOF_SimHits_RE2_Ring1_C_Plus);   TOF_SimHits_RE2_Ring1_Plus.push_back(TOF_SimHits_RE2_Ring1_D_Plus);   TOF_SimHits_RE2_Ring1_Plus.push_back(TOF_SimHits_RE2_Ring1_E_Plus);
  TOF_SimHits_RE3_Ring1_Plus.push_back(TOF_SimHits_RE3_Ring1_A_Plus);   TOF_SimHits_RE3_Ring1_Plus.push_back(TOF_SimHits_RE3_Ring1_B_Plus); 
  TOF_SimHits_RE3_Ring1_Plus.push_back(TOF_SimHits_RE3_Ring1_C_Plus);   TOF_SimHits_RE3_Ring1_Plus.push_back(TOF_SimHits_RE3_Ring1_D_Plus);   TOF_SimHits_RE3_Ring1_Plus.push_back(TOF_SimHits_RE3_Ring1_E_Plus);
  TOF_SimHits_RE4_Ring1_Plus.push_back(TOF_SimHits_RE4_Ring1_A_Plus);   TOF_SimHits_RE4_Ring1_Plus.push_back(TOF_SimHits_RE4_Ring1_B_Plus); 
  TOF_SimHits_RE4_Ring1_Plus.push_back(TOF_SimHits_RE4_Ring1_C_Plus);   TOF_SimHits_RE4_Ring1_Plus.push_back(TOF_SimHits_RE4_Ring1_D_Plus);   TOF_SimHits_RE4_Ring1_Plus.push_back(TOF_SimHits_RE4_Ring1_E_Plus);
  TOF_SimHits_RE_Ring1_Plus.push_back(TOF_SimHits_RE1_Ring1_Plus);  TOF_SimHits_RE_Ring1_Plus.push_back(TOF_SimHits_RE2_Ring1_Plus);
  TOF_SimHits_RE_Ring1_Plus.push_back(TOF_SimHits_RE3_Ring1_Plus);  TOF_SimHits_RE_Ring1_Plus.push_back(TOF_SimHits_RE4_Ring1_Plus);
  TOF_SimHits_RE1_Ring1_A_Minus  = new TH1F("TOF_SimHits_RE1_Ring1_A_Minus",  "TOF_SimHits_RE1_Ring1_A_Minus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_B_Minus  = new TH1F("TOF_SimHits_RE1_Ring1_B_Minus",  "TOF_SimHits_RE1_Ring1_B_Minus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_C_Minus  = new TH1F("TOF_SimHits_RE1_Ring1_C_Minus",  "TOF_SimHits_RE1_Ring1_C_Minus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_D_Minus  = new TH1F("TOF_SimHits_RE1_Ring1_D_Minus",  "TOF_SimHits_RE1_Ring1_D_Minus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE1_Ring1_E_Minus  = new TH1F("TOF_SimHits_RE1_Ring1_E_Minus",  "TOF_SimHits_RE1_Ring1_E_Minus", n_tof_re11, n1_tof_re11, n2_tof_re11);
  TOF_SimHits_RE2_Ring1_A_Minus  = new TH1F("TOF_SimHits_RE2_Ring1_A_Minus",  "TOF_SimHits_RE2_Ring1_A_Minus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_B_Minus  = new TH1F("TOF_SimHits_RE2_Ring1_B_Minus",  "TOF_SimHits_RE2_Ring1_B_Minus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_C_Minus  = new TH1F("TOF_SimHits_RE2_Ring1_C_Minus",  "TOF_SimHits_RE2_Ring1_C_Minus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_D_Minus  = new TH1F("TOF_SimHits_RE2_Ring1_D_Minus",  "TOF_SimHits_RE2_Ring1_D_Minus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE2_Ring1_E_Minus  = new TH1F("TOF_SimHits_RE2_Ring1_E_Minus",  "TOF_SimHits_RE2_Ring1_E_Minus", n_tof_re21, n1_tof_re21, n2_tof_re21);
  TOF_SimHits_RE3_Ring1_A_Minus  = new TH1F("TOF_SimHits_RE3_Ring1_A_Minus",  "TOF_SimHits_RE3_Ring1_A_Minus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_B_Minus  = new TH1F("TOF_SimHits_RE3_Ring1_B_Minus",  "TOF_SimHits_RE3_Ring1_B_Minus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_C_Minus  = new TH1F("TOF_SimHits_RE3_Ring1_C_Minus",  "TOF_SimHits_RE3_Ring1_C_Minus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_D_Minus  = new TH1F("TOF_SimHits_RE3_Ring1_D_Minus",  "TOF_SimHits_RE3_Ring1_D_Minus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE3_Ring1_E_Minus  = new TH1F("TOF_SimHits_RE3_Ring1_E_Minus",  "TOF_SimHits_RE3_Ring1_E_Minus", n_tof_re31, n1_tof_re31, n2_tof_re31);
  TOF_SimHits_RE4_Ring1_A_Minus  = new TH1F("TOF_SimHits_RE4_Ring1_A_Minus",  "TOF_SimHits_RE4_Ring1_A_Minus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_B_Minus  = new TH1F("TOF_SimHits_RE4_Ring1_B_Minus",  "TOF_SimHits_RE4_Ring1_B_Minus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_C_Minus  = new TH1F("TOF_SimHits_RE4_Ring1_C_Minus",  "TOF_SimHits_RE4_Ring1_C_Minus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_D_Minus  = new TH1F("TOF_SimHits_RE4_Ring1_D_Minus",  "TOF_SimHits_RE4_Ring1_D_Minus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE4_Ring1_E_Minus  = new TH1F("TOF_SimHits_RE4_Ring1_E_Minus",  "TOF_SimHits_RE4_Ring1_E_Minus", n_tof_re41, n1_tof_re41, n2_tof_re41);
  TOF_SimHits_RE1_Ring1_Minus.push_back(TOF_SimHits_RE1_Ring1_A_Minus);   TOF_SimHits_RE1_Ring1_Minus.push_back(TOF_SimHits_RE1_Ring1_B_Minus); 
  TOF_SimHits_RE1_Ring1_Minus.push_back(TOF_SimHits_RE1_Ring1_C_Minus);   TOF_SimHits_RE1_Ring1_Minus.push_back(TOF_SimHits_RE1_Ring1_D_Minus);   TOF_SimHits_RE1_Ring1_Minus.push_back(TOF_SimHits_RE1_Ring1_E_Minus);
  TOF_SimHits_RE2_Ring1_Minus.push_back(TOF_SimHits_RE2_Ring1_A_Minus);   TOF_SimHits_RE2_Ring1_Minus.push_back(TOF_SimHits_RE2_Ring1_B_Minus); 
  TOF_SimHits_RE2_Ring1_Minus.push_back(TOF_SimHits_RE2_Ring1_C_Minus);   TOF_SimHits_RE2_Ring1_Minus.push_back(TOF_SimHits_RE2_Ring1_D_Minus);   TOF_SimHits_RE2_Ring1_Minus.push_back(TOF_SimHits_RE2_Ring1_E_Minus);
  TOF_SimHits_RE3_Ring1_Minus.push_back(TOF_SimHits_RE3_Ring1_A_Minus);   TOF_SimHits_RE3_Ring1_Minus.push_back(TOF_SimHits_RE3_Ring1_B_Minus); 
  TOF_SimHits_RE3_Ring1_Minus.push_back(TOF_SimHits_RE3_Ring1_C_Minus);   TOF_SimHits_RE3_Ring1_Minus.push_back(TOF_SimHits_RE3_Ring1_D_Minus);   TOF_SimHits_RE3_Ring1_Minus.push_back(TOF_SimHits_RE3_Ring1_E_Minus);
  TOF_SimHits_RE4_Ring1_Minus.push_back(TOF_SimHits_RE4_Ring1_A_Minus);   TOF_SimHits_RE4_Ring1_Minus.push_back(TOF_SimHits_RE4_Ring1_B_Minus); 
  TOF_SimHits_RE4_Ring1_Minus.push_back(TOF_SimHits_RE4_Ring1_C_Minus);   TOF_SimHits_RE4_Ring1_Minus.push_back(TOF_SimHits_RE4_Ring1_D_Minus);   TOF_SimHits_RE4_Ring1_Minus.push_back(TOF_SimHits_RE4_Ring1_E_Minus);
  TOF_SimHits_RE_Ring1_Minus.push_back(TOF_SimHits_RE1_Ring1_Minus);  TOF_SimHits_RE_Ring1_Minus.push_back(TOF_SimHits_RE2_Ring1_Minus);
  TOF_SimHits_RE_Ring1_Minus.push_back(TOF_SimHits_RE3_Ring1_Minus);  TOF_SimHits_RE_Ring1_Minus.push_back(TOF_SimHits_RE4_Ring1_Minus);
  // ======================================================================================================================================

  // TOF All Hits inside the different rolls  :: Fits =====================================================================================
  Fit_SimHits_RE1_Ring1_A_Plus  = new TF1("Fit_SimHits_RE1_Ring1_A_Plus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_B_Plus  = new TF1("Fit_SimHits_RE1_Ring1_B_Plus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_C_Plus  = new TF1("Fit_SimHits_RE1_Ring1_C_Plus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_D_Plus  = new TF1("Fit_SimHits_RE1_Ring1_D_Plus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_E_Plus  = new TF1("Fit_SimHits_RE1_Ring1_E_Plus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE2_Ring1_A_Plus  = new TF1("Fit_SimHits_RE2_Ring1_A_Plus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_B_Plus  = new TF1("Fit_SimHits_RE2_Ring1_B_Plus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_C_Plus  = new TF1("Fit_SimHits_RE2_Ring1_C_Plus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_D_Plus  = new TF1("Fit_SimHits_RE2_Ring1_D_Plus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_E_Plus  = new TF1("Fit_SimHits_RE2_Ring1_E_Plus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE3_Ring1_A_Plus  = new TF1("Fit_SimHits_RE3_Ring1_A_Plus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_B_Plus  = new TF1("Fit_SimHits_RE3_Ring1_B_Plus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_C_Plus  = new TF1("Fit_SimHits_RE3_Ring1_C_Plus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_D_Plus  = new TF1("Fit_SimHits_RE3_Ring1_D_Plus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_E_Plus  = new TF1("Fit_SimHits_RE3_Ring1_E_Plus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE4_Ring1_A_Plus  = new TF1("Fit_SimHits_RE4_Ring1_A_Plus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_B_Plus  = new TF1("Fit_SimHits_RE4_Ring1_B_Plus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_C_Plus  = new TF1("Fit_SimHits_RE4_Ring1_C_Plus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_D_Plus  = new TF1("Fit_SimHits_RE4_Ring1_D_Plus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_E_Plus  = new TF1("Fit_SimHits_RE4_Ring1_E_Plus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE1_Ring1_Plus.push_back(Fit_SimHits_RE1_Ring1_A_Plus);   Fit_SimHits_RE1_Ring1_Plus.push_back(Fit_SimHits_RE1_Ring1_B_Plus); 
  Fit_SimHits_RE1_Ring1_Plus.push_back(Fit_SimHits_RE1_Ring1_C_Plus);   Fit_SimHits_RE1_Ring1_Plus.push_back(Fit_SimHits_RE1_Ring1_D_Plus);   Fit_SimHits_RE1_Ring1_Plus.push_back(Fit_SimHits_RE1_Ring1_E_Plus);
  Fit_SimHits_RE2_Ring1_Plus.push_back(Fit_SimHits_RE2_Ring1_A_Plus);   Fit_SimHits_RE2_Ring1_Plus.push_back(Fit_SimHits_RE2_Ring1_B_Plus); 
  Fit_SimHits_RE2_Ring1_Plus.push_back(Fit_SimHits_RE2_Ring1_C_Plus);   Fit_SimHits_RE2_Ring1_Plus.push_back(Fit_SimHits_RE2_Ring1_D_Plus);   Fit_SimHits_RE2_Ring1_Plus.push_back(Fit_SimHits_RE2_Ring1_E_Plus);
  Fit_SimHits_RE3_Ring1_Plus.push_back(Fit_SimHits_RE3_Ring1_A_Plus);   Fit_SimHits_RE3_Ring1_Plus.push_back(Fit_SimHits_RE3_Ring1_B_Plus); 
  Fit_SimHits_RE3_Ring1_Plus.push_back(Fit_SimHits_RE3_Ring1_C_Plus);   Fit_SimHits_RE3_Ring1_Plus.push_back(Fit_SimHits_RE3_Ring1_D_Plus);   Fit_SimHits_RE3_Ring1_Plus.push_back(Fit_SimHits_RE3_Ring1_E_Plus);
  Fit_SimHits_RE4_Ring1_Plus.push_back(Fit_SimHits_RE4_Ring1_A_Plus);   Fit_SimHits_RE4_Ring1_Plus.push_back(Fit_SimHits_RE4_Ring1_B_Plus); 
  Fit_SimHits_RE4_Ring1_Plus.push_back(Fit_SimHits_RE4_Ring1_C_Plus);   Fit_SimHits_RE4_Ring1_Plus.push_back(Fit_SimHits_RE4_Ring1_D_Plus);   Fit_SimHits_RE4_Ring1_Plus.push_back(Fit_SimHits_RE4_Ring1_E_Plus);
  Fit_SimHits_RE_Ring1_Plus.push_back(Fit_SimHits_RE1_Ring1_Plus);  Fit_SimHits_RE_Ring1_Plus.push_back(Fit_SimHits_RE2_Ring1_Plus);
  Fit_SimHits_RE_Ring1_Plus.push_back(Fit_SimHits_RE3_Ring1_Plus);  Fit_SimHits_RE_Ring1_Plus.push_back(Fit_SimHits_RE4_Ring1_Plus);
  Fit_SimHits_RE1_Ring1_A_Minus  = new TF1("Fit_SimHits_RE1_Ring1_A_Minus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_B_Minus  = new TF1("Fit_SimHits_RE1_Ring1_B_Minus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_C_Minus  = new TF1("Fit_SimHits_RE1_Ring1_C_Minus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_D_Minus  = new TF1("Fit_SimHits_RE1_Ring1_D_Minus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE1_Ring1_E_Minus  = new TF1("Fit_SimHits_RE1_Ring1_E_Minus",  "gaus", n1_tof_re11, n2_tof_re11);
  Fit_SimHits_RE2_Ring1_A_Minus  = new TF1("Fit_SimHits_RE2_Ring1_A_Minus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_B_Minus  = new TF1("Fit_SimHits_RE2_Ring1_B_Minus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_C_Minus  = new TF1("Fit_SimHits_RE2_Ring1_C_Minus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_D_Minus  = new TF1("Fit_SimHits_RE2_Ring1_D_Minus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE2_Ring1_E_Minus  = new TF1("Fit_SimHits_RE2_Ring1_E_Minus",  "gaus", n1_tof_re21, n2_tof_re21);
  Fit_SimHits_RE3_Ring1_A_Minus  = new TF1("Fit_SimHits_RE3_Ring1_A_Minus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_B_Minus  = new TF1("Fit_SimHits_RE3_Ring1_B_Minus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_C_Minus  = new TF1("Fit_SimHits_RE3_Ring1_C_Minus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_D_Minus  = new TF1("Fit_SimHits_RE3_Ring1_D_Minus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE3_Ring1_E_Minus  = new TF1("Fit_SimHits_RE3_Ring1_E_Minus",  "gaus", n1_tof_re31, n2_tof_re31);
  Fit_SimHits_RE4_Ring1_A_Minus  = new TF1("Fit_SimHits_RE4_Ring1_A_Minus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_B_Minus  = new TF1("Fit_SimHits_RE4_Ring1_B_Minus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_C_Minus  = new TF1("Fit_SimHits_RE4_Ring1_C_Minus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_D_Minus  = new TF1("Fit_SimHits_RE4_Ring1_D_Minus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE4_Ring1_E_Minus  = new TF1("Fit_SimHits_RE4_Ring1_E_Minus",  "gaus", n1_tof_re41, n2_tof_re41);
  Fit_SimHits_RE1_Ring1_Minus.push_back(Fit_SimHits_RE1_Ring1_A_Minus);   Fit_SimHits_RE1_Ring1_Minus.push_back(Fit_SimHits_RE1_Ring1_B_Minus); 
  Fit_SimHits_RE1_Ring1_Minus.push_back(Fit_SimHits_RE1_Ring1_C_Minus);   Fit_SimHits_RE1_Ring1_Minus.push_back(Fit_SimHits_RE1_Ring1_D_Minus);   Fit_SimHits_RE1_Ring1_Minus.push_back(Fit_SimHits_RE1_Ring1_E_Minus);
  Fit_SimHits_RE2_Ring1_Minus.push_back(Fit_SimHits_RE2_Ring1_A_Minus);   Fit_SimHits_RE2_Ring1_Minus.push_back(Fit_SimHits_RE2_Ring1_B_Minus); 
  Fit_SimHits_RE2_Ring1_Minus.push_back(Fit_SimHits_RE2_Ring1_C_Minus);   Fit_SimHits_RE2_Ring1_Minus.push_back(Fit_SimHits_RE2_Ring1_D_Minus);   Fit_SimHits_RE2_Ring1_Minus.push_back(Fit_SimHits_RE2_Ring1_E_Minus);
  Fit_SimHits_RE3_Ring1_Minus.push_back(Fit_SimHits_RE3_Ring1_A_Minus);   Fit_SimHits_RE3_Ring1_Minus.push_back(Fit_SimHits_RE3_Ring1_B_Minus); 
  Fit_SimHits_RE3_Ring1_Minus.push_back(Fit_SimHits_RE3_Ring1_C_Minus);   Fit_SimHits_RE3_Ring1_Minus.push_back(Fit_SimHits_RE3_Ring1_D_Minus);   Fit_SimHits_RE3_Ring1_Minus.push_back(Fit_SimHits_RE3_Ring1_E_Minus);
  Fit_SimHits_RE4_Ring1_Minus.push_back(Fit_SimHits_RE4_Ring1_A_Minus);   Fit_SimHits_RE4_Ring1_Minus.push_back(Fit_SimHits_RE4_Ring1_B_Minus); 
  Fit_SimHits_RE4_Ring1_Minus.push_back(Fit_SimHits_RE4_Ring1_C_Minus);   Fit_SimHits_RE4_Ring1_Minus.push_back(Fit_SimHits_RE4_Ring1_D_Minus);   Fit_SimHits_RE4_Ring1_Minus.push_back(Fit_SimHits_RE4_Ring1_E_Minus);
  Fit_SimHits_RE_Ring1_Minus.push_back(Fit_SimHits_RE1_Ring1_Minus);  Fit_SimHits_RE_Ring1_Minus.push_back(Fit_SimHits_RE2_Ring1_Minus);
  Fit_SimHits_RE_Ring1_Minus.push_back(Fit_SimHits_RE3_Ring1_Minus);  Fit_SimHits_RE_Ring1_Minus.push_back(Fit_SimHits_RE4_Ring1_Minus);
  // ======================================================================================================================================

  // TOF All Hits inside the different rolls  :: Statboxes  ===================================================================================
  Box_SimHits_RE1_Ring1_A_Plus  = new TPaveStats(); 
  Box_SimHits_RE1_Ring1_B_Plus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_C_Plus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_D_Plus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_E_Plus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_A_Plus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_B_Plus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_C_Plus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_D_Plus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_E_Plus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_A_Plus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_B_Plus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_C_Plus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_D_Plus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_E_Plus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_A_Plus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_B_Plus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_C_Plus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_D_Plus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_E_Plus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_Plus.push_back(Box_SimHits_RE1_Ring1_A_Plus);   Box_SimHits_RE1_Ring1_Plus.push_back(Box_SimHits_RE1_Ring1_B_Plus); 
  Box_SimHits_RE1_Ring1_Plus.push_back(Box_SimHits_RE1_Ring1_C_Plus);   Box_SimHits_RE1_Ring1_Plus.push_back(Box_SimHits_RE1_Ring1_D_Plus);   Box_SimHits_RE1_Ring1_Plus.push_back(Box_SimHits_RE1_Ring1_E_Plus);
  Box_SimHits_RE2_Ring1_Plus.push_back(Box_SimHits_RE2_Ring1_A_Plus);   Box_SimHits_RE2_Ring1_Plus.push_back(Box_SimHits_RE2_Ring1_B_Plus); 
  Box_SimHits_RE2_Ring1_Plus.push_back(Box_SimHits_RE2_Ring1_C_Plus);   Box_SimHits_RE2_Ring1_Plus.push_back(Box_SimHits_RE2_Ring1_D_Plus);   Box_SimHits_RE2_Ring1_Plus.push_back(Box_SimHits_RE2_Ring1_E_Plus);
  Box_SimHits_RE3_Ring1_Plus.push_back(Box_SimHits_RE3_Ring1_A_Plus);   Box_SimHits_RE3_Ring1_Plus.push_back(Box_SimHits_RE3_Ring1_B_Plus); 
  Box_SimHits_RE3_Ring1_Plus.push_back(Box_SimHits_RE3_Ring1_C_Plus);   Box_SimHits_RE3_Ring1_Plus.push_back(Box_SimHits_RE3_Ring1_D_Plus);   Box_SimHits_RE3_Ring1_Plus.push_back(Box_SimHits_RE3_Ring1_E_Plus);
  Box_SimHits_RE4_Ring1_Plus.push_back(Box_SimHits_RE4_Ring1_A_Plus);   Box_SimHits_RE4_Ring1_Plus.push_back(Box_SimHits_RE4_Ring1_B_Plus); 
  Box_SimHits_RE4_Ring1_Plus.push_back(Box_SimHits_RE4_Ring1_C_Plus);   Box_SimHits_RE4_Ring1_Plus.push_back(Box_SimHits_RE4_Ring1_D_Plus);   Box_SimHits_RE4_Ring1_Plus.push_back(Box_SimHits_RE4_Ring1_E_Plus);
  Box_SimHits_RE_Ring1_Plus.push_back(Box_SimHits_RE1_Ring1_Plus);  Box_SimHits_RE_Ring1_Plus.push_back(Box_SimHits_RE2_Ring1_Plus);
  Box_SimHits_RE_Ring1_Plus.push_back(Box_SimHits_RE3_Ring1_Plus);  Box_SimHits_RE_Ring1_Plus.push_back(Box_SimHits_RE4_Ring1_Plus);
  Box_SimHits_RE1_Ring1_A_Minus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_B_Minus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_C_Minus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_D_Minus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_E_Minus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_A_Minus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_B_Minus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_C_Minus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_D_Minus  = new TPaveStats();
  Box_SimHits_RE2_Ring1_E_Minus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_A_Minus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_B_Minus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_C_Minus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_D_Minus  = new TPaveStats();
  Box_SimHits_RE3_Ring1_E_Minus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_A_Minus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_B_Minus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_C_Minus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_D_Minus  = new TPaveStats();
  Box_SimHits_RE4_Ring1_E_Minus  = new TPaveStats();
  Box_SimHits_RE1_Ring1_Minus.push_back(Box_SimHits_RE1_Ring1_A_Minus);   Box_SimHits_RE1_Ring1_Minus.push_back(Box_SimHits_RE1_Ring1_B_Minus); 
  Box_SimHits_RE1_Ring1_Minus.push_back(Box_SimHits_RE1_Ring1_C_Minus);   Box_SimHits_RE1_Ring1_Minus.push_back(Box_SimHits_RE1_Ring1_D_Minus);   Box_SimHits_RE1_Ring1_Minus.push_back(Box_SimHits_RE1_Ring1_E_Minus);
  Box_SimHits_RE2_Ring1_Minus.push_back(Box_SimHits_RE2_Ring1_A_Minus);   Box_SimHits_RE2_Ring1_Minus.push_back(Box_SimHits_RE2_Ring1_B_Minus); 
  Box_SimHits_RE2_Ring1_Minus.push_back(Box_SimHits_RE2_Ring1_C_Minus);   Box_SimHits_RE2_Ring1_Minus.push_back(Box_SimHits_RE2_Ring1_D_Minus);   Box_SimHits_RE2_Ring1_Minus.push_back(Box_SimHits_RE2_Ring1_E_Minus);
  Box_SimHits_RE3_Ring1_Minus.push_back(Box_SimHits_RE3_Ring1_A_Minus);   Box_SimHits_RE3_Ring1_Minus.push_back(Box_SimHits_RE3_Ring1_B_Minus); 
  Box_SimHits_RE3_Ring1_Minus.push_back(Box_SimHits_RE3_Ring1_C_Minus);   Box_SimHits_RE3_Ring1_Minus.push_back(Box_SimHits_RE3_Ring1_D_Minus);   Box_SimHits_RE3_Ring1_Minus.push_back(Box_SimHits_RE3_Ring1_E_Minus);
  Box_SimHits_RE4_Ring1_Minus.push_back(Box_SimHits_RE4_Ring1_A_Minus);   Box_SimHits_RE4_Ring1_Minus.push_back(Box_SimHits_RE4_Ring1_B_Minus); 
  Box_SimHits_RE4_Ring1_Minus.push_back(Box_SimHits_RE4_Ring1_C_Minus);   Box_SimHits_RE4_Ring1_Minus.push_back(Box_SimHits_RE4_Ring1_D_Minus);   Box_SimHits_RE4_Ring1_Minus.push_back(Box_SimHits_RE4_Ring1_E_Minus);
  Box_SimHits_RE_Ring1_Minus.push_back(Box_SimHits_RE1_Ring1_Minus);  Box_SimHits_RE_Ring1_Minus.push_back(Box_SimHits_RE2_Ring1_Minus);
  Box_SimHits_RE_Ring1_Minus.push_back(Box_SimHits_RE3_Ring1_Minus);  Box_SimHits_RE_Ring1_Minus.push_back(Box_SimHits_RE4_Ring1_Minus);
  // ======================================================================================================================================


  VTX_PX = new TH1F("VTX_PX", "x position of Primary Vertex", n_v_x, n1_v_x, n2_v_x);
  VTX_PY = new TH1F("VTX_PY", "y position of Primary Vertex", n_v_x, n1_v_x, n2_v_x);
  VTX_PZ = new TH1F("VTX_PZ", "z position of Primary Vertex", n_v_z, n1_v_z, n2_v_z);

  VTX_AX = new TH1F("VTX_AX", "x position of All Vertices", p_v_x, p1_v_x, p2_v_x);
  VTX_AY = new TH1F("VTX_AY", "y position of All vertices", p_v_x, p1_v_x, p2_v_x);
  VTX_AZ = new TH1F("VTX_AZ", "z position of All vertices", p_v_z, p1_v_z, p2_v_z);

  Histo_XY_RE11_P = new TH2F("Histo_XY_RE11_P", "XY distribution of SimHits in RE+1/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE21_P = new TH2F("Histo_XY_RE21_P", "XY distribution of SimHits in RE+2/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE31_P = new TH2F("Histo_XY_RE31_P", "XY distribution of SimHits in RE+3/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE41_P = new TH2F("Histo_XY_RE41_P", "XY distribution of SimHits in RE+4/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE11_N = new TH2F("Histo_XY_RE11_N", "XY distribution of SimHits in RE-1/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE21_N = new TH2F("Histo_XY_RE21_N", "XY distribution of SimHits in RE-2/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE31_N = new TH2F("Histo_XY_RE31_N", "XY distribution of SimHits in RE-3/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_XY_RE41_N = new TH2F("Histo_XY_RE41_N", "XY distribution of SimHits in RE-4/1", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Histo_RZ_RE     = new TH2F("Histo_RZ_RE",     "RZ distribution of SimHits in RE",     n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Histo_XY_RE_P.push_back(Histo_XY_RE11_P); Histo_XY_RE_P.push_back(Histo_XY_RE21_P); Histo_XY_RE_P.push_back(Histo_XY_RE31_P); Histo_XY_RE_P.push_back(Histo_XY_RE41_P);
  Histo_XY_RE_N.push_back(Histo_XY_RE11_N); Histo_XY_RE_N.push_back(Histo_XY_RE21_N); Histo_XY_RE_N.push_back(Histo_XY_RE31_N); Histo_XY_RE_N.push_back(Histo_XY_RE41_N);


  for(int i=0; i<n_stations; ++i) {
    std::vector< std::vector < double > > dummy_xp;
    std::vector< std::vector < double > > dummy_yp;
    std::vector< std::vector < double > > dummy_zp;
    std::vector< std::vector < double > > dummy_rp;
    std::vector< std::vector < double > > dummy_xn;
    std::vector< std::vector < double > > dummy_yn;
    std::vector< std::vector < double > > dummy_zn;
    std::vector< std::vector < double > > dummy_rn;
    for(int j=0; j<n_rolls; ++j) {
      std::vector<double> dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4;
      std::vector<double> dummy_yp_1, dummy_yp_2, dummy_yp_3, dummy_yp_4;
      std::vector<double> dummy_zp_1, dummy_zp_2, dummy_zp_3, dummy_zp_4;
      std::vector<double> dummy_rp_1, dummy_rp_2, dummy_rp_3, dummy_rp_4;
      std::vector<double> dummy_xn_1, dummy_xn_2, dummy_xn_3, dummy_xn_4;
      std::vector<double> dummy_yn_1, dummy_yn_2, dummy_yn_3, dummy_yn_4;
      std::vector<double> dummy_zn_1, dummy_zn_2, dummy_zn_3, dummy_zn_4;
      std::vector<double> dummy_rn_1, dummy_rn_2, dummy_rn_3, dummy_rn_4;
      dummy_xp.push_back(dummy_xp_1);       dummy_xp.push_back(dummy_xp_2);      dummy_xp.push_back(dummy_xp_3);       dummy_xp.push_back(dummy_xp_4);
      dummy_yp.push_back(dummy_yp_1);       dummy_zp.push_back(dummy_yp_2);      dummy_yp.push_back(dummy_yp_3);       dummy_yp.push_back(dummy_yp_4);
      dummy_zp.push_back(dummy_zp_1);       dummy_yp.push_back(dummy_zp_2);      dummy_zp.push_back(dummy_zp_3);       dummy_zp.push_back(dummy_zp_4);
      dummy_rp.push_back(dummy_rp_1);       dummy_rp.push_back(dummy_rp_2);      dummy_rp.push_back(dummy_rp_3);       dummy_rp.push_back(dummy_rp_4);
      dummy_xn.push_back(dummy_xn_1);       dummy_xn.push_back(dummy_xn_2);      dummy_xn.push_back(dummy_xn_3);       dummy_xn.push_back(dummy_xn_4);
      dummy_yn.push_back(dummy_yn_1);       dummy_yn.push_back(dummy_yn_2);      dummy_yn.push_back(dummy_yn_3);       dummy_yn.push_back(dummy_yn_4);
      dummy_zn.push_back(dummy_zn_1);       dummy_zn.push_back(dummy_zn_2);      dummy_zn.push_back(dummy_zn_3);       dummy_zn.push_back(dummy_zn_4);
      dummy_rn.push_back(dummy_rn_1);       dummy_rn.push_back(dummy_rn_2);      dummy_rn.push_back(dummy_rn_3);       dummy_rn.push_back(dummy_rn_4);
    }
    x_p.push_back(dummy_xp);
    y_p.push_back(dummy_yp);
    z_p.push_back(dummy_zp);
    r_p.push_back(dummy_rp);
    x_n.push_back(dummy_xn);
    y_n.push_back(dummy_yn);
    z_n.push_back(dummy_zn);
    r_n.push_back(dummy_rn);
  }
  for(int i=0; i<n_rolls; ++i) {
    std::vector<double> dummy_rr, dummy_zr;
    r_r.push_back(dummy_rr);
    z_r.push_back(dummy_zr);
  }
}


MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer()
{
  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer DESTRUCTOR]"<<std::endl;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  ParticleIDs->Write();
  MuonHitEta->Write();
  MuonHitELoss->Write();
  MuonHitPhi->Write();

  EleHitELoss->Write();
  OtherHitELoss->Write();

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: Colours and Legends ]"<<std::endl;
  for(int i=0; i<n_stations; ++i) {
    TOF_SimHits_REn_Ring1_All_Plus[i]->SetLineColor(1);
    TOF_SimHits_REn_Ring1_All_Minus[i]->SetLineColor(1);
    for(int j=0; j<n_rolls; ++j) {
      TOF_SimHits_RE_Ring1_Plus[i][j]->SetLineColor(colours[j]);
      TOF_SimHits_RE_Ring1_Minus[i][j]->SetLineColor(colours[j]);
      TOF_SimHits_RE_Ring1_Plus[i][j]->SetMarkerStyle(5);
      TOF_SimHits_RE_Ring1_Minus[i][j]->SetMarkerStyle(5);
    }
  }

  TLegend *l1 = new TLegend(0.125,0.725,0.200,0.875,NULL,"brNDC");
  l1->SetLineColor(1); l1->SetLineStyle(1); l1->SetLineWidth(2); l1->SetFillColor(4000); l1->SetBorderSize(1);
  l1->AddEntry(TOF_SimHits_RE1_Ring1_All_Plus, "Roll","");
  l1->AddEntry(TOF_SimHits_RE1_Ring1_All_Plus, "All","f");
  l1->AddEntry(TOF_SimHits_RE1_Ring1_A_Plus, "A","f");
  l1->AddEntry(TOF_SimHits_RE1_Ring1_B_Plus, "B","f");
  l1->AddEntry(TOF_SimHits_RE1_Ring1_C_Plus, "C","f");
  l1->AddEntry(TOF_SimHits_RE1_Ring1_D_Plus, "D","f");
  l1->AddEntry(TOF_SimHits_RE1_Ring1_E_Plus, "E","f");

  TLegend *l2 = new TLegend(0.125,0.725,0.200,0.875,NULL,"brNDC");
  l2->SetLineColor(1); l2->SetLineStyle(1); l2->SetLineWidth(2); l2->SetFillColor(4000); l2->SetBorderSize(1);
  l2->AddEntry(TOF_SimHits_RE1_Ring1_All_Plus, "Roll","p");
  l2->AddEntry(TOF_SimHits_RE1_Ring1_All_Plus, "All","p");
  l2->AddEntry(TOF_SimHits_RE1_Ring1_A_Plus, "A","p");
  l2->AddEntry(TOF_SimHits_RE1_Ring1_B_Plus, "B","p");
  l2->AddEntry(TOF_SimHits_RE1_Ring1_C_Plus, "C","p");
  l2->AddEntry(TOF_SimHits_RE1_Ring1_D_Plus, "D","p");

  TLegend *l3 = new TLegend(0.125,0.725,0.200,0.875,NULL,"brNDC");
  l3->SetLineColor(1); l3->SetLineStyle(1); l3->SetLineWidth(2); l3->SetFillColor(4000); l3->SetBorderSize(1);
  l3->AddEntry(TOF_SimHits_RE1_Ring1_All_Plus, "Roll","p");
  l3->AddEntry(TOF_SimHits_RE1_Ring1_All_Plus, "All","p");
  l3->AddEntry(TOF_SimHits_RE1_Ring1_A_Plus, "A","p");
  l3->AddEntry(TOF_SimHits_RE1_Ring1_B_Plus, "B","p");


  // double DetBoundaries_Radii[4][5] = {{0.0}};
  // DetBoundaries_Radii[2][0] = 241.9;   DetBoundaries_Radii[2][1] = 284.6; DetBoundaries_Radii[2][2] = 330.0; 
  // DetBoundaries_Radii[3][0] = 264.4;   DetBoundaries_Radii[3][1] = 285.1; DetBoundaries_Radii[3][2] = 330.0; 

  // max and min values of the sensitive volumes
  for(int i=0; i<n_stations; ++i) {
    if(i<2) continue;
    for(int j=0; j<n_rolls; ++j) {
      // double r_min = *std::min_element(mybegin(r_r[i][j]), myend(r_r[i][j]));
      // double r_max = *std::max_element(mybegin(r_r[i][j]), myend(r_r[i][j]));
      // double z_min = *std::min_element(mybegin(z_r[i][j]), myend(z_r[i][j]));
      // double z_max = *std::max_element(mybegin(z_r[i][j]), myend(z_r[i][j]))
      double rp_min = +10000.0, rp_max = -10000.0, zp_min = +10000.0, zp_max = -10000.0, rn_min = +10000.0, rn_max = -10000.0, zn_min = +10000.0, zn_max = -10000.0;
      for (unsigned int k = 0; k < r_p[i][j].size(); ++k) { if (r_p[i][j][k] < rp_min) {rp_min = r_p[i][j][k];} if (r_p[i][j][k] > rp_max) {rp_max = r_p[i][j][k];}}
      for (unsigned int k = 0; k < r_n[i][j].size(); ++k) { if (r_n[i][j][k] < rn_min) {rn_min = r_n[i][j][k];} if (r_n[i][j][k] > rn_max) {rn_max = r_n[i][j][k];}}
      for (unsigned int k = 0; k < z_p[i][j].size(); ++k) { if (z_p[i][j][k] < zp_min) {zp_min = z_p[i][j][k];} if (z_p[i][j][k] > zp_max) {zp_max = z_p[i][j][k];}}
      for (unsigned int k = 0; k < z_n[i][j].size(); ++k) { if (z_n[i][j][k] < zn_min) {zn_min = z_n[i][j][k];} if (z_n[i][j][k] > zn_max) {zn_max = z_n[i][j][k];}}

      std::string letter = "";
      switch(j) {
      case 0 : letter = "A"; break;
      case 1 : letter = "B"; break;
      case 2 : letter = "C"; break;
      case 3 : letter = "D"; break;
      case 4 : letter = "E"; break;
      default : letter = "X";
      }
      std::cout<<"RE+"<<i+1<<"/1"<<letter<<" R min = "<<std::setw(9)<<rp_min<<" R max = "<<std::setw(9)<<rp_max<<" Z min = "<<std::setw(9)<<zp_min<<" Z max = "<<std::setw(9)<<zp_max<<std::endl;
      std::cout<<"RE-"<<i+1<<"/1"<<letter<<" R min = "<<std::setw(9)<<rn_min<<" R max = "<<std::setw(9)<<rn_max<<" Z min = "<<std::setw(9)<<zn_min<<" Z max = "<<std::setw(9)<<zn_max<<std::endl;
    }
  }
  // std::cout << *std::max_element(mybegin(cloud), myend(cloud)) << '\n';
  // std::cout << *std::min_element(mybegin(cloud), myend(cloud)) << '\n';
  // x_p[i][j], y_p[i][j], x_n[i][j], y_n[i][j], r_r[i][k]; z_r[i][k];


  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: XY and YZ Canvases :: 1]"<<std::endl;
  for(int i=0; i<n_stations; ++i) {
    std::vector<TGraph*> RE_Plus_XY_dummy, RE_Minus_XY_dummy;
    // std::vector<TEllipse*> DetBound_dummy;
    for(int j=0; j<n_rolls; ++j) {
      const int n_p = x_p[i][j].size();  double x_ap[n_p]; double y_ap[n_p];  // double z_ap[n_p];
      const int n_n = x_n[i][j].size();  double x_an[n_n]; double y_an[n_n];  // double z_an[n_n];

      // std::cout<<"[Station = "<<i+1<<"][Roll = "<<j+1<<"] Positive :: n_p = "<<n_p<<" := "<<x_p[i][j].size()<<" Negative :: n_n = "<<n_n<<" := "<<x_n[i][j].size()<<std::endl;
      // std::cout<<"x_p["<<i<<"]["<<j<<"].size() = "<<x_p[i][j].size()<<" and "<<y_p[i][j].size()<<std::endl;
      // std::cout<<"x_n["<<i<<"]["<<j<<"].size() = "<<x_n[i][j].size()<<" and "<<y_n[i][j].size()<<std::endl;
      for(int k=0; k<n_p; ++k) { x_ap[k] = x_p[i][j][k]; y_ap[k] = y_p[i][j][k]; /*z_ap[i] = z_p[i];*/ }
      for(int k=0; k<n_n; ++k) { x_an[k] = x_n[i][j][k]; y_an[k] = y_n[i][j][k]; /*z_an[i] = z_n[i];*/ }
      // std::cout<<"for loops ok, now will make the graphs"<<std::endl;
      TGraph * RE_Plus_XY_dd =  new TGraph(n_p, x_ap, y_ap);
      TGraph * RE_Minus_XY_dd =  new TGraph(n_n, x_an, y_an);
      RE_Plus_XY_dummy.push_back(RE_Plus_XY_dd);
      RE_Minus_XY_dummy.push_back(RE_Minus_XY_dd);
      // RE_Plus_XY[i][j]  = new TGraph(n_p, x_ap, y_ap);
      // RE_Minus_XY[i][j] = new TGraph(n_p, x_an, y_an);

      // if(DetBoundaries_Radii[i][j] !=0 ) {
      // TEllipse * ell = new TEllipse(0.0,0.0,DetBoundaries_Radii[i][j],DetBoundaries_Radii[i][j]);
      // DetBound_dummy.push_back(ell);
      // }
    }
    // outer one
    // TEllipse * ell = new TEllipse(0.0,0.0,DetBoundaries_Radii[i][4],DetBoundaries_Radii[i][4]);
    // DetBound_dummy.push_back(ell);

    RE_Plus_XY.push_back(RE_Plus_XY_dummy);
    RE_Minus_XY.push_back(RE_Minus_XY_dummy);
    // DetBound.push_back(DetBound_dummy);
  }
  for(int i=0; i<n_rolls; ++i) {
    const int n_r = r_r[i].size();  double r_ar[n_r]; double z_ar[n_r];
    for(int k=0; k<n_r; ++k) { r_ar[k] = r_r[i][k]; z_ar[k] = z_r[i][k]; }
    TGraph * RE_YZ_dummy = new TGraph(n_r, z_ar, r_ar);
    RE_YZ.push_back(RE_YZ_dummy); 
  }

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: XY and YZ Canvases :: 2]"<<std::endl;

  TCanvas * Canvas_RE1_Plus_XY, * Canvas_RE1_Minus_XY, * Canvas_RE1_Plus_TOF,  * Canvas_RE1_Minus_TOF; 

  Canvas_RE1_Plus_XY  = new TCanvas("Canvas_RE1_Plus_XY",  "Canvas_RE1_Plus_XY", 600, 600);
  Canvas_RE1_Minus_XY = new TCanvas("Canvas_RE1_Minus_XY", "Canvas_RE1_Minus_XY", 600, 600);
  Canvas_RE2_Plus_XY  = new TCanvas("Canvas_RE2_Plus_XY",  "Canvas_RE2_Plus_XY", 600, 600);
  Canvas_RE2_Minus_XY = new TCanvas("Canvas_RE2_Minus_XY", "Canvas_RE2_Minus_XY", 600, 600);
  Canvas_RE3_Plus_XY  = new TCanvas("Canvas_RE3_Plus_XY",  "Canvas_RE3_Plus_XY", 600, 600);
  Canvas_RE3_Minus_XY = new TCanvas("Canvas_RE3_Minus_XY", "Canvas_RE3_Minus_XY", 600, 600);
  Canvas_RE4_Plus_XY  = new TCanvas("Canvas_RE4_Plus_XY",  "Canvas_RE4_Plus_XY", 600, 600);
  Canvas_RE4_Minus_XY = new TCanvas("Canvas_RE4_Minus_XY", "Canvas_RE4_Minus_XY", 600, 600);
  // std::vector<TCanvas*> Canvas_RE_Plus_XY, Canvas_RE_Minus_XY;
  Canvas_RE_Plus_XY.push_back(Canvas_RE1_Plus_XY);   Canvas_RE_Plus_XY.push_back(Canvas_RE2_Plus_XY);   Canvas_RE_Plus_XY.push_back(Canvas_RE3_Plus_XY);   Canvas_RE_Plus_XY.push_back(Canvas_RE4_Plus_XY);
  Canvas_RE_Minus_XY.push_back(Canvas_RE1_Minus_XY); Canvas_RE_Minus_XY.push_back(Canvas_RE2_Minus_XY); Canvas_RE_Minus_XY.push_back(Canvas_RE3_Minus_XY); Canvas_RE_Minus_XY.push_back(Canvas_RE4_Minus_XY);
  Canvas_RE_YZ        = new TCanvas("Canvas_RE_YZ",        "Canvas_RE_YZ", 600, 600);

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: XY and YZ Canvases :: 3]"<<std::endl;

  // XY and RZ Graphs
  for(int i=0; i<n_stations; ++i) {
    Canvas_RE_Plus_XY[i]->cd();  
    for(int j=0; j<n_rolls; ++j) {
      RE_Plus_XY[i][j]->SetMarkerStyle(5);  
      RE_Plus_XY[i][j]->SetMarkerColor(colours[j]);
      if(j==0) { 
	RE_Plus_XY[i][j]->Draw("AP");
	RE_Plus_XY[i][j]->GetXaxis()->SetTitle("X [cm]"); RE_Plus_XY[i][j]->GetYaxis()->SetTitle("Y [cm]"); RE_Plus_XY[i][j]->SetTitle(pos_stations[i].c_str());
      }
      else RE_Plus_XY[i][j]->Draw("Psame");
      // DetBound[i][j]->Draw();
    }
    // DetBound[i][4]->Draw();
    l1->Draw();
  }

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: XY and YZ Canvases :: 4]"<<std::endl;

  for(int i=0; i<n_stations; ++i) {
    Canvas_RE_Minus_XY[i]->cd();  
    for(int j=0; j<n_rolls; ++j) {
      RE_Minus_XY[i][j]->SetMarkerStyle(5);  
      RE_Minus_XY[i][j]->SetMarkerColor(colours[j]);
      if(j==0) { 
	RE_Minus_XY[i][j]->Draw("AP");
	RE_Minus_XY[i][j]->GetXaxis()->SetTitle("X [cm]"); RE_Minus_XY[i][j]->GetYaxis()->SetTitle("Y [cm]"); RE_Minus_XY[i][j]->SetTitle(neg_stations[i].c_str());
      }
      else RE_Minus_XY[i][j]->Draw("Psame");
      // DetBound[i][j]->Draw();
    }
    // DetBound[i][4]->Draw();
    l1->Draw();
  }

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: XY and YZ Canvases :: 5]"<<std::endl;

  Canvas_RE_YZ->cd();
  for(int i=0; i<n_rolls; ++i) {
    RE_YZ[i]->SetMarkerStyle(5);
    RE_YZ[i]->SetMarkerColor(colours[i]);
    if(i==0) {
      RE_YZ[i]->Draw("AP");        
      RE_YZ[i]->GetXaxis()->SetTitle("Z [cm]"); RE_YZ[i]->GetYaxis()->SetTitle("R [cm]"); RE_YZ[i]->SetTitle("RE SimHits");
      RE_YZ[i]->GetXaxis()->SetRangeUser(-1500,+1500);
      RE_YZ[i]->GetYaxis()->SetRangeUser(100,350);
    }
    else RE_YZ[i]->Draw("Psame");
  }
  l1->Draw();
  
  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: XY and YZ Canvases :: 6]"<<std::endl;

  // VTX Canvas
  std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: VTX Canvas ]"<<std::endl;
  
  Canvas_VTX = new TCanvas("Canvas_VTX",          "Canvas_VTX", 800, 600);
  
  Canvas_VTX->cd();  Canvas_VTX->Divide(3,2);
  Canvas_VTX->cd(1); VTX_PX->Draw(); VTX_PX->GetXaxis()->SetTitle("x [cm]"); VTX_PX->GetYaxis()->SetTitle("entries [-]"); VTX_PX->SetTitle("Primary Vertex X Position");
  Canvas_VTX->cd(2); VTX_PY->Draw(); VTX_PY->GetXaxis()->SetTitle("y [cm]"); VTX_PY->GetYaxis()->SetTitle("entries [-]"); VTX_PY->SetTitle("Primary Vertex Y Position");
  Canvas_VTX->cd(3); VTX_PZ->Draw(); VTX_PZ->GetXaxis()->SetTitle("z [cm]"); VTX_PZ->GetYaxis()->SetTitle("entries [-]"); VTX_PZ->SetTitle("Primary Vertex Z Position");
  Canvas_VTX->cd(4); VTX_AX->Draw(); VTX_AX->GetXaxis()->SetTitle("x [cm]"); VTX_AX->GetYaxis()->SetTitle("entries [-]"); VTX_AX->SetTitle("SimVertices in RPC :: X");
  Canvas_VTX->cd(5); VTX_AY->Draw(); VTX_AY->GetXaxis()->SetTitle("y [cm]"); VTX_AY->GetYaxis()->SetTitle("entries [-]"); VTX_AY->SetTitle("SimVertices in RPC :: Y");
  Canvas_VTX->cd(6); VTX_AZ->Draw(); VTX_AZ->GetXaxis()->SetTitle("z [cm]"); VTX_AZ->GetYaxis()->SetTitle("entries [-]"); VTX_AZ->SetTitle("SimVertices in RPC :: Z");
  
  
  // TOF Canvas
  std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: TOF Canvases ]"<<std::endl;
  
  Canvas_RE1_Plus_TOF            = new TCanvas("Canvas_RE1_Plus_TOF",      "Canvas_RE1_Plus_TOF", 600, 600);
  Canvas_RE2_Plus_TOF            = new TCanvas("Canvas_RE2_Plus_TOF",      "Canvas_RE2_Plus_TOF", 600, 600);
  Canvas_RE3_Plus_TOF            = new TCanvas("Canvas_RE3_Plus_TOF",      "Canvas_RE3_Plus_TOF", 600, 600);
  Canvas_RE4_Plus_TOF            = new TCanvas("Canvas_RE4_Plus_TOF",      "Canvas_RE4_Plus_TOF", 600, 600);
  // std::vector<TCanvas*> Canvas_RE_Plus_TOF;
  Canvas_RE_Plus_TOF.push_back(Canvas_RE1_Plus_TOF); Canvas_RE_Plus_TOF.push_back(Canvas_RE2_Plus_TOF); Canvas_RE_Plus_TOF.push_back(Canvas_RE3_Plus_TOF); Canvas_RE_Plus_TOF.push_back(Canvas_RE4_Plus_TOF);

  Canvas_RE1_Minus_TOF            = new TCanvas("Canvas_RE1_Minus_TOF",      "Canvas_RE1_Minus_TOF", 600, 600);
  Canvas_RE2_Minus_TOF            = new TCanvas("Canvas_RE2_Minus_TOF",      "Canvas_RE2_Minus_TOF", 600, 600);
  Canvas_RE3_Minus_TOF            = new TCanvas("Canvas_RE3_Minus_TOF",      "Canvas_RE3_Minus_TOF", 600, 600);
  Canvas_RE4_Minus_TOF            = new TCanvas("Canvas_RE4_Minus_TOF",      "Canvas_RE4_Minus_TOF", 600, 600);
  // std::vector<TCanvas*> Canvas_RE_Minus_TOF;
  Canvas_RE_Minus_TOF.push_back(Canvas_RE1_Minus_TOF); Canvas_RE_Minus_TOF.push_back(Canvas_RE2_Minus_TOF); Canvas_RE_Minus_TOF.push_back(Canvas_RE3_Minus_TOF); Canvas_RE_Minus_TOF.push_back(Canvas_RE4_Minus_TOF);

  std::vector< TF1* > pos_all_fits, neg_all_fits; 
  std::vector< std::vector< TF1* > > pos_fits, neg_fits; 
  for(int i=0; i<n_stations; ++i) {
    TF1 * g1 = new TF1("g1","gaus",n1_tof,n2_tof);
    TF1 * h1 = new TF1("h1","gaus",n1_tof,n2_tof);
    pos_all_fits.push_back(g1);
    neg_all_fits.push_back(h1);
  }
  for(int i=0; i<n_stations; ++i) {
    std::vector< TF1* > pos_dummy, neg_dummy;
    for(int j=0; j<n_rolls; ++j) {
      TF1 * g1 = new TF1("g1","gaus",n1_tof,n2_tof);
      TF1 * h1 = new TF1("h1","gaus",n1_tof,n2_tof);
      pos_dummy.push_back(g1);
      neg_dummy.push_back(h1);
    }
    pos_fits.push_back(pos_dummy); 
    neg_fits.push_back(neg_dummy); 
  }
  // TF1 *g1p    = new TF1("g1p","gaus",n1_tof,n2_tof);  TF1 *g2p    = new TF1("g2p","gaus",n1_tof,n2_tof);  TF1 *g3p    = new TF1("g3p","gaus",n1_tof,n2_tof);  TF1 *g4p    = new TF1("g4p","gaus",n1_tof,n2_tof);
  // TF1 *g1n    = new TF1("g1p","gaus",n1_tof,n2_tof);  TF1 *g2n    = new TF1("g2p","gaus",n1_tof,n2_tof);  TF1 *g3n    = new TF1("g3p","gaus",n1_tof,n2_tof);  TF1 *g4n    = new TF1("g4p","gaus",n1_tof,n2_tof);

  TStyle *plain  = new TStyle("Plain","Plain Style (no colours/fill areas)");
  plain->cd();
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(1011);

  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::~MyREn1SimHitAnalyzer :: Draw & Save ]"<<std::endl;

  // Draw Simhits Positive Endcap
  for(int i=0; i<n_stations; ++i) {
    Canvas_RE_Plus_TOF[i]->cd();
    TOF_SimHits_REn_Ring1_All_Plus[i]->Draw("H");
    TOF_SimHits_REn_Ring1_All_Plus[i]->GetXaxis()->SetTitle("tof [ns]");
    TOF_SimHits_REn_Ring1_All_Plus[i]->GetYaxis()->SetTitle("entries [-]");
    TOF_SimHits_REn_Ring1_All_Plus[i]->SetTitle(pos_stations[i].c_str());
    pos_all_fits[i]->SetLineColor(1);
    pos_all_fits[i]->SetLineWidth(2);
    pos_all_fits[i]->SetLineStyle(2);
    TOF_SimHits_REn_Ring1_All_Plus[i]->Fit(pos_all_fits[i],"R");

    for(int j=0; j<n_rolls; ++j) {
      TOF_SimHits_RE_Ring1_Plus[i][j]->Draw("HsameS");
      Fit_SimHits_RE_Ring1_Plus[i][j]->SetLineColor(colours[j]);
      Fit_SimHits_RE_Ring1_Plus[i][j]->SetLineWidth(2);
      Fit_SimHits_RE_Ring1_Plus[i][j]->SetLineStyle(2);
      TOF_SimHits_RE_Ring1_Plus[i][j]->Fit(Fit_SimHits_RE_Ring1_Plus[i][j], "R+", "SameS");
    }
    l1->Draw();

    Canvas_RE_Plus_TOF[i]->Update();
    Box_SimHits_RE1_Ring1_Plus[i] = (TPaveStats*) TOF_SimHits_REn_Ring1_All_Plus[i]->GetListOfFunctions()->FindObject("stats");
    std::cout<<"ECCCO...."<<Box_SimHits_RE1_Ring1_Plus[i]<<std::endl;
    Box_SimHits_RE1_Ring1_Plus[i]->SetY1NDC(0.85); Box_SimHits_RE1_Ring1_Plus[i]->SetY2NDC(1.0);
    Box_SimHits_RE1_Ring1_Plus[i]->SetX1NDC(0.75); //Box_SimHits_RE1_Ring1_Plus[i]->SetX2NDC(1.0);                                                                                                                                  
    Canvas_RE_Plus_TOF[i]->Modified();

    for(int j=0; j<n_rolls; ++j) {
      Canvas_RE_Plus_TOF[i]->Update();
      Box_SimHits_RE_Ring1_Plus[i][j] = (TPaveStats*) TOF_SimHits_RE_Ring1_Plus[i][j]->GetListOfFunctions()->FindObject("stats");
      std::cout<<"Got statbox :: "<<Box_SimHits_RE_Ring1_Plus[i][j]<<std::endl;
      Box_SimHits_RE_Ring1_Plus[i][j]->SetTextColor(colours[j]);
      Box_SimHits_RE_Ring1_Plus[i][j]->SetX1NDC(0.75); //Box_SimHits_RE_Ring1_Plus[i][j]->SetX2NDC(0.95);                                                                                                                          
      if(j==0) {Box_SimHits_RE_Ring1_Plus[i][j]->SetY1NDC(0.70); Box_SimHits_RE_Ring1_Plus[i][j]->SetY2NDC(0.85);}
      if(j==1) {Box_SimHits_RE_Ring1_Plus[i][j]->SetY1NDC(0.55); Box_SimHits_RE_Ring1_Plus[i][j]->SetY2NDC(0.70);}
      if(j==2) {Box_SimHits_RE_Ring1_Plus[i][j]->SetY1NDC(0.40); Box_SimHits_RE_Ring1_Plus[i][j]->SetY2NDC(0.55);}
      if(j==3) {Box_SimHits_RE_Ring1_Plus[i][j]->SetY1NDC(0.25); Box_SimHits_RE_Ring1_Plus[i][j]->SetY2NDC(0.40);}
      if(j==4) {Box_SimHits_RE_Ring1_Plus[i][j]->SetY1NDC(0.10); Box_SimHits_RE_Ring1_Plus[i][j]->SetY2NDC(0.25);}
      Canvas_RE_Plus_TOF[i]->Modified();
    }
    Canvas_RE_Plus_TOF[i]->Write();
  }

  // Draw Simhits Negative Endcap
  for(int i=0; i<n_stations; ++i) {
    Canvas_RE_Minus_TOF[i]->cd();
    TOF_SimHits_REn_Ring1_All_Minus[i]->Draw("H");
    TOF_SimHits_REn_Ring1_All_Minus[i]->GetXaxis()->SetTitle("tof [ns]");
    TOF_SimHits_REn_Ring1_All_Minus[i]->GetYaxis()->SetTitle("entries [-]");
    TOF_SimHits_REn_Ring1_All_Minus[i]->SetTitle(neg_stations[i].c_str());
    neg_all_fits[i]->SetLineColor(1);
    neg_all_fits[i]->SetLineWidth(2);
    neg_all_fits[i]->SetLineStyle(2);
    TOF_SimHits_REn_Ring1_All_Minus[i]->Fit(neg_all_fits[i],"R");
    for(int j=0; j<n_rolls; ++j) {
      TOF_SimHits_RE_Ring1_Minus[i][j]->Draw("HsameS");
      Fit_SimHits_RE_Ring1_Minus[i][j]->SetLineColor(colours[j]);
      Fit_SimHits_RE_Ring1_Minus[i][j]->SetLineWidth(2);
      Fit_SimHits_RE_Ring1_Minus[i][j]->SetLineStyle(2);
      TOF_SimHits_RE_Ring1_Minus[i][j]->Fit(Fit_SimHits_RE_Ring1_Minus[i][j], "R+", "SameS");
    }
    l1->Draw();
    Canvas_RE_Minus_TOF[i]->Update();
    Box_SimHits_RE1_Ring1_Minus[i] = (TPaveStats*) TOF_SimHits_REn_Ring1_All_Minus[i]->GetListOfFunctions()->FindObject("stats");
    std::cout<<"ECCCO...."<<Box_SimHits_RE1_Ring1_Minus[i]<<std::endl;
    Box_SimHits_RE1_Ring1_Minus[i]->SetY1NDC(0.85); Box_SimHits_RE1_Ring1_Minus[i]->SetY2NDC(1.0);
    Box_SimHits_RE1_Ring1_Minus[i]->SetX1NDC(0.75); //Box_SimHits_RE1_Ring1_Minus[i]->SetX2NDC(1.0);                                                                                                                                
    Canvas_RE_Minus_TOF[i]->Modified();

    for(int j=0; j<n_rolls; ++j) {
      Canvas_RE_Minus_TOF[i]->Update();
      Box_SimHits_RE_Ring1_Minus[i][j] = (TPaveStats*) TOF_SimHits_RE_Ring1_Minus[i][j]->GetListOfFunctions()->FindObject("stats");
      std::cout<<"Got statbox :: "<<Box_SimHits_RE_Ring1_Minus[i][j]<<std::endl;
      Box_SimHits_RE_Ring1_Minus[i][j]->SetTextColor(colours[j]);
      Box_SimHits_RE_Ring1_Minus[i][j]->SetX1NDC(0.75); //Box_SimHits_RE_Ring1_Minus[i][j]->SetX2NDC(0.95);                                                                                                                        
      if(j==0) {Box_SimHits_RE_Ring1_Minus[i][j]->SetY1NDC(0.70); Box_SimHits_RE_Ring1_Minus[i][j]->SetY2NDC(0.85);}
      if(j==1) {Box_SimHits_RE_Ring1_Minus[i][j]->SetY1NDC(0.55); Box_SimHits_RE_Ring1_Minus[i][j]->SetY2NDC(0.70);}
      if(j==2) {Box_SimHits_RE_Ring1_Minus[i][j]->SetY1NDC(0.40); Box_SimHits_RE_Ring1_Minus[i][j]->SetY2NDC(0.55);}
      if(j==3) {Box_SimHits_RE_Ring1_Minus[i][j]->SetY1NDC(0.25); Box_SimHits_RE_Ring1_Minus[i][j]->SetY2NDC(0.40);}
      if(j==4) {Box_SimHits_RE_Ring1_Minus[i][j]->SetY1NDC(0.10); Box_SimHits_RE_Ring1_Minus[i][j]->SetY2NDC(0.25);}
      Canvas_RE_Minus_TOF[i]->Modified();
    }
    Canvas_RE_Minus_TOF[i]->Write();
  }




  for(int i=0; i<n_stations; ++i) { Canvas_RE_Plus_XY[i]->Write(); Canvas_RE_Minus_XY[i]->Write(); }
  Canvas_RE_YZ->Write();
  Canvas_VTX->Write();
  // for(int i=0; i<4; ++i) { Canvas_RE_Plus_TOF[i]->Write(); Canvas_RE_Minus_TOF[i]->Write(); }

  for(int i=0; i<n_stations; ++i) {
    Histo_XY_RE_P[i]->Write();
    Histo_XY_RE_N[i]->Write();
  }
  Histo_RZ_RE->Write();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyREn1SimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::analyze ]"<<std::endl;

  // SimHits
  if(debug)  std::cout << " Getting the SimHits " <<std::endl;
  std::vector<edm::Handle<edm::PSimHitContainer> > theSimHitContainers;
  iEvent.getManyByType(theSimHitContainers);
  if(debug)  std::cout << " The number of SimHit Containers is  " << theSimHitContainers.size() <<std::endl;
  std::vector<PSimHit> theSimHits;
  for (int i = 0; i < int(theSimHitContainers.size()); ++i) {
    theSimHits.insert(theSimHits.end(),theSimHitContainers.at(i)->begin(),theSimHitContainers.at(i)->end());
  }
  // SimTracks
  std::vector<SimTrack> theSimTracks;
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel("g4SimHits",SimTk);
  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  if(debug)  std::cout << "This Event has " <<  theSimTracks.size() << " sim tracks " << std::endl;
  // SimVertices
  std::vector<SimVertex> theSimVertices; 
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel("g4SimHits",SimVtx);
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  if(debug)  std::cout << "This Event has " <<  theSimVertices.size() << " sim vertices " << std::endl;
  // SimHits
  std::vector<PSimHit> theRPCSimHits;
  edm::Handle<edm::PSimHitContainer> RPCHits;
  iEvent.getByToken(RPCSimHit_Token, RPCHits);
  theRPCSimHits.insert(theRPCSimHits.end(),RPCHits->begin(),RPCHits->end());
  if(debug)  std::cout << "This Event has " <<  theRPCSimHits.size() << " rpc sim hits " << std::endl;


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
      // if(debug)  std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<simtrack.momentum().pt()<<" GeV/c from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      // double Px = simtrack.momentum().x(); double Py = simtrack.momentum().y(); double mu_pt = sqrt(pow(Px,2)+pow(Py,2));
      // std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<mu_pt<<" GeV/c from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      // std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      if(mu_pt > 0.95*simPtCut) {
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
      // if(debug)  std::cout<<"Muon SimVertex Found: x = "<<vtx_x<<" y = "<<vtx_y<<" z = "<<vtx_z<<std::endl;
    }
    else if(fabs(vtx_z)<25 && fabs(vtx_x)<5 && fabs(vtx_y)<5) {
      VTX_AX->Fill(simvertex.position().x()); VTX_AY->Fill(simvertex.position().y()); VTX_AZ->Fill(simvertex.position().z());
      // if(debug)  std::cout<<"PileUp SimVertex Found: x = "<<simvertex.position().x()<<" y = "<<simvertex.position().y()<<" z = "<<simvertex.position().z()<<std::endl;
    }
    else {
      // VTX_AX->Fill(simvertex.position().x()); VTX_AY->Fill(simvertex.position().y()); VTX_AZ->Fill(simvertex.position().z());
      // if(debug)  std::cout<<"RPC SimVertex Found: x = "<<simvertex.position().x()<<" y = "<<simvertex.position().y()<<" z = "<<simvertex.position().z()<<std::endl;
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

      ParticleIDs->Fill((*iHit).particleType());

      if (fabs((*iHit).particleType())==11)                                     { EleHitELoss->Fill((*iHit).energyLoss());   /*std::cout<<"Electron hit loss = "<<(*iHit).energyLoss()<<std::endl;*/}
      if (fabs((*iHit).particleType())!=11 && fabs((*iHit).particleType())!=13) { OtherHitELoss->Fill((*iHit).energyLoss()); /*std::cout<<"PDG ID "<<(*iHit).particleType()<<" hit loss = "<<(*iHit).energyLoss()<<std::endl;*/}
      if (fabs((*iHit).particleType())==13) {

	MuonHitEta->Fill(RPCGlobalPoint.eta());
	MuonHitPhi->Fill(RPCGlobalPoint.phi());
	MuonHitELoss->Fill((*iHit).energyLoss());
	// std::cout<<"Muon hit loss = "<<(*iHit).energyLoss()<<std::endl;

	// Only REn/1
	if ((rollId.ring()==1) && rollId.region() != 0) {
	  if(debug) {
	    std::cout<<"RPC SimHit in "<<std::setw(12)<<(int)rollId<<" a.k.a. "<<std::setw(24)<<rpcsrv.name()<<" details: "<<std::setw(24)<<rollId<<" | time t = "<<std::setw(8)<<(*iHit).timeOfFlight();
	    std::cout<<" | z = "<<std::setw(8)<<RPCGlobalPoint.z()<<" | r = "<<std::setw(8)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(8)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(8)<<RPCGlobalPoint.eta();
	    std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
	    // std::cout<<"              "<<" | id = "<<std::setw(4)<<(*iHit).particleType()<<" | tof = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | theta = "<<std::setw(12)<<(*iHit).thetaAtEntry()<<std::endl;
	  }
	  // Timing as calculated from the global point
	  // all distances are in cm and all time is in ns, therefore c is in cm/ns
	  // GlobalPoint gp = RPCSurface.toGlobal(LocalPoint(0,0,0)); // global point of the middle of the rol
	  // float c = 29.9792458;
	  // float time0 = gp.mag()/c; // time w.r.t. the middle of the roll
	  // float time1 = RPCGlobalPoint.mag()/c;
	  // float time2 = sqrt(pow(RPCGlobalPoint.x()-VtxGlobalPoint.x(),2)+pow(RPCGlobalPoint.y()-VtxGlobalPoint.y(),2)+pow(RPCGlobalPoint.z()-VtxGlobalPoint.z(),2))/c;
	  if (fabs((*iHit).particleType())==13) {
	    if(rollId.region() == 1) {
	      // if(debug) std::cout<<"Positive Endcap"<<std::endl;
	      x_p[rollId.station()-1][rollId.roll()-1].push_back(RPCGlobalPoint.x()); 
	      y_p[rollId.station()-1][rollId.roll()-1].push_back(RPCGlobalPoint.y()); 
	      // std::cout<<"Positive Endcap :: before pushing back in z-vector"<<std::endl;
	      z_p[rollId.station()-1][rollId.roll()-1].push_back(RPCGlobalPoint.z());
	      // std::cout<<"Positive Endcap :: before pushing back in r-vector"<<std::endl;
	      r_p[rollId.station()-1][rollId.roll()-1].push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2)));
	      // std::cout<<"Positive Endcap :: after  pushing back in r-vector"<<std::endl;
	      // TOF_SimHits_RE_Ring1_Plus[rollId.station()-1]->Fill((*iHit).timeOfFlight());
	      TOF_SimHits_REn_Ring1_All_Plus[rollId.station()-1]->Fill((*iHit).timeOfFlight());
	      // if(debug) std::cout<<"Before TOF Fill"<<std::endl;
	      TOF_SimHits_RE_Ring1_Plus[rollId.station()-1][rollId.roll()-1]->Fill((*iHit).timeOfFlight());
	      // if(debug) std::cout<<"After TOF Fill"<<std::endl;
	      // std::cout<<"TOF_SimHits_RE_Ring1_Plus["<<rollId.station()-1<<"]["<<rollId.roll()-1<<"]->Fill("<<(*iHit).timeOfFlight()<<")"<<std::endl;
	      Histo_XY_RE_P[rollId.station()-1]->Fill(RPCGlobalPoint.x(),RPCGlobalPoint.y());
	    }
	    if(rollId.region() == -1) {
	      // if(debug) std::cout<<"Negative Endcap"<<std::endl;
	      x_n[rollId.station()-1][rollId.roll()-1].push_back(RPCGlobalPoint.x()); 
	      y_n[rollId.station()-1][rollId.roll()-1].push_back(RPCGlobalPoint.y()); 
	      // std::cout<<"Negative Endcap :: before pushing back in z-vector"<<std::endl;
	      z_n[rollId.station()-1][rollId.roll()-1].push_back(RPCGlobalPoint.z());
	      // std::cout<<"Negative Endcap :: before pushing back in r-vector"<<std::endl;
	      r_n[rollId.station()-1][rollId.roll()-1].push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2)));
	      // std::cout<<"Negative Endcap :: after  pushing back in r-vector"<<std::endl;
	      // TOF_SimHits_RE_Ring1_Minus[rollId.station()-1]->Fill((*iHit).timeOfFlight());
	      TOF_SimHits_REn_Ring1_All_Minus[rollId.station()-1]->Fill((*iHit).timeOfFlight());
	      // if(debug) std::cout<<"Before TOF Fill"<<std::endl;
	      TOF_SimHits_RE_Ring1_Minus[rollId.station()-1][rollId.roll()-1]->Fill((*iHit).timeOfFlight());
	      // if(debug) std::cout<<"After TOF Fill"<<std::endl;
	      // std::cout<<"TOF_SimHits_RE_Ring1_Minus["<<rollId.station()-1<<"]["<<rollId.roll()-1<<"]->Fill("<<(*iHit).timeOfFlight()<<")"<<std::endl;
	      Histo_XY_RE_N[rollId.station()-1]->Fill(RPCGlobalPoint.x(),RPCGlobalPoint.y());
	    }
	    // RZ graph endcap
	    // std::cout<<"Before Filling R-Z plot"<<std::endl;
	    // std::cout<<"Filling R-Z plot["<<rollId.roll()-1<<"] with R = "<<sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2))<<" and Z = "<<RPCGlobalPoint.z()<<std::endl;
	    r_r[rollId.roll()-1].push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2))); 
	    z_r[rollId.roll()-1].push_back(RPCGlobalPoint.z());
	    // std::cout<<"After Filling R-Z plot"<<std::endl;
	    Histo_RZ_RE->Fill(RPCGlobalPoint.z(),sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2)));
	  }
	}
      }
      /*
      if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::CSC){ // Only CSCs
        CSCDetId rollId(theDetUnitId);
        // CSCGeomServ rpcsrv(rollId);
        GlobalPoint CSCGlobalPoint = cscGeom->idToDet(rollId)->toGlobal((*iHit).localPosition());
        if(debug) {
	  // std::cout<<"CSC SimHit in "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
	  std::cout<<"CSC SimHit in "<<std::setw(24)<<rollId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
	  std::cout<<" | r = "<<std::setw(12)<<CSCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<CSCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<CSCGlobalPoint.eta();
	  std::cout<<" | global position = "<<CSCGlobalPoint<<std::endl;
        }
      }
      */
    }
  }
  // std::cout<<"End of Analyze loop"<<std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyREn1SimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyREn1SimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyREn1SimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  if(debug) std::cout<<"[ MyREn1SimHitAnalyzer::beginRun ]"<<std::endl;
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
}


// ------------ method called when ending the processing of a run  ------------
void 
MyREn1SimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyREn1SimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyREn1SimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyREn1SimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyREn1SimHitAnalyzer);
