// -*- C++ -*-
//
// Package:    MyME11SimHitAnalyzer
// Class:      MyME11SimHitAnalyzer
// 
/**\class MyME11SimHitAnalyzer MyME11SimHitAnalyzer.cc MyAnalyzers/MyME11SimHitAnalyzer/src/MyME11SimHitAnalyzer.cc

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
#include <math.h>


// root include files
#include <TRandom.h>
#include "TROOT.h"
#include "TStyle.h"
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
#include "TLatex.h"



// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
// OLD :: #include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// why is this one not included by default
#include "FWCore/Framework/interface/ESHandle.h"


// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
// not in 71X release
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>

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
// #include "FastSimulation/Tracking/test/FastTrackAnalyzer.h"
// DetIds
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include "DataFormats/MuonDetId/interface/DTWireId.h"
// not in 71X release
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/MuonDetId/interface/ME0DetId.h>

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

// OLD :: class MyME11SimHitAnalyzer : public edm::EDAnalyzer {
class MyME11SimHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit MyME11SimHitAnalyzer(const edm::ParameterSet&);
      ~MyME11SimHitAnalyzer();
  edm::ESHandle <RPCGeometry> rpcGeom;
  edm::ESHandle <CSCGeometry> cscGeom;
  edm::ESHandle <DTGeometry>  dtGeom;
  // not in 71X release
  edm::ESHandle <GEMGeometry> gemGeom;
  edm::ESHandle <ME0Geometry> me0Geom;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  double getLumi(int, int, int);
  double getPU(double, int, int);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  std::string pdfFileNameBase, pdfFileName, rootFileName;
  double bunchspacing;
  double comenergy;
  double maxsimtime;
  double varEDepCuteV;
  bool phys_debug, tech_debug;
  TFile * outputfile;


  // required for 7XY: registration of the data access
  // see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent
  edm::EDGetTokenT<edm::SimTrackContainer>  SIMTrack_Token;
  edm::EDGetTokenT<edm::SimVertexContainer> SIMVertex_Token;

  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeom_Token;
  edm::ESGetToken<DTGeometry,  MuonGeometryRecord> dtGeom_Token;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeom_Token;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeom_Token;

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_deposits;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH2F * RPCb_el_hits, * RPCb_mu_hits, * RPCb_pi_hits, * RPCb_ka_hits, * RPCb_p_hits, * RPCb_n_hits, * RPCb_g_hits, * RPCb_N_hits, * RPCb_OH_hits, * RPCb_All_hits, * RPCb_HIP_hits;
  TH2F * RPCf_el_hits, * RPCf_mu_hits, * RPCf_pi_hits, * RPCf_ka_hits, * RPCf_p_hits, * RPCf_n_hits, * RPCf_g_hits, * RPCf_N_hits, * RPCf_OH_hits, * RPCf_All_hits, * RPCf_HIP_hits;
  TH2F * CSC_el_hits,  * CSC_mu_hits,  * CSC_pi_hits,  * CSC_ka_hits,  * CSC_p_hits,  * CSC_n_hits,  * CSC_g_hits,  * CSC_N_hits,  * CSC_OH_hits,  * CSC_All_hits,  * CSC_HIP_hits;
  TH2F * DT_el_hits,   * DT_mu_hits,   * DT_pi_hits,   * DT_ka_hits,   * DT_p_hits,   * DT_n_hits,   * DT_g_hits,   * DT_N_hits,   * DT_OH_hits,   * DT_All_hits,   * DT_HIP_hits;
  TH2F * GEM_el_hits,  * GEM_mu_hits,  * GEM_pi_hits,  * GEM_ka_hits,  * GEM_p_hits,  * GEM_n_hits,  * GEM_g_hits,  * GEM_N_hits,  * GEM_OH_hits,  * GEM_All_hits,  * GEM_HIP_hits;
  TH2F * ME0_el_hits,  * ME0_mu_hits,  * ME0_pi_hits,  * ME0_ka_hits,  * ME0_p_hits,  * ME0_n_hits,  * ME0_g_hits,  * ME0_N_hits,  * ME0_OH_hits,  * ME0_All_hits,  * ME0_HIP_hits;

  TH1F * RPCb_el_kins, * RPCb_mu_kins, * RPCb_ha_kins;
  TH1F * RPCf_el_kins, * RPCf_mu_kins, * RPCf_ha_kins;
  TH1F * CSC_el_kins,  * CSC_mu_kins,  * CSC_ha_kins;
  TH1F * DT_el_kins,   * DT_mu_kins,   * DT_ha_kins;
  TH1F * GEM_el_kins,  * GEM_mu_kins,  * GEM_ha_kins;
  TH1F * ME0_el_kins,  * ME0_mu_kins,  * ME0_ha_kins;

  TH2F * RPCb_el_deposits, * RPCb_mu_deposits, * RPCb_pi_deposits, * RPCb_ka_deposits, * RPCb_p_deposits, * RPCb_n_deposits, * RPCb_g_deposits, * RPCb_N_deposits, * RPCb_OH_deposits, * RPCb_All_deposits, * RPCb_HIP_deposits;
  TH2F * RPCf_el_deposits, * RPCf_mu_deposits, * RPCf_pi_deposits, * RPCf_ka_deposits, * RPCf_p_deposits, * RPCf_n_deposits, * RPCf_g_deposits, * RPCf_N_deposits, * RPCf_OH_deposits, * RPCf_All_deposits, * RPCf_HIP_deposits;
  TH2F * CSC_el_deposits,  * CSC_mu_deposits,  * CSC_pi_deposits,  * CSC_ka_deposits,  * CSC_p_deposits,  * CSC_n_deposits,  * CSC_g_deposits,  * CSC_N_deposits,  * CSC_OH_deposits,  * CSC_All_deposits,  * CSC_HIP_deposits;
  TH2F * DT_el_deposits,   * DT_mu_deposits,   * DT_pi_deposits,   * DT_ka_deposits,   * DT_p_deposits,   * DT_n_deposits,   * DT_g_deposits,   * DT_N_deposits,   * DT_OH_deposits,   * DT_All_deposits,   * DT_HIP_deposits;
  TH2F * GEM_el_deposits,  * GEM_mu_deposits,  * GEM_pi_deposits,  * GEM_ka_deposits,  * GEM_p_deposits,  * GEM_n_deposits,  * GEM_g_deposits,  * GEM_N_deposits,  * GEM_OH_deposits,  * GEM_All_deposits,  * GEM_HIP_deposits;
  TH2F * ME0_el_deposits,  * ME0_mu_deposits,  * ME0_pi_deposits,  * ME0_ka_deposits,  * ME0_p_deposits,  * ME0_n_deposits,  * ME0_g_deposits,  * ME0_N_deposits,  * ME0_OH_deposits,  * ME0_All_deposits,  * ME0_HIP_deposits;

  TH1F * RPCb_el_deps, * RPCb_mu_deps, * RPCb_ha_deps, * RPCb_el_tof, * RPCb_mu_tof, * RPCb_ha_tof;
  TH1F * RPCf_el_deps, * RPCf_mu_deps, * RPCf_ha_deps, * RPCf_el_tof, * RPCf_mu_tof, * RPCf_ha_tof;
  TH1F * CSC_el_deps,  * CSC_mu_deps,  * CSC_ha_deps,  * CSC_el_tof,  * CSC_mu_tof,  * CSC_ha_tof;
  TH1F * DT_el_deps,   * DT_mu_deps,   * DT_ha_deps,   * DT_el_tof,   * DT_mu_tof,   * DT_ha_tof;
  TH1F * GEM_el_deps,  * GEM_mu_deps,  * GEM_ha_deps,  * GEM_el_tof,  * GEM_mu_tof,  * GEM_ha_tof;
  TH1F * ME0_el_deps,  * ME0_mu_deps,  * ME0_ha_deps,  * ME0_el_tof,  * ME0_mu_tof,  * ME0_ha_tof;

  TH1F * RPCb_hits_tof, * RPCf_hits_tof, * CSC_hits_tof, * DT_hits_tof, * GEM_hits_tof, * ME0_hits_tof;
  TH1F * RPCb_hits_eta, * RPCf_hits_eta, * CSC_hits_eta, * DT_hits_eta, * GEM_hits_eta, * ME0_hits_eta;
  TH1F * RPCb_hits_phi, * RPCf_hits_phi, * CSC_hits_phi, * DT_hits_phi, * GEM_hits_phi, * ME0_hits_phi; 
  TH1F * RPCb_hits_lin, * RPCf_hits_lin, * CSC_hits_lin, * DT_hits_lin, * GEM_hits_lin, * ME0_hits_lin;

  TH1F * MB1_hits_phi, * MB2_hits_phi, *MB3_hits_phi, * MB4_hits_phi;
  TH1F * RB1_hits_phi, * RB2_hits_phi, *RB3_hits_phi, * RB4_hits_phi;
  TH1F * ME11_hits_phi, * ME12_hits_phi, * ME13_hits_phi, * ME21_hits_phi, * ME22_hits_phi, * ME31_hits_phi,  *ME32_hits_phi, * ME41_hits_phi, * ME42_hits_phi;
  TH1F * RE12_hits_phi, * RE13_hits_phi, * RE22_hits_phi, * RE23_hits_phi, * RE32_hits_phi,  *RE33_hits_phi, * RE42_hits_phi, * RE43_hits_phi;
  // TH1F * RB4_hits_tof, * RE4_hits_tof, * ME4_hits_tof, * MB4_hits_tof;
  // TH1F * RB4_hits_phi, * RE4_hits_phi, * ME4_hits_phi, * MB4_hits_phi;

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_Nuclei;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH2F * Muon_Nuclei_A_Z, * RPCb_Nuclei_A_Z, * RPCf_Nuclei_A_Z, * CSC_Nuclei_A_Z, * DT_Nuclei_A_Z, * GEM_Nuclei_A_Z, * ME0_Nuclei_A_Z;
  TH1F * Muon_Nuclei_A, * RPCb_Nuclei_A, * RPCf_Nuclei_A, * CSC_Nuclei_A, * DT_Nuclei_A, * GEM_Nuclei_A, * ME0_Nuclei_A;
  TH1F * Muon_Nuclei_Z, * RPCb_Nuclei_Z, * RPCf_Nuclei_Z, * CSC_Nuclei_Z, * DT_Nuclei_Z, * GEM_Nuclei_Z, * ME0_Nuclei_Z;
  TH1F * Muon_Nuclei_List, * RPCb_Nuclei_List, * RPCf_Nuclei_List, * CSC_Nuclei_List, * DT_Nuclei_List, * GEM_Nuclei_List, * ME0_Nuclei_List;

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_radius;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * ME0_el_hits_R,  * ME0_mu_hits_R,  * ME0_pi_hits_R,  * ME0_ka_hits_R,  * ME0_p_hits_R,  * ME0_n_hits_R,  * ME0_g_hits_R,  * ME0_N_hits_R, * ME0_OH_hits_R,  * ME0_All_hits_R,  * ME0_HIP_hits_R;
  TH1F * G11_el_hits_R,  * G11_mu_hits_R,  * G11_pi_hits_R,  * G11_ka_hits_R,  * G11_p_hits_R,  * G11_n_hits_R,  * G11_g_hits_R,  * G11_N_hits_R, * G11_OH_hits_R,  * G11_All_hits_R,  * G11_HIP_hits_R;
  TH1F * G21_el_hits_R,  * G21_mu_hits_R,  * G21_pi_hits_R,  * G21_ka_hits_R,  * G21_p_hits_R,  * G21_n_hits_R,  * G21_g_hits_R,  * G21_N_hits_R, * G21_OH_hits_R,  * G21_All_hits_R,  * G21_HIP_hits_R;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * G11_L1_el_hits_R,  * G11_L1_mu_hits_R,  * G11_L1_pi_hits_R,  * G11_L1_ka_hits_R,  * G11_L1_p_hits_R,  * G11_L1_n_hits_R,  * G11_L1_g_hits_R,  * G11_L1_N_hits_R, * G11_L1_OH_hits_R,  * G11_L1_All_hits_R,  * G11_L1_HIP_hits_R;
  TH1F * G11_L2_el_hits_R,  * G11_L2_mu_hits_R,  * G11_L2_pi_hits_R,  * G11_L2_ka_hits_R,  * G11_L2_p_hits_R,  * G11_L2_n_hits_R,  * G11_L2_g_hits_R,  * G11_L2_N_hits_R, * G11_L2_OH_hits_R,  * G11_L2_All_hits_R,  * G11_L2_HIP_hits_R;
  TH1F * G11_Od_el_hits_R,  * G11_Od_mu_hits_R,  * G11_Od_pi_hits_R,  * G11_Od_ka_hits_R,  * G11_Od_p_hits_R,  * G11_Od_n_hits_R,  * G11_Od_g_hits_R,  * G11_Od_N_hits_R, * G11_Od_OH_hits_R,  * G11_Od_All_hits_R,  * G11_Od_HIP_hits_R;
  TH1F * G11_Ev_el_hits_R,  * G11_Ev_mu_hits_R,  * G11_Ev_pi_hits_R,  * G11_Ev_ka_hits_R,  * G11_Ev_p_hits_R,  * G11_Ev_n_hits_R,  * G11_Ev_g_hits_R,  * G11_Ev_N_hits_R, * G11_Ev_OH_hits_R,  * G11_Ev_All_hits_R,  * G11_Ev_HIP_hits_R;
  TH1F * M11_Od_el_hits_R,  * M11_Od_mu_hits_R,  * M11_Od_pi_hits_R,  * M11_Od_ka_hits_R,  * M11_Od_p_hits_R,  * M11_Od_n_hits_R,  * M11_Od_g_hits_R,  * M11_Od_N_hits_R, * M11_Od_OH_hits_R,  * M11_Od_All_hits_R,  * M11_Od_HIP_hits_R;
  TH1F * M11_Ev_el_hits_R,  * M11_Ev_mu_hits_R,  * M11_Ev_pi_hits_R,  * M11_Ev_ka_hits_R,  * M11_Ev_p_hits_R,  * M11_Ev_n_hits_R,  * M11_Ev_g_hits_R,  * M11_Ev_N_hits_R, * M11_Ev_OH_hits_R,  * M11_Ev_All_hits_R,  * M11_Ev_HIP_hits_R;

  // 2B Further implemented ....
  TH1F * M11_Od_L1_All_hits_R, * M11_Od_L2_All_hits_R, * M11_Od_L3_All_hits_R, * M11_Od_L4_All_hits_R, * M11_Od_L5_All_hits_R, * M11_Od_L6_All_hits_R; 
  TH1F * M11_Ev_L1_All_hits_R, * M11_Ev_L2_All_hits_R, * M11_Ev_L3_All_hits_R, * M11_Ev_L4_All_hits_R, * M11_Ev_L5_All_hits_R, * M11_Ev_L6_All_hits_R; 


  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_etapart;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * ME0_el_hits_E,  * ME0_mu_hits_E,  * ME0_pi_hits_E,  * ME0_ka_hits_E,  * ME0_p_hits_E,  * ME0_n_hits_E,  * ME0_g_hits_E,  * ME0_N_hits_E, * ME0_OH_hits_E,  * ME0_All_hits_E,  * ME0_HIP_hits_E;
  TH1F * G11_el_hits_E,  * G11_mu_hits_E,  * G11_pi_hits_E,  * G11_ka_hits_E,  * G11_p_hits_E,  * G11_n_hits_E,  * G11_g_hits_E,  * G11_N_hits_E, * G11_OH_hits_E,  * G11_All_hits_E,  * G11_HIP_hits_E;
  TH1F * G21_el_hits_E,  * G21_mu_hits_E,  * G21_pi_hits_E,  * G21_ka_hits_E,  * G21_p_hits_E,  * G21_n_hits_E,  * G21_g_hits_E,  * G21_N_hits_E, * G21_OH_hits_E,  * G21_All_hits_E,  * G21_HIP_hits_E;

};

//
// constants, enums and typedefs
//

int n_tof = 1100,  n1_tof = 1,  n2_tof = 12;
int m_tof = 110,   m1_tof = 1,  m2_tof = 12;
int m_eta = 50;   double m1_eta =  0.0,  m2_eta = 2.5;
int m_phi = 36;   double m1_phi = -3.1415, m2_phi = 3.1415;
int m_lin = 1000, m1_lin = 0,  m2_lin = 100000000;
int n_hits = 20;  double n1_hits = 0.5, n2_hits = 20.5;
int n_lay  = 12;  double n1_lay  = 0.5, n2_lay  = 12.5;

// logaritmic bins ...
int n_D   = 900,  n1_D   = -3, n2_D   = 6;
int n_E   = 900,  n1_E   = -4, n2_E   = 5;
int n_F   = 90,   n1_F   = -3, n2_F   = 6;

// linear bins
int n_FL = 2500,  n1_FL   = 0, n2FL = 250000; // 1-250keV in bins of 100 eV

int n_xy_x =  800; double n_xy_x1 = -800;  double n_xy_x2 = +800;
int n_xy_y =  800; double n_xy_y1 = -800;  double n_xy_y2 = +800;
int n_zr_z =  550; double n_zr_z1 =    0;  double n_zr_z2 = 1100;
int n_zr_r =  400; double n_zr_r1 =    0;  double n_zr_r2 =  800;

int n_dz = 100; double n1_dz = 0.0; double n2_dz = 2.0;  double n2_dz2 = n2_dz*1.0/20;
int n_dR = 100; double n1_dR = 0.0; double n2_dR = 20.0; double n2_dR2 = n2_dR*1.0/20;
int n_dG = 500; double n1_dG = 0.0; double n2_dG = 1.0; 

int n_zrc_z = 1350; double n_zrc_z1 = 0; double n_zrc_z2 = 2700;
int n_zrc_r =  650; double n_zrc_r1 = 0; double n_zrc_r2 = 1300;

int n_pv_z = 50, n1_pv_z = -25, n2_pv_z = 25;
int n_pv_r = 10, n1_pv_r = 0,   n2_pv_r = 2;
int n_cat  = 27; double n1_cat  = 0.5, n2_cat  = 27.5;
int n_pro  = 15; double n1_pro  = 0.5, n2_pro  = 15.5;

std::string allcat[27] = {"All", "Barrel", "W0_RB1", "W0_RB2", "W0_RB3", "W0_RB4", "W1_RB1", "W1_RB2", "W1_RB3", "W1_RB4", "W2_RB1", "W2_RB2", "W2_RB3", "W2_RB4", 
                                 "Endcap", "RE11", "RE12", "RE13", "RE21", "RE22", "RE23", "RE31", "RE32", "RE33", "RE41", "RE42", "RE43"};
std::string endcat[24] = {"RE12C", "RE12B", "RE12A", "RE13C", "RE13B", "RE13A", "RE22C", "RE22B", "RE22A", "RE23C", "RE23B", "RE23A", 
                          "RE32C", "RE32B", "RE32A", "RE33C", "RE33B", "RE33A", "RE42C", "RE42B", "RE42A", "RE43C", "RE43B", "RE43A"};
std::string proc[15] = {"notDefined", "fCoulombScattering", "fIonisation", "fBremsstrahlung", "fPairProdByCharged", "fAnnihilation", "fAnnihilationToMuMu", "fAnnihilationToHadrons", 
                        "fNuclearStopping", "notDefined", "fMultipleScattering", "fRayleigh", "fPhotoElectricEffect", "fComptonScattering", "fGammaConversion"};
std::string nuclei[21] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Other"};

int r_CSC = 75, r1_CSC = 0, r2_CSC = 750; // we should do this with a precision of 5 cm ...
int n_layers_csc = 6;

//
// static data member definitions
//

//
// constructors and destructor
//
MyME11SimHitAnalyzer::MyME11SimHitAnalyzer(const edm::ParameterSet& iConfig)

{
  // now do what ever initialization is needed
  pdfFileNameBase = iConfig.getUntrackedParameter<std::string>("PdfFileNameBase");
  rootFileName    = iConfig.getUntrackedParameter<std::string>("RootFileName");
  bunchspacing    = iConfig.getUntrackedParameter<double>("BunchSpacing");
  comenergy       = iConfig.getUntrackedParameter<double>("COMEnergy");
  maxsimtime      = iConfig.getUntrackedParameter<double>("MaxSimTime");
  varEDepCuteV    = iConfig.getUntrackedParameter<double>("VarEDepCuteV");
  phys_debug      = iConfig.getUntrackedParameter<bool>("PhysicsDebug");
  tech_debug      = iConfig.getUntrackedParameter<bool>("TechnicDebug");

  outputfile      = new TFile(rootFileName.c_str(), "RECREATE" );

  // required for 7XY: registration of the data access
  consumesMany<edm::PSimHitContainer>();
  SIMTrack_Token  = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
  SIMVertex_Token = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits")); // or consumes<std::vector<SimVertex>>

  cscGeom_Token = esConsumes<CSCGeometry, MuonGeometryRecord>();
  dtGeom_Token  = esConsumes<DTGeometry,  MuonGeometryRecord>();
  gemGeom_Token = esConsumes<GEMGeometry, MuonGeometryRecord>();
  rpcGeom_Token = esConsumes<RPCGeometry, MuonGeometryRecord>();

  if(tech_debug) std::cout<<"[MyME11SimHitAnalyzer :: Constructor]"<<std::endl;
  
  TDir_Muon_hits_deposits = new TDirectoryFile("Muon_hits_deposits", "Muon_hits_deposits");
  TDir_Muon_hits_radius  = new TDirectoryFile("Muon_hits_radius",  "Muon_hits_radius");
  TDir_Muon_hits_etapart = new TDirectoryFile("Muon_hits_etapart", "Muon_hits_etapart");
  TDir_Muon_hits_Nuclei  = new TDirectoryFile("Muon_hits_Nuclei", "Muon_hits_Nuclei");



  ME0_el_hits_R  = new TH1F("ME0_el_hits_R",  "Simhit Radius :: ME0 :: Electrons",       90, 60, 150);
  ME0_mu_hits_R  = new TH1F("ME0_mu_hits_R",  "Simhit Radius :: ME0 :: Muons",           90, 60, 150);
  ME0_pi_hits_R  = new TH1F("ME0_pi_hits_R",  "Simhit Radius :: ME0 :: Pions",           90, 60, 150);
  ME0_ka_hits_R  = new TH1F("ME0_ka_hits_R",  "Simhit Radius :: ME0 :: Kaons",           90, 60, 150);
  ME0_p_hits_R   = new TH1F("ME0_p_hits_R",   "Simhit Radius :: ME0 :: Protons",         90, 60, 150);
  ME0_n_hits_R   = new TH1F("ME0_n_hits_R",   "Simhit Radius :: ME0 :: Neutrons",        90, 60, 150);
  ME0_g_hits_R   = new TH1F("ME0_g_hits_R",   "Simhit Radius :: ME0 :: Photons",         90, 60, 150);
  ME0_N_hits_R   = new TH1F("ME0_N_hits_R",   "Simhit Radius :: ME0 :: Nuclei",          90, 60, 150);
  ME0_OH_hits_R  = new TH1F("ME0_OH_hits_R",  "Simhit Radius :: ME0 :: Other Hadrons",   90, 60, 150);
  ME0_All_hits_R = new TH1F("ME0_All_hits_R", "Simhit Radius :: ME0 :: All Particles",   90, 60, 150);
  ME0_HIP_hits_R = new TH1F("ME0_HIP_hits_R", "Simhit Radius :: ME0 :: Highly Ionizing", 90, 60, 150); 

  G11_el_hits_R  = new TH1F("G11_el_hits_R",  "Simhit Radius :: GE11 :: Electrons",       121, 130, 251);
  G11_mu_hits_R  = new TH1F("G11_mu_hits_R",  "Simhit Radius :: GE11 :: Muons",           121, 130, 251);
  G11_pi_hits_R  = new TH1F("G11_pi_hits_R",  "Simhit Radius :: GE11 :: Pions",           121, 130, 251);
  G11_ka_hits_R  = new TH1F("G11_ka_hits_R",  "Simhit Radius :: GE11 :: Kaons",           121, 130, 251);
  G11_p_hits_R   = new TH1F("G11_p_hits_R",   "Simhit Radius :: GE11 :: Protons",         121, 130, 251);
  G11_n_hits_R   = new TH1F("G11_n_hits_R",   "Simhit Radius :: GE11 :: Neutrons",        121, 130, 251);
  G11_g_hits_R   = new TH1F("G11_g_hits_R",   "Simhit Radius :: GE11 :: Photons",         121, 130, 251);
  G11_N_hits_R   = new TH1F("G11_N_hits_R",   "Simhit Radius :: GE11 :: Nuclei",          121, 130, 251);
  G11_OH_hits_R  = new TH1F("G11_OH_hits_R",  "Simhit Radius :: GE11 :: Other Hadrons",   121, 130, 251);
  G11_All_hits_R = new TH1F("G11_All_hits_R", "Simhit Radius :: GE11 :: All Particles",   121, 130, 251);
  G11_HIP_hits_R = new TH1F("G11_HIP_hits_R", "Simhit Radius :: GE11 :: Highly Ionizing", 121, 130, 251); 

  G21_el_hits_R  = new TH1F("G21_el_hits_R",  "Simhit Radius :: GE21 :: Electrons",       184, 136, 320);
  G21_mu_hits_R  = new TH1F("G21_mu_hits_R",  "Simhit Radius :: GE21 :: Muons",           184, 136, 320);
  G21_pi_hits_R  = new TH1F("G21_pi_hits_R",  "Simhit Radius :: GE21 :: Pions",           184, 136, 320);
  G21_ka_hits_R  = new TH1F("G21_ka_hits_R",  "Simhit Radius :: GE21 :: Kaons",           184, 136, 320);
  G21_p_hits_R   = new TH1F("G21_p_hits_R",   "Simhit Radius :: GE21 :: Protons",         184, 136, 320);
  G21_n_hits_R   = new TH1F("G21_n_hits_R",   "Simhit Radius :: GE21 :: Neutrons",        184, 136, 320);
  G21_g_hits_R   = new TH1F("G21_g_hits_R",   "Simhit Radius :: GE21 :: Photons",         184, 136, 320);
  G21_N_hits_R   = new TH1F("G21_N_hits_R",   "Simhit Radius :: GE21 :: Nuclei",          184, 136, 320);
  G21_OH_hits_R  = new TH1F("G21_OH_hits_R",  "Simhit Radius :: GE21 :: Other Hadrons",   184, 136, 320);
  G21_All_hits_R = new TH1F("G21_All_hits_R", "Simhit Radius :: GE21 :: All Particles",   184, 136, 320);
  G21_HIP_hits_R = new TH1F("G21_HIP_hits_R", "Simhit Radius :: GE21 :: Highly Ionizing", 184, 136, 320); 

  G11_L1_el_hits_R  = new TH1F("G11_L1_el_hits_R",  "Simhit Radius :: G1E1 L1 :: Electrons",       121, 130, 251);
  G11_L1_mu_hits_R  = new TH1F("G11_L1_mu_hits_R",  "Simhit Radius :: GE11 L1 :: Muons",           121, 130, 251);
  G11_L1_pi_hits_R  = new TH1F("G11_L1_pi_hits_R",  "Simhit Radius :: GE11 L1 :: Pions",           121, 130, 251);
  G11_L1_ka_hits_R  = new TH1F("G11_L1_ka_hits_R",  "Simhit Radius :: GE11 L1 :: Kaons",           121, 130, 251);
  G11_L1_p_hits_R   = new TH1F("G11_L1_p_hits_R",   "Simhit Radius :: GE11 L1 :: Protons",         121, 130, 251);
  G11_L1_n_hits_R   = new TH1F("G11_L1_n_hits_R",   "Simhit Radius :: GE11 L1 :: Neutrons",        121, 130, 251);
  G11_L1_g_hits_R   = new TH1F("G11_L1_g_hits_R",   "Simhit Radius :: GE11 L1 :: Photons",         121, 130, 251);
  G11_L1_N_hits_R   = new TH1F("G11_L1_N_hits_R",   "Simhit Radius :: GE11 L1 :: Nuclei",          121, 130, 251);
  G11_L1_OH_hits_R  = new TH1F("G11_L1_OH_hits_R",  "Simhit Radius :: GE11 L1 :: Other Hadrons",   121, 130, 251);
  G11_L1_All_hits_R = new TH1F("G11_L1_All_hits_R", "Simhit Radius :: GE11 L1 :: All Particles",   121, 130, 251);
  G11_L1_HIP_hits_R = new TH1F("G11_L1_HIP_hits_R", "Simhit Radius :: GE11 L1 :: Highly Ionizing", 121, 130, 251); 

  G11_L2_el_hits_R  = new TH1F("G11_L2_el_hits_R",  "Simhit Radius :: GE11 L2 :: Electrons",       121, 130, 251);
  G11_L2_mu_hits_R  = new TH1F("G11_L2_mu_hits_R",  "Simhit Radius :: GE11 L2 :: Muons",           121, 130, 251);
  G11_L2_pi_hits_R  = new TH1F("G11_L2_pi_hits_R",  "Simhit Radius :: GE11 L2 :: Pions",           121, 130, 251);
  G11_L2_ka_hits_R  = new TH1F("G11_L2_ka_hits_R",  "Simhit Radius :: GE11 L2 :: Kaons",           121, 130, 251);
  G11_L2_p_hits_R   = new TH1F("G11_L2_p_hits_R",   "Simhit Radius :: GE11 L2 :: Protons",         121, 130, 251);
  G11_L2_n_hits_R   = new TH1F("G11_L2_n_hits_R",   "Simhit Radius :: GE11 L2 :: Neutrons",        121, 130, 251);
  G11_L2_g_hits_R   = new TH1F("G11_L2_g_hits_R",   "Simhit Radius :: GE11 L2 :: Photons",         121, 130, 251);
  G11_L2_N_hits_R   = new TH1F("G11_L2_N_hits_R",   "Simhit Radius :: GE11 L2 :: Nuclei",          121, 130, 251);
  G11_L2_OH_hits_R  = new TH1F("G11_L2_OH_hits_R",  "Simhit Radius :: GE11 L2 :: Other Hadrons",   121, 130, 251);
  G11_L2_All_hits_R = new TH1F("G11_L2_All_hits_R", "Simhit Radius :: GE11 L2 :: All Particles",   121, 130, 251);
  G11_L2_HIP_hits_R = new TH1F("G11_L2_HIP_hits_R", "Simhit Radius :: GE11 L2 :: Highly Ionizing", 121, 130, 251); 

  G11_Od_el_hits_R  = new TH1F("G11_Od_el_hits_R",  "Simhit Radius :: GE11 Odd :: Electrons",       121, 130, 251);
  G11_Od_mu_hits_R  = new TH1F("G11_Od_mu_hits_R",  "Simhit Radius :: GE11 Odd :: Muons",           121, 130, 251);
  G11_Od_pi_hits_R  = new TH1F("G11_Od_pi_hits_R",  "Simhit Radius :: GE11 Odd :: Pions",           121, 130, 251);
  G11_Od_ka_hits_R  = new TH1F("G11_Od_ka_hits_R",  "Simhit Radius :: GE11 Odd :: Kaons",           121, 130, 251);
  G11_Od_p_hits_R   = new TH1F("G11_Od_p_hits_R",   "Simhit Radius :: GE11 Odd :: Protons",         121, 130, 251);
  G11_Od_n_hits_R   = new TH1F("G11_Od_n_hits_R",   "Simhit Radius :: GE11 Odd :: Neutrons",        121, 130, 251);
  G11_Od_g_hits_R   = new TH1F("G11_Od_g_hits_R",   "Simhit Radius :: GE11 Odd :: Photons",         121, 130, 251);
  G11_Od_N_hits_R   = new TH1F("G11_Od_N_hits_R",   "Simhit Radius :: GE11 Odd :: Nuclei",          121, 130, 251);
  G11_Od_OH_hits_R  = new TH1F("G11_Od_OH_hits_R",  "Simhit Radius :: GE11 Odd :: Other Hadrons",   121, 130, 251);
  G11_Od_All_hits_R = new TH1F("G11_Od_All_hits_R", "Simhit Radius :: GE11 Odd :: All Particles",   121, 130, 251);
  G11_Od_HIP_hits_R = new TH1F("G11_Od_HIP_hits_R", "Simhit Radius :: GE11 Odd :: Highly Ionizing", 121, 130, 251); 

  G11_Ev_el_hits_R  = new TH1F("G11_Ev_el_hits_R",  "Simhit Radius :: GE11 Even :: Electrons",       121, 130, 251);
  G11_Ev_mu_hits_R  = new TH1F("G11_Ev_mu_hits_R",  "Simhit Radius :: GE11 Even :: Muons",           121, 130, 251);
  G11_Ev_pi_hits_R  = new TH1F("G11_Ev_pi_hits_R",  "Simhit Radius :: GE11 Even :: Pions",           121, 130, 251);
  G11_Ev_ka_hits_R  = new TH1F("G11_Ev_ka_hits_R",  "Simhit Radius :: GE11 Even :: Kaons",           121, 130, 251);
  G11_Ev_p_hits_R   = new TH1F("G11_Ev_p_hits_R",   "Simhit Radius :: GE11 Even :: Protons",         121, 130, 251);
  G11_Ev_n_hits_R   = new TH1F("G11_Ev_n_hits_R",   "Simhit Radius :: GE11 Even :: Neutrons",        121, 130, 251);
  G11_Ev_g_hits_R   = new TH1F("G11_Ev_g_hits_R",   "Simhit Radius :: GE11 Even :: Photons",         121, 130, 251);
  G11_Ev_N_hits_R   = new TH1F("G11_Ev_N_hits_R",   "Simhit Radius :: GE11 Even :: Nuclei",          121, 130, 251);
  G11_Ev_OH_hits_R  = new TH1F("G11_Ev_OH_hits_R",  "Simhit Radius :: GE11 Even :: Other Hadrons",   121, 130, 251);
  G11_Ev_All_hits_R = new TH1F("G11_Ev_All_hits_R", "Simhit Radius :: GE11 Even :: All Particles",   121, 130, 251);
  G11_Ev_HIP_hits_R = new TH1F("G11_Ev_HIP_hits_R", "Simhit Radius :: GE11 Even :: Highly Ionizing", 121, 130, 251); 

  M11_Od_el_hits_R  = new TH1F("M11_Od_el_hits_R",  "Simhit Radius :: ME11 Odd :: Electrons",       121, 130, 251);
  M11_Od_mu_hits_R  = new TH1F("M11_Od_mu_hits_R",  "Simhit Radius :: ME11 Odd :: Muons",           121, 130, 251);
  M11_Od_pi_hits_R  = new TH1F("M11_Od_pi_hits_R",  "Simhit Radius :: ME11 Odd :: Pions",           121, 130, 251);
  M11_Od_ka_hits_R  = new TH1F("M11_Od_ka_hits_R",  "Simhit Radius :: ME11 Odd :: Kaons",           121, 130, 251);
  M11_Od_p_hits_R   = new TH1F("M11_Od_p_hits_R",   "Simhit Radius :: ME11 Odd :: Protons",         121, 130, 251);
  M11_Od_n_hits_R   = new TH1F("M11_Od_n_hits_R",   "Simhit Radius :: ME11 Odd :: Neutrons",        121, 130, 251);
  M11_Od_g_hits_R   = new TH1F("M11_Od_g_hits_R",   "Simhit Radius :: ME11 Odd :: Photons",         121, 130, 251);
  M11_Od_N_hits_R   = new TH1F("M11_Od_N_hits_R",   "Simhit Radius :: ME11 Odd :: Nuclei",          121, 130, 251);
  M11_Od_OH_hits_R  = new TH1F("M11_Od_OH_hits_R",  "Simhit Radius :: ME11 Odd :: Other Hadrons",   121, 130, 251);
  M11_Od_All_hits_R = new TH1F("M11_Od_All_hits_R", "Simhit Radius :: ME11 Odd :: All Particles",   121, 130, 251);
  M11_Od_HIP_hits_R = new TH1F("M11_Od_HIP_hits_R", "Simhit Radius :: ME11 Odd :: Highly Ionizing", 121, 130, 251); 

  M11_Ev_el_hits_R  = new TH1F("M11_Ev_el_hits_R",  "Simhit Radius :: ME11 Even :: Electrons",       121, 130, 251);
  M11_Ev_mu_hits_R  = new TH1F("M11_Ev_mu_hits_R",  "Simhit Radius :: ME11 Even :: Muons",           121, 130, 251);
  M11_Ev_pi_hits_R  = new TH1F("M11_Ev_pi_hits_R",  "Simhit Radius :: ME11 Even :: Pions",           121, 130, 251);
  M11_Ev_ka_hits_R  = new TH1F("M11_Ev_ka_hits_R",  "Simhit Radius :: ME11 Even :: Kaons",           121, 130, 251);
  M11_Ev_p_hits_R   = new TH1F("M11_Ev_p_hits_R",   "Simhit Radius :: ME11 Even :: Protons",         121, 130, 251);
  M11_Ev_n_hits_R   = new TH1F("M11_Ev_n_hits_R",   "Simhit Radius :: ME11 Even :: Neutrons",        121, 130, 251);
  M11_Ev_g_hits_R   = new TH1F("M11_Ev_g_hits_R",   "Simhit Radius :: ME11 Even :: Photons",         121, 130, 251);
  M11_Ev_N_hits_R   = new TH1F("M11_Ev_N_hits_R",   "Simhit Radius :: ME11 Even :: Nuclei",          121, 130, 251);
  M11_Ev_OH_hits_R  = new TH1F("M11_Ev_OH_hits_R",  "Simhit Radius :: ME11 Even :: Other Hadrons",   121, 130, 251);
  M11_Ev_All_hits_R = new TH1F("M11_Ev_All_hits_R", "Simhit Radius :: ME11 Even :: All Particles",   121, 130, 251);
  M11_Ev_HIP_hits_R = new TH1F("M11_Ev_HIP_hits_R", "Simhit Radius :: ME11 Even :: Highly Ionizing", 121, 130, 251); 






  ME0_el_hits_E  = new TH1F("ME0_el_hits_E",  "Simhit EtaPart :: ME0 :: Electrons",       8, 0.5, 8.5);
  ME0_mu_hits_E  = new TH1F("ME0_mu_hits_E",  "Simhit EtaPart :: ME0 :: Muons",           8, 0.5, 8.5);
  ME0_pi_hits_E  = new TH1F("ME0_pi_hits_E",  "Simhit EtaPart :: ME0 :: Pions",           8, 0.5, 8.5);
  ME0_ka_hits_E  = new TH1F("ME0_ka_hits_E",  "Simhit EtaPart :: ME0 :: Kaons",           8, 0.5, 8.5);
  ME0_p_hits_E   = new TH1F("ME0_p_hits_E",   "Simhit EtaPart :: ME0 :: Protons",         8, 0.5, 8.5);
  ME0_n_hits_E   = new TH1F("ME0_n_hits_E",   "Simhit EtaPart :: ME0 :: Neutrons",        8, 0.5, 8.5);
  ME0_g_hits_E   = new TH1F("ME0_g_hits_E",   "Simhit EtaPart :: ME0 :: Photons",         8, 0.5, 8.5);
  ME0_N_hits_E   = new TH1F("ME0_N_hits_E",   "Simhit EtaPart :: ME0 :: Nuclei",          8, 0.5, 8.5);
  ME0_OH_hits_E  = new TH1F("ME0_OH_hits_E",  "Simhit EtaPart :: ME0 :: Other Hadrons",   8, 0.5, 8.5);
  ME0_All_hits_E = new TH1F("ME0_All_hits_E", "Simhit EtaPart :: ME0 :: All Particles",   8, 0.5, 8.5);
  ME0_HIP_hits_E = new TH1F("ME0_HIP_hits_E", "Simhit EtaPart :: ME0 :: Highly Ionizing", 8, 0.5, 8.5); 

  G11_el_hits_E  = new TH1F("G11_el_hits_E",  "Simhit EtaPart :: GE11 :: Electrons",       8, 0.5, 8.5);
  G11_mu_hits_E  = new TH1F("G11_mu_hits_E",  "Simhit EtaPart :: GE11 :: Muons",           8, 0.5, 8.5);
  G11_pi_hits_E  = new TH1F("G11_pi_hits_E",  "Simhit EtaPart :: GE11 :: Pions",           8, 0.5, 8.5);
  G11_ka_hits_E  = new TH1F("G11_ka_hits_E",  "Simhit EtaPart :: GE11 :: Kaons",           8, 0.5, 8.5);
  G11_p_hits_E   = new TH1F("G11_p_hits_E",   "Simhit EtaPart :: GE11 :: Protons",         8, 0.5, 8.5);
  G11_n_hits_E   = new TH1F("G11_n_hits_E",   "Simhit EtaPart :: GE11 :: Neutrons",        8, 0.5, 8.5);
  G11_g_hits_E   = new TH1F("G11_g_hits_E",   "Simhit EtaPart :: GE11 :: Photons",         8, 0.5, 8.5);
  G11_N_hits_E   = new TH1F("G11_N_hits_E",   "Simhit EtaPart :: GE11 :: Nuclei",          8, 0.5, 8.5);
  G11_OH_hits_E  = new TH1F("G11_OH_hits_E",  "Simhit EtaPart :: GE11 :: Other Hadrons",   8, 0.5, 8.5);
  G11_All_hits_E = new TH1F("G11_All_hits_E", "Simhit EtaPart :: GE11 :: All Particles",   8, 0.5, 8.5);
  G11_HIP_hits_E = new TH1F("G11_HIP_hits_E", "Simhit EtaPart :: GE11 :: Highly Ionizing", 8, 0.5, 8.5); 

  G21_el_hits_E  = new TH1F("G21_el_hits_E",  "Simhit EtaPart :: GE21 :: Electrons",       16, 0.5, 16.5);
  G21_mu_hits_E  = new TH1F("G21_mu_hits_E",  "Simhit EtaPart :: GE21 :: Muons",           16, 0.5, 16.5);
  G21_pi_hits_E  = new TH1F("G21_pi_hits_E",  "Simhit EtaPart :: GE21 :: Pions",           16, 0.5, 16.5);
  G21_ka_hits_E  = new TH1F("G21_ka_hits_E",  "Simhit EtaPart :: GE21 :: Kaons",           16, 0.5, 16.5);
  G21_p_hits_E   = new TH1F("G21_p_hits_E",   "Simhit EtaPart :: GE21 :: Protons",         16, 0.5, 16.5);
  G21_n_hits_E   = new TH1F("G21_n_hits_E",   "Simhit EtaPart :: GE21 :: Neutrons",        16, 0.5, 16.5);
  G21_g_hits_E   = new TH1F("G21_g_hits_E",   "Simhit EtaPart :: GE21 :: Photons",         16, 0.5, 16.5);
  G21_N_hits_E   = new TH1F("G21_N_hits_E",   "Simhit EtaPart :: GE21 :: Nuclei",          16, 0.5, 16.5);
  G21_OH_hits_E  = new TH1F("G21_OH_hits_E",  "Simhit EtaPart :: GE21 :: Other Hadrons",   16, 0.5, 16.5);
  G21_All_hits_E = new TH1F("G21_All_hits_E", "Simhit EtaPart :: GE21 :: All Particles",   16, 0.5, 16.5);
  G21_HIP_hits_E = new TH1F("G21_HIP_hits_E", "Simhit EtaPart :: GE21 :: Highly Ionizing", 16, 0.5, 16.5);

  Muon_Nuclei_A_Z = new TH2F("Muon_Nuclei_A_Z", "Nuclei A vs Z :: All",  100, 0.5, 100.5, 200, 0.5, 200.5); 
  RPCb_Nuclei_A_Z = new TH2F("RPCb_Nuclei_A_Z", "Nuclei A vs Z :: RPCb", 100, 0.5, 100.5, 200, 0.5, 200.5); 
  RPCf_Nuclei_A_Z = new TH2F("RPCf_Nuclei_A_Z", "Nuclei A vs Z :: RPCf", 100, 0.5, 100.5, 200, 0.5, 200.5); 
  CSC_Nuclei_A_Z  = new TH2F("CSC_Nuclei_A_Z",  "Nuclei A vs Z :: CSC",  100, 0.5, 100.5, 200, 0.5, 200.5); 
  DT_Nuclei_A_Z   = new TH2F("DT_Nuclei_A_Z",   "Nuclei A vs Z :: DT",   100, 0.5, 100.5, 200, 0.5, 200.5); 
  GEM_Nuclei_A_Z  = new TH2F("GEM_Nuclei_A_Z",  "Nuclei A vs Z :: GEM",  100, 0.5, 100.5, 200, 0.5, 200.5); 
  ME0_Nuclei_A_Z  = new TH2F("ME0_Nuclei_A_Z",  "Nuclei A vs Z :: ME0",  100, 0.5, 100.5, 200, 0.5, 200.5); 

  Muon_Nuclei_A = new TH1F("Muon_Nuclei_A", "Nuclei A :: All",  200, 0.5, 200.5);
  RPCb_Nuclei_A = new TH1F("RPCb_Nuclei_A", "Nuclei A :: RPCb", 200, 0.5, 200.5);
  RPCf_Nuclei_A = new TH1F("RPCf_Nuclei_A", "Nuclei A :: RPCf", 200, 0.5, 200.5);
  CSC_Nuclei_A  = new TH1F("CSC_Nuclei_A",  "Nuclei A :: CSC",  200, 0.5, 200.5);
  DT_Nuclei_A   = new TH1F("DT_Nuclei_A",   "Nuclei A :: DT",   200, 0.5, 200.5); 
  GEM_Nuclei_A  = new TH1F("GEM_Nuclei_A",  "Nuclei A :: GEM",  200, 0.5, 200.5); 
  ME0_Nuclei_A  = new TH1F("ME0_Nuclei_A",  "Nuclei A :: ME0",  200, 0.5, 200.5); 

  Muon_Nuclei_Z = new TH1F("Muon_Nuclei_Z", "Nuclei Z :: All",  100, 0.5, 100.5);
  RPCb_Nuclei_Z = new TH1F("RPCb_Nuclei_Z", "Nuclei Z :: RPCb", 100, 0.5, 100.5);
  RPCf_Nuclei_Z = new TH1F("RPCf_Nuclei_Z", "Nuclei Z :: RPCf", 100, 0.5, 100.5);
  CSC_Nuclei_Z  = new TH1F("CSC_Nuclei_Z",  "Nuclei Z :: CSC",  100, 0.5, 100.5);
  DT_Nuclei_Z   = new TH1F("DT_Nuclei_Z",   "Nuclei Z :: DT",   100, 0.5, 100.5); 
  GEM_Nuclei_Z  = new TH1F("GEM_Nuclei_Z",  "Nuclei Z :: GEM",  100, 0.5, 100.5); 
  ME0_Nuclei_Z  = new TH1F("ME0_Nuclei_Z",  "Nuclei Z :: ME0",  100, 0.5, 100.5); 

  std::vector<TH1F*> NucListVec; 
  Muon_Nuclei_List = new TH1F("Muon_Nuclei_List", "Nuclei List :: All",  21, 0.5, 21.5); NucListVec.push_back(Muon_Nuclei_List);
  RPCb_Nuclei_List = new TH1F("RPCb_Nuclei_List", "Nuclei List :: RPCb", 21, 0.5, 21.5); NucListVec.push_back(RPCb_Nuclei_List);
  RPCf_Nuclei_List = new TH1F("RPCf_Nuclei_List", "Nuclei List :: RPCf", 21, 0.5, 21.5); NucListVec.push_back(RPCf_Nuclei_List);
  CSC_Nuclei_List  = new TH1F("CSC_Nuclei_List",  "Nuclei List :: CSC",  21, 0.5, 21.5); NucListVec.push_back(CSC_Nuclei_List);
  DT_Nuclei_List   = new TH1F("DT_Nuclei_List",   "Nuclei List :: DT",   21, 0.5, 21.5); NucListVec.push_back(DT_Nuclei_List);
  GEM_Nuclei_List  = new TH1F("GEM_Nuclei_List",  "Nuclei List :: GEM",  21, 0.5, 21.5); NucListVec.push_back(GEM_Nuclei_List);
  ME0_Nuclei_List  = new TH1F("ME0_Nuclei_List",  "Nuclei List :: ME0",  21, 0.5, 21.5); NucListVec.push_back(ME0_Nuclei_List);
  
  for(int i=0; i<7; ++i) {
    for(int j=0; j<21; ++j) {
      NucListVec[i]->GetXaxis()->SetBinLabel(j+1,nuclei[j].c_str());
    }
  }


  // Simhit Time vs Ekin
  RPCb_el_hits = new TH2F("RPCb_el_hits", "Simhit time vs E_{kin} :: RPCb :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_mu_hits = new TH2F("RPCb_mu_hits", "Simhit time vs E_{kin} :: RPCb :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_pi_hits = new TH2F("RPCb_pi_hits", "Simhit time vs E_{kin} :: RPCb :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_ka_hits = new TH2F("RPCb_ka_hits", "Simhit time vs E_{kin} :: RPCb :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_p_hits  = new TH2F("RPCb_p_hits",  "Simhit time vs E_{kin} :: RPCb :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_n_hits  = new TH2F("RPCb_n_hits",  "Simhit time vs E_{kin} :: RPCb :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_g_hits  = new TH2F("RPCb_g_hits",  "Simhit time vs E_{kin} :: RPCb :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_N_hits  = new TH2F("RPCb_N_hits",  "Simhit time vs E_{kin} :: RPCb :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_OH_hits = new TH2F("RPCb_OH_hits", "Simhit time vs E_{kin} :: RPCb :: Other Hadrons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_All_hits= new TH2F("RPCb_All_hits","Simhit time vs E_{kin} :: RPCb :: All particles",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_HIP_hits= new TH2F("RPCb_HIP_hits","Simhit time vs E_{kin} :: RPCb :: Highly Ionizing", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);


  RPCf_el_hits = new TH2F("RPCf_el_hits", "Simhit time vs E_{kin} :: RPCf :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_mu_hits = new TH2F("RPCf_mu_hits", "Simhit time vs E_{kin} :: RPCf :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_pi_hits = new TH2F("RPCf_pi_hits", "Simhit time vs E_{kin} :: RPCf :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_ka_hits = new TH2F("RPCf_ka_hits", "Simhit time vs E_{kin} :: RPCf :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_p_hits  = new TH2F("RPCf_p_hits",  "Simhit time vs E_{kin} :: RPCf :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_n_hits  = new TH2F("RPCf_n_hits",  "Simhit time vs E_{kin} :: RPCf :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_g_hits  = new TH2F("RPCf_g_hits",  "Simhit time vs E_{kin} :: RPCf :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_N_hits  = new TH2F("RPCf_N_hits",  "Simhit time vs E_{kin} :: RPCf :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_OH_hits = new TH2F("RPCf_OH_hits", "Simhit time vs E_{kin} :: RPCf :: Other Hadrons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_All_hits= new TH2F("RPCf_All_hits","Simhit time vs E_{kin} :: RPCf :: All particles",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_HIP_hits= new TH2F("RPCf_HIP_hits","Simhit time vs E_{kin} :: RPCf :: Highly Ionizing", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  CSC_el_hits = new TH2F("CSC_el_hits", "Simhit time vs E_{kin} :: CSC :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_mu_hits = new TH2F("CSC_mu_hits", "Simhit time vs E_{kin} :: CSC :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_pi_hits = new TH2F("CSC_pi_hits", "Simhit time vs E_{kin} :: CSC :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_ka_hits = new TH2F("CSC_ka_hits", "Simhit time vs E_{kin} :: CSC :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_p_hits  = new TH2F("CSC_p_hits",  "Simhit time vs E_{kin} :: CSC :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_n_hits  = new TH2F("CSC_n_hits",  "Simhit time vs E_{kin} :: CSC :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_g_hits  = new TH2F("CSC_g_hits",  "Simhit time vs E_{kin} :: CSC :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_N_hits  = new TH2F("CSC_N_hits",  "Simhit time vs E_{kin} :: CSC :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_OH_hits = new TH2F("CSC_OH_hits", "Simhit time vs E_{kin} :: CSC :: Other Hadrons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_All_hits= new TH2F("CSC_All_hits","Simhit time vs E_{kin} :: CSC :: All particles",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_HIP_hits= new TH2F("CSC_HIP_hits","Simhit time vs E_{kin} :: CSC :: Highly Ionizing", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  DT_el_hits = new TH2F("DT_el_hits", "Simhit time vs E_{kin} :: DT :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_mu_hits = new TH2F("DT_mu_hits", "Simhit time vs E_{kin} :: DT :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_pi_hits = new TH2F("DT_pi_hits", "Simhit time vs E_{kin} :: DT :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_ka_hits = new TH2F("DT_ka_hits", "Simhit time vs E_{kin} :: DT :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_p_hits  = new TH2F("DT_p_hits",  "Simhit time vs E_{kin} :: DT :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_n_hits  = new TH2F("DT_n_hits",  "Simhit time vs E_{kin} :: DT :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_g_hits  = new TH2F("DT_g_hits",  "Simhit time vs E_{kin} :: DT :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_N_hits  = new TH2F("DT_N_hits",  "Simhit time vs E_{kin} :: DT :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_OH_hits = new TH2F("DT_OH_hits", "Simhit time vs E_{kin} :: DT :: Other Hadrons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_All_hits= new TH2F("DT_All_hits","Simhit time vs E_{kin} :: DT :: All particles",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_HIP_hits= new TH2F("DT_HIP_hits","Simhit time vs E_{kin} :: DT :: Highly Ionizing", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  GEM_el_hits = new TH2F("GEM_el_hits", "Simhit time vs E_{kin} :: GEM :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_mu_hits = new TH2F("GEM_mu_hits", "Simhit time vs E_{kin} :: GEM :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_pi_hits = new TH2F("GEM_pi_hits", "Simhit time vs E_{kin} :: GEM :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_ka_hits = new TH2F("GEM_ka_hits", "Simhit time vs E_{kin} :: GEM :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_p_hits  = new TH2F("GEM_p_hits",  "Simhit time vs E_{kin} :: GEM :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_n_hits  = new TH2F("GEM_n_hits",  "Simhit time vs E_{kin} :: GEM :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_g_hits  = new TH2F("GEM_g_hits",  "Simhit time vs E_{kin} :: GEM :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_N_hits  = new TH2F("GEM_N_hits",  "Simhit time vs E_{kin} :: GEM :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_OH_hits = new TH2F("GEM_OH_hits", "Simhit time vs E_{kin} :: GEM :: Other Hadrons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_All_hits= new TH2F("GEM_All_hits","Simhit time vs E_{kin} :: GEM :: All particles",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_HIP_hits= new TH2F("GEM_HIP_hits","Simhit time vs E_{kin} :: GEM :: Highly Ionizing", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  ME0_el_hits = new TH2F("ME0_el_hits", "Simhit time vs E_{kin} :: ME0 :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_mu_hits = new TH2F("ME0_mu_hits", "Simhit time vs E_{kin} :: ME0 :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_pi_hits = new TH2F("ME0_pi_hits", "Simhit time vs E_{kin} :: ME0 :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_ka_hits = new TH2F("ME0_ka_hits", "Simhit time vs E_{kin} :: ME0 :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_p_hits  = new TH2F("ME0_p_hits",  "Simhit time vs E_{kin} :: ME0 :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_n_hits  = new TH2F("ME0_n_hits",  "Simhit time vs E_{kin} :: ME0 :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_g_hits  = new TH2F("ME0_g_hits",  "Simhit time vs E_{kin} :: ME0 :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_N_hits  = new TH2F("ME0_N_hits",  "Simhit time vs E_{kin} :: ME0 :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_OH_hits = new TH2F("ME0_OH_hits", "Simhit time vs E_{kin} :: ME0 :: Other Hadrons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_All_hits= new TH2F("ME0_All_hits","Simhit time vs E_{kin} :: ME0 :: All particles",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_HIP_hits= new TH2F("ME0_HIP_hits","Simhit time vs E_{kin} :: ME0 :: Highly Ionizing", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  RPCb_el_kins = new TH1F("RPCb_el_kins", "E_{deposit} :: RPCb :: Electrons", n_E, n1_E, n2_E);
  RPCb_mu_kins = new TH1F("RPCb_mu_kins", "E_{deposit} :: RPCb :: Muons",     n_E, n1_E, n2_E);
  RPCb_ha_kins = new TH1F("RPCb_ha_kins", "E_{deposit} :: RPCb :: Hadrons",   n_E, n1_E, n2_E);
  RPCf_el_kins = new TH1F("RPCf_el_kins", "E_{deposit} :: RPCf :: Electrons", n_E, n1_E, n2_E);
  RPCf_mu_kins = new TH1F("RPCf_mu_kins", "E_{deposit} :: RPCf :: Muons",     n_E, n1_E, n2_E);
  RPCf_ha_kins = new TH1F("RPCf_ha_kins", "E_{deposit} :: RPCf :: Hadrons",   n_E, n1_E, n2_E);
  CSC_el_kins  = new TH1F("CSC_el_kins",  "E_{deposit} :: CSC :: Electrons",  n_E, n1_E, n2_E);
  CSC_mu_kins  = new TH1F("CSC_mu_kins",  "E_{deposit} :: CSC :: Muons",      n_E, n1_E, n2_E);
  CSC_ha_kins  = new TH1F("CSC_ha_kins",  "E_{deposit} :: CSC :: Hadrons",    n_E, n1_E, n2_E);
  DT_el_kins   = new TH1F("DT_el_kins",   "E_{deposit} :: DT :: Electrons",   n_E, n1_E, n2_E);
  DT_mu_kins   = new TH1F("DT_mu_kins",   "E_{deposit} :: DT :: Muons",       n_E, n1_E, n2_E);
  DT_ha_kins   = new TH1F("DT_ha_kins",   "E_{deposit} :: DT :: Hadrons",     n_E, n1_E, n2_E);
  GEM_el_kins  = new TH1F("GEM_el_kins",  "E_{deposit} :: GEM :: Electrons",  n_E, n1_E, n2_E);
  GEM_mu_kins  = new TH1F("GEM_mu_kins",  "E_{deposit} :: GEM :: Muons",      n_E, n1_E, n2_E);
  GEM_ha_kins  = new TH1F("GEM_ha_kins",  "E_{deposit} :: GEM :: Hadrons",    n_E, n1_E, n2_E);
  ME0_el_kins  = new TH1F("ME0_el_kins",  "E_{deposit} :: ME0 :: Electrons",  n_E, n1_E, n2_E);
  ME0_mu_kins  = new TH1F("ME0_mu_kins",  "E_{deposit} :: ME0 :: Muons",      n_E, n1_E, n2_E);
  ME0_ha_kins  = new TH1F("ME0_ha_kins",  "E_{deposit} :: ME0 :: Hadrons",    n_E, n1_E, n2_E);

  // Simhit Time vs E deposit
  RPCb_el_deposits = new TH2F("RPCb_el_deposits", "Simhit time vs E_{deposit} :: RPCb :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_mu_deposits = new TH2F("RPCb_mu_deposits", "Simhit time vs E_{deposit} :: RPCb :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_pi_deposits = new TH2F("RPCb_pi_deposits", "Simhit time vs E_{deposit} :: RPCb :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_ka_deposits = new TH2F("RPCb_ka_deposits", "Simhit time vs E_{deposit} :: RPCb :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_p_deposits  = new TH2F("RPCb_p_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_n_deposits  = new TH2F("RPCb_n_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_g_deposits  = new TH2F("RPCb_g_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_N_deposits  = new TH2F("RPCb_N_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_OH_deposits = new TH2F("RPCb_OH_deposits", "Simhit time vs E_{deposit} :: RPCb :: Other Hadrons",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_All_deposits= new TH2F("RPCb_All_deposits","Simhit time vs E_{deposit} :: RPCb :: All Particles",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_HIP_deposits= new TH2F("RPCb_HIP_deposits","Simhit time vs E_{deposit} :: RPCb :: Highly Ionizing",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  RPCf_el_deposits = new TH2F("RPCf_el_deposits", "Simhit time vs E_{deposit} :: RPCf :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_mu_deposits = new TH2F("RPCf_mu_deposits", "Simhit time vs E_{deposit} :: RPCf :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_pi_deposits = new TH2F("RPCf_pi_deposits", "Simhit time vs E_{deposit} :: RPCf :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_ka_deposits = new TH2F("RPCf_ka_deposits", "Simhit time vs E_{deposit} :: RPCf :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_p_deposits  = new TH2F("RPCf_p_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_n_deposits  = new TH2F("RPCf_n_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_g_deposits  = new TH2F("RPCf_g_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_N_deposits  = new TH2F("RPCf_N_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_OH_deposits = new TH2F("RPCf_OH_deposits", "Simhit time vs E_{deposit} :: RPCf :: Other Hadrons",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_All_deposits= new TH2F("RPCf_All_deposits","Simhit time vs E_{deposit} :: RPCf :: All Particles",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_HIP_deposits= new TH2F("RPCf_HIP_deposits","Simhit time vs E_{deposit} :: RPCf :: Highly Ionizing",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  CSC_el_deposits = new TH2F("CSC_el_deposits", "Simhit time vs E_{deposit} :: CSC :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_mu_deposits = new TH2F("CSC_mu_deposits", "Simhit time vs E_{deposit} :: CSC :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_pi_deposits = new TH2F("CSC_pi_deposits", "Simhit time vs E_{deposit} :: CSC :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_ka_deposits = new TH2F("CSC_ka_deposits", "Simhit time vs E_{deposit} :: CSC :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_p_deposits  = new TH2F("CSC_p_deposits",  "Simhit time vs E_{deposit} :: CSC :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_n_deposits  = new TH2F("CSC_n_deposits",  "Simhit time vs E_{deposit} :: CSC :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_g_deposits  = new TH2F("CSC_g_deposits",  "Simhit time vs E_{deposit} :: CSC :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_N_deposits  = new TH2F("CSC_N_deposits",  "Simhit time vs E_{deposit} :: CSC :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_OH_deposits = new TH2F("CSC_OH_deposits", "Simhit time vs E_{deposit} :: CSC :: Other Hadrons",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_All_deposits= new TH2F("CSC_All_deposits","Simhit time vs E_{deposit} :: CSC :: All Particles",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_HIP_deposits= new TH2F("CSC_HIP_deposits","Simhit time vs E_{deposit} :: CSC :: Highly Ionizing",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  DT_el_deposits = new TH2F("DT_el_deposits", "Simhit time vs E_{deposit} :: DT :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_mu_deposits = new TH2F("DT_mu_deposits", "Simhit time vs E_{deposit} :: DT :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_pi_deposits = new TH2F("DT_pi_deposits", "Simhit time vs E_{deposit} :: DT :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_ka_deposits = new TH2F("DT_ka_deposits", "Simhit time vs E_{deposit} :: DT :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_p_deposits  = new TH2F("DT_p_deposits",  "Simhit time vs E_{deposit} :: DT :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_n_deposits  = new TH2F("DT_n_deposits",  "Simhit time vs E_{deposit} :: DT :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_g_deposits  = new TH2F("DT_g_deposits",  "Simhit time vs E_{deposit} :: DT :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_N_deposits  = new TH2F("DT_N_deposits",  "Simhit time vs E_{deposit} :: DT :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_OH_deposits = new TH2F("DT_OH_deposits", "Simhit time vs E_{deposit} :: DT :: Other Hadrons",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_All_deposits= new TH2F("DT_All_deposits","Simhit time vs E_{deposit} :: DT :: All Particles",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_HIP_deposits= new TH2F("DT_HIP_deposits","Simhit time vs E_{deposit} :: DT :: Highly Ionizing",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  GEM_el_deposits = new TH2F("GEM_el_deposits", "Simhit time vs E_{deposit} :: GEM :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_mu_deposits = new TH2F("GEM_mu_deposits", "Simhit time vs E_{deposit} :: GEM :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_pi_deposits = new TH2F("GEM_pi_deposits", "Simhit time vs E_{deposit} :: GEM :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_ka_deposits = new TH2F("GEM_ka_deposits", "Simhit time vs E_{deposit} :: GEM :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_p_deposits  = new TH2F("GEM_p_deposits",  "Simhit time vs E_{deposit} :: GEM :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_n_deposits  = new TH2F("GEM_n_deposits",  "Simhit time vs E_{deposit} :: GEM :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_g_deposits  = new TH2F("GEM_g_deposits",  "Simhit time vs E_{deposit} :: GEM :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_N_deposits  = new TH2F("GEM_N_deposits",  "Simhit time vs E_{deposit} :: GEM :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_OH_deposits = new TH2F("GEM_OH_deposits", "Simhit time vs E_{deposit} :: GEM :: Other Hadrons",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_All_deposits= new TH2F("GEM_All_deposits","Simhit time vs E_{deposit} :: GEM :: All Particles",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_HIP_deposits= new TH2F("GEM_HIP_deposits","Simhit time vs E_{deposit} :: GEM :: Highly Ionizing",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  ME0_el_deposits = new TH2F("ME0_el_deposits", "Simhit time vs E_{deposit} :: ME0 :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_mu_deposits = new TH2F("ME0_mu_deposits", "Simhit time vs E_{deposit} :: ME0 :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_pi_deposits = new TH2F("ME0_pi_deposits", "Simhit time vs E_{deposit} :: ME0 :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_ka_deposits = new TH2F("ME0_ka_deposits", "Simhit time vs E_{deposit} :: ME0 :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_p_deposits  = new TH2F("ME0_p_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_n_deposits  = new TH2F("ME0_n_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_g_deposits  = new TH2F("ME0_g_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_N_deposits  = new TH2F("ME0_N_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_OH_deposits = new TH2F("ME0_OH_deposits", "Simhit time vs E_{deposit} :: ME0 :: Other Hadrons",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_All_deposits= new TH2F("ME0_All_deposits","Simhit time vs E_{deposit} :: ME0 :: All Particles",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_HIP_deposits= new TH2F("ME0_HIP_deposits","Simhit time vs E_{deposit} :: ME0 :: Highly Ionizing",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  // 1D Linear Energy Deposits
  // Probably not necessary --- is implemented in MyME0SimHitAnalyzer
  /*
  GEM_el_lindep = new TH1F("GEM_el_lindep", "E_{deposit} :: GEM :: Electrons", n_FL, n1_FL, n2_FL);
  GEM_mu_lindep = new TH1F("GEM_mu_lindep", "E_{deposit} :: GEM :: Muons",     n_FL, n1_FL, n2_FL);
  GEM_pi_lindep = new TH1F("GEM_pi_lindep", "E_{deposit} :: GEM :: Pions",     n_FL, n1_FL, n2_FL);
  GEM_ka_lindep = new TH1F("GEM_ka_lindep", "E_{deposit} :: GEM :: Kaons",     n_FL, n1_FL, n2_FL);
  GEM_p_lindep  = new TH1F("GEM_p_lindep",  "E_{deposit} :: GEM :: Protons",   n_FL, n1_FL, n2_FL);
  GEM_n_lindep  = new TH1F("GEM_n_lindep",  "E_{deposit} :: GEM :: Neutrons",  n_FL, n1_FL, n2_FL);
  GEM_g_lindep  = new TH1F("GEM_g_lindep",  "E_{deposit} :: GEM :: Photons",   n_FL, n1_FL, n2_FL);
  GEM_N_lindep  = new TH1F("GEM_N_lindep",  "E_{deposit} :: GEM :: Nuclei",    n_FL, n1_FL, n2_FL,);
  GEM_OH_lindep = new TH1F("GEM_OH_lindep", "E_{deposit} :: GEM :: Other Hadrons",    n_FL, n1_FL, n2_FL);
  GEM_All_lindep= new TH1F("GEM_All_lindep","E_{deposit} :: GEM :: All Particles",    n_FL, n1_FL, n2_FL);
  GEM_HIP_lindep= new TH1F("GEM_HIP_lindep","E_{deposit} :: GEM :: Highly Ionizing",  n_FL, n1_FL, n2_FL);

  ME0_el_lindep = new TH1F("ME0_el_lindep", "E_{deposit} :: ME0 :: Electrons", n_FL, n1_FL, n2_FL);
  ME0_mu_lindep = new TH1F("ME0_mu_lindep", "E_{deposit} :: ME0 :: Muons",     n_FL, n1_FL, n2_FL);
  ME0_pi_lindep = new TH1F("ME0_pi_lindep", "E_{deposit} :: ME0 :: Pions",     n_FL, n1_FL, n2_FL);
  ME0_ka_lindep = new TH1F("ME0_ka_lindep", "E_{deposit} :: ME0 :: Kaons",     n_FL, n1_FL, n2_FL);
  ME0_p_lindep  = new TH1F("ME0_p_lindep",  "E_{deposit} :: ME0 :: Protons",   n_FL, n1_FL, n2_FL);
  ME0_n_lindep  = new TH1F("ME0_n_lindep",  "E_{deposit} :: ME0 :: Neutrons",  n_FL, n1_FL, n2_FL);
  ME0_g_lindep  = new TH1F("ME0_g_lindep",  "E_{deposit} :: ME0 :: Photons",   n_FL, n1_FL, n2_FL);
  ME0_N_lindep  = new TH1F("ME0_N_lindep",  "E_{deposit} :: ME0 :: Nuclei",    n_FL, n1_FL, n2_FL,);
  ME0_OH_lindep = new TH1F("ME0_OH_lindep", "E_{deposit} :: ME0 :: Other Hadrons",    n_FL, n1_FL, n2_FL);
  ME0_All_lindep= new TH1F("ME0_All_lindep","E_{deposit} :: ME0 :: All Particles",    n_FL, n1_FL, n2_FL);
  ME0_HIP_lindep= new TH1F("ME0_HIP_lindep","E_{deposit} :: ME0 :: Highly Ionizing",  n_FL, n1_FL, n2_FL);
  */ 

  // Old 1D Logaritmic Energy Deposits
  RPCb_el_deps = new TH1F("RPCb_el_deps", "E_{deposit} :: RPCb :: Electrons", n_F, n1_F, n2_F);
  RPCb_mu_deps = new TH1F("RPCb_mu_deps", "E_{deposit} :: RPCb :: Muons",     n_F, n1_F, n2_F);
  RPCb_ha_deps = new TH1F("RPCb_ha_deps", "E_{deposit} :: RPCb :: Hadrons",   n_F, n1_F, n2_F);
  RPCf_el_deps = new TH1F("RPCf_el_deps", "E_{deposit} :: RPCf :: Electrons", n_F, n1_F, n2_F);
  RPCf_mu_deps = new TH1F("RPCf_mu_deps", "E_{deposit} :: RPCf :: Muons",     n_F, n1_F, n2_F);
  RPCf_ha_deps = new TH1F("RPCf_ha_deps", "E_{deposit} :: RPCf :: Hadrons",   n_F, n1_F, n2_F);
  CSC_el_deps  = new TH1F("CSC_el_deps",  "E_{deposit} :: CSC :: Electrons",  n_F, n1_F, n2_F);
  CSC_mu_deps  = new TH1F("CSC_mu_deps",  "E_{deposit} :: CSC :: Muons",      n_F, n1_F, n2_F);
  CSC_ha_deps  = new TH1F("CSC_ha_deps",  "E_{deposit} :: CSC :: Hadrons",    n_F, n1_F, n2_F);
  DT_el_deps   = new TH1F("DT_el_deps",   "E_{deposit} :: DT :: Electrons",   n_F, n1_F, n2_F);
  DT_mu_deps   = new TH1F("DT_mu_deps",   "E_{deposit} :: DT :: Muons",       n_F, n1_F, n2_F);
  DT_ha_deps   = new TH1F("DT_ha_deps",   "E_{deposit} :: DT :: Hadrons",     n_F, n1_F, n2_F);
  GEM_el_deps  = new TH1F("GEM_el_deps",  "E_{deposit} :: GEM :: Electrons",  n_F, n1_F, n2_F);
  GEM_mu_deps  = new TH1F("GEM_mu_deps",  "E_{deposit} :: GEM :: Muons",      n_F, n1_F, n2_F);
  GEM_ha_deps  = new TH1F("GEM_ha_deps",  "E_{deposit} :: GEM :: Hadrons",    n_F, n1_F, n2_F);
  ME0_el_deps  = new TH1F("ME0_el_deps",  "E_{deposit} :: ME0 :: Electrons",  n_F, n1_F, n2_F);
  ME0_mu_deps  = new TH1F("ME0_mu_deps",  "E_{deposit} :: ME0 :: Muons",      n_F, n1_F, n2_F);
  ME0_ha_deps  = new TH1F("ME0_ha_deps",  "E_{deposit} :: ME0 :: Hadrons",    n_F, n1_F, n2_F);

  RPCb_el_tof = new TH1F("RPCb_el_tof", "Time Of Flight :: RPCb :: Electrons", n_tof, n1_tof, n2_tof);
  RPCb_mu_tof = new TH1F("RPCb_mu_tof", "Time Of Flight :: RPCb :: Muons",     n_tof, n1_tof, n2_tof);
  RPCb_ha_tof = new TH1F("RPCb_ha_tof", "Time Of Flight :: RPCb :: Hadrons",   n_tof, n1_tof, n2_tof);
  RPCf_el_tof = new TH1F("RPCf_el_tof", "Time Of Flight :: RPCf :: Electrons", n_tof, n1_tof, n2_tof);
  RPCf_mu_tof = new TH1F("RPCf_mu_tof", "Time Of Flight :: RPCf :: Muons",     n_tof, n1_tof, n2_tof);
  RPCf_ha_tof = new TH1F("RPCf_ha_tof", "Time Of Flight :: RPCf :: Hadrons",   n_tof, n1_tof, n2_tof);
  CSC_el_tof  = new TH1F("CSC_el_tof",  "Time Of Flight :: CSC :: Electrons",  n_tof, n1_tof, n2_tof);
  CSC_mu_tof  = new TH1F("CSC_mu_tof",  "Time Of Flight :: CSC :: Muons",      n_tof, n1_tof, n2_tof);
  CSC_ha_tof  = new TH1F("CSC_ha_tof",  "Time Of Flight :: CSC :: Hadrons",    n_tof, n1_tof, n2_tof);
  DT_el_tof   = new TH1F("DT_el_tof",   "Time Of Flight :: DT :: Electrons",   n_tof, n1_tof, n2_tof);
  DT_mu_tof   = new TH1F("DT_mu_tof",   "Time Of Flight :: DT :: Muons",       n_tof, n1_tof, n2_tof);
  DT_ha_tof   = new TH1F("DT_ha_tof",   "Time Of Flight :: DT :: Hadrons",     n_tof, n1_tof, n2_tof);
  GEM_el_tof  = new TH1F("GEM_el_tof",  "Time Of Flight :: GEM :: Electrons",  n_tof, n1_tof, n2_tof);
  GEM_mu_tof  = new TH1F("GEM_mu_tof",  "Time Of Flight :: GEM :: Muons",      n_tof, n1_tof, n2_tof);
  GEM_ha_tof  = new TH1F("GEM_ha_tof",  "Time Of Flight :: GEM :: Hadrons",    n_tof, n1_tof, n2_tof);
  ME0_el_tof  = new TH1F("ME0_el_tof",  "Time Of Flight :: ME0 :: Electrons",  n_tof, n1_tof, n2_tof);
  ME0_mu_tof  = new TH1F("ME0_mu_tof",  "Time Of Flight :: ME0 :: Muons",      n_tof, n1_tof, n2_tof);
  ME0_ha_tof  = new TH1F("ME0_ha_tof",  "Time Of Flight :: ME0 :: Hadrons",    n_tof, n1_tof, n2_tof);


}


MyME11SimHitAnalyzer::~MyME11SimHitAnalyzer(){

  if(tech_debug) std::cout<<"[MyME11SimHitAnalyzer :: Destructor]"<<std::endl; 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  RPCf_el_hits->Write(); 
  RPCf_mu_hits->Write(); 
  RPCf_pi_hits->Write(); 
  RPCf_ka_hits->Write(); 
  RPCf_p_hits ->Write(); 
  RPCf_n_hits ->Write(); 
  RPCf_g_hits ->Write(); 
  RPCf_N_hits ->Write(); 
  RPCf_OH_hits ->Write(); 
  RPCf_All_hits->Write(); 
  RPCf_HIP_hits->Write(); 

  RPCb_el_hits->Write(); 
  RPCb_mu_hits->Write(); 
  RPCb_pi_hits->Write(); 
  RPCb_ka_hits->Write(); 
  RPCb_p_hits ->Write(); 
  RPCb_n_hits ->Write(); 
  RPCb_g_hits ->Write(); 
  RPCb_N_hits ->Write(); 
  RPCb_OH_hits ->Write(); 
  RPCb_All_hits->Write(); 
  RPCb_HIP_hits->Write(); 

  CSC_el_hits->Write(); 
  CSC_mu_hits->Write(); 
  CSC_pi_hits->Write(); 
  CSC_ka_hits->Write(); 
  CSC_p_hits ->Write(); 
  CSC_n_hits ->Write(); 
  CSC_g_hits ->Write(); 
  CSC_N_hits ->Write(); 
  CSC_OH_hits ->Write(); 
  CSC_All_hits->Write(); 
  CSC_HIP_hits->Write(); 

  DT_el_hits->Write(); 
  DT_mu_hits->Write(); 
  DT_pi_hits->Write(); 
  DT_ka_hits->Write(); 
  DT_p_hits ->Write(); 
  DT_n_hits ->Write(); 
  DT_g_hits ->Write(); 
  DT_N_hits ->Write(); 
  DT_OH_hits ->Write(); 
  DT_All_hits->Write(); 
  DT_HIP_hits->Write(); 

  GEM_el_hits->Write(); 
  GEM_mu_hits->Write(); 
  GEM_pi_hits->Write(); 
  GEM_ka_hits->Write(); 
  GEM_p_hits ->Write(); 
  GEM_n_hits ->Write(); 
  GEM_g_hits ->Write(); 
  GEM_N_hits ->Write(); 
  GEM_OH_hits ->Write(); 
  GEM_All_hits->Write(); 
  GEM_HIP_hits->Write(); 

  ME0_el_hits->Write(); 
  ME0_mu_hits->Write(); 
  ME0_pi_hits->Write(); 
  ME0_ka_hits->Write(); 
  ME0_p_hits ->Write(); 
  ME0_n_hits ->Write(); 
  ME0_g_hits ->Write(); 
  ME0_N_hits ->Write(); 
  ME0_OH_hits ->Write(); 
  ME0_All_hits->Write(); 
  ME0_HIP_hits->Write(); 
  // --------------------------- 
  outputfile->cd();

  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  RPCf_el_deposits->Write(); 
  RPCf_mu_deposits->Write(); 
  RPCf_pi_deposits->Write(); 
  RPCf_ka_deposits->Write(); 
  RPCf_p_deposits ->Write(); 
  RPCf_n_deposits ->Write(); 
  RPCf_g_deposits ->Write(); 
  RPCf_N_deposits ->Write(); 
  RPCf_OH_deposits ->Write(); 
  RPCf_All_deposits->Write(); 
  RPCf_HIP_deposits->Write(); 

  RPCb_el_deposits->Write(); 
  RPCb_mu_deposits->Write(); 
  RPCb_pi_deposits->Write(); 
  RPCb_ka_deposits->Write(); 
  RPCb_p_deposits ->Write(); 
  RPCb_n_deposits ->Write(); 
  RPCb_g_deposits ->Write(); 
  RPCb_N_deposits ->Write(); 
  RPCb_OH_deposits ->Write(); 
  RPCb_All_deposits->Write(); 
  RPCb_HIP_deposits->Write(); 

  CSC_el_deposits->Write(); 
  CSC_mu_deposits->Write(); 
  CSC_pi_deposits->Write(); 
  CSC_ka_deposits->Write(); 
  CSC_p_deposits ->Write(); 
  CSC_n_deposits ->Write(); 
  CSC_g_deposits ->Write(); 
  CSC_N_deposits ->Write(); 
  CSC_OH_deposits ->Write(); 
  CSC_All_deposits->Write(); 
  CSC_HIP_deposits->Write(); 

  DT_el_deposits->Write(); 
  DT_mu_deposits->Write(); 
  DT_pi_deposits->Write(); 
  DT_ka_deposits->Write(); 
  DT_p_deposits ->Write(); 
  DT_n_deposits ->Write(); 
  DT_g_deposits ->Write(); 
  DT_N_deposits ->Write(); 
  DT_OH_deposits ->Write(); 
  DT_All_deposits->Write(); 
  DT_HIP_deposits->Write(); 

  GEM_el_deposits->Write(); 
  GEM_mu_deposits->Write(); 
  GEM_pi_deposits->Write(); 
  GEM_ka_deposits->Write(); 
  GEM_p_deposits ->Write(); 
  GEM_n_deposits ->Write(); 
  GEM_g_deposits ->Write(); 
  GEM_N_deposits ->Write(); 
  GEM_OH_deposits ->Write(); 
  GEM_All_deposits->Write(); 
  GEM_HIP_deposits->Write(); 

  ME0_el_deposits->Write(); 
  ME0_mu_deposits->Write(); 
  ME0_pi_deposits->Write(); 
  ME0_ka_deposits->Write(); 
  ME0_p_deposits ->Write(); 
  ME0_n_deposits ->Write(); 
  ME0_g_deposits ->Write(); 
  ME0_N_deposits ->Write(); 
  ME0_OH_deposits ->Write(); 
  ME0_All_deposits->Write(); 
  ME0_HIP_deposits->Write(); 
  // --------------------------- 
  outputfile->cd();

  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  RPCb_el_deps->Write();
  RPCb_mu_deps->Write();
  RPCb_ha_deps->Write();
  RPCf_el_deps->Write();
  RPCf_mu_deps->Write();
  RPCf_ha_deps->Write();
  CSC_el_deps->Write();
  CSC_mu_deps->Write();
  CSC_ha_deps->Write();
  DT_el_deps->Write();
  DT_mu_deps->Write();
  DT_ha_deps->Write();
  GEM_el_deps->Write();
  GEM_mu_deps->Write();
  GEM_ha_deps->Write();
  ME0_el_deps->Write();
  ME0_mu_deps->Write();
  ME0_ha_deps->Write();
  // --------------------------- 


  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  RPCb_el_tof->Write();
  RPCb_mu_tof->Write();
  RPCb_ha_tof->Write();
  RPCf_el_tof->Write();
  RPCf_mu_tof->Write();
  RPCf_ha_tof->Write();
  CSC_el_tof->Write();
  CSC_mu_tof->Write();
  CSC_ha_tof->Write();
  DT_el_tof->Write();
  DT_mu_tof->Write();
  DT_ha_tof->Write();
  GEM_el_tof->Write();
  GEM_mu_tof->Write();
  GEM_ha_tof->Write();
  ME0_el_tof->Write();
  ME0_mu_tof->Write();
  ME0_ha_tof->Write();
  // --------------------------- 
  outputfile->cd();


  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  Muon_000ns_el_RZ->Write();
  Muon_000ns_mu_RZ->Write();
  Muon_000ns_ha_RZ->Write();
  Muon_250ns_el_RZ->Write();
  Muon_250ns_mu_RZ->Write();
  Muon_250ns_ha_RZ->Write();
  // --------------------------- 
  outputfile->cd();


  TDir_Muon_hits_radius->cd();
  // ---------------------------
  ME0_el_hits_R->Write(); ME0_mu_hits_R->Write(); ME0_pi_hits_R->Write(); ME0_ka_hits_R->Write(); ME0_p_hits_R->Write();  
  ME0_n_hits_R->Write();  ME0_g_hits_R->Write();  ME0_N_hits_R->Write();  ME0_OH_hits_R->Write(); ME0_All_hits_R->Write(); ME0_HIP_hits_R->Write();
  G11_el_hits_R->Write(); G11_mu_hits_R->Write(); G11_pi_hits_R->Write(); G11_ka_hits_R->Write(); G11_p_hits_R->Write(); 
  G11_n_hits_R->Write();  G11_g_hits_R->Write();  G11_N_hits_R->Write();  G11_OH_hits_R->Write(); G11_All_hits_R->Write(); G11_HIP_hits_R->Write();
  G21_el_hits_R->Write(); G21_mu_hits_R->Write(); G21_pi_hits_R->Write(); G21_ka_hits_R->Write(); G21_p_hits_R->Write(); 
  G21_n_hits_R->Write();  G21_g_hits_R->Write();  G21_N_hits_R->Write();  G21_OH_hits_R->Write(); G21_All_hits_R->Write(); G21_HIP_hits_R->Write();

  G11_L1_el_hits_R->Write(); G11_L1_mu_hits_R->Write(); G11_L1_pi_hits_R->Write(); G11_L1_ka_hits_R->Write(); G11_L1_p_hits_R->Write(); 
  G11_L1_n_hits_R->Write(); G11_L1_g_hits_R->Write(); G11_L1_N_hits_R->Write(); G11_L1_OH_hits_R->Write(); G11_L1_All_hits_R->Write(); G11_L1_HIP_hits_R->Write();
  G11_L2_el_hits_R->Write(); G11_L2_mu_hits_R->Write(); G11_L2_pi_hits_R->Write(); G11_L2_ka_hits_R->Write(); G11_L2_p_hits_R->Write(); 
  G11_L2_n_hits_R->Write(); G11_L2_g_hits_R->Write(); G11_L2_N_hits_R->Write(); G11_L2_OH_hits_R->Write(); G11_L2_All_hits_R->Write(); G11_L2_HIP_hits_R->Write();
  G11_Od_el_hits_R->Write(); G11_Od_mu_hits_R->Write(); G11_Od_pi_hits_R->Write(); G11_Od_ka_hits_R->Write(); G11_Od_p_hits_R->Write(); 
  G11_Od_n_hits_R->Write(); G11_Od_g_hits_R->Write(); G11_Od_N_hits_R->Write(); G11_Od_OH_hits_R->Write(); G11_Od_All_hits_R->Write(); G11_Od_HIP_hits_R->Write();
  G11_Ev_el_hits_R->Write(); G11_Ev_mu_hits_R->Write(); G11_Ev_pi_hits_R->Write(); G11_Ev_ka_hits_R->Write(); G11_Ev_p_hits_R->Write(); 
  G11_Ev_n_hits_R->Write(); G11_Ev_g_hits_R->Write(); G11_Ev_N_hits_R->Write(); G11_Ev_OH_hits_R->Write(); G11_Ev_All_hits_R->Write(); G11_Ev_HIP_hits_R->Write();
  M11_Od_el_hits_R->Write(); M11_Od_mu_hits_R->Write(); M11_Od_pi_hits_R->Write(); M11_Od_ka_hits_R->Write(); M11_Od_p_hits_R->Write(); 
  M11_Od_n_hits_R->Write(); M11_Od_g_hits_R->Write(); M11_Od_N_hits_R->Write(); M11_Od_OH_hits_R->Write(); M11_Od_All_hits_R->Write(); M11_Od_HIP_hits_R->Write();
  M11_Ev_el_hits_R->Write(); M11_Ev_mu_hits_R->Write(); M11_Ev_pi_hits_R->Write(); M11_Ev_ka_hits_R->Write(); M11_Ev_p_hits_R->Write(); 
  M11_Ev_n_hits_R->Write(); M11_Ev_g_hits_R->Write(); M11_Ev_N_hits_R->Write(); M11_Ev_OH_hits_R->Write(); M11_Ev_All_hits_R->Write(); M11_Ev_HIP_hits_R->Write();



  // ---------------------------
  TDir_Muon_hits_etapart->cd();
  // ---------------------------
  ME0_el_hits_E->Write(); ME0_mu_hits_E->Write(); ME0_pi_hits_E->Write(); ME0_ka_hits_E->Write(); ME0_p_hits_E->Write(); 
  ME0_n_hits_E->Write();  ME0_g_hits_E->Write();  ME0_N_hits_E->Write();  ME0_OH_hits_E->Write(); ME0_All_hits_E->Write(); ME0_HIP_hits_E->Write();
  G11_el_hits_E->Write(); G11_mu_hits_E->Write(); G11_pi_hits_E->Write(); G11_ka_hits_E->Write(); G11_p_hits_E->Write(); 
  G11_n_hits_E->Write();  G11_g_hits_E->Write();  G11_N_hits_E->Write();  G11_OH_hits_E->Write(); G11_All_hits_E->Write(); G11_HIP_hits_E->Write();
  G21_el_hits_E->Write(); G21_mu_hits_E->Write(); G21_pi_hits_E->Write(); G21_ka_hits_E->Write(); G21_p_hits_E->Write();   
  G21_n_hits_E->Write();  G21_g_hits_E->Write();  G21_N_hits_E->Write();  G21_OH_hits_E->Write(); G21_All_hits_E->Write(); G21_HIP_hits_E->Write();


  TDir_Muon_hits_Nuclei->cd();
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  Muon_Nuclei_A_Z->Write(); RPCb_Nuclei_A_Z->Write(); RPCf_Nuclei_A_Z->Write(); CSC_Nuclei_A_Z->Write(); DT_Nuclei_A_Z->Write(); GEM_Nuclei_A_Z->Write(); ME0_Nuclei_A_Z->Write();
  Muon_Nuclei_A->Write(); RPCb_Nuclei_A->Write(); RPCf_Nuclei_A->Write(); CSC_Nuclei_A->Write(); DT_Nuclei_A->Write(); GEM_Nuclei_A->Write(); ME0_Nuclei_A->Write();
  Muon_Nuclei_Z->Write(); RPCb_Nuclei_Z->Write(); RPCf_Nuclei_Z->Write(); CSC_Nuclei_Z->Write(); DT_Nuclei_Z->Write(); GEM_Nuclei_Z->Write(); ME0_Nuclei_Z->Write();
  Muon_Nuclei_List->Write(); RPCb_Nuclei_List->Write(); RPCf_Nuclei_List->Write(); CSC_Nuclei_List->Write(); DT_Nuclei_List->Write(); GEM_Nuclei_List->Write(); ME0_Nuclei_List->Write();

  outputfile->cd();

  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  RPCb_hits_tof->Write();
  RPCb_hits_eta->Write();
  RPCb_hits_phi->Write();
  RPCb_hits_lin->Write();

  RPCf_hits_tof->Write();
  RPCf_hits_eta->Write();
  RPCf_hits_phi->Write();
  RPCf_hits_lin->Write();

  CSC_hits_tof->Write();
  CSC_hits_eta->Write();
  CSC_hits_phi->Write();
  CSC_hits_lin->Write();

  DT_hits_tof->Write();
  DT_hits_eta->Write();
  DT_hits_phi->Write();
  DT_hits_lin->Write();

  GEM_hits_tof->Write();
  GEM_hits_eta->Write();
  GEM_hits_phi->Write();
  GEM_hits_lin->Write();

  ME0_hits_tof->Write();
  ME0_hits_eta->Write();
  ME0_hits_phi->Write();
  ME0_hits_lin->Write();
  /*
  RB4_hits_tof->Write();
  RE4_hits_tof->Write();
  MB4_hits_tof->Write();
  ME4_hits_tof->Write();
  */
  /*
  RB4_hits_phi->Write();
  RE4_hits_phi->Write();
  MB4_hits_phi->Write();
  ME4_hits_phi->Write();
  */

  RB1_hits_phi->Write();
  RB2_hits_phi->Write();
  RB3_hits_phi->Write();
  RB4_hits_phi->Write();
  MB1_hits_phi->Write();
  MB2_hits_phi->Write();
  MB3_hits_phi->Write();
  MB4_hits_phi->Write();
  RE12_hits_phi->Write();
  RE13_hits_phi->Write();
  RE22_hits_phi->Write();
  RE23_hits_phi->Write();
  RE32_hits_phi->Write();
  RE33_hits_phi->Write();
  RE42_hits_phi->Write();
  RE43_hits_phi->Write();
  ME11_hits_phi->Write();
  ME12_hits_phi->Write();
  ME13_hits_phi->Write();
  ME21_hits_phi->Write();
  ME22_hits_phi->Write();
  ME31_hits_phi->Write();
  ME32_hits_phi->Write();
  ME41_hits_phi->Write();
  ME42_hits_phi->Write();

  // ??? Probably I sacced twice ???
  /*
  RPCb_el_deps->Write();
  RPCb_mu_deps->Write();
  RPCb_ha_deps->Write();
  RPCf_el_deps->Write();
  RPCf_mu_deps->Write();
  RPCf_ha_deps->Write();
  CSC_el_deps->Write();
  CSC_mu_deps->Write();
  CSC_ha_deps->Write();
  DT_el_deps->Write();
  DT_mu_deps->Write();
  DT_ha_deps->Write();
  GEM_el_deps->Write();
  GEM_mu_deps->Write();
  GEM_ha_deps->Write();
  ME0_el_deps->Write();
  ME0_mu_deps->Write();
  ME0_ha_deps->Write();
  */
  RPCb_el_tof->Write();
  RPCb_mu_tof->Write();
  RPCb_ha_tof->Write();
  RPCf_el_tof->Write();
  RPCf_mu_tof->Write();
  RPCf_ha_tof->Write();
  CSC_el_tof->Write();
  CSC_mu_tof->Write();
  CSC_ha_tof->Write();
  DT_el_tof->Write();
  DT_mu_tof->Write();
  DT_ha_tof->Write();
  GEM_el_tof->Write();
  GEM_mu_tof->Write();
  GEM_ha_tof->Write();
  ME0_el_tof->Write();
  ME0_mu_tof->Write();
  ME0_ha_tof->Write();
  // --------------------------- 
  outputfile->cd();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyME11SimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(tech_debug) std::cout<<"[MyME11SimHitAnalyzer :: Analyze]"<<std::endl;

  const edm::ESTransientHandle<CSCGeometry> cscGeom = iSetup.getTransientHandle(cscGeom_Token);
  const edm::ESTransientHandle<DTGeometry>  dtGeom  = iSetup.getTransientHandle(dtGeom_Token);
  const edm::ESTransientHandle<GEMGeometry> gemGeom = iSetup.getTransientHandle(gemGeom_Token);
  const edm::ESTransientHandle<RPCGeometry> rpcGeom = iSetup.getTransientHandle(rpcGeom_Token);

  // ===============================
  //      General Collections
  // ===============================
  
  // SimHits
  // std::cout << " Getting the SimHits " <<std::endl;
  std::vector<edm::Handle<edm::PSimHitContainer> > theSimHitContainers;
  iEvent.getManyByType(theSimHitContainers);  // both 6XY and 7XY
  if(phys_debug) std::cout << " The number of SimHit Containers is  " << theSimHitContainers.size() <<std::endl;
  std::vector<PSimHit> theSimHits;
  for (int i = 0; i < int(theSimHitContainers.size()); ++i) {
    theSimHits.insert(theSimHits.end(),theSimHitContainers.at(i)->begin(),theSimHitContainers.at(i)->end());
  }
  // SimTracks
  std::vector<SimTrack> theSimTracks;
  edm::Handle<edm::SimTrackContainer> SimTk;
  // iEvent.getByLabel("g4SimHits",SimTk); // 6XY   
  iEvent.getByToken(SIMTrack_Token,SimTk); // 7XY 
  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  if(phys_debug) std::cout << "This Event has " <<  theSimTracks.size() << " sim tracks " << std::endl;
  // SimVertices
  std::vector<SimVertex> theSimVertices; 
  edm::Handle<edm::SimVertexContainer> SimVtx;
  // iEvent.getByLabel("g4SimHits",SimVtx);  // 6XY
  iEvent.getByToken(SIMVertex_Token,SimVtx); // 7XY 
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  if(phys_debug) std::cout << "This Event has " <<  theSimVertices.size() << " sim vertices " << std::endl;

  double vtx_r = 0.0, vtx_z = 0.0;
  for (std::vector<SimVertex>::const_iterator iVertex = theSimVertices.begin(); iVertex != theSimVertices.end(); ++iVertex) {
    SimVertex simvertex = (*iVertex);    
    // vtx_x = simvertex.position().x();
    // vtx_y = simvertex.position().y();
    vtx_r = sqrt(pow(simvertex.position().x(),2)+pow(simvertex.position().y(),2));
    vtx_z = simvertex.position().z();
    SimVertices_RZ->Fill(fabs(vtx_z),vtx_r);
    if( vtx_r < 2 && vtx_z < 25 ) {
      PrimVertices_Z->Fill(vtx_z);
      PrimVertices_R->Fill(vtx_r);
    }
    // if(vtx_r > 400 || vtx_z > 500 ) {               // muon + forward
    //   SimVertices_Muon_RZ->Fill(fabs(vtx_z),vtx_r);
    // }
    if(vtx_r > 400 || (vtx_z > 500 && vtx_r > 100) ) { // muon only
      SimVertices_Muon_RZ->Fill(fabs(vtx_z),vtx_r);
    }
  }

  
  // create a map associating geant particle id and position in the 
  // event SimTrack vector
  std::map<unsigned, unsigned> geantToIndex;
  for( unsigned it=0; it<theSimTracks.size(); ++it ) {
    geantToIndex[ theSimTracks[it].trackId() ] = it;
    // std::cout<<"geantToIndex[ "<<theSimTracks[it].trackId()<<", "<<it<<"]"<<std::endl;
  }  



  // ===============================
  //      Loop over the SimHits
  // ===============================

  for (std::vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); ++iHit) {
    DetId theDetUnitId((*iHit).detUnitId());
    DetId simdetid= DetId((*iHit).detUnitId());


    int pid            = (*iHit).particleType();
    int process        = (*iHit).processType();
    double time        = (*iHit).timeOfFlight();
    double log_time    = log10((*iHit).timeOfFlight());
    double log_energy  = log10((*iHit).momentumAtEntry().perp()*1000); // MeV
    double log_deposit = log10((*iHit).energyLoss()*1000000);          // keV
    double e_deposit   = (*iHit).energyLoss()*1000000; // keV

    // Highly Ionizing Particles
    bool HIP = false;
    if(e_deposit < varEDepCuteV*1E-3) continue;  // VarEDepCut is in eV
    if(e_deposit > 30.) HIP = true;              // HIP definition is in keV


    // Extract Nuclei (A,Z) :: Nuclear codes are given as 10-digit numbers: ±10LZZZAAAI. 
    int I=0, A=0, Z=0;
    if(abs(pid) > 999999999) {
      I = (int) ( abs(pid) %10);
      A = (int) ( ((abs(pid)-I) %10000) /10) ;
      Z = (int) ( ((abs(pid)-A-I) %10000000) / 10000) ;
      if(tech_debug) std::cout<<" SimHit with PDG ID = "<<pid<<" Extracted Isomer level I = "<<I<<" Atomic Number A = "<<A<<" Charge Number Z = "<<Z<<std::endl;
      if(Z<21) Muon_Nuclei_List->Fill(Z);
      else Muon_Nuclei_List->Fill(21);
    }

    if(tech_debug) std::cout<<"RPC Hits"<<std::endl;
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){ // Only RPCs

      // RPC Geometry
      // ============
      RPCDetId rollId(theDetUnitId);
      RPCGeomServ rpcsrv(rollId);
      const RPCRoll* rollasociated = rpcGeom->roll(rollId);
      const BoundPlane & RPCSurface = rollasociated->surface();
      GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal((*iHit).localPosition());
      GlobalPoint RPCGlobalEntry = RPCSurface.toGlobal((*iHit).entryPoint());
      GlobalPoint RPCGlobalExit  = RPCSurface.toGlobal((*iHit).exitPoint());
      double RPCGlobalEntryExitDZ = fabs(RPCGlobalEntry.z()-RPCGlobalExit.z());
      double RPCGlobalEntryExitDR = fabs(sqrt(pow(RPCGlobalEntry.x(),2)+pow(RPCGlobalEntry.y(),2))-sqrt(pow(RPCGlobalExit.x(),2)+pow(RPCGlobalExit.y(),2)));
      double RPCLocalEntryExitDZ  = fabs((*iHit).entryPoint().z()-(*iHit).exitPoint().z());

      if(phys_debug) {
	std::cout<<"RPC SimHit in "<<std::setw(12)<<(int)rollId<< /*" a.k.a. "<<std::setw(24)<<rpcsrv.name()<<" details: "<<*/ std::setw(24)<<rollId;
	std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
	std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta();
	std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
      }
      double RPC_GlobalPoint_R = sqrt(pow(RPCGlobalPoint.x(),2)+pow(RPCGlobalPoint.y(),2));

      // SIMHIT Investigations
      // =====================
      /*
      // get track
      int trackId = (*iHit).trackId();
      const SimTrack& track = theSimTracks[trackId];
      // get vertex
      int vertexId = track.vertIndex();
      const SimVertex& vertex = theSimVertices[vertexId];
      // get mother
      int motherId = -1;
      if( !vertex.noParent() ) { // there is a parent to this vertex
	// geant id of the mother
	unsigned motherGeantId =   vertex.parentIndex(); 
	std::map<unsigned, unsigned >::iterator association = geantToIndex.find( motherGeantId );
	if(association != geantToIndex.end() ) {
	  motherId = association->second;
	}
      }
      // int originId = motherId == - 1 ? -1 : myTracks[motherId];
      std::cout<<"TrackId = "<<trackId<<" VertexId = "<<vertexId<<" OriginId = "<<motherId<<std::endl;
      */

      // count hits
      // ----------
      int region = abs(rollId.region());
      int station = abs(rollId.station())-1;
      int ring = (region==0)? abs(rollId.ring()) : abs(rollId.ring())-1;
      ++RPC_hits_array[region][station][ring];


      // RPCb
      // ----
      if(rollId.region() == 0) {
	if(abs(pid)==11)      { RPCb_Electrons_SHPT->Fill(process);        RPCb_el_deps->Fill(log_deposit);  RPCb_el_kins->Fill(log_energy);  RPCb_el_tof->Fill(log_time);}
	else if(abs(pid)==13) { RPCb_Muons_SHPT->Fill(process);            RPCb_mu_deps->Fill(log_deposit);  RPCb_mu_kins->Fill(log_energy);  RPCb_mu_tof->Fill(log_time);}
	else                  { RPCb_Hadrons_SHPT->Fill(process);          RPCb_ha_deps->Fill(log_deposit);  RPCb_ha_kins->Fill(log_energy);  RPCb_ha_tof->Fill(log_time);}
	RPCb_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());             
	RPCb_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));  
	Muon_Barrel_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	Muon_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));  
	RPCb_EntryExit_All_Glob_dr->Fill(RPCGlobalEntryExitDR);
	RPCb_EntryExit_All_Loc_dz->Fill(RPCLocalEntryExitDZ);

	if((*iHit).timeOfFlight()<250) {
	  if(abs(pid)==11)      { Muon_000ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_000ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  else if(abs(pid)==13) { Muon_000ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_000ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  else                  { Muon_000ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_000ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  if((*iHit).timeOfFlight()<50) {
	    if(abs(pid)==11)      { Muon_00ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_00ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y()); RPCb_Electrons_000ns_SHPT->Fill(process);}
	    else if(abs(pid)==13) { Muon_00ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_00ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	    else                  { Muon_00ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_00ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  }
	  if((*iHit).timeOfFlight()>50) {
	    if(abs(pid)==11)      { Muon_50ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_50ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());RPCb_Electrons_050ns_SHPT->Fill(process);}
	    else if(abs(pid)==13) { Muon_50ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_50ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	    else                  { Muon_50ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_50ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  }
	  RPCb_000ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  RPCb_000ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));
	  Muon_Barrel_000ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  Muon_000ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); 
	}
	if((*iHit).timeOfFlight()>250) {
	  if(abs(pid)==11)      { Muon_250ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_250ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y()); RPCb_Electrons_250ns_SHPT->Fill(process);}
	  else if(abs(pid)==13) { Muon_250ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_250ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  else {                  Muon_250ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Barrel_250ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  RPCb_250ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  RPCb_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));
	  Muon_Barrel_250ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  Muon_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); 
	}

	// All Hits and HIP Hits
	RPCb_All_hits->Fill(log_energy,log_time); RPCb_All_deposits->Fill(log_deposit,log_time);
	if(HIP) {
	  RPCb_HIP_hits->Fill(log_energy,log_time); RPCb_HIP_deposits->Fill(log_deposit,log_time);
	}
	// Catalog Nuclei
	if(abs(pid) > 999999999) {
	  if(Z<21) RPCb_Nuclei_List->Fill(Z);
	  else RPCb_Nuclei_List->Fill(21);
	}

	// std::cout<<"RPCb :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and deposit deposit "<<(*iHit).energyLoss()<< " [GeV]";
	// std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
	if(abs(pid)==2212)      {RPCb_p_hits->Fill(log_energy,log_time); RPCb_p_deposits->Fill(log_deposit,log_time);}
	else if(abs(pid)==2112) {RPCb_n_hits->Fill(log_energy,log_time); RPCb_n_deposits->Fill(log_deposit,log_time);}
	else { 
	  switch (abs(pid)%1000) {
	    // leptopdfns
	  case 11:  RPCb_el_hits->Fill(log_energy,log_time); RPCb_el_deposits->Fill(log_deposit,log_time); break;
	  case 13:  RPCb_mu_hits->Fill(log_energy,log_time); RPCb_mu_deposits->Fill(log_deposit,log_time); break;
	    // Pions
	  case 111: RPCb_pi_hits->Fill(log_energy,log_time); RPCb_pi_deposits->Fill(log_deposit,log_time); break;
	  case 211: RPCb_pi_hits->Fill(log_energy,log_time); RPCb_pi_deposits->Fill(log_deposit,log_time); break;
	  case 130: RPCb_pi_hits->Fill(log_energy,log_time); RPCb_pi_deposits->Fill(log_deposit,log_time); break;
	    // Kaons
	  case 310: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 311: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 321: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 313: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 323: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 315: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 325: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 317: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 327: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 319: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	  case 329: RPCb_ka_hits->Fill(log_energy,log_time); RPCb_ka_deposits->Fill(log_deposit,log_time); break;
	    // Protons
	    // case 2212: RPCb_p_hits->Fill(log_energy,log_time); RPCb_p_deposits->Fill(log_deposit,log_time); break;
	    // Neutrons
	    // case 2112: RPCb_n_hits->Fill(log_energy,log_time); RPCb_n_deposits->Fill(log_deposit,log_time); break;
	    // Photons
	  case 22:   RPCb_g_hits->Fill(log_energy,log_time); RPCb_g_deposits->Fill(log_deposit,log_time); break;
	    // Nucleons
	    // Examples:
	    // 1000060120 => A = 012, Z = 006 => Carbon
	    // 1000080160 => A = 016, Z = 008 => Oxygen
	    // 1000170360 => A = 036, Z = 017 => Chlorine
	    // 1000010020 => A = 002, Z = 001 => Deuterium
	    // 1000010030 => A = 003, Z = 001 => Tritium
            // 0999999999
	  default:   {
	    if(abs(pid) > 999999999) { // Nucleons (10-digit numbers) 
	      RPCb_N_hits->Fill(log_energy,log_time); RPCb_N_deposits->Fill(log_deposit,log_time);
	      std::cout<<"RPCb :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and deposit deposit "<<(*iHit).energyLoss()<< " [GeV]";
	      std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
	    }
	    else { // non-Nucleon ... catalogue as "Other Hadron" 
	      RPCb_OH_hits->Fill(log_energy,log_time); RPCb_OH_deposits->Fill(log_deposit,log_time);
	      std::cout<<"RPCf :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<<" [GeV]";
	      std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;
	    }
	    break;
	  }
	  }
	}
	RPCb_hits_tof->Fill(log_time);
	RPCb_hits_eta->Fill(fabs(RPCGlobalPoint.eta()));
	RPCb_hits_phi->Fill(RPCGlobalPoint.phi());
	RPCb_hits_lin->Fill(time);
	if(rollId.station() == 1) { RB1_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 2) { RB2_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 3) { RB3_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 4) { RB4_hits_phi->Fill(RPCGlobalPoint.phi()); /*RB4_hits_tof->Fill(log_time);*/}
      }
      if(rollId.region() != 0) {
	if(abs(pid)==11)      { RPCf_Electrons_SHPT->Fill(process);        RPCf_el_deps->Fill(log_deposit);  RPCf_el_kins->Fill(log_energy);  RPCf_el_tof->Fill(log_time);}
	else if(abs(pid)==13) { RPCf_Muons_SHPT->Fill(process);            RPCf_mu_deps->Fill(log_deposit);  RPCf_mu_kins->Fill(log_energy);  RPCf_mu_tof->Fill(log_time);}
	else                  { RPCf_Hadrons_SHPT->Fill(process);          RPCf_ha_deps->Fill(log_deposit);  RPCf_ha_kins->Fill(log_energy);  RPCf_ha_tof->Fill(log_time);}
	RPCf_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());             
	RPCf_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));  
	Muon_Endcap_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	Muon_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));  
	RPCf_EntryExit_All_Glob_dz->Fill(RPCGlobalEntryExitDZ);
	RPCf_EntryExit_All_Loc_dz->Fill(RPCLocalEntryExitDZ);

	if((*iHit).timeOfFlight()<250) {
	  if(abs(pid)==11)      { Muon_000ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_000ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  else if(abs(pid)==13) { Muon_000ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_000ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  else                  { Muon_000ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_000ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  if((*iHit).timeOfFlight()<50) {
	    if(abs(pid)==11)      { Muon_00ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_00ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y()); RPCf_Electrons_000ns_SHPT->Fill(process);}
	    else if(abs(pid)==13) { Muon_00ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_00ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	    else                  { Muon_00ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_00ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  }
	  if((*iHit).timeOfFlight()>50) {
	    if(abs(pid)==11)      { Muon_50ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_50ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y()); RPCf_Electrons_050ns_SHPT->Fill(process);}
	    else if(abs(pid)==13) { Muon_50ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_50ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	    else                  { Muon_50ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_50ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  }
	  RPCf_000ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  RPCf_000ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));
	  Muon_Endcap_000ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  Muon_000ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); 
	}
	if((*iHit).timeOfFlight()>250) {
	  if(abs(pid)==11)      { Muon_250ns_el_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_250ns_el_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y()); RPCf_Electrons_250ns_SHPT->Fill(process);}
	  else if(abs(pid)==13) { Muon_250ns_mu_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_250ns_mu_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  else                  { Muon_250ns_ha_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); Muon_Endcap_250ns_ha_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	  RPCf_250ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  RPCf_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));
	  Muon_Endcap_250ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());
	  Muon_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); 
	}

	// All Hits and HIP Hits
	RPCf_All_hits->Fill(log_energy,log_time); RPCf_All_deposits->Fill(log_deposit,log_time);
	if(HIP) {
	  RPCf_HIP_hits->Fill(log_energy,log_time); RPCf_HIP_deposits->Fill(log_deposit,log_time);
	}
	// Catalog Nuclei
	if(abs(pid) > 999999999) {
	  if(Z<21) RPCf_Nuclei_List->Fill(Z);
	  else RPCf_Nuclei_List->Fill(21);
	}

	// std::cout<<"RPCf :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<<" [GeV]";
	// std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
	if(abs(pid)==2212)      {RPCf_p_hits->Fill(log_energy,log_time); RPCf_p_deposits->Fill(log_deposit,log_time);}
	else if(abs(pid)==2112) {RPCf_n_hits->Fill(log_energy,log_time); RPCf_n_deposits->Fill(log_deposit,log_time);}
	else { 
	  switch (abs(pid)%1000) {
	    // leptons
	  case 11:  RPCf_el_hits->Fill(log_energy,log_time); RPCf_el_deposits->Fill(log_deposit,log_time); break;
	  case 13:  RPCf_mu_hits->Fill(log_energy,log_time); RPCf_mu_deposits->Fill(log_deposit,log_time); break;
	    // Pions
	  case 111: RPCf_pi_hits->Fill(log_energy,log_time); RPCf_pi_deposits->Fill(log_deposit,log_time); break;
	  case 211: RPCf_pi_hits->Fill(log_energy,log_time); RPCf_pi_deposits->Fill(log_deposit,log_time); break;
	  case 130: RPCf_pi_hits->Fill(log_energy,log_time); RPCf_pi_deposits->Fill(log_deposit,log_time); break;
	    // Kaons
	  case 310: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 311: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 321: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 313: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 323: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 315: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 325: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 317: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 327: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 319: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	  case 329: RPCf_ka_hits->Fill(log_energy,log_time); RPCf_ka_deposits->Fill(log_deposit,log_time); break;
	    // Protons
	    // case 2212: RPCf_p_hits->Fill(log_energy,log_time); RPCf_p_deposits->Fill(log_deposit,log_time); break;
	    // Neutrons
	    // case 2112: RPCf_n_hits->Fill(log_energy,log_time); RPCf_n_deposits->Fill(log_deposit,log_time); break;
	    // Photons
	  case 22:   RPCf_g_hits->Fill(log_energy,log_time); RPCf_g_deposits->Fill(log_deposit,log_time); break;
	    // Nucleons
	  default:   {
	    if(abs(pid) > 999999999) { // Nucleons (10-digit numbers) 
	      RPCf_N_hits->Fill(log_energy,log_time); RPCf_N_deposits->Fill(log_deposit,log_time);
	      std::cout<<"RPCf :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and deposit deposit "<<(*iHit).energyLoss()<< " [GeV]";
	      std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
	    }
	    else { // non-Nucleon ... catalogue as "Other Hadron" 
	      RPCf_OH_hits->Fill(log_energy,log_time); RPCf_OH_deposits->Fill(log_deposit,log_time);
	      std::cout<<"RPCf :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<<" [GeV]";
	      std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;
	    }
	    break;
	  }
	  }
	}
	RPCf_hits_tof->Fill(log_time);
	RPCf_hits_eta->Fill(fabs(RPCGlobalPoint.eta()));
	RPCf_hits_phi->Fill(RPCGlobalPoint.phi());
	RPCf_hits_lin->Fill(time);
	if(rollId.station() == 1 && rollId.ring()==2) { RE12_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 1 && rollId.ring()==3) { RE13_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 2 && rollId.ring()==2) { RE22_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 2 && rollId.ring()==3) { RE23_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 3 && rollId.ring()==2) { RE32_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 3 && rollId.ring()==3) { RE33_hits_phi->Fill(RPCGlobalPoint.phi()); }
	if(rollId.station() == 4 && rollId.ring()==2) { RE42_hits_phi->Fill(RPCGlobalPoint.phi()); /*RE4_hits_tof->Fill(log_time);*/ }
	if(rollId.station() == 4 && rollId.ring()==3) { RE43_hits_phi->Fill(RPCGlobalPoint.phi()); /*RE4_hits_tof->Fill(log_time);*/ }
      }
    }

    if(tech_debug) std::cout<<"CSC Hits"<<std::endl;
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::CSC){ // Only CSCs
      CSCDetId rollId(theDetUnitId);
      GlobalPoint CSCGlobalPoint = cscGeom->idToDet(rollId)->toGlobal((*iHit).localPosition());
      GlobalPoint CSCGlobalEntry = cscGeom->idToDet(rollId)->toGlobal((*iHit).entryPoint());
      GlobalPoint CSCGlobalExit  = cscGeom->idToDet(rollId)->toGlobal((*iHit).exitPoint());
      double CSCGlobalEntryExitDZ = fabs(CSCGlobalEntry.z()-CSCGlobalExit.z());
      double CSCLocalEntryExitDZ = fabs((*iHit).entryPoint().z()-(*iHit).exitPoint().z());
      double CSCGlobalEntryExitDR = fabs(sqrt(pow(CSCGlobalEntry.x(),2)+pow(CSCGlobalEntry.y(),2)) - sqrt(pow(CSCGlobalExit.x(),2)+pow(CSCGlobalExit.y(),2)));
      LocalVector GapVector = ((*iHit).exitPoint() - (*iHit).entryPoint());
      double GapLength= GapVector.mag();
      int layer = rollId.layer();

      bool isME11 = false, isME11odd = false, isME11even = false;
      if(rollId.station() == 1 && (rollId.ring() == 1 || rollId.ring() == 4) ) {
	isME11 = true;
	if(rollId.chamber()%2==1) isME11odd = true; else isME11even = true;
      }


      if (phys_debug) {
       std::cout<<"CSC SimHit in "<<std::setw(24)<<rollId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
       std::cout<<" | r = "<<std::setw(12)<<CSCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<CSCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<CSCGlobalPoint.eta();
       std::cout<<" | global position = "<<CSCGlobalPoint<<std::endl;
      }
      /*
      if(fabs(CSCGlobalEntry.z()-CSCGlobalExit.z()) < 4) {
	std::cout<<"CSC Problematic Simhit :: dz = "<<fabs(CSCGlobalEntry.z()-CSCGlobalExit.z())<<" unchanged during coord transfo :: dz = "<<fabs(CSCLocalEntry.z()-CSCLocalExit.z())<<std::endl;
	std::cout<<" Local Entry Point :: "<<CSCLocalEntry<<" Local Exit Point :: "<<CSCLocalExit<<std::endl;
	std::cout<<" Global Entry Point :: "<<CSCGlobalEntry<<" Global Exit Point :: "<<CSCGlobalExit<<std::endl;
      }
      if(fabs(CSCGlobalEntry.z()-CSCGlobalExit.z()) > 4) {
	std::cout<<"CSC Good Simhit :: dz = "<<fabs(CSCGlobalEntry.z()-CSCGlobalExit.z())<<std::endl;
	std::cout<<" Local Entry Point :: "<<CSCLocalEntry<<" Local Exit Point :: "<<CSCLocalExit<<std::endl;
	std::cout<<" Global Entry Point :: "<<CSCGlobalEntry<<" Global Exit Point :: "<<CSCGlobalExit<<std::endl;
      }
      */
      double CSC_GlobalPoint_R = sqrt(pow(CSCGlobalPoint.x(),2)+pow(CSCGlobalPoint.y(),2));
      if(abs(pid)==11)      { 
	CSC_Electrons_SHPT->Fill(process);        CSC_el_deps->Fill(log_deposit);  CSC_el_kins->Fill(log_energy);  CSC_el_tof->Fill(log_time);
	CSC_EntryExit_Electrons_Glob_dz->Fill(CSCGlobalEntryExitDZ); 
	CSC_EntryExit_Electrons_Loc_dz->Fill(CSCLocalEntryExitDZ);
	// CSC_EntryExit_Electrons_Glob_dR->Fill(CSCGlobalEntryExitDR);
	// CSC_EntryExit_Electrons_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ);
	CSC_EntryExit_Electrons_Glob_dGap->Fill(GapLength); 
      }
      else if(abs(pid)==13) { 
	CSC_Muons_SHPT->Fill(process);            CSC_mu_deps->Fill(log_deposit); CSC_mu_kins->Fill(log_energy);  CSC_mu_tof->Fill(log_time);
	CSC_EntryExit_Muons_Glob_dz->Fill(CSCGlobalEntryExitDZ);     
	CSC_EntryExit_Muons_Loc_dz->Fill(CSCLocalEntryExitDZ);
	// CSC_EntryExit_Muons_Glob_dR->Fill(CSCGlobalEntryExitDR);
	// CSC_EntryExit_Muons_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ);
	CSC_EntryExit_Muons_Glob_dGap->Fill(GapLength); 
      }
      else                  { 
	CSC_Hadrons_SHPT->Fill(process);          CSC_ha_deps->Fill(log_deposit); CSC_ha_kins->Fill(log_energy);  CSC_ha_tof->Fill(log_time);
	CSC_EntryExit_Hadrons_Glob_dz->Fill(CSCGlobalEntryExitDZ);   
	CSC_EntryExit_Hadrons_Loc_dz->Fill(CSCLocalEntryExitDZ);
	// CSC_EntryExit_Hadrons_Glob_dR->Fill(CSCGlobalEntryExitDR);
	// CSC_EntryExit_Hadrons_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ);
	CSC_EntryExit_Hadrons_Glob_dGap->Fill(GapLength); 
      }
      CSC_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());             
      CSC_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));  
      Muon_Endcap_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());
      Muon_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));  
      CSC_EntryExit_All_Glob_dz->Fill(CSCGlobalEntryExitDZ);
      CSC_EntryExit_All_Loc_dz->Fill(CSCLocalEntryExitDZ);
      CSC_EntryExit_All_Glob_dGap->Fill(GapLength);
      // CSC_EntryExit_All_Glob_dR->Fill(CSCGlobalEntryExitDR);
      // CSC_EntryExit_All_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ);

      if(rollId.station() == 1) { CSC_ME1_all_hits->Fill(fabs(CSC_GlobalPoint_R)); CSC_Geom_ME1_all_hits->Fill(fabs(CSC_GlobalPoint_R));}
      if(rollId.station() == 2) { CSC_ME2_all_hits->Fill(fabs(CSC_GlobalPoint_R)); CSC_Geom_ME2_all_hits->Fill(fabs(CSC_GlobalPoint_R));}
      if(rollId.station() == 3) { CSC_ME3_all_hits->Fill(fabs(CSC_GlobalPoint_R)); CSC_Geom_ME3_all_hits->Fill(fabs(CSC_GlobalPoint_R));}
      if(rollId.station() == 4) { CSC_ME4_all_hits->Fill(fabs(CSC_GlobalPoint_R)); CSC_Geom_ME4_all_hits->Fill(fabs(CSC_GlobalPoint_R));}

      if((*iHit).timeOfFlight()<250) {
	if(abs(pid)==11)      { Muon_000ns_el_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_000ns_el_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	else if(abs(pid)==13) { Muon_000ns_mu_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_000ns_mu_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	else                  { Muon_000ns_ha_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_000ns_ha_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	// CSC hit rate plots
	if(rollId.station() == 1) { CSC_ME1_000ns_hits->Fill(fabs(CSC_GlobalPoint_R)); }

	if((*iHit).timeOfFlight()<50) {
	  if(abs(pid)==11)      { Muon_00ns_el_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_00ns_el_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y()); CSC_Electrons_000ns_SHPT->Fill(process);}
	  else if(abs(pid)==13) { Muon_00ns_mu_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_00ns_mu_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	  else                  { Muon_00ns_ha_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_00ns_ha_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	  // CSC hit rate plots
	  if(rollId.station() == 1) { CSC_ME1_00ns_hits->Fill(fabs(CSC_GlobalPoint_R)); }
	}

	if((*iHit).timeOfFlight()>50) {
	  if(abs(pid)==11)      { Muon_50ns_el_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_50ns_el_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y()); CSC_Electrons_050ns_SHPT->Fill(process);}
	  else if(abs(pid)==13) { Muon_50ns_mu_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_50ns_mu_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	  else                  { Muon_50ns_ha_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_50ns_ha_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	  // CSC hit rate plots
	  if(rollId.station() == 1) { CSC_ME1_50ns_hits->Fill(fabs(CSC_GlobalPoint_R)); }
	}
	CSC_000ns_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());
	CSC_000ns_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));
	Muon_Endcap_000ns_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());
	Muon_000ns_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); 
      }
      if((*iHit).timeOfFlight()>250) {
	if(abs(pid)==11)      { Muon_250ns_el_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_000ns_el_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y()); CSC_Electrons_250ns_SHPT->Fill(process);}
	else if(abs(pid)==13) { Muon_250ns_mu_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_250ns_mu_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	else                  { Muon_250ns_ha_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); Muon_Endcap_250ns_ha_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
	CSC_250ns_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());
	CSC_250ns_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));
	Muon_Endcap_250ns_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());
	Muon_250ns_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); 
	// CSC hit rate plots
	if(rollId.station() == 1) { CSC_ME1_250ns_hits->Fill(fabs(CSC_GlobalPoint_R)); }
      }
      
      // int pid           = (*iHit).particleType();
      // int process       = (*iHit).processType();
      if(phys_debug) std::cout<<"SimHit from particle = "<<pid<<" created in process = "<<process<<" aka "<<proc[process]<<std::endl;
      // There is a problem from here on ... 
      /*
      // get track
      int trackIdG = (*iHit).trackId();                     // GEANT TRACK ID
      int trackIdC = geantToIndex.find(trackIdG)->second;   // CMSSW VECTOR ID
      const SimTrack& track = theSimTracks[trackIdC];
      std::cout<<"TrackId Geant = "<<std::setw(7)<<trackIdG<<" TrackId CMSSW = "<<std::setw(5)<<trackIdC<<" track = "<<&track;
      std::cout<<" SimHit Particle Type = "<<std::setw(5)<<pid<<" SimTrack Particle Type = "<<std::setw(5)<<track.type()<<std::endl;
      */
      /*
      // get vertex        
      int vertexId = track.vertIndex();
      std::cout<<" VertexId = "<<vertexId<<std::endl;
      const SimVertex& vertex = theSimVertices[vertexId];
      process = vertex.processType();
      std::cout<<" VertexId = "<<vertexId<<" Vertex = "<<&vertex<<" Process = "<<process<<" aka "<<proc[process]<<std::endl;
      */
      // get mother
      /*
      int motherId = -1;
      if( !vertex.noParent() ) { // there is a parent to this vertex
	// geant id of the mother
	unsigned motherGeantId =   vertex.parentIndex(); 
	std::map<unsigned, unsigned >::iterator association = geantToIndex.find( motherGeantId );
	if(association != geantToIndex.end() ) {
	  motherId = association->second;
	}
      }
      // int originId = motherId == - 1 ? -1 : myTracks[motherId];
      std::cout<<"TrackId = "<<trackIdG<<" VertexId = "<<vertexId<<" OriginId = "<<motherId<<std::endl;
      const SimTrack& mother = theSimTracks[motherId];
      int motherGenId =  mother.genpartIndex();
      std::cout<<"TrackId = "<<trackIdG<<" VertexId = "<<vertexId<<" OriginId = "<<motherId<<" GenParicleId = "<<motherGenId<<std::endl;
      */

      // All Hits and HIP Hits
      CSC_All_hits->Fill(log_energy,log_time); CSC_All_deposits->Fill(log_deposit,log_time);
      if(isME11) { if(isME11odd) {M11_Od_All_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_All_hits_R->Fill(CSC_GlobalPoint_R);} }
      if(HIP) {
	CSC_HIP_hits->Fill(log_energy,log_time); CSC_HIP_deposits->Fill(log_deposit,log_time);
	if(isME11) { if(isME11odd) {M11_Od_HIP_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_HIP_hits_R->Fill(CSC_GlobalPoint_R);} }
      }
      // Catalog Nuclei
      if(abs(pid) > 999999999) {
	if(Z<21) CSC_Nuclei_List->Fill(Z);
	else CSC_Nuclei_List->Fill(21);
      }

      // std::cout<<"CSC :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
      // std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
      if(abs(pid)==2212)      {
	CSC_p_hits->Fill(log_energy,log_time); CSC_p_deposits->Fill(log_deposit,log_time); CSC_EntryExit_p_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); 
        if(isME11) { if(isME11odd) {M11_Od_p_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_p_hits_R->Fill(CSC_GlobalPoint_R);} }
	CSC_EntryExit_p_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_p_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR);
	CSC_EntryExit_p_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_p_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); 
	CSC_EntryExit_p_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);  CSC_EntryExit_p_GapLength_Deposit->Fill(log_deposit,GapLength);
	if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_p_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_p_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));

      }
      else if(abs(pid)==2112) {
	CSC_n_hits->Fill(log_energy,log_time); CSC_n_deposits->Fill(log_deposit,log_time); CSC_EntryExit_n_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); 
	if(isME11) { if(isME11odd) {M11_Od_n_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_n_hits_R->Fill(CSC_GlobalPoint_R);} }
	CSC_EntryExit_n_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_n_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR);
	CSC_EntryExit_n_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_n_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); 
	CSC_EntryExit_n_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);  CSC_EntryExit_n_GapLength_Deposit->Fill(log_deposit,GapLength);
	if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_n_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_n_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
      }
      else {
	switch (abs(pid)%1000) {
	  // leptons
	case 11:  
	  CSC_el_hits->Fill(log_energy,log_time); CSC_el_deposits->Fill(log_deposit,log_time); CSC_EntryExit_el_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); CSC_EntryExit_el_Time_dz->Fill(CSCGlobalEntryExitDZ,log_time);
	  if(isME11) { if(isME11odd) {M11_Od_el_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_el_hits_R->Fill(CSC_GlobalPoint_R);} }
	  CSC_EntryExit_el_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_el_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR); 
	  CSC_EntryExit_el_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_el_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); 
	  CSC_EntryExit_el_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);  CSC_EntryExit_el_GapLength_Deposit->Fill(log_deposit,GapLength);
	  if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_el_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	  if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_el_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
	  break;
	case 13:  
	  CSC_mu_hits->Fill(log_energy,log_time); CSC_mu_deposits->Fill(log_deposit,log_time); CSC_EntryExit_mu_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); CSC_EntryExit_mu_Time_dz->Fill(CSCGlobalEntryExitDZ,log_time);
	  if(isME11) { if(isME11odd) {M11_Od_mu_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_mu_hits_R->Fill(CSC_GlobalPoint_R);} }
	  CSC_EntryExit_mu_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_mu_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR); 
	  CSC_EntryExit_mu_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_mu_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);
	  CSC_EntryExit_mu_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); CSC_EntryExit_mu_GapLength_Deposit->Fill(log_deposit,GapLength);
	  if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_mu_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	  if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_mu_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
	  break;
	  // Photons
	case 22:   CSC_g_hits->Fill(log_energy,log_time); CSC_g_deposits->Fill(log_deposit,log_time); CSC_EntryExit_g_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); CSC_EntryExit_g_Time_dz->Fill(CSCGlobalEntryExitDZ,log_time);
	  if(isME11) { if(isME11odd) {M11_Od_g_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_g_hits_R->Fill(CSC_GlobalPoint_R);} }
	  CSC_EntryExit_g_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_g_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR); 
	  CSC_EntryExit_g_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_g_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);
	  CSC_EntryExit_g_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); CSC_EntryExit_g_GapLength_Deposit->Fill(log_deposit,GapLength);
	  if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_g_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	  if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_g_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
	  break;
	  // Pions
	case 111:
	case 211:
        case 130:
	  CSC_pi_hits->Fill(log_energy,log_time); CSC_pi_deposits->Fill(log_deposit,log_time); CSC_EntryExit_pi_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); CSC_EntryExit_pi_Time_dz->Fill(CSCGlobalEntryExitDZ,log_time);
	  if(isME11) { if(isME11odd) {M11_Od_pi_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_pi_hits_R->Fill(CSC_GlobalPoint_R);} }
	  CSC_EntryExit_pi_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_pi_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR); 
	  CSC_EntryExit_pi_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_pi_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);
	  CSC_EntryExit_pi_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); CSC_EntryExit_pi_GapLength_Deposit->Fill(log_deposit,GapLength);
	  if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_pi_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	  if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_pi_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
	  break;
	  // Kaons
	case 310: 
	case 311:
	case 313:
	case 315:
	case 317:
	case 319:
	case 321:
	case 323:
	case 325:
	case 327:
	case 329:
	  CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); CSC_EntryExit_ka_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); CSC_EntryExit_ka_Time_dz->Fill(CSCGlobalEntryExitDZ,log_time);
	  if(isME11) { if(isME11odd) {M11_Od_ka_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_ka_hits_R->Fill(CSC_GlobalPoint_R);} }
	  CSC_EntryExit_ka_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_ka_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR); 
	  CSC_EntryExit_ka_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR);  CSC_EntryExit_ka_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);
	  CSC_EntryExit_ka_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); CSC_EntryExit_ka_GapLength_Deposit->Fill(log_deposit,GapLength);
	  if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_ka_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	  if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) CSC_EntryExit_ka_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
	  break;
	  // Protons
	  // case 2212: CSC_p_hits->Fill(log_energy,log_time); CSC_p_deposits->Fill(log_deposit,log_time); break;
	  // Neutrons
	  // case 2112: CSC_n_hits->Fill(log_energy,log_time); CSC_n_deposits->Fill(log_deposit,log_time); break;

	  // Nucleons
	default:   {
	    if(abs(pid) > 999999999) { // Nucleons (10-digit numbers) 
	      CSC_N_hits->Fill(log_energy,log_time); CSC_N_deposits->Fill(log_deposit,log_time); 
	      if(isME11) { if(isME11odd) {M11_Od_N_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_N_hits_R->Fill(CSC_GlobalPoint_R);} }
	      CSC_EntryExit_N_Deposit_dz->Fill(log_deposit,CSCGlobalEntryExitDZ); CSC_EntryExit_N_Time_dz->Fill(CSCGlobalEntryExitDZ,log_time);
	      CSC_EntryExit_N_KinEn_dz->Fill(log_energy,CSCGlobalEntryExitDZ); CSC_EntryExit_N_Deposit_dR->Fill(log_deposit,CSCGlobalEntryExitDR); 
	      CSC_EntryExit_N_KinEn_dR->Fill(log_energy,CSCGlobalEntryExitDR); CSC_EntryExit_N_dz_dR->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR);
	      CSC_EntryExit_N_dz_dR_detail->Fill(CSCGlobalEntryExitDZ,CSCGlobalEntryExitDR); CSC_EntryExit_N_GapLength_Deposit->Fill(log_deposit,GapLength);
	      if(CSCGlobalEntryExitDZ < 0.05) {CSC_EntryExit_N_Deposit_pidR2->Fill(log_deposit,3.141592*pow(CSCGlobalEntryExitDR*1.0/2,2)); }
	      if((CSCGlobalEntryExitDZ > 0.65 && CSCGlobalEntryExitDZ < 0.75) || (CSCGlobalEntryExitDZ > 0.95 && CSCGlobalEntryExitDZ < 1.05)) 
		CSC_EntryExit_N_Deposit_dRdz->Fill(log_deposit,pow(pow(CSCGlobalEntryExitDR,2)+pow(CSCGlobalEntryExitDZ,2),0.5));
	      std::cout<<"CSC :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	      std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
	    }
	    else { // non-Nucleon ... catalogue as "Other Hadron" 
	      CSC_OH_hits->Fill(log_energy,log_time); CSC_OH_deposits->Fill(log_deposit,log_time); 
	      if(isME11) { if(isME11odd) {M11_Od_OH_hits_R->Fill(CSC_GlobalPoint_R);} if(isME11even) {M11_Ev_OH_hits_R->Fill(CSC_GlobalPoint_R);} }
	      std::cout<<"CSC :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	      std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;
	    }
	  break;
	}
	}
      }
      CSC_hits_tof->Fill(log_time);
      CSC_hits_eta->Fill(fabs(CSCGlobalPoint.eta()));
      CSC_hits_phi->Fill(CSCGlobalPoint.phi());
      CSC_hits_lin->Fill(time);
      if(rollId.station() == 1 && rollId.ring() == 1) { ME11_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 1 && rollId.ring() == 2) { ME12_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 1 && rollId.ring() == 3) { ME13_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 1 && rollId.ring() == 4) { ME11_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 2 && rollId.ring() == 1) { ME21_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 2 && rollId.ring() == 2) { ME22_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 3 && rollId.ring() == 1) { ME31_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 3 && rollId.ring() == 2) { ME32_hits_phi->Fill(CSCGlobalPoint.phi()); }
      if(rollId.station() == 4 && rollId.ring() == 1) { ME41_hits_phi->Fill(CSCGlobalPoint.phi()); /*ME4_hits_tof->Fill(log_time);*/}
      if(rollId.station() == 4 && rollId.ring() == 2) { ME42_hits_phi->Fill(CSCGlobalPoint.phi()); /*ME4_hits_tof->Fill(log_time);*/}

      // Count the amount of layers within a single chamber that have a hit
      // ------------------------------------------------------------------
      // Approach :: count the number of hits in a single chamber for each simtrack
      // Assumption :: each simtrack gives only rise to one hit / chamber
      // Cannot work because electron (pid = 11) simhits are due to hadron tracks (321, 2212, ...)
      // ------------------------------------------------------------------
      // Approach :: use ordering of SimHits in CSC SimHit containter
      // Cannot work, it is possible that 1 muon radiates an electron in each gas layer of CSC
      // ------------------------------------------------------------------
      // have to ask Tim Cox for smart idea
      // ------------------------------------------------------------------
      // int trackIdG = (*iHit).trackId();                     // GEANT TRACK ID
      // int trackIdC = geantToIndex.find(trackIdG)->second;   // CMSSW VECTOR ID
      // const SimTrack& track = theSimTracks[trackIdC]; 
      // ------------------------------------------------------------------
      if(abs(pid)==11) {
	CSC_el_HPL->Fill(layer);
	if(rollId.station() == 1) ME1_el_HPL->Fill(layer);      
	if(rollId.station() == 2) ME2_el_HPL->Fill(layer);      
	if(rollId.station() == 3) ME3_el_HPL->Fill(layer);      
	if(rollId.station() == 4) ME4_el_HPL->Fill(layer);      
      }
      else if(abs(pid)==13) {
	CSC_mu_HPL->Fill(layer);
	if(rollId.station() == 1) ME1_mu_HPL->Fill(layer);      
	if(rollId.station() == 2) ME2_mu_HPL->Fill(layer);      
	if(rollId.station() == 3) ME3_mu_HPL->Fill(layer);      
	if(rollId.station() == 4) ME4_mu_HPL->Fill(layer);      
      }
      else {
	CSC_ha_HPL->Fill(layer);
	if(rollId.station() == 1) ME1_ha_HPL->Fill(layer);      
	if(rollId.station() == 2) ME2_ha_HPL->Fill(layer);      
	if(rollId.station() == 3) ME3_ha_HPL->Fill(layer);      
	if(rollId.station() == 4) ME4_ha_HPL->Fill(layer);      
      } 
    }

    if(tech_debug) std::cout<<"DT Hits"<<std::endl;
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::DT){
      DTWireId wireId(theDetUnitId);
      GlobalPoint DTGlobalPoint = dtGeom->idToDet(wireId)->toGlobal((*iHit).localPosition());
      GlobalPoint DTGlobalEntry = dtGeom->idToDet(wireId)->toGlobal((*iHit).entryPoint());
      GlobalPoint DTGlobalExit  = dtGeom->idToDet(wireId)->toGlobal((*iHit).exitPoint());
      // double DTGlobalEntryExitDZ = fabs(DTGlobalEntry.z()-DTGlobalExit.z());
      double DTGlobalEntryExitDR = fabs(sqrt(pow(DTGlobalEntry.x(),2)+pow(DTGlobalEntry.y(),2))-sqrt(pow(DTGlobalExit.x(),2)+pow(DTGlobalExit.y(),2)));
      double DTLocalEntryExitDZ  = fabs((*iHit).entryPoint().z()-(*iHit).exitPoint().z());
      int layer = (wireId.superlayer()-1)*4+wireId.layer();

      // dtGeometry->idToDet(id)->surface().toGlobal(LocalPoint(point->displacement.x(), point->displacement.y(), point->displacement.z()));
      if (phys_debug) {
       std::cout<<"DT SimHit in "<<std::setw(24)<<wireId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<DTGlobalPoint.z();
       std::cout<<" | r = "<<std::setw(12)<<DTGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<DTGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<DTGlobalPoint.eta();
       std::cout<<" | global position = "<<DTGlobalPoint<<std::endl;
      }
      double DT_GlobalPoint_R = sqrt(pow(DTGlobalPoint.x(),2)+pow(DTGlobalPoint.y(),2));
      if(abs(pid)==11)      { DT_Electrons_SHPT->Fill(process);        DT_el_deps->Fill(log_deposit);  DT_el_kins->Fill(log_energy);  DT_el_tof->Fill(log_time);}
      else if(abs(pid)==13) { DT_Muons_SHPT->Fill(process);            DT_mu_deps->Fill(log_deposit);  DT_mu_kins->Fill(log_energy);  DT_mu_tof->Fill(log_time);}
      else                  { DT_Hadrons_SHPT->Fill(process);          DT_ha_deps->Fill(log_deposit);  DT_ha_kins->Fill(log_energy);  DT_ha_tof->Fill(log_time);}
      DT_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());             
      DT_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));  
      Muon_Barrel_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());
      Muon_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));  
      DT_EntryExit_All_Glob_dr->Fill(DTGlobalEntryExitDR);
      DT_EntryExit_All_Loc_dz->Fill(DTLocalEntryExitDZ);
      
      if((*iHit).timeOfFlight()<250) {
	if(abs(pid)==11)      { Muon_000ns_el_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_000ns_el_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); DT_Electrons_000ns_SHPT->Fill(process);}
	else if(abs(pid)==13) { Muon_000ns_mu_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_000ns_mu_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); }
	else                  { Muon_000ns_ha_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_000ns_ha_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); }
	if((*iHit).timeOfFlight()<50) {
	  if(abs(pid)==11)      { Muon_00ns_el_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_00ns_el_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());}
	  else if(abs(pid)==13) { Muon_00ns_mu_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_00ns_mu_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());}
	  else                  { Muon_00ns_ha_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_00ns_ha_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());}
	}
	if((*iHit).timeOfFlight()>50) {
	  if(abs(pid)==11)      { Muon_50ns_el_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_50ns_el_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); DT_Electrons_050ns_SHPT->Fill(process);}
	  else if(abs(pid)==13) { Muon_50ns_mu_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_50ns_mu_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());}
	  else                  { Muon_50ns_ha_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_50ns_ha_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());}
	}
	DT_000ns_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());
	DT_000ns_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));
	Muon_Barrel_000ns_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());
	Muon_000ns_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); 
      }
      if((*iHit).timeOfFlight()>250) {
	if(abs(pid)==11)      { Muon_250ns_el_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_250ns_el_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); DT_Electrons_250ns_SHPT->Fill(process);}
	else if(abs(pid)==13) { Muon_250ns_mu_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_250ns_mu_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); }
	else                  { Muon_250ns_ha_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); Muon_Barrel_250ns_ha_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y()); }
	DT_250ns_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());
	DT_250ns_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));
	Muon_Barrel_250ns_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());
	Muon_250ns_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R)); 
      }

      // int pid             = (*iHit).particleType();
      // int trackid         = (*iHit).trackId();
      // std::cout<<"TrackId = "<<trackid<<std::endl;

      // All Hits and HIP Hits
      DT_All_hits->Fill(log_energy,log_time); DT_All_deposits->Fill(log_deposit,log_time);
      if(HIP) {
	DT_HIP_hits->Fill(log_energy,log_time); DT_HIP_deposits->Fill(log_deposit,log_time);
      }
      // Catalog Nuclei
      if(abs(pid) > 999999999) {
	if(Z<21) DT_Nuclei_List->Fill(Z);
	else DT_Nuclei_List->Fill(21);
      }
      
      // std::cout<<"DT :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
      // std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
      if(abs(pid)==2212)      {DT_p_hits->Fill(log_energy,log_time); DT_p_deposits->Fill(log_deposit,log_time);}
      else if(abs(pid)==2112) {DT_n_hits->Fill(log_energy,log_time); DT_n_deposits->Fill(log_deposit,log_time);}
      else {
	switch (abs(pid)%1000) {
	  // leptons
	case 11:  DT_el_hits->Fill(log_energy,log_time); DT_el_deposits->Fill(log_deposit,log_time); break;
	case 13:  DT_mu_hits->Fill(log_energy,log_time); DT_mu_deposits->Fill(log_deposit,log_time); break;
	  // Pions
	case 111: DT_pi_hits->Fill(log_energy,log_time); DT_pi_deposits->Fill(log_deposit,log_time); break;
	case 211: DT_pi_hits->Fill(log_energy,log_time); DT_pi_deposits->Fill(log_deposit,log_time); break;
	case 130: DT_pi_hits->Fill(log_energy,log_time); DT_pi_deposits->Fill(log_deposit,log_time); break;
	  // Kaons
	case 310: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 311: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 321: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 313: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 323: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 315: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 325: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 317: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 327: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 319: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	case 329: DT_ka_hits->Fill(log_energy,log_time); DT_ka_deposits->Fill(log_deposit,log_time); break;
	  // Protons
	  // case 2212: DT_p_hits->Fill(log_energy,log_time); DT_p_deposits->Fill(log_deposit,log_time); break;
	  // Neutrons
	  // case 2112: DT_n_hits->Fill(log_energy,log_time); DT_n_deposits->Fill(log_deposit,log_time); break;
	  // Photons
	case 22:   DT_g_hits->Fill(log_energy,log_time); DT_g_deposits->Fill(log_deposit,log_time); break;
	  // Nucleons
	default:   {
	  if(abs(pid) > 999999999) { // Nucleons (10-digit numbers) 
	    DT_OH_hits->Fill(log_energy,log_time); DT_OH_deposits->Fill(log_deposit,log_time); 
	    std::cout<<"DT :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	    std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
	  }
	  else { // non-Nucleon ... catalogue as "Other Hadron" 
	    DT_OH_hits->Fill(log_energy,log_time); DT_OH_deposits->Fill(log_deposit,log_time); 
	    std::cout<<"DT :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	    std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;  
	  }
	  break;
	}
	}
      }
      DT_hits_tof->Fill(log_time);
      DT_hits_eta->Fill(fabs(DTGlobalPoint.eta()));
      DT_hits_phi->Fill(DTGlobalPoint.phi());
      DT_hits_lin->Fill(time);
      if(wireId.station() == 1) { MB1_hits_phi->Fill(DTGlobalPoint.phi()); }
      if(wireId.station() == 2) { MB2_hits_phi->Fill(DTGlobalPoint.phi()); }
      if(wireId.station() == 3) { MB3_hits_phi->Fill(DTGlobalPoint.phi()); }
      if(wireId.station() == 4) { MB4_hits_phi->Fill(DTGlobalPoint.phi()); /*MB4_hits_tof->Fill(log_time);*/ }

      if(abs(pid)==11) {
	DT_el_HPL->Fill(layer);
	if(wireId.station() == 1) MB1_el_HPL->Fill(layer);      
	if(wireId.station() == 2) MB2_el_HPL->Fill(layer);      
	if(wireId.station() == 3) MB3_el_HPL->Fill(layer);      
	if(wireId.station() == 4) MB4_el_HPL->Fill(layer);      
      }
      else if(abs(pid)==13) {
	DT_mu_HPL->Fill(layer);
	if(wireId.station() == 1) MB1_mu_HPL->Fill(layer);      
	if(wireId.station() == 2) MB2_mu_HPL->Fill(layer);      
	if(wireId.station() == 3) MB3_mu_HPL->Fill(layer);      
	if(wireId.station() == 4) MB4_mu_HPL->Fill(layer);      
      }
      else {
	DT_ha_HPL->Fill(layer);
	if(wireId.station() == 1) MB1_ha_HPL->Fill(layer);      
	if(wireId.station() == 2) MB2_ha_HPL->Fill(layer);      
	if(wireId.station() == 3) MB3_ha_HPL->Fill(layer);      
	if(wireId.station() == 4) MB4_ha_HPL->Fill(layer);      
      } 
    }
    // not in 71X
    // /*
    // if(tech_debug) std::cout<<"GEM Hits"<<std::endl;
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::GEM){ // Only GEMs


      GEMDetId rollId(theDetUnitId);
      const int station = rollId.station();
      const int etapart = rollId.roll();

      if(station!=0)
        {
	  if(tech_debug) std::cout<<"SimHit is a GEM SimHit"<<std::endl;

	  // GEM Geometry
	  // ============
	  const GEMEtaPartition* rollasociated = gemGeom->etaPartition(rollId);
	  const BoundPlane & GEMSurface = rollasociated->surface(); 
	  GlobalPoint GEMGlobalPoint = GEMSurface.toGlobal((*iHit).localPosition());
	  GlobalPoint GEMGlobalEntry = GEMSurface.toGlobal((*iHit).entryPoint());
	  GlobalPoint GEMGlobalExit  = GEMSurface.toGlobal((*iHit).exitPoint());
	  double GEMGlobalEntryExitDZ = fabs(GEMGlobalEntry.z()-GEMGlobalExit.z());
	  double GEMLocalEntryExitDZ  = fabs((*iHit).entryPoint().z()-(*iHit).exitPoint().z());
	  
	  bool isGE11L1 = false, isGE11L2 = false, isGE11odd = false, isGE11even = false;
	  if(rollId.station() == 1) {
	    if(rollId.layer() == 1) isGE11L1 = true;
	    if(rollId.layer() == 2) isGE11L2 = true;
	    if(rollId.chamber()%2 == 1) isGE11odd = true;
	    if(rollId.chamber()%2 == 0) isGE11even = true;
	  }

	  if(phys_debug) {
	    std::cout<<"GEM SimHit in "<<std::setw(12)<<(int)rollId<<std::setw(24)<<rollId;
	    std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<GEMGlobalPoint.z();
	    std::cout<<" | r = "<<std::setw(12)<<GEMGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<GEMGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<GEMGlobalPoint.eta();
	    std::cout<<" | global position = "<<GEMGlobalPoint<<std::endl;
	  }
	  double GEM_GlobalPoint_R = sqrt(pow(GEMGlobalPoint.x(),2)+pow(GEMGlobalPoint.y(),2));
	  if(abs(pid)==11)      { GEM_Electrons_SHPT->Fill(process);        GEM_el_deps->Fill(log_deposit);}
	  else if(abs(pid)==13) { GEM_Muons_SHPT->Fill(process);            GEM_mu_deps->Fill(log_deposit);}
	  else                  { GEM_Hadrons_SHPT->Fill(process);          GEM_ha_deps->Fill(log_deposit);}
	  GEM_XY->Fill(GEMGlobalPoint.x(), GEMGlobalPoint.y());             
	  GEM_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));  
	  Muon_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));  
	  GEM_EntryExit_All_Glob_dz->Fill(GEMGlobalEntryExitDZ);
	  GEM_EntryExit_All_Loc_dz->Fill(GEMLocalEntryExitDZ);
	  
	  if((*iHit).timeOfFlight()<250) {
	    if(abs(pid)==11)      { Muon_000ns_el_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    else if(abs(pid)==13) { Muon_000ns_mu_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    else                  { Muon_000ns_ha_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    
	    if((*iHit).timeOfFlight()<50) {
	      if(abs(pid)==11)      { Muon_00ns_el_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R)); GEM_Electrons_000ns_SHPT->Fill(process);}
	      else if(abs(pid)==13) { Muon_00ns_mu_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	      else                  { Muon_00ns_ha_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    }
	    
	    if((*iHit).timeOfFlight()>50) {
	      if(abs(pid)==11)      { Muon_50ns_el_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R)); GEM_Electrons_050ns_SHPT->Fill(process);}
	      else if(abs(pid)==13) { Muon_50ns_mu_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	      else                  { Muon_50ns_ha_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    }
	    GEM_000ns_XY->Fill(GEMGlobalPoint.x(), GEMGlobalPoint.y());
	    GEM_000ns_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));
	    Muon_000ns_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R)); 
	  }
	  if((*iHit).timeOfFlight()>250) {
	    if(abs(pid)==11)      { Muon_250ns_el_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R)); GEM_Electrons_250ns_SHPT->Fill(process);}
	    else if(abs(pid)==13) { Muon_250ns_mu_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    else                  { Muon_250ns_ha_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));}
	    GEM_250ns_XY->Fill(GEMGlobalPoint.x(), GEMGlobalPoint.y());
	    GEM_250ns_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));
	    Muon_250ns_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R)); 
	  }
	  
	  if(phys_debug) std::cout<<"SimHit from particle = "<<pid<<" created in process = "<<process<<" aka "<<proc[process]<<std::endl;
	  
	  // std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	  // std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
	  if(station==1) {

	    // All Hits and HIP Hits
	    GEM_All_hits->Fill(log_energy,log_time); GEM_All_deposits->Fill(log_deposit,log_time);
	    G11_All_hits_R->Fill(GEM_GlobalPoint_R); G11_All_hits_E->Fill(etapart);
	    if(isGE11odd) {G11_Od_All_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_All_hits_R->Fill(GEM_GlobalPoint_R);}
	    if(isGE11L1)  {G11_L1_All_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_All_hits_R->Fill(GEM_GlobalPoint_R);}
	    if(HIP) {
	      GEM_HIP_hits->Fill(log_energy,log_time); GEM_HIP_deposits->Fill(log_deposit,log_time);
	      G11_HIP_hits_R->Fill(GEM_GlobalPoint_R); G11_HIP_hits_E->Fill(etapart);
	      if(isGE11odd) {G11_Od_HIP_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_HIP_hits_R->Fill(GEM_GlobalPoint_R);}
	      if(isGE11L1)  {G11_L1_HIP_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_HIP_hits_R->Fill(GEM_GlobalPoint_R);}
	    }
	    // Catalog Nuclei
	    if(abs(pid) > 999999999) {
	      if(Z<21) GEM_Nuclei_List->Fill(Z);
	      else GEM_Nuclei_List->Fill(21);
	    }

	    if(abs(pid)==2212)      {
	      GEM_p_hits->Fill(log_energy,log_time); GEM_p_deposits->Fill(log_deposit,log_time); G11_p_hits_R->Fill(GEM_GlobalPoint_R); G11_p_hits_E->Fill(etapart);
	      if(isGE11odd) {G11_Od_p_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_p_hits_R->Fill(GEM_GlobalPoint_R);}
	      if(isGE11L1)  {G11_L1_p_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_p_hits_R->Fill(GEM_GlobalPoint_R);}
	    }
	    else if(abs(pid)==2112) {
	      GEM_n_hits->Fill(log_energy,log_time); GEM_n_deposits->Fill(log_deposit,log_time); G11_n_hits_R->Fill(GEM_GlobalPoint_R); G11_n_hits_E->Fill(etapart);
	      if(isGE11odd) {G11_Od_n_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_n_hits_R->Fill(GEM_GlobalPoint_R);}
	      if(isGE11L1)  {G11_L1_n_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_n_hits_R->Fill(GEM_GlobalPoint_R);}
	    }
	    else {
	      switch (abs(pid)%1000) {
		// leptons
	      case 11:  GEM_el_hits->Fill(log_energy,log_time); GEM_el_deposits->Fill(log_deposit,log_time); G11_el_hits_R->Fill(GEM_GlobalPoint_R); G11_el_hits_E->Fill(etapart);
 		if(isGE11odd) {G11_Od_el_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_el_hits_R->Fill(GEM_GlobalPoint_R);} 
		if(isGE11L1)  {G11_L1_el_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_el_hits_R->Fill(GEM_GlobalPoint_R);}
		break;
	      case 13:  GEM_mu_hits->Fill(log_energy,log_time); GEM_mu_deposits->Fill(log_deposit,log_time); G11_mu_hits_R->Fill(GEM_GlobalPoint_R); G11_mu_hits_E->Fill(etapart);
 		if(isGE11odd) {G11_Od_mu_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_mu_hits_R->Fill(GEM_GlobalPoint_R);} 
		if(isGE11L1)  {G11_L1_mu_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_mu_hits_R->Fill(GEM_GlobalPoint_R);}
		break;
		// Pions
	      case 111: 
	      case 211: 
	      case 130: 
		GEM_pi_hits->Fill(log_energy,log_time); GEM_pi_deposits->Fill(log_deposit,log_time); G11_pi_hits_R->Fill(GEM_GlobalPoint_R); G11_pi_hits_E->Fill(etapart); 
 		if(isGE11odd) {G11_Od_pi_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_pi_hits_R->Fill(GEM_GlobalPoint_R);} 
		if(isGE11L1)  {G11_L1_pi_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_pi_hits_R->Fill(GEM_GlobalPoint_R);}
		break;
		// Kaons
	      case 310:
	      case 311:
	      case 321:
	      case 313:
	      case 323:
	      case 315:
	      case 325:
	      case 317:
	      case 327:
	      case 319:
	      case 329: 
		GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G11_ka_hits_R->Fill(GEM_GlobalPoint_R); G11_ka_hits_E->Fill(etapart);
 		if(isGE11odd) {G11_Od_ka_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_ka_hits_R->Fill(GEM_GlobalPoint_R);} 
		if(isGE11L1)  {G11_L1_ka_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_ka_hits_R->Fill(GEM_GlobalPoint_R);}
		break;
		// Protons
		// case 2212: GEM_p_hits->Fill(log_energy,log_time); GEM_p_deposits->Fill(log_deposit,log_time); break;
		// Neutrons
		// case 2112: GEM_n_hits->Fill(log_energy,log_time); GEM_n_deposits->Fill(log_deposit,log_time); break;
		// Photons
	      case 22:   
		GEM_g_hits->Fill(log_energy,log_time); GEM_g_deposits->Fill(log_deposit,log_time); G11_g_hits_R->Fill(GEM_GlobalPoint_R); G11_g_hits_E->Fill(etapart); 
 		if(isGE11odd) {G11_Od_g_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_g_hits_R->Fill(GEM_GlobalPoint_R);} 
		if(isGE11L1)  {G11_L1_g_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_g_hits_R->Fill(GEM_GlobalPoint_R);}
		break;
		// Nucleons
	      default:   {
		if(abs(pid) > 999999999) { // Nucleons (10-digit numbers)
		  GEM_N_hits->Fill(log_energy,log_time); GEM_N_deposits->Fill(log_deposit,log_time); G11_N_hits_R->Fill(GEM_GlobalPoint_R); G11_N_hits_E->Fill(etapart);
		  if(isGE11odd) {G11_Od_N_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_N_hits_R->Fill(GEM_GlobalPoint_R);} 
		  if(isGE11L1)  {G11_L1_N_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_N_hits_R->Fill(GEM_GlobalPoint_R);}
		  std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
		  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
		}
		else { // non-Nucleon ... catalogue as "Other Hadron"
		  GEM_OH_hits->Fill(log_energy,log_time); GEM_OH_deposits->Fill(log_deposit,log_time); G11_OH_hits_R->Fill(GEM_GlobalPoint_R); G11_OH_hits_E->Fill(etapart);
		  if(isGE11odd) {G11_Od_OH_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11even) {G11_Ev_OH_hits_R->Fill(GEM_GlobalPoint_R);} 
		  if(isGE11L1)  {G11_L1_OH_hits_R->Fill(GEM_GlobalPoint_R);} if(isGE11L2)   {G11_L2_OH_hits_R->Fill(GEM_GlobalPoint_R);}
		  std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
		  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;
		}
		break;
	      }
	      } // end switch
	    }
	  } // end if station==1
	  if(station==2) {

	    // All Hits and HIP Hits
	    GEM_All_hits->Fill(log_energy,log_time); GEM_All_deposits->Fill(log_deposit,log_time);
	    G21_All_hits_R->Fill(GEM_GlobalPoint_R); G21_All_hits_E->Fill(etapart);
	    if(HIP) {
	      GEM_HIP_hits->Fill(log_energy,log_time); GEM_HIP_deposits->Fill(log_deposit,log_time);
	      G21_HIP_hits_R->Fill(GEM_GlobalPoint_R); G21_HIP_hits_E->Fill(etapart);
	    }

	    if(abs(pid)==2212)      {GEM_p_hits->Fill(log_energy,log_time); GEM_p_deposits->Fill(log_deposit,log_time); G21_p_hits_R->Fill(GEM_GlobalPoint_R); G21_p_hits_E->Fill(etapart);}
	    else if(abs(pid)==2112) {GEM_n_hits->Fill(log_energy,log_time); GEM_n_deposits->Fill(log_deposit,log_time); G21_n_hits_R->Fill(GEM_GlobalPoint_R); G21_n_hits_E->Fill(etapart);}
	    else {
	      switch (abs(pid)%1000) {
		// leptons
	      case 11:  GEM_el_hits->Fill(log_energy,log_time); GEM_el_deposits->Fill(log_deposit,log_time); G21_el_hits_R->Fill(GEM_GlobalPoint_R); G21_el_hits_E->Fill(etapart); break;
	      case 13:  GEM_mu_hits->Fill(log_energy,log_time); GEM_mu_deposits->Fill(log_deposit,log_time); G21_mu_hits_R->Fill(GEM_GlobalPoint_R); G21_mu_hits_E->Fill(etapart);break;
		// Pions
	      case 111: GEM_pi_hits->Fill(log_energy,log_time); GEM_pi_deposits->Fill(log_deposit,log_time); G21_pi_hits_R->Fill(GEM_GlobalPoint_R); G21_pi_hits_E->Fill(etapart); break;
	      case 211: GEM_pi_hits->Fill(log_energy,log_time); GEM_pi_deposits->Fill(log_deposit,log_time); G21_pi_hits_R->Fill(GEM_GlobalPoint_R); G21_pi_hits_E->Fill(etapart); break;
	      case 130: GEM_pi_hits->Fill(log_energy,log_time); GEM_pi_deposits->Fill(log_deposit,log_time); G21_pi_hits_R->Fill(GEM_GlobalPoint_R); G21_pi_hits_E->Fill(etapart); break;
		// Kaons
	      case 310: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 311: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 321: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 313: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 323: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 315: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 325: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 317: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 327: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 319: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
	      case 329: GEM_ka_hits->Fill(log_energy,log_time); GEM_ka_deposits->Fill(log_deposit,log_time); G21_ka_hits_R->Fill(GEM_GlobalPoint_R); G21_ka_hits_E->Fill(etapart); break;
		// Protons
		// case 2212: GEM_p_hits->Fill(log_energy,log_time); GEM_p_deposits->Fill(log_deposit,log_time); break;
		// Neutrons
		// case 2112: GEM_n_hits->Fill(log_energy,log_time); GEM_n_deposits->Fill(log_deposit,log_time); break;
		// Photons
	      case 22:   GEM_g_hits->Fill(log_energy,log_time); GEM_g_deposits->Fill(log_deposit,log_time); G21_g_hits_R->Fill(GEM_GlobalPoint_R); G21_g_hits_E->Fill(etapart); break;
		// Nucleons
	      default:   {
		if(abs(pid) > 999999999) { // Nucleons (10-digit numbers)
		  GEM_N_hits->Fill(log_energy,log_time); GEM_N_deposits->Fill(log_deposit,log_time); G21_N_hits_R->Fill(GEM_GlobalPoint_R); G21_N_hits_E->Fill(etapart); 
		  std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
		  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
		}
		else { // non-Nucleon ... catalogue as "Other Hadron"
		  GEM_OH_hits->Fill(log_energy,log_time); GEM_OH_deposits->Fill(log_deposit,log_time); G21_OH_hits_R->Fill(GEM_GlobalPoint_R); G21_OH_hits_E->Fill(etapart); 
		  std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
		  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;
		}
		break;
	      }
	      } // end switch
	    }
	  } // end if station==2

	  GEM_hits_tof->Fill(log_time);
	  GEM_hits_eta->Fill(fabs(GEMGlobalPoint.eta()));
	  GEM_hits_phi->Fill(GEMGlobalPoint.phi());
	  GEM_hits_lin->Fill(time);     
	} // end station != 0

    // */
    // not in 71X 
    // /*
    // if(tech_debug) std::cout<<"ME0 Hits"<<std::endl;
    // if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::ME0){ // Only ME0

      if(station==0)
        {
          if(tech_debug) std::cout<<"SimHit is a ME0 SimHit"<<std::endl;

	  // ME0 Geometry
	  // ============
	  // ME0DetId rollId(theDetUnitId);
	  // const ME0EtaPartition* rollasociated = me0Geom->etaPartition(rollId);
	  const GEMEtaPartition* rollasociated = gemGeom->etaPartition(rollId);
	  const int etapart = rollId.roll();
	  const BoundPlane & ME0Surface = rollasociated->surface(); 
	  GlobalPoint ME0GlobalPoint = ME0Surface.toGlobal((*iHit).localPosition());
	  GlobalPoint ME0GlobalEntry = ME0Surface.toGlobal((*iHit).entryPoint());
	  GlobalPoint ME0GlobalExit  = ME0Surface.toGlobal((*iHit).exitPoint());
	  double ME0GlobalEntryExitDZ = fabs(ME0GlobalEntry.z()-ME0GlobalExit.z());
	  double ME0LocalEntryExitDZ  = fabs((*iHit).entryPoint().z()-(*iHit).exitPoint().z());

	  if(phys_debug) {
	    std::cout<<"ME0 SimHit in "<<std::setw(12)<<(int)rollId<<std::setw(24)<<rollId;
	    std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<ME0GlobalPoint.z();
	    std::cout<<" | r = "<<std::setw(12)<<ME0GlobalPoint.mag()<<" | phi = "<<std::setw(12)<<ME0GlobalPoint.phi()<<" | eta = "<<std::setw(12)<<ME0GlobalPoint.eta();
	    std::cout<<" | global position = "<<ME0GlobalPoint<<std::endl;
	  }
	  double ME0_GlobalPoint_R = sqrt(pow(ME0GlobalPoint.x(),2)+pow(ME0GlobalPoint.y(),2));
	  if(abs(pid)==11)      { ME0_Electrons_SHPT->Fill(process);        ME0_el_deps->Fill(log_deposit);}
	  else if(abs(pid)==13) { ME0_Muons_SHPT->Fill(process);            ME0_mu_deps->Fill(log_deposit);}
	  else                  { ME0_Hadrons_SHPT->Fill(process);          ME0_ha_deps->Fill(log_deposit);}
	  ME0_XY->Fill(ME0GlobalPoint.x(), ME0GlobalPoint.y());             
	  ME0_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));  
	  Muon_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));  
	  ME0_EntryExit_All_Glob_dz->Fill(ME0GlobalEntryExitDZ);
	  ME0_EntryExit_All_Loc_dz->Fill(ME0LocalEntryExitDZ);
	  
	  if((*iHit).timeOfFlight()<250) {
	    if(abs(pid)==11)      { Muon_000ns_el_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    else if(abs(pid)==13) { Muon_000ns_mu_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    else                  { Muon_000ns_ha_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    
	    if((*iHit).timeOfFlight()<50) {
	      if(abs(pid)==11)      { Muon_00ns_el_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R)); ME0_Electrons_000ns_SHPT->Fill(process);}
	      else if(abs(pid)==13) { Muon_00ns_mu_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	      else                  { Muon_00ns_ha_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    }
	    
	    if((*iHit).timeOfFlight()>50) {
	      if(abs(pid)==11)      { Muon_50ns_el_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R)); ME0_Electrons_050ns_SHPT->Fill(process);}
	      else if(abs(pid)==13) { Muon_50ns_mu_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	      else                  { Muon_50ns_ha_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    }
	    ME0_000ns_XY->Fill(ME0GlobalPoint.x(), ME0GlobalPoint.y());
	    ME0_000ns_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));
	    Muon_000ns_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R)); 
	  }
	  if((*iHit).timeOfFlight()>250) {
	    if(abs(pid)==11)      { Muon_250ns_el_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R)); ME0_Electrons_250ns_SHPT->Fill(process);}
	    else if(abs(pid)==13) { Muon_250ns_mu_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    else                  { Muon_250ns_ha_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));}
	    ME0_250ns_XY->Fill(ME0GlobalPoint.x(), ME0GlobalPoint.y());
	    ME0_250ns_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));
	    Muon_250ns_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R)); 
	  }
	  
	  // All Hits and HIP Hits
	  ME0_All_hits->Fill(log_energy,log_time); ME0_All_deposits->Fill(log_deposit,log_time);
	  ME0_All_hits_R->Fill(ME0_GlobalPoint_R); ME0_All_hits_E->Fill(etapart);
	  if(HIP) {
	    ME0_HIP_hits->Fill(log_energy,log_time); ME0_HIP_deposits->Fill(log_deposit,log_time);
	    ME0_HIP_hits_R->Fill(ME0_GlobalPoint_R); ME0_HIP_hits_E->Fill(etapart);
	  }
	  // Catalog Nuclei
	  if(abs(pid) > 999999999) {
	    if(Z<21) ME0_Nuclei_List->Fill(Z);
	    else ME0_Nuclei_List->Fill(21);
	  }

	  if(phys_debug) std::cout<<"SimHit from particle = "<<pid<<" created in process = "<<process<<" aka "<<proc[process]<<std::endl;
	  
	  // std::cout<<"ME0 :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	  // std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
	  if(abs(pid)==2212)      {ME0_p_hits->Fill(log_energy,log_time); ME0_p_deposits->Fill(log_deposit,log_time); ME0_p_hits_R->Fill(ME0_GlobalPoint_R); ME0_p_hits_E->Fill(etapart);}
	  else if(abs(pid)==2112) {ME0_n_hits->Fill(log_energy,log_time); ME0_n_deposits->Fill(log_deposit,log_time); ME0_n_hits_R->Fill(ME0_GlobalPoint_R); ME0_n_hits_E->Fill(etapart);}
	  else {
	    switch (abs(pid)%1000) {
	      // leptons
	    case 11:  ME0_el_hits->Fill(log_energy,log_time); ME0_el_deposits->Fill(log_deposit,log_time); ME0_el_hits_R->Fill(ME0_GlobalPoint_R); ME0_el_hits_E->Fill(etapart); break;
	    case 13:  ME0_mu_hits->Fill(log_energy,log_time); ME0_mu_deposits->Fill(log_deposit,log_time); ME0_mu_hits_R->Fill(ME0_GlobalPoint_R); ME0_mu_hits_E->Fill(etapart); break;
	      // Pions
	    case 111: ME0_pi_hits->Fill(log_energy,log_time); ME0_pi_deposits->Fill(log_deposit,log_time); ME0_pi_hits_R->Fill(ME0_GlobalPoint_R); ME0_pi_hits_E->Fill(etapart); break;
	    case 211: ME0_pi_hits->Fill(log_energy,log_time); ME0_pi_deposits->Fill(log_deposit,log_time); ME0_pi_hits_R->Fill(ME0_GlobalPoint_R); ME0_pi_hits_E->Fill(etapart); break;
	    case 130: ME0_pi_hits->Fill(log_energy,log_time); ME0_pi_deposits->Fill(log_deposit,log_time); ME0_pi_hits_R->Fill(ME0_GlobalPoint_R); ME0_pi_hits_E->Fill(etapart); break;
	      // Kaons
	    case 310: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 311: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 321: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 313: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 323: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 315: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 325: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 317: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 327: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 319: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	    case 329: ME0_ka_hits->Fill(log_energy,log_time); ME0_ka_deposits->Fill(log_deposit,log_time); ME0_ka_hits_R->Fill(ME0_GlobalPoint_R); ME0_ka_hits_E->Fill(etapart); break;
	      // Protons
	      // case 2212: ME0_p_hits->Fill(log_energy,log_time); ME0_p_deposits->Fill(log_deposit,log_time); break;
	      // Neutrons
	      // case 2112: ME0_n_hits->Fill(log_energy,log_time); ME0_n_deposits->Fill(log_deposit,log_time); break;
	      // Photons
	    case 22:   ME0_g_hits->Fill(log_energy,log_time); ME0_g_deposits->Fill(log_deposit,log_time); ME0_g_hits_R->Fill(ME0_GlobalPoint_R); ME0_g_hits_E->Fill(etapart); break;
	      // Nucleons
	    default:   {
		if(abs(pid) > 999999999) { // Nucleons (10-digit numbers)
		  ME0_N_hits->Fill(log_energy,log_time); ME0_N_deposits->Fill(log_deposit,log_time); ME0_N_hits_R->Fill(ME0_GlobalPoint_R); ME0_N_hits_E->Fill(etapart);
		  std::cout<<"ME0 :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
		  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as NUCLEI"<<std::endl;
		}
		else { // non-Nucleon ... catalogue as "Other Hadron"
		  ME0_OH_hits->Fill(log_energy,log_time); ME0_OH_deposits->Fill(log_deposit,log_time); ME0_OH_hits_R->Fill(ME0_GlobalPoint_R); ME0_OH_hits_E->Fill(etapart);
		  std::cout<<"ME0 :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
		  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV] catalogued as OTHER HADRON"<<std::endl;
		}
	      break;
	    }
	    }
	  }
	  ME0_hits_tof->Fill(log_time);
	  ME0_hits_eta->Fill(fabs(ME0GlobalPoint.eta()));
	  ME0_hits_phi->Fill(ME0GlobalPoint.phi());
	  ME0_hits_lin->Fill(time);     
	} // end if station==0
    } // end GEMHits
      // */
    if(tech_debug) std::cout<<"END LOOP"<<std::endl;
  } // end SimHit loop
  if(tech_debug) std::cout<<"End loop over SimHit collections"<<std::endl;
}

double
MyME11SimHitAnalyzer::getLumi(int pu, int space, int com)
{
  double rev   = 11245.0;
  double scale = 1.0e34;
  std::map<int,double> cross;
  std::map<int,double> bunches;

  cross[7]=7.31e-26;  // TOTEM Measurement
  cross[8]=7.47e-26;  // TOTEM Measurement
  cross[13]=8.0e-26;

  bunches[25]=2808.0;
  bunches[50]=1380.0;
 
  double instlumi = pu*bunches[space]*rev/cross[com]/scale; 
  return instlumi;
}

double MyME11SimHitAnalyzer::getPU(double lumi, int space, int com)
{
  double rev   = 11245.0;
  double scale = 1.0e34;
  std::map<int,double> cross;
  std::map<int,double> bunches;

  cross[7]=7.31e-26;  // TOTEM Measurement
  cross[8]=7.47e-26;  // TOTEM Measurement
  cross[13]=8.0e-26;

  bunches[25]=2808.0;
  bunches[50]=1380.0;

  double pu = lumi/bunches[space]/rev*cross[com]*scale;
  return pu;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyME11SimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyME11SimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyME11SimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  if(tech_debug) std::cout<<"[MyME11SimHitAnalyzer :: BeginRun]"<<std::endl;

  const edm::ESTransientHandle<RPCGeometry> rpcGeom = iSetup.getTransientHandle(rpcGeom_Token);

  // iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  // iSetup.get<MuonGeometryRecord>().get(cscGeom);
  // iSetup.get<MuonGeometryRecord>().get(dtGeom);
  // not in 71X release
  // iSetup.get<MuonGeometryRecord>().get(gemGeom);
  // iSetup.get<MuonGeometryRecord>().get(me0Geom);

  // Loop over RPC Geometry
  // std::cout <<"Analyze RPC Geometry :: Loop over RPC Chambers"<<std::endl;
  // for (TrackingGeometry::DetContainer::const_iterator it=rpcGeom->dets().begin(); it<rpcGeom->dets().end(); ++it) {
  //    if( dynamic_cast< RPCChamber* const >( *it ) != 0 ){
  //      RPCChamber* ch = dynamic_cast< RPCChamber* const >( *it );
  //      std::vector< const RPCRoll*> rolls = (ch->rolls());
  //      for(std::vector<const RPCRoll*>::const_iterator r = rolls.begin();r != rolls.end(); ++r) {

  // 71X :: use const
  for(std::vector<const RPCRoll*>::const_iterator it = rpcGeom->rolls().begin(); it != rpcGeom->rolls().end(); ++it) {
    const RPCRoll* roll = (*it);
  // 62X, 70X :: do not use const
  // for(std::vector<RPCRoll*>::const_iterator it = rpcGeom->rolls().begin(); it != rpcGeom->rolls().end(); ++it) {
  // RPCRoll* roll = (*it);
    RPCDetId rpcId = roll->id();
    int n_strips   = roll->nstrips();

    int region = abs(rpcId.region());
    int station = abs(rpcId.station())-1;
    int ring = (region==0)? abs(rpcId.ring()) : abs(rpcId.ring())-1; // region == 0 ? Barrel : Endcap
    
    if (region == 0) {
      const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(roll->topology()));
      float stripl = top_->stripLength();
      float stripw = top_->pitch();
      rpc_barrel_area += stripw*stripl*n_strips;
      RPC_area_array[region][station][ring] += stripw*stripl*n_strips;
      // std::cout<<(int) rpcId<<" aka "<<rpcId<<" = "<< stripw*stripl*n_strips <<" cm2";
      // std::cout<<"==> RPC_area_array[region="<<region<<"][station="<<station+1<<"][ring="<<ring<<"] = "<<RPC_area_array[region][station][ring]<<" cm2"<<std::endl;
    }
    if(region == 1) {
      const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&(roll->topology()));
      float stripl = top_->stripLength();
      float stripw = top_->pitch();
      rpc_endcap_area += stripw*stripl*n_strips;
      RPC_area_array[region][station][ring] += stripw*stripl*n_strips;
      // std::cout<<(int) rpcId<<" aka "<<rpcId<<" = "<< stripw*stripl*n_strips <<" cm2";
      // std::cout<<" ==> RPC_area_array[region="<<region<<"][station="<<station+1<<"][ring"<<ring+1<<"] = "<<RPC_area_array[region][station][ring]<<" cm2"<<std::endl;
    }
  }
  // std::cout<<"RPC Barrel Active area = "<<rpc_barrel_area<<" cm2 || RPC Endcap Active area = "<<rpc_endcap_area<<" cm2"<<std::endl;

  // for (TrackingGeometry::DetContainer::const_iterator it=cscGeom->dets().begin(); it<cscGeom->dets().end(); ++it) {
  // if( dynamic_cast< CSCChamber* >( *it ) != 0 ){
  // CSCChamber* ch = dynamic_cast< CSCChamber* >( *it );
  // std::vector< const CSCRoll*> rolls = (ch->rolls());
  // for(std::vector<const RPCRoll*>::const_iterator r = rolls.begin();r != rolls.end(); ++r) {
  // RPCDetId rpcId = (*r)->id();
  // }
  // }

  /*
  const std::vector<CSCChamber*> ChamberContainer = cscGeom->chambers();
  for(unsigned int nCh=0; nCh<ChamberContainer.size(); ++nCh){
    const CSCChamber *cscchamber = ChamberContainer[nCh];
    std::cout<<"CSC ID = "<<cscchamber->id()<<" in Endcap "<<cscchamber->id().endcap()<<" in Station "<<cscchamber->id().station()<<" and Ring "<<cscchamber->id().ring()<<std::endl;
    const std::vector< const CSCLayer* > LayerContainer = cscchamber->layers();
  }
  */


  // CSC AREA First TRY
  // ------------------

  double area_me11 = 0.0, area_me12 = 0.0, area_me13 = 0.0;
  double area_me21 = 0.0, area_me31 = 0.0, area_me41 = 0.0;
  double area_me22 = 0.0, area_me32 = 0.0, area_me42 = 0.0;

  double dy = 10.0;
  int j_me11 = 0, j_me12 = 0, j_me13 = 0; // local counter 
  int j_me21 = 0, j_me31 = 0, j_me41 = 0; // local counter 
  int j_meX2 = 0; 

  for(int i=0; i<r_CSC; ++i) {

    // From mf.xml :: <Tubs name="ME11Space" rMin="990*mm" rMax="2725*mm" 
    // From mf.xml :: <Tubs name="ME12Space" rMin="2.734*m  " rMax="4.685*m
    // From mf.xml :: <Tubs name="ME13Space" rMin="5.055*m  " rMax="6.95*m

    // From data above I would assume/round ::
    // ME11 :: from 1200 mm to 2700 mm 
    // ME12 :: from 2800 mm to 4500 mm
    // ME13 :: from 5000 mm to 6600 mm
    // Instead from CSC RecHits I get ::
    // ME11 :: from 1000 mm to 2600 mm  ==> 16 bins of 10 cm
    // ME12 :: from 2800 mm to 4600 mm  ==> 18 bins of 10 cm
    // ME13 :: from 5000 mm to 6600 mm  ==> 16 bins of 10 cm


    // ME1/1
    // 150 cm long chamber split in bins of 10cm so 15 bins
    // From mf.xml :: ME1a_ActiveGasVol :: dz="21.56*cm" dx1="9.35*cm" dx2="13.365*cm"
    // From mf.xml :: ME11_ActiveGasVol :: dz="52.81*cm" dx1="13.365*cm" dx2="23.31*cm"
    // ==> lower trapezoid  = 2 *  9.35 = 18.70
    // ==> higher trapezoid = 2 * 23.31 = 46.61  
    // chamberlength = 2*52.81 + 2 * 21.56 = 109.62 + 43.12 = 152.74
    // --> check ME1/1 in CMSSW71X
    // double dx_me11 = (46.61-18.70)/16; double x_me11 = 18.70;
    // Tim's numbers ??? why are they different ??? --> new gas volume of ME11, made a little bit wider
    // double dx = (47.40-19.14)/15; double x_me11 = 19.14;
    double dx_me11 = (47.40-19.14)/16; double x_me11 = 19.14;
    if(r1_CSC+i*10 > 99.99 && r1_CSC+(i+1)*10 < 260.01) {
      double top = x_me11 + j_me11*dx_me11;
      double bottom = x_me11 + (j_me11+1)*dx_me11;
      double area = 2*36*0.5*(top + bottom)*dy; area_me11 += area;
      // std::cout<<"CSC ME1/1 Area :: i = "<<i<<" j = "<<j_me11<<" R_min = "<<r1_CSC+i*10<<" cm and R_max = "<<r1_CSC+(i+1)*10<<" cm ";
      // std::cout<<" || x1 = "<<top<<"cm | x2 = "<<bottom<<" cm | Area of little trapezoid = "<<area<<" cm^2"<<std::endl;
      CSC_ME1_area->SetBinContent(i+1, area);
      ++j_me11;
    }

    // ME1/2
    // From mf.xml :: ME12_ActiveGasVol :: dz="87.245*cm dx1="25.5*cm " dx2="41.87*cm
    // chamberlength = 174.490
    double dx_me12 = (83.74-51.00)/18; double x_me12 = 51.00;
    if(r1_CSC+i*10 > 279.99 && r1_CSC+(i+1)*10 < 460.01) {
      double top = x_me12 + j_me12*dx_me12;
      double bottom = x_me12 + (j_me12+1)*dx_me12;
      double area = 2*36*0.5*(top + bottom)*dy; area_me12 += area;
      CSC_ME1_area->SetBinContent(i+1, area);
      ++j_me12;
    }
    // ME1/3
    // from mf.xml :: ME13_ActiveGasVol :: dz="82.08*cm dx1="31.70*cm " dx2="46.05*cm
    // chamberlength = 164.16
    double dx_me13 = (92.10-63.40)/16; double x_me13 = 63.40;
    if(r1_CSC+i*10 > 499.99 && r1_CSC+(i+1)*10 < 660.01) {
      double top = x_me13 + j_me13*dx_me13;
      double bottom = x_me13 + (j_me13+1)*dx_me13;
      double area = 2*36*0.5*(top + bottom)*dy; area_me13 += area;
      CSC_ME1_area->SetBinContent(i+1, area);
      ++j_me13;
    }
  }


  for(int i=0; i<r_CSC; ++i) {

    // ME2/1 + ME3/1 + ME4/1
    // <Trd1 name="ME21_ActiveGasVol" dz="94.83*cm " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="27.00*cm " dx2="62.855*cm "/>
    // <Trd1 name="ME31_ActiveGasVol" dz="84.85*cm " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="30.70*cm " dx2="62.855*cm "/>
    // <Trd1 name="ME41_ActiveGasVol" dz="74.710*cm " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="34.505*cm " dx2="62.825*cm "/>  ??? typo should be 62.855 ???
    
    // chamberlength = 189.66 --- 169.70 --- 149.42 ==> 19 - 17 - 15 bins
    double dx_me21 = (125.71-54.00)/21; double x_me21 = 54.00;
    if(r1_CSC+i*10 > 139.99 && r1_CSC+(i+1)*10 < 340.01) {
      double top = x_me21 + j_me21*dx_me21;
      double bottom = x_me21 + (j_me21+1)*dx_me21;
      double area = 2*18*0.5*(top + bottom)*dy; area_me21 += area;
      CSC_ME2_area->SetBinContent(i+1, area);
      ++j_me21;
    }
    double dx_me31 = (125.71-61.40)/19; double x_me31 = 61.40;
    if(r1_CSC+i*10 > 159.99 && r1_CSC+(i+1)*10 < 340.01) {
      double top = x_me31 + j_me31*dx_me31;
      double bottom = x_me31 + (j_me31+1)*dx_me31;
      double area = 2*18*0.5*(top + bottom)*dy; area_me31 += area;
      CSC_ME3_area->SetBinContent(i+1, area);
      ++j_me31;
    }
    double dx_me41 = (125.71-69.01)/17; double x_me41 = 69.01;
    if(r1_CSC+i*10 > 179.99 && r1_CSC+(i+1)*10 < 340.01) {
      double top = x_me41 + j_me41*dx_me41;
      double bottom = x_me41 + (j_me41+1)*dx_me41;
      double area = 2*18*0.5*(top + bottom)*dy; area_me41 += area;
      CSC_ME4_area->SetBinContent(i+1, area);
      ++j_me41;
    }
    
    // ME2/2 + ME3/2 + ME4/2
    // <Trd1 name="ME22_ActiveGasVol" dz="1.6153*m  " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="33.230*cm" dx2="63.575*cm"/>
    // <Trd1 name="ME32_ActiveGasVol" dz="1.6153*m  " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="33.230*cm " dx2="63.575*cm "/>
    // <Trd1 name="ME42_ActiveGasVol" dz="1.6153*m  " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="33.230*cm " dx2="63.575*cm "/>
    
    // chamberlength = 323.06 ==> 32 bins
    double dx_meX2 = (127.15-66.46)/32; double x_meX2 = 66.46;
    if(r1_CSC+i*10 > 359.99 && r1_CSC+(i+1)*10 < 680.01) {
      double top = x_meX2 + j_meX2*dx_meX2;
      double bottom = x_meX2 + (j_meX2+1)*dx_meX2;
      double area = 2*36*0.5*(top + bottom)*dy; area_me22 += area; area_me32 += area; area_me42 += area;
      CSC_ME2_area->SetBinContent(i+1, area);
      CSC_ME3_area->SetBinContent(i+1, area);
      CSC_ME4_area->SetBinContent(i+1, area);
      ++j_meX2;
    }
  }
  // std::cout<<"--- AREA of CSC system ---"<<std::endl;
  // std::cout<<"ME1/1 :: "<<area_me11/36/2<<" cm^2 / chamber ==> "<<area_me11<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME1/2 :: "<<area_me12/36/2<<" cm^2 / chamber ==> "<<area_me12<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME1/3 :: "<<area_me13/36/2<<" cm^2 / chamber ==> "<<area_me13<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME2/1 :: "<<area_me21/36/2<<" cm^2 / chamber ==> "<<area_me21<<" cm^2 in total [2 endcaps x 18 chambers]"<<std::endl;
  // std::cout<<"ME2/2 :: "<<area_me22/36/2<<" cm^2 / chamber ==> "<<area_me22<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME3/1 :: "<<area_me31/36/2<<" cm^2 / chamber ==> "<<area_me31<<" cm^2 in total [2 endcaps x 18 chambers]"<<std::endl;
  // std::cout<<"ME3/2 :: "<<area_me32/36/2<<" cm^2 / chamber ==> "<<area_me32<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME4/1 :: "<<area_me41/36/2<<" cm^2 / chamber ==> "<<area_me41<<" cm^2 in total [2 endcaps x 18 chambers]"<<std::endl;
  // std::cout<<"ME4/2 :: "<<area_me42/36/2<<" cm^2 / chamber ==> "<<area_me42<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;

  // CSC AREA SECOND TRY
  // -------------------
  // |      |  ME11 |  ME12 | ME13  |  ME21 |  ME31 |  ME41 | ME2x  |
  // |------|-------|-------|-------|-------|-------|-------|-------|
  // | R_1  | 105.5 | 281.5 | 512.0 | 146.9 | 166.9 | 187.0 | 364.0 |
  // | R_2  | 257.5 | 456.0 | 676.2 | 336.6 | 336.6 | 336.4 | 687.1 |
  // |------|-------|-------|-------|-------|-------|-------|-------|
  // | Diff | 152.0 | 174.5 | 164.2 | 189.7 | 169.7 | 149.6 | 323.1 |
  // |------|-------|-------|-------|-------|-------|-------|-------|
  // | Bins |    15 |    17 |    16 |    19 |    17 |    15 |    32 | 
  // |   No | 02-16 | 18-34 | 36-51 | 02-20 | 02-20 | 02-20 | 22-53 |
  // |------|-------|-------|-------|-------|-------|-------|-------|
  // | count| 01-15 | 17-33 | 35-50 | 01-19 | 01-19 | 01-19 | 21-52 |
  // |------|-------|-------|-------|-------|-------|-------|-------|

  // reset values
  area_me11 = 0.0; area_me12 = 0.0; area_me13 = 0.0;
  area_me21 = 0.0; area_me31 = 0.0; area_me41 = 0.0;
  area_me22 = 0.0; area_me32 = 0.0; area_me42 = 0.0;

  // the area of the last strip of each chamber is much larger
  // because of the longer length ... this effect is even magnified
  // because of the multiplication with the number of chambers ...
  // redoing the calculations it actually is fine
  // reducing the length to standard 10 cm instead of 15 shows value
  // which is in line with the value of the one before last bin 

  std::cout<<"----------------------------------------------------------------------"<<std::endl;
  for(int i=0; i<51; ++i) {
    if(i==0) {
      CSC_Geom_ME1_area->SetBinContent(i+1, 0); // Fill Bin 1 
    }
    // ME1/1 :: Bins 02--16 :: Counter 1--15
    // 
    // Correspondence :: 19.14 at R_1 = 105.5
    // Correspondence :: 47.40 at R_2 = 257.5
    double dx_me11 = (47.40-19.14)/152.0; double x_me11 = 19.14;
    if(i > 0 && i< 16) {
      double top = x_me11 + (r_ME1_geom[i]  -r_ME1_geom[1])*dx_me11;
      double bottom = top + (r_ME1_geom[i+1]-r_ME1_geom[i])*dx_me11;
      double area = 2*36*0.5*(top + bottom)*(r_ME1_geom[i+1]-r_ME1_geom[i]); area_me11 += area;
      // std::cout<<"CSC ME1/1 Area :: i = "<<i<<" R_min = "<<r_ME1_geom[i]<<" cm and R_max = "<<r_ME1_geom[i+1]<<" cm ==> dR = "<<(r_ME1_geom[i+1]-r_ME1_geom[i])<<" cm";
      // std::cout<<" || x1 = "<<top<<"cm | x2 = "<<bottom<<" cm | Area of little trapezoid = "<<area<<" cm^2"<<std::endl;
      CSC_Geom_ME1_area->SetBinContent(i+1, area);
    }
    if(i==16) {
      CSC_Geom_ME1_area->SetBinContent(i+1, 0); // Fill Bin 17
    }
    // ME1/2 :: Bins 18--34 :: Counter 17--33
    // From mf.xml :: ME12_ActiveGasVol :: dz="87.245*cm dx1="25.5*cm " dx2="41.87*cm
    // chamberlength = 174.490
    double dx_me12 = (83.74-51.00)/174.5; double x_me12 = 51.00;
    if(i > 16 && i< 34) {
      double top = x_me12 + (r_ME1_geom[i]  -r_ME1_geom[17])*dx_me12;
      double bottom = top + (r_ME1_geom[i+1]-r_ME1_geom[i])*dx_me12;
      double area = 2*36*0.5*(top + bottom)*(r_ME1_geom[i+1]-r_ME1_geom[i]); area_me12 += area;
      // std::cout<<"CSC ME1/2 Area :: i = "<<i<<" R_min = "<<r_ME1_geom[i]<<" cm and R_max = "<<r_ME1_geom[i+1]<<" cm ==> dR = "<<(r_ME1_geom[i+1]-r_ME1_geom[i])<<" cm";
      // std::cout<<" || x1 = "<<top<<"cm | x2 = "<<bottom<<" cm | Area of little trapezoid = "<<area<<" cm^2"<<std::endl;
      CSC_Geom_ME1_area->SetBinContent(i+1, area);
    }
    if(i==34) {
      CSC_Geom_ME1_area->SetBinContent(i+1, 0); // Fill Bin 17
    }
    // ME1/3 :: Bins 36--51 :: Counter 35--50
    // from mf.xml :: ME13_ActiveGasVol :: dz="82.08*cm dx1="31.70*cm " dx2="46.05*cm
    // chamberlength = 164.16
    double dx_me13 = (92.10-63.40)/164.2; double x_me13 = 63.40;
    if(i>34 && i<51) {
      double top = x_me13 + (r_ME1_geom[i]  -r_ME1_geom[35])*dx_me13;
      double bottom = top + (r_ME1_geom[i+1]-r_ME1_geom[i])*dx_me13;
      double area = 2*36*0.5*(top + bottom)*(r_ME1_geom[i+1]-r_ME1_geom[i]); area_me13 += area;
      CSC_Geom_ME1_area->SetBinContent(i+1, area);
    }
  }

  // std::cout<<"----------------------------------------------------------------------"<<std::endl;
  for(int i=0; i<53; ++i) {
    if(i==0) {
      CSC_Geom_ME2_area->SetBinContent(i+1, 0); // Fill Bin 1
    }
    if(i < 3) {
      CSC_Geom_ME3_area->SetBinContent(i+1, 0); // Fill Bin 1-3
    }
    if(i < 5) {
      CSC_Geom_ME4_area->SetBinContent(i+1, 0); // Fill Bin 1-5
    }
    // ME2/1 + ME3/1 + ME4/1 :: Bins 02--20
    // <Trd1 name="ME21_ActiveGasVol" dz="94.83*cm " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="27.00*cm " dx2="62.855*cm "/>
    // <Trd1 name="ME31_ActiveGasVol" dz="84.85*cm " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="30.70*cm " dx2="62.855*cm "/>
    // <Trd1 name="ME41_ActiveGasVol" dz="74.710*cm " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="34.505*cm " dx2="62.825*cm "/>  ??? typo should be 62.855 ???
    // chamberlength = 189.66 --- 169.70 --- 149.42 ==> 19 - 17 - 15 bins
    double dx_me21 = (125.71-54.00)/189.7; double x_me21 = 54.00;
    if(i > 0 && i< 20) {
      double top = x_me21 + (r_MEX_geom[i]  -r_MEX_geom[1])*dx_me21;
      double bottom = top + (r_MEX_geom[i+1]-r_MEX_geom[i])*dx_me21;
      double area = 2*18*0.5*(top + bottom)*(r_MEX_geom[i+1]-r_MEX_geom[i]); area_me21 += area;
      // std::cout<<"CSC ME2/1 Area :: i = "<<i<<" R_min = "<<r_MEX_geom[i]<<" cm and R_max = "<<r_MEX_geom[i+1]<<" cm ==> dR = "<<(r_MEX_geom[i+1]-r_MEX_geom[i])<<" cm";
      // std::cout<<" || x1 = "<<top<<"cm | x2 = "<<bottom<<" cm | Area of little trapezoid = "<<area<<" cm^2"<<std::endl;
      CSC_Geom_ME2_area->SetBinContent(i+1, area);
    }
    double dx_me31 = (125.71-61.40)/169.7; double x_me31 = 61.40;
    if(i > 2 && i < 20) {
      double top = x_me31 + (r_MEX_geom[i]  -r_MEX_geom[3])*dx_me31;
      double bottom = top + (r_MEX_geom[i+1]-r_MEX_geom[i])*dx_me31;
      double area = 2*18*0.5*(top + bottom)*(r_MEX_geom[i+1]-r_MEX_geom[i]); area_me31 += area;
      CSC_Geom_ME3_area->SetBinContent(i+1, area);
    }
    double dx_me41 = (125.71-69.01)/149.6; double x_me41 = 69.01;
    if(i > 4 && i < 20) {
      double top = x_me41 + (r_MEX_geom[i]  -r_MEX_geom[5])*dx_me41;
      double bottom = top + (r_MEX_geom[i+1]-r_MEX_geom[i])*dx_me41;
      double area = 2*18*0.5*(top + bottom)*(r_MEX_geom[i+1]-r_MEX_geom[i]); area_me41 += area;
      CSC_Geom_ME4_area->SetBinContent(i+1, area);
    }

    if(i==20) {
      CSC_Geom_ME2_area->SetBinContent(i+1, 0); // Fill Bin 21
      CSC_Geom_ME3_area->SetBinContent(i+1, 0); // Fill Bin 21
      CSC_Geom_ME4_area->SetBinContent(i+1, 0); // Fill Bin 21
    }
    // ME2/2 + ME3/2 + ME4/2
    // <Trd1 name="ME22_ActiveGasVol" dz="1.6153*m  " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="33.230*cm" dx2="63.575*cm"/>
    // <Trd1 name="ME32_ActiveGasVol" dz="1.6153*m  " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="33.230*cm " dx2="63.575*cm "/>
    // <Trd1 name="ME42_ActiveGasVol" dz="1.6153*m  " dy1="[ActiveGasVolHalfThick]" dy2="[ActiveGasVolHalfThick]" dx1="33.230*cm " dx2="63.575*cm "/>
    // chamberlength = 323.06 ==> 32 bins
    double dx_meX2 = (127.15-66.46)/323.1; double x_meX2 = 66.46;
    if(i>20 && i<53) {
      double top = x_meX2 + (r_MEX_geom[i]  -r_MEX_geom[21])*dx_meX2;
      double bottom = top + (r_MEX_geom[i+1]-r_MEX_geom[i])*dx_meX2;
      double area = 2*36*0.5*(top + bottom)*(r_MEX_geom[i+1]-r_MEX_geom[i]); area_me22 += area; area_me32 += area; area_me42 += area;
      // std::cout<<"CSC MEX/2 Area :: i = "<<i<<" R_min = "<<r_MEX_geom[i]<<" cm and R_max = "<<r_MEX_geom[i+1]<<" cm ==> dR = "<<(r_MEX_geom[i+1]-r_MEX_geom[i])<<" cm";
      // std::cout<<" || x1 = "<<top<<"cm | x2 = "<<bottom<<" cm | Area of little trapezoid = "<<area<<" cm^2"<<std::endl;
      CSC_Geom_ME2_area->SetBinContent(i+1, area);
      CSC_Geom_ME3_area->SetBinContent(i+1, area);
      CSC_Geom_ME4_area->SetBinContent(i+1, area);
    }
  }
  // std::cout<<"----------------------------------------------------------------------"<<std::endl;
  // std::cout<<"--- AREA of CSC system ---"<<std::endl;
  // std::cout<<"ME1/1 :: "<<area_me11/36/2<<" cm^2 / chamber ==> "<<area_me11<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME1/2 :: "<<area_me12/36/2<<" cm^2 / chamber ==> "<<area_me12<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME1/3 :: "<<area_me13/36/2<<" cm^2 / chamber ==> "<<area_me13<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME2/1 :: "<<area_me21/36/2<<" cm^2 / chamber ==> "<<area_me21<<" cm^2 in total [2 endcaps x 18 chambers]"<<std::endl;
  // std::cout<<"ME2/2 :: "<<area_me22/36/2<<" cm^2 / chamber ==> "<<area_me22<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME3/1 :: "<<area_me31/36/2<<" cm^2 / chamber ==> "<<area_me31<<" cm^2 in total [2 endcaps x 18 chambers]"<<std::endl;
  // std::cout<<"ME3/2 :: "<<area_me32/36/2<<" cm^2 / chamber ==> "<<area_me32<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;
  // std::cout<<"ME4/1 :: "<<area_me41/36/2<<" cm^2 / chamber ==> "<<area_me41<<" cm^2 in total [2 endcaps x 18 chambers]"<<std::endl;
  // std::cout<<"ME4/2 :: "<<area_me42/36/2<<" cm^2 / chamber ==> "<<area_me42<<" cm^2 in total [2 endcaps x 36 chambers]"<<std::endl;

  // RPC
  // <Trd1 name="REB1" dz="26.75*cm " dy1="1*mm " dy2="1*mm " dx1="35.721*cm " dx2="40.40*cm "/>
  // <Trd1 name="REB2" dz="24.5*cm " dy1="1*mm " dy2="1*mm " dx1="31.347*cm " dx2="35.633*cm "/>
  // <Trd1 name="REB3" dz="39.7*cm " dy1="1*mm " dy2="1*mm " dx1="24.315*cm " dx2="31.26*cm "/>
  // <Trd1 name="REC1" dz="20.15*cm " dy1="1*mm " dy2="1*mm " dx1="56.097*cm " dx2="59.623*cm "/>
  // <Trd1 name="REC2" dz="37.95*cm " dy1="1*mm " dy2="1*mm " dx1="49.369*cm " dx2="56.009*cm "/>
  // <Trd1 name="REC3" dz="28*cm " dy1="1*mm " dy2="1*mm " dx1="44.382*cm " dx2="49.281*cm "/>
  // <Trd1 name="REE1" dz="30.85*cm " dy1="1*mm " dy2="1*mm " dx1="38.118*cm " dx2="43.516*cm "/>
  // <Trd1 name="REE2" dz="26.95*cm " dy1="1*mm " dy2="1*mm " dx1="33.315*cm " dx2="38.03*cm "/>
  // <Trd1 name="REE3" dz="23.95*cm " dy1="1*mm " dy2="1*mm " dx1="29.036*cm " dx2="33.227*cm "/>
  // <Trd1 name="REF1" dz="25.95*cm " dy1="1*mm " dy2="1*mm " dx1="56.14*cm " dx2="60.681*cm "/>
  // <Trd1 name="REF2" dz="37.95*cm " dy1="1*mm " dy2="1*mm " dx1="49.413*cm " dx2="56.053*cm "/>
  // <Trd1 name="REF3" dz="30.2*cm " dy1="1*mm " dy2="1*mm " dx1="44.041*cm " dx2="49.325*cm "/>

  // RPC Upgrade
  // <Trd1 name="REG1" dz="16.05*cm " dy1="1*mm " dy2="1*mm " dx1="50.68*cm " dx2="56.34*cm " />
  // <Trd1 name="REG2" dz="16.05*cm " dy1="1*mm " dy2="1*mm " dx1="44.84*cm " dx2="50.5*cm " />
  // <Trd1 name="REG3" dz="16.05*cm " dy1="1*mm " dy2="1*mm " dx1="39.*cm " dx2="44.46*cm " />
  // <Trd1 name="REG4" dz="16.05*cm " dy1="1*mm " dy2="1*mm " dx1="33.17*cm " dx2="38.83*cm " />
  // <Trd1 name="REG5" dz="16.05*cm " dy1="1*mm " dy2="1*mm " dx1="27.33*cm " dx2="32.99*cm " />
  // <Trd1 name="REH1" dz="[dzS1]" dy1="[dyGap]" dy2="[dyGap]" dx1="[dx11]" dx2="[dx12]" />
  // <Trd1 name="REH2" dz="[dzS2]" dy1="[dyGap]" dy2="[dyGap]" dx1="[dx21]" dx2="[dx22]" />
  // <Trd1 name="REH3" dz="[dzS3]" dy1="[dyGap]" dy2="[dyGap]" dx1="[dx31]" dx2="[dx32]" />
  // <Trd1 name="REH4" dz="[dzS4]" dy1="[dyGap]" dy2="[dyGap]" dx1="[dx41]" dx2="[dx42]" />
  // <Trd1 name="REH5" dz="[dzS5]" dy1="[dyGap]" dy2="[dyGap]" dx1="[dx51]" dx2="[dx52]" />

  // GEM Upgrade



}


// ------------ method called when ending the processing of a run  ------------
void 
MyME11SimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{

}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyME11SimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyME11SimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyME11SimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyME11SimHitAnalyzer);

//  LocalWords:  str
