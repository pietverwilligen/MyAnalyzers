// -*- C++ -*-
//
// Package:    MyNeutronSimHitAnalyzer
// Class:      MyNeutronSimHitAnalyzer
// 
/**\class MyNeutronSimHitAnalyzer MyNeutronSimHitAnalyzer.cc MyAnalyzers/MyNeutronSimHitAnalyzer/src/MyNeutronSimHitAnalyzer.cc

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

// OLD :: class MyNeutronSimHitAnalyzer : public edm::EDAnalyzer {
class MyNeutronSimHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit MyNeutronSimHitAnalyzer(const edm::ParameterSet&);
      ~MyNeutronSimHitAnalyzer();
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
  bool phys_debug, tech_debug, pdf_output;
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

  // Linear Histograms - 1ms bin width - 1000 bins 0-1s
  TH1F * DT_MB1_all_hits_tof, * DT_MB1_W2_hits_tof, * DT_MB1_W1_hits_tof, * DT_MB1_W0_hits_tof, * DT_MB2_all_hits_tof, * DT_MB3_all_hits_tof, * DT_MB4_all_hits_tof, * DT_MBX_all_hits_tof;
  TH1F * CSC_ME11_hits_tof, * CSC_ME12_hits_tof, * CSC_ME13_hits_tof, * CSC_ME21_hits_tof, * CSC_ME22_hits_tof, * CSC_ME31_hits_tof, * CSC_ME32_hits_tof, * CSC_ME41_hits_tof, * CSC_ME42_hits_tof;
  TH1F * CSC_MEX1_hits_tof, * CSC_MEX2_hits_tof, * CSC_MEXX_hits_tof;



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


  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TCanvas * Canvas_RPCb_hits,           * Canvas_RPCf_hits,           * Canvas_CSC_hits,           * Canvas_DT_hits,           * Canvas_GEM_hits,           * Canvas_ME0_hits;
  TCanvas * Canvas_RPCb_deposits,       * Canvas_RPCf_deposits,       * Canvas_CSC_deposits,       * Canvas_DT_deposits,       * Canvas_GEM_deposits,       * Canvas_ME0_deposits;
  TCanvas * Canvas_RPCb_hits_fancy,     * Canvas_RPCf_hits_fancy,     * Canvas_CSC_hits_fancy,     * Canvas_DT_hits_fancy,     * Canvas_GEM_hits_fancy,     * Canvas_ME0_hits_fancy;
  TCanvas * Canvas_RPCb_deposits_fancy, * Canvas_RPCf_deposits_fancy, * Canvas_CSC_deposits_fancy, * Canvas_DT_deposits_fancy, * Canvas_GEM_deposits_fancy, * Canvas_ME0_deposits_fancy;
  TCanvas * Canvas_RPCb_1D_tof,         * Canvas_RPCf_1D_tof,         * Canvas_CSC_1D_tof,         * Canvas_DT_1D_tof,         * Canvas_GEM_1D_tof,         * Canvas_ME0_1D_tof;
  TCanvas * Canvas_RPCb_1D_deps,        * Canvas_RPCf_1D_deps,        * Canvas_CSC_1D_deps,        * Canvas_DT_1D_deps,        * Canvas_GEM_1D_deps,        * Canvas_ME0_1D_deps;
  TCanvas * Canvas_RPCb_hits_phi,       * Canvas_RPCf_hits_phi,       * Canvas_CSC_hits_phi,       * Canvas_DT_hits_phi;
  TCanvas * Canvas_CSC_EntryExit_vs_deposits, * Canvas_CSC_EntryExit_vs_kinenergy, * Canvas_CSC_EntryExit_vs_time; 
  TCanvas * Canvas_CSC_PathLength_vs_deposits, * Canvas_CSC_PathLength_vs_kinenergy;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 



  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_XY_RZ_views;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH2F    * RPCb_XY, * RPCb_RZ, * RPCf_XY, * RPCf_RZ, *CSC_XY, * CSC_RZ, * DT_XY, * DT_RZ, * Muon_Barrel_XY, * Muon_Endcap_XY, * Muon_RZ, * GEM_XY, *GEM_RZ, * ME0_XY, *ME0_RZ;
  TH2F    * RPCb_000ns_RZ, * RPCf_000ns_RZ, * CSC_000ns_RZ, * DT_000ns_RZ, * Muon_000ns_RZ, * GEM_000ns_RZ, * ME0_000ns_RZ;
  TH2F    * RPCb_250ns_RZ, * RPCf_250ns_RZ, * CSC_250ns_RZ, * DT_250ns_RZ, * Muon_250ns_RZ, * GEM_250ns_RZ, * ME0_250ns_RZ;
  TH2F    * RPCb_000ns_XY, * RPCf_000ns_XY, * CSC_000ns_XY, * DT_000ns_XY, * GEM_000ns_XY, * ME0_000ns_XY, * Muon_Barrel_000ns_XY, * Muon_Endcap_000ns_XY;
  TH2F    * RPCb_250ns_XY, * RPCf_250ns_XY, * CSC_250ns_XY, * DT_250ns_XY, * GEM_250ns_XY, * ME0_250ns_XY, * Muon_Barrel_250ns_XY, * Muon_Endcap_250ns_XY;

  TH2F * Muon_000ns_el_RZ, * Muon_000ns_mu_RZ, * Muon_000ns_ha_RZ; // hits with 000 < tof < 250 ns
  TH2F * Muon_250ns_el_RZ, * Muon_250ns_mu_RZ, * Muon_250ns_ha_RZ; // hits with 250 < tof < 10^8 ns
  TH2F * Muon_00ns_el_RZ,  * Muon_00ns_mu_RZ,  * Muon_00ns_ha_RZ;  // hits with 00 < tof < 050 ns
  TH2F * Muon_50ns_el_RZ,  * Muon_50ns_mu_RZ,  * Muon_50ns_ha_RZ;  // hirs with 50 < tof < 250 ns

  TH2F * Muon_Barrel_000ns_el_XY, * Muon_Barrel_000ns_mu_XY, * Muon_Barrel_000ns_ha_XY; // hits with 000 < tof < 250 ns
  TH2F * Muon_Barrel_250ns_el_XY, * Muon_Barrel_250ns_mu_XY, * Muon_Barrel_250ns_ha_XY; // hits with 250 < tof < 10^8 ns
  TH2F * Muon_Barrel_00ns_el_XY,  * Muon_Barrel_00ns_mu_XY,  * Muon_Barrel_00ns_ha_XY;  // hits with 00 < tof < 050 ns
  TH2F * Muon_Barrel_50ns_el_XY,  * Muon_Barrel_50ns_mu_XY,  * Muon_Barrel_50ns_ha_XY;  // hirs with 50 < tof < 250 ns

  TH2F * Muon_Endcap_000ns_el_XY, * Muon_Endcap_000ns_mu_XY, * Muon_Endcap_000ns_ha_XY; // hits with 000 < tof < 250 ns
  TH2F * Muon_Endcap_250ns_el_XY, * Muon_Endcap_250ns_mu_XY, * Muon_Endcap_250ns_ha_XY; // hits with 250 < tof < 10^8 ns
  TH2F * Muon_Endcap_00ns_el_XY,  * Muon_Endcap_00ns_mu_XY,  * Muon_Endcap_00ns_ha_XY;  // hits with 00 < tof < 050 ns
  TH2F * Muon_Endcap_50ns_el_XY,  * Muon_Endcap_50ns_mu_XY,  * Muon_Endcap_50ns_ha_XY;  // hirs with 50 < tof < 250 ns

  TH2F * SimVertices_RZ, * SimVertices_Muon_RZ;
  TH1F * PrimVertices_Z, * PrimVertices_R;

  TH1F * RPC_hits, * RPC_area, * RPC_rates;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TCanvas * Canvas_Muon_RZ,        * Canvas_Muon_000ns_RZ,        * Canvas_Muon_250ns_RZ,        * Canvas_Muon_000ns_Cont_RZ,        * Canvas_Muon_250ns_Cont_RZ; 
  TCanvas                          * Canvas_Muon_00ns_RZ,         * Canvas_Muon_50ns_RZ,         * Canvas_Muon_00ns_Cont_RZ,         * Canvas_Muon_50ns_Cont_RZ; 
  TCanvas * Canvas_Muon_Barrel_XY, * Canvas_Muon_Barrel_000ns_XY, * Canvas_Muon_Barrel_250ns_XY, * Canvas_Muon_Barrel_000ns_Cont_XY, * Canvas_Muon_Barrel_250ns_Cont_XY; 
  TCanvas                          * Canvas_Muon_Barrel_00ns_XY,  * Canvas_Muon_Barrel_50ns_XY,  * Canvas_Muon_Barrel_00ns_Cont_XY,  * Canvas_Muon_Barrel_50ns_Cont_XY; 
  TCanvas * Canvas_Muon_Endcap_XY, * Canvas_Muon_Endcap_000ns_XY, * Canvas_Muon_Endcap_250ns_XY, * Canvas_Muon_Endcap_000ns_Cont_XY, * Canvas_Muon_Endcap_250ns_Cont_XY; 
  TCanvas                          * Canvas_Muon_Endcap_00ns_XY,  * Canvas_Muon_Endcap_50ns_XY,  * Canvas_Muon_Endcap_00ns_Cont_XY,  * Canvas_Muon_Endcap_50ns_Cont_XY; 
  TCanvas * Canvas_SimVertices_RZ, * Canvas_SimVertices_Muon_RZ;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 



  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hit_info;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * RPCb_Muons_SHPT,           * RPCf_Muons_SHPT,           * CSC_Muons_SHPT,           *DT_Muons_SHPT,           *GEM_Muons_SHPT,           *ME0_Muons_SHPT;           
  TH1F * RPCb_Hadrons_SHPT,         * RPCf_Hadrons_SHPT,         * CSC_Hadrons_SHPT,         *DT_Hadrons_SHPT,         *GEM_Hadrons_SHPT,         *ME0_Hadrons_SHPT;
  TH1F * RPCb_Electrons_SHPT,       * RPCf_Electrons_SHPT,       * CSC_Electrons_SHPT,       *DT_Electrons_SHPT,       *GEM_Electrons_SHPT,       *ME0_Electrons_SHPT;
  TH1F * RPCb_Electrons_000ns_SHPT, * RPCf_Electrons_000ns_SHPT, * CSC_Electrons_000ns_SHPT, *DT_Electrons_000ns_SHPT, *GEM_Electrons_000ns_SHPT, *ME0_Electrons_000ns_SHPT;
  TH1F * RPCb_Electrons_050ns_SHPT, * RPCf_Electrons_050ns_SHPT, * CSC_Electrons_050ns_SHPT, *DT_Electrons_050ns_SHPT, *GEM_Electrons_050ns_SHPT, *ME0_Electrons_050ns_SHPT;
  TH1F * RPCb_Electrons_250ns_SHPT, * RPCf_Electrons_250ns_SHPT, * CSC_Electrons_250ns_SHPT, *DT_Electrons_250ns_SHPT, *GEM_Electrons_250ns_SHPT, *ME0_Electrons_250ns_SHPT;

  // TH1F * RPCb_Electrons_SVPT,       * RPCf_Electrons_SVPT,       * CSC_Electrons_SVPT,       *DT_Electrons_SVPT;
  // hits per chamber (HPC)
  TH1F * RPCb_HPC, * RPCb_el_HPC;
  TH1F * RPCf_HPC, * RPCf_el_HPC;
  TH1F * CSC_HPC,  * CSC_el_HPC;
  TH1F * DT_HPC,   * DT_el_HPC;
  // hits per layer (HPL)
  TH1F * RPCb_el_HPL;  TH1F * RPCb_mu_HPL; TH1F * RPCb_ha_HPL;
  TH1F * RPCf_el_HPL;  TH1F * RPCf_mu_HPL; TH1F * RPCf_ha_HPL;
  TH1F * CSC_el_HPL;   TH1F * CSC_mu_HPL;  TH1F * CSC_ha_HPL;
  TH1F * DT_el_HPL;    TH1F * DT_mu_HPL;   TH1F * DT_ha_HPL;
  TH1F * MB1_el_HPL, * MB1_mu_HPL, * MB1_ha_HPL, * MB2_el_HPL, * MB2_mu_HPL, * MB2_ha_HPL, * MB3_el_HPL, * MB3_mu_HPL, * MB3_ha_HPL, * MB4_el_HPL, * MB4_mu_HPL, * MB4_ha_HPL;
  TH1F * ME1_el_HPL, * ME1_mu_HPL, * ME1_ha_HPL, * ME2_el_HPL, * ME2_mu_HPL, * ME2_ha_HPL, * ME3_el_HPL, * ME3_mu_HPL, * ME3_ha_HPL, * ME4_el_HPL, * ME4_mu_HPL, * ME4_ha_HPL;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TCanvas * Canvas_RPCb_Layers, * Canvas_RPCf_Layers, * Canvas_CSC_Layers, * Canvas_DT_Layers; 
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 



  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_entry_info;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F *RPCb_EntryExit_All_Glob_dr,      *DT_EntryExit_All_Glob_dr;                                                                    // Barrel Detectors
  TH1F *RPCf_EntryExit_All_Glob_dz,      *CSC_EntryExit_All_Glob_dz,   *GEM_EntryExit_All_Glob_dz,      *ME0_EntryExit_All_Glob_dz;    // Endcap Detectors
  TH1F *CSC_EntryExit_Electrons_Glob_dz, *CSC_EntryExit_Muons_Glob_dz, *CSC_EntryExit_Hadrons_Glob_dz;
  TH1F *RPCb_EntryExit_All_Loc_dz,       *RPCf_EntryExit_All_Loc_dz,   *CSC_EntryExit_All_Loc_dz,       *DT_EntryExit_All_Loc_dz, *GEM_EntryExit_All_Loc_dz, *ME0_EntryExit_All_Loc_dz;
  TH1F                                   *CSC_EntryExit_Electrons_Loc_dz,    *CSC_EntryExit_Muons_Loc_dz,    *CSC_EntryExit_Hadrons_Loc_dz;
  TH1F *CSC_EntryExit_All_Glob_dGap,     *CSC_EntryExit_Electrons_Glob_dGap, *CSC_EntryExit_Muons_Glob_dGap, *CSC_EntryExit_Hadrons_Glob_dGap;

  // TH2F *CSC_EntryExit_All_Deposit_dz, *CSC_EntryExit_Electrons_Deposit_dz, *CSC_EntryExit_Muons_Deposit_dz, *CSC_EntryExit_Hadrons_Deposit_dz;
  TH2F *CSC_EntryExit_el_Deposit_dz, *CSC_EntryExit_mu_Deposit_dz, *CSC_EntryExit_pi_Deposit_dz, *CSC_EntryExit_ka_Deposit_dz;
  TH2F *CSC_EntryExit_p_Deposit_dz,  *CSC_EntryExit_n_Deposit_dz,  *CSC_EntryExit_g_Deposit_dz,  *CSC_EntryExit_N_Deposit_dz;
  TH2F *CSC_EntryExit_el_KinEn_dz,   *CSC_EntryExit_mu_KinEn_dz,   *CSC_EntryExit_pi_KinEn_dz,   *CSC_EntryExit_ka_KinEn_dz;
  TH2F *CSC_EntryExit_p_KinEn_dz,    *CSC_EntryExit_n_KinEn_dz,    *CSC_EntryExit_g_KinEn_dz,    *CSC_EntryExit_N_KinEn_dz;
  TH2F *CSC_EntryExit_el_Time_dz,    *CSC_EntryExit_mu_Time_dz,    *CSC_EntryExit_pi_Time_dz,    *CSC_EntryExit_ka_Time_dz;
  TH2F *CSC_EntryExit_p_Time_dz,     *CSC_EntryExit_n_Time_dz,     *CSC_EntryExit_g_Time_dz,     *CSC_EntryExit_N_Time_dz;

  TH2F *CSC_EntryExit_el_Deposit_dR, *CSC_EntryExit_mu_Deposit_dR, *CSC_EntryExit_pi_Deposit_dR, *CSC_EntryExit_ka_Deposit_dR;
  TH2F *CSC_EntryExit_p_Deposit_dR,  *CSC_EntryExit_n_Deposit_dR,  *CSC_EntryExit_g_Deposit_dR,  *CSC_EntryExit_N_Deposit_dR;
  TH2F *CSC_EntryExit_el_KinEn_dR,   *CSC_EntryExit_mu_KinEn_dR,   *CSC_EntryExit_pi_KinEn_dR,   *CSC_EntryExit_ka_KinEn_dR;
  TH2F *CSC_EntryExit_p_KinEn_dR,    *CSC_EntryExit_n_KinEn_dR,    *CSC_EntryExit_g_KinEn_dR,    *CSC_EntryExit_N_KinEn_dR;

  TH2F *CSC_EntryExit_el_dz_dR, *CSC_EntryExit_mu_dz_dR, *CSC_EntryExit_pi_dz_dR, *CSC_EntryExit_ka_dz_dR;
  TH2F *CSC_EntryExit_p_dz_dR,  *CSC_EntryExit_n_dz_dR,  *CSC_EntryExit_g_dz_dR,  *CSC_EntryExit_N_dz_dR;

  TH2F *CSC_EntryExit_el_dz_dR_detail, *CSC_EntryExit_mu_dz_dR_detail, *CSC_EntryExit_pi_dz_dR_detail, *CSC_EntryExit_ka_dz_dR_detail;
  TH2F *CSC_EntryExit_p_dz_dR_detail,  *CSC_EntryExit_n_dz_dR_detail,  *CSC_EntryExit_g_dz_dR_detail,  *CSC_EntryExit_N_dz_dR_detail;

  TH2F *CSC_EntryExit_el_GapLength_Deposit, *CSC_EntryExit_mu_GapLength_Deposit, *CSC_EntryExit_pi_GapLength_Deposit, *CSC_EntryExit_ka_GapLength_Deposit;
  TH2F *CSC_EntryExit_p_GapLength_Deposit,  *CSC_EntryExit_n_GapLength_Deposit,  *CSC_EntryExit_g_GapLength_Deposit,  *CSC_EntryExit_N_GapLength_Deposit;

  TH2F *CSC_EntryExit_el_Deposit_pidR2, *CSC_EntryExit_mu_Deposit_pidR2, *CSC_EntryExit_pi_Deposit_pidR2, *CSC_EntryExit_ka_Deposit_pidR2;
  TH2F *CSC_EntryExit_p_Deposit_pidR2,  *CSC_EntryExit_n_Deposit_pidR2,  *CSC_EntryExit_g_Deposit_pidR2,  *CSC_EntryExit_N_Deposit_pidR2;
  TH2F *CSC_EntryExit_el_Deposit_dRdz,  *CSC_EntryExit_mu_Deposit_dRdz,  *CSC_EntryExit_pi_Deposit_dRdz,  *CSC_EntryExit_ka_Deposit_dRdz;
  TH2F *CSC_EntryExit_p_Deposit_dRdz,   *CSC_EntryExit_n_Deposit_dRdz,   *CSC_EntryExit_g_Deposit_dRdz,   *CSC_EntryExit_N_Deposit_dRdz;

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TCanvas * Canvas_CSC_dz_vs_dR, * Canvas_CSC_dz_vs_dR_detail, * Canvas_CSC_deposits_vs_pidR2, * Canvas_CSC_deposits_vs_dRdz, *Canvas_CSC_deposits_vs_GapLength;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 



  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_rates;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TGraphErrors * gr_RPC_Rates_All, * gr_RPC_Rates_Barrel, * gr_RPC_Rates_Endcap;
  TGraphErrors * gr_RPC_Rates_RE12, * gr_RPC_Rates_RE13, * gr_RPC_Rates_RE22, * gr_RPC_Rates_RE23, * gr_RPC_Rates_RE32, * gr_RPC_Rates_RE33, * gr_RPC_Rates_RE42, * gr_RPC_Rates_RE43;
  // TGraphErrors * gr_CSC_Rates_ME1;
  TH1F * CSC_ME1_all_hits, * CSC_ME1_area, * CSC_ME1_all_rates;
  TH1F * CSC_ME2_all_hits, * CSC_ME2_area, * CSC_ME2_all_rates;
  TH1F * CSC_ME3_all_hits, * CSC_ME3_area, * CSC_ME3_all_rates;
  TH1F * CSC_ME4_all_hits, * CSC_ME4_area, * CSC_ME4_all_rates;

  TH1F * CSC_ME1_000ns_hits, * CSC_ME1_000ns_rates, * CSC_ME1_250ns_hits, * CSC_ME1_250ns_rates, * CSC_ME1_00ns_hits, * CSC_ME1_00ns_rates, * CSC_ME1_50ns_hits, * CSC_ME1_50ns_rates;

  TH1F * CSC_Geom_ME1_all_hits, * CSC_Geom_ME1_area, * CSC_Geom_ME1_all_rates;
  TH1F * CSC_Geom_ME2_all_hits, * CSC_Geom_ME2_area, * CSC_Geom_ME2_all_rates;
  TH1F * CSC_Geom_ME3_all_hits, * CSC_Geom_ME3_area, * CSC_Geom_ME3_all_rates;
  TH1F * CSC_Geom_ME4_all_hits, * CSC_Geom_ME4_area, * CSC_Geom_ME4_all_rates;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TCanvas * Canvas_RPC_Rates, * Canvas_RPC_Rates_RE1, * Canvas_RPC_Rates_RE2, * Canvas_RPC_Rates_RE3, * Canvas_RPC_Rates_RE4;
  int    RPC_hits_array[2][4][3];
  double RPC_area_array[2][4][3], RPC_rates_array[2][4][3], InstLumi[50], RPC_rates_Summary[15][50], RPC_uncer_Rate[15][50], RPC_uncer_Lumi[15][50]; 
  double rpc_barrel_area, rpc_endcap_area = 0.0;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 



};

//
// constants, enums and typedefs
//

int n_tof = 1100,  n1_tof = 1,  n2_tof = 12;
int m_tof = 110,   m1_tof = 1,  m2_tof = 12;
int n_toflin = 1000, n1_toflin = 0, n2_toflin = 50; // 0-1000ms in 1ms | 0-100ms in 100us | 0-50ms in steps of 50us
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

// CSC Areas
// ---------
// |      |  ME11 |  ME12 | ME13  |  ME21 |  ME31 |  ME41 | ME2x  |
// |------|-------|-------|-------|-------|-------|-------|-------|
// | R_1  | 105.5 | 281.5 | 512.0 | 146.9 | 166.9 | 187.0 | 364.0 |
// | R_2  | 257.5 | 456.0 | 676.2 | 336.6 | 336.6 | 336.4 | 687.1 |
// |------|-------|-------|-------|-------|-------|-------|-------|
// | Diff | 152.0 | 174.5 | 164.2 | 189.7 | 169.7 | 149.6 | 323.1 |
// |------|-------|-------|-------|-------|-------|-------|-------|
// | Bins |    15 |    17 |    16 |    19 |    17 |    15 |    32 |
// |------|-------|-------|-------|-------|-------|-------|-------|

double r_ME1_geom[] = {000.0, 
		       105.5, 115.5, 125.5, 135.5, 145.5, 155.5, 165.5, 175.5, 185.5, 195.5, 205.5, 215.5, 225.5, 235.5, 245.5, 257.5,                                                          // ME11 15 bins
		       281.5, 291.5, 301.5, 311.5, 321.5, 331.5, 341.5, 351.5, 361.5, 371.5, 381.5, 391.5, 401.5, 411.5, 421.5, 431.5, 441.5, 456.0,                                            // ME21 17 bins
                       512.0, 522.0, 532.0, 542.0, 552.0, 562.0, 572.0, 582.0, 592.0, 602.0, 612.0, 622.0, 632.0, 642.0, 652.0, 662.0, 676.2};                                                  // ME31 16 bins ==> tot = 51 bins
/*
double r_ME2_geom[] = {000.0, 
                       146.9, 156.9, 166.9, 176.9, 186.9, 196.9, 206.9, 216.9, 226.9, 236.9, 246.9, 256.9, 266.9, 276.9, 286.9, 296.9, 306.9, 316.9, 326.9, 336.6,                              // ME21 19 bins
		       364.0, 374.0, 384.0, 394.0, 404.0, 414.0, 424.0, 434.0, 444.4, 454.0, 464.0, 474.0, 484.0, 494.0, 504.0, 514.0, 524.0, 534.0, 544.4, 554.0, 564.0, 574.0, 584.0, 594.0,  
                       604.0, 614.0, 624.0, 634.0, 644.4, 654.0, 664.0, 674.0, 687.1};                                                                                                          // ME22 32 bins ==> tot = 53 bins
double r_ME3_geom[] = {000.0, 
		       166.9, 176.9, 186.9, 196.9, 206.9, 216.9, 226.9, 236.9, 246.9, 256.9, 266.9, 276.9, 286.9, 296.9, 306.9, 316.9, 326.9, 336.6,                                            // ME31 17 bins 
                       364.0, 374.0, 384.0, 394.0, 404.0, 414.0, 424.0,434.0, 444.4, 454.0, 464.0, 474.0, 484.0, 494.0, 504.0, 514.0, 524.0, 534.0, 544.4, 554.0, 564.0, 574.0, 584.0, 594.0,  
                       604.0, 614.0, 624.0, 634.0, 644.4, 654.0, 664.0, 674.0, 687.1};                                                                                                          // ME32 32 bins ==> tot = 51 bins
double r_ME4_geom[] = {000.0,
                       187.0, 197.0, 207.0, 217.0, 227.0, 237.0, 247.0, 257.0, 267.0, 277.0, 287.0, 297.0, 307.0, 317.0, 327.0, 336.4,                                                          // ME41 15 bins
                       364.0, 374.0, 384.0, 394.0, 404.0, 414.0, 424.0,434.0, 444.4, 454.0, 464.0, 474.0, 484.0, 494.0, 504.0, 514.0, 524.0, 534.0, 544.4, 554.0, 564.0, 574.0, 584.0, 594.0,
		       604.0, 614.0, 624.0, 634.0, 644.4, 654.0, 664.0, 674.0, 687.1};                                                                                                          // ME42 32 bins ==> tot = 49 bins
*/
double r_MEX_geom[] = {000.0, 
                       146.9, 156.9, 166.9, 176.9, 186.9, 196.9, 206.9, 216.9, 226.9, 236.9, 246.9, 256.9, 266.9, 276.9, 286.9, 296.9, 306.9, 316.9, 326.9, 336.6,                              // ME21 19 bins
		       364.0, 374.0, 384.0, 394.0, 404.0, 414.0, 424.0, 434.0, 444.0, 454.0, 464.0, 474.0, 484.0, 494.0, 504.0, 514.0, 524.0, 534.0, 544.0, 554.0, 564.0, 574.0, 584.0, 594.0,  
                       604.0, 614.0, 624.0, 634.0, 644.0, 654.0, 664.0, 674.0, 687.1};                                                                                                          // ME22 32 bins ==> tot = 53 bins
		       

//
// static data member definitions
//

//
// constructors and destructor
//
MyNeutronSimHitAnalyzer::MyNeutronSimHitAnalyzer(const edm::ParameterSet& iConfig)

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
  pdf_output      = iConfig.getUntrackedParameter<bool>("PDFOutput");
  outputfile      = new TFile(rootFileName.c_str(), "RECREATE" );

  // required for 7XY: registration of the data access
  consumesMany<edm::PSimHitContainer>();
  SIMTrack_Token  = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
  SIMVertex_Token = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits")); // or consumes<std::vector<SimVertex>>

  cscGeom_Token = esConsumes<CSCGeometry, MuonGeometryRecord>();
  dtGeom_Token  = esConsumes<DTGeometry,  MuonGeometryRecord>();
  gemGeom_Token = esConsumes<GEMGeometry, MuonGeometryRecord>();
  rpcGeom_Token = esConsumes<RPCGeometry, MuonGeometryRecord>();

  if(tech_debug) std::cout<<"[MyNeutronSimHitAnalyzer :: Constructor]"<<std::endl;
  
  TDir_Muon_hits_deposits = new TDirectoryFile("Muon_hits_deposits", "Muon_hits_deposits");
  TDir_Muon_XY_RZ_views   = new TDirectoryFile("Muon_XY_RZ_views",   "Muon_XY_RZ_views");
  TDir_Muon_hit_info      = new TDirectoryFile("Muon_hit_info",      "Muon_hit_info");
  TDir_Muon_entry_info    = new TDirectoryFile("Muon_entry_info",    "Muon_entry_info");
  TDir_Muon_rates         = new TDirectoryFile("Muon_rates",         "Muon_rates");
  // TDir_Muon_Barrel_asymm  = new TDirectoryFile("Muon_Barrel_asymm",  "Muon_Barrel_assym");
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

  DT_MB1_all_hits_tof = new TH1F("DT_MB1_all_hits_tof", "Time of Flight :: DT MB1 :: all hits", n_toflin, n1_toflin, n2_toflin);
  DT_MB2_all_hits_tof = new TH1F("DT_MB2_all_hits_tof", "Time of Flight :: DT MB1 :: all hits", n_toflin, n1_toflin, n2_toflin);
  DT_MB3_all_hits_tof = new TH1F("DT_MB3_all_hits_tof", "Time of Flight :: DT MB1 :: all hits", n_toflin, n1_toflin, n2_toflin);
  DT_MB4_all_hits_tof = new TH1F("DT_MB4_all_hits_tof", "Time of Flight :: DT MB1 :: all hits", n_toflin, n1_toflin, n2_toflin);
  DT_MB1_W2_hits_tof  = new TH1F("DT_MB1_W2_hits_tof",  "Time of Flight :: DT MB1 Wheel+/-2 :: all hits", n_toflin, n1_toflin, n2_toflin);
  DT_MB1_W1_hits_tof  = new TH1F("DT_MB1_W1_hits_tof",  "Time of Flight :: DT MB1 Wheel+/-1 :: all hits", n_toflin, n1_toflin, n2_toflin);
  DT_MB1_W0_hits_tof  = new TH1F("DT_MB1_W0_hits_tof",  "Time of Flight :: DT MB1 Wheel 0 :: all hits",   n_toflin, n1_toflin, n2_toflin);
  DT_MBX_all_hits_tof = new TH1F("DT_MBX_all_hits_tof", "Time of Flight :: DT MBX :: all hits", n_toflin, n1_toflin, n2_toflin);
  CSC_ME11_hits_tof= new TH1F("CSC_ME11_hits_tof", "Time of Flight :: CSC ME11 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME12_hits_tof= new TH1F("CSC_ME12_hits_tof", "Time of Flight :: CSC ME12 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME13_hits_tof= new TH1F("CSC_ME13_hits_tof", "Time of Flight :: CSC ME13 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME21_hits_tof= new TH1F("CSC_ME21_hits_tof", "Time of Flight :: CSC ME21 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME22_hits_tof= new TH1F("CSC_ME22_hits_tof", "Time of Flight :: CSC ME22 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME31_hits_tof= new TH1F("CSC_ME31_hits_tof", "Time of Flight :: CSC ME31 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME32_hits_tof= new TH1F("CSC_ME32_hits_tof", "Time of Flight :: CSC ME32 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME41_hits_tof= new TH1F("CSC_ME41_hits_tof", "Time of Flight :: CSC ME41 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_ME42_hits_tof= new TH1F("CSC_ME42_hits_tof", "Time of Flight :: CSC ME42 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_MEX1_hits_tof= new TH1F("CSC_MEX1_hits_tof", "Time of Flight :: CSC MEX1 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_MEX2_hits_tof= new TH1F("CSC_MEX2_hits_tof", "Time of Flight :: CSC MEX2 :: all hits" , n_toflin, n1_toflin, n2_toflin);
  CSC_MEXX_hits_tof= new TH1F("CSC_MEXX_hits_tof", "Time of Flight :: CSC MEXX :: all hits" , n_toflin, n1_toflin, n2_toflin);




  RPCb_XY = new TH2F("RPCb_XY", "Simhits in XY :: RPCb", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCb_RZ = new TH2F("RPCb_RZ", "Simhits in RZ :: RPCb", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  RPCf_XY = new TH2F("RPCf_XY", "Simhits in XY :: RPCf", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCf_RZ = new TH2F("RPCf_RZ", "Simhits in RZ :: RPCf", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  CSC_XY  = new TH2F("CSC_XY",  "Simhits in XY :: CSC", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  CSC_RZ  = new TH2F("CSC_RZ",  "Simhits in RZ :: CSC", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  DT_XY   = new TH2F("DT_XY",   "Simhits in XY :: DT", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  DT_RZ   = new TH2F("DT_RZ",   "Simhits in RZ :: DT", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  GEM_XY  = new TH2F("GEM_XY",  "Simhits in XY :: GEM", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  GEM_RZ  = new TH2F("GEM_RZ",  "Simhits in RZ :: GEM", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  ME0_XY  = new TH2F("ME0_XY",  "Simhits in XY :: ME0", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  ME0_RZ  = new TH2F("ME0_RZ",  "Simhits in RZ :: ME0", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_RZ = new TH2F("Muon_RZ", "Simhits in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_Barrel_XY = new TH2F("Muon_Barrel_XY", "Simhits in XY :: Muon Barrel", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_XY = new TH2F("Muon_Endcap_XY", "Simhits in XY :: Muon Endcap", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);

  RPCb_000ns_XY = new TH2F("RPCb_000ns_XY", "Simhits with tof < 250ns in XY :: RPCb", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCb_000ns_RZ = new TH2F("RPCb_000ns_RZ", "Simhits with tof < 250ns in RZ :: RPCb", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  RPCf_000ns_XY = new TH2F("RPCf_000ns_XY", "Simhits with tof < 250ns in XY :: RPCf", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCf_000ns_RZ = new TH2F("RPCf_000ns_RZ", "Simhits with tof < 250ns in RZ :: RPCf", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  CSC_000ns_XY  = new TH2F("CSC_000ns_XY",  "Simhits with tof < 250ns in XY :: CSC", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  CSC_000ns_RZ  = new TH2F("CSC_000ns_RZ",  "Simhits with tof < 250ns in RZ :: CSC", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  DT_000ns_XY   = new TH2F("DT_000ns_XY",   "Simhits with tof < 250ns in XY :: DT", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  DT_000ns_RZ   = new TH2F("DT_000ns_RZ",   "Simhits with tof < 250ns in RZ :: DT", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  GEM_000ns_XY  = new TH2F("GEM_000ns_XY",  "Simhits with tof < 250ns in XY :: GEM", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  GEM_000ns_RZ  = new TH2F("GEM_000ns_RZ",  "Simhits with tof < 250ns in RZ :: GEM", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  ME0_000ns_XY  = new TH2F("ME0_000ns_XY",  "Simhits with tof < 250ns in XY :: ME0", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  ME0_000ns_RZ  = new TH2F("ME0_000ns_RZ",  "Simhits with tof < 250ns in RZ :: ME0", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_000ns_RZ = new TH2F("Muon_000ns_RZ", "Simhits with tof < 250ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_000ns_el_RZ = new TH2F("Muon_000ns_el_RZ", "Electron Simhits with tof < 250ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_000ns_mu_RZ = new TH2F("Muon_000ns_mu_RZ", "Muon Simhits with tof < 250ns in RZ :: Muon",     n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_000ns_ha_RZ = new TH2F("Muon_000ns_ha_RZ", "Hadron Simhits with tof < 250ns in RZ :: Muon",   n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_Barrel_000ns_XY    = new TH2F("Muon_Barrel_000ns_XY", "Simhits with tof < 250ns in XY :: Muon Barrel", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_000ns_el_XY = new TH2F("Muon_Barrel_000ns_el_XY", "Electron Simhits with tof < 250ns in XY :: Muon Barrel", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_000ns_mu_XY = new TH2F("Muon_Barrel_000ns_mu_XY", "Muon Simhits with tof < 250ns in XY :: Muon Barrel",     n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_000ns_ha_XY = new TH2F("Muon_Barrel_000ns_ha_XY", "Hadron Simhits with tof < 250ns in XY :: Muon Barrel",   n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_000ns_XY    = new TH2F("Muon_Endcap_000ns_XY", "Simhits with tof < 250ns in XY :: Muon Endcap", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_000ns_el_XY = new TH2F("Muon_Endcap_000ns_el_XY", "Electron Simhits with tof < 250ns in XY :: Muon Endcap", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_000ns_mu_XY = new TH2F("Muon_Endcap_000ns_mu_XY", "Muon Simhits with tof < 250ns in XY :: Muon Endcap",     n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_000ns_ha_XY = new TH2F("Muon_Endcap_000ns_ha_XY", "Hadron Simhits with tof < 250ns in XY :: Muon Endcap",   n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);

  RPCb_250ns_XY = new TH2F("RPCb_250ns_XY", "Simhits with tof > 250ns in XY :: RPCb", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCb_250ns_RZ = new TH2F("RPCb_250ns_RZ", "Simhits with tof > 250ns in RZ :: RPCb", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  RPCf_250ns_XY = new TH2F("RPCf_250ns_XY", "Simhits with tof > 250ns in XY :: RPCf", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCf_250ns_RZ = new TH2F("RPCf_250ns_RZ", "Simhits with tof > 250ns in RZ :: RPCf", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  CSC_250ns_XY  = new TH2F("CSC_250ns_XY",  "Simhits with tof > 250ns in XY :: CSC", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  CSC_250ns_RZ  = new TH2F("CSC_250ns_RZ",  "Simhits with tof > 250ns in RZ :: CSC", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  DT_250ns_XY   = new TH2F("DT_250ns_XY",   "Simhits with tof > 250ns in XY :: DT", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  DT_250ns_RZ   = new TH2F("DT_250ns_RZ",   "Simhits with tof > 250ns in RZ :: DT", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  GEM_250ns_XY  = new TH2F("GEM_250ns_XY",  "Simhits with tof > 250ns in XY :: GEM", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  GEM_250ns_RZ  = new TH2F("GEM_250ns_RZ",  "Simhits with tof > 250ns in RZ :: GEM", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  ME0_250ns_XY  = new TH2F("ME0_250ns_XY",  "Simhits with tof > 250ns in XY :: ME0", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  ME0_250ns_RZ  = new TH2F("ME0_250ns_RZ",  "Simhits with tof > 250ns in RZ :: ME0", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_250ns_RZ = new TH2F("Muon_250ns_RZ", "Simhits with tof > 250ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_250ns_el_RZ = new TH2F("Muon_250ns_el_RZ", "Electron Simhits with tof > 250ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_250ns_mu_RZ = new TH2F("Muon_250ns_mu_RZ", "Muon Simhits with tof > 250ns in RZ :: Muon",     n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_250ns_ha_RZ = new TH2F("Muon_250ns_ha_RZ", "Hadron Simhits with tof > 250ns in RZ :: Muon",   n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_Barrel_250ns_XY = new TH2F("Muon_Barrel_250ns_XY", "Simhits with tof > 250ns in XY :: Muon Barrel",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_250ns_el_XY = new TH2F("Muon_Barrel_250ns_el_XY", "Electron Simhits with tof > 250ns in XY :: Muon Barrel",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_250ns_mu_XY = new TH2F("Muon_Barrel_250ns_mu_XY", "Muon_Barrel Simhits with tof > 250ns in XY :: Muon Barrel",    n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_250ns_ha_XY = new TH2F("Muon_Barrel_250ns_ha_XY", "Hadron Simhits with tof > 250ns in XY :: Muon Barrel",  n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_250ns_XY = new TH2F("Muon_Endcap_250ns_XY", "Simhits with tof > 250ns in XY :: Muon Endcap",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_250ns_el_XY = new TH2F("Muon_Endcap_250ns_el_XY", "Electron Simhits with tof > 250ns in XY :: Muon Endcap",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_250ns_mu_XY = new TH2F("Muon_Endcap_250ns_mu_XY", "Muon_Endcap Simhits with tof > 250ns in XY :: Muon Endcap",    n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_250ns_ha_XY = new TH2F("Muon_Endcap_250ns_ha_XY", "Hadron Simhits with tof > 250ns in XY :: Muon Endcap",  n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);


  Muon_00ns_el_RZ = new TH2F("Muon_00ns_el_RZ", "Electron Simhits with tof < 50ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_00ns_mu_RZ = new TH2F("Muon_00ns_mu_RZ", "Muon Simhits with tof < 50ns in RZ :: Muon",     n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_00ns_ha_RZ = new TH2F("Muon_00ns_ha_RZ", "Hadron Simhits with tof < 50ns in RZ :: Muon",   n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_50ns_el_RZ = new TH2F("Muon_50ns_el_RZ", "Electron Simhits with 50 < tof < 250ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_50ns_mu_RZ = new TH2F("Muon_50ns_mu_RZ", "Muon Simhits with 50 < tof < 250ns in RZ :: Muon",     n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_50ns_ha_RZ = new TH2F("Muon_50ns_ha_RZ", "Hadron Simhits with 50 < tof < 250ns in RZ :: Muon",   n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_Barrel_00ns_el_XY = new TH2F("Muon_Barrel_00ns_el_XY", "Electron Simhits with tof < 50ns in XY :: Muon Barrel",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_00ns_mu_XY = new TH2F("Muon_Barrel_00ns_mu_XY", "Muon_Barrel Simhits with tof < 50ns in XY :: Muon Barrel",    n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_00ns_ha_XY = new TH2F("Muon_Barrel_00ns_ha_XY", "Hadron Simhits with tof < 50ns in XY :: Muon Barrel",  n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_50ns_el_XY = new TH2F("Muon_Barrel_50ns_el_XY", "Electron Simhits with 50 < tof < 250ns in XY :: Muon Barrel",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_50ns_mu_XY = new TH2F("Muon_Barrel_50ns_mu_XY", "Muon_Barrel Simhits with 50 < tof < 250ns in XY :: Muon Barrel",    n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Barrel_50ns_ha_XY = new TH2F("Muon_Barrel_50ns_ha_XY", "Hadron Simhits with 50 < tof < 250ns in XY :: Muon Barrel",  n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_00ns_el_XY = new TH2F("Muon_Endcap_00ns_el_XY", "Electron Simhits with tof < 50ns in XY :: Muon Endcap",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_00ns_mu_XY = new TH2F("Muon_Endcap_00ns_mu_XY", "Muon_Endcap Simhits with tof < 50ns in XY :: Muon Endcap",    n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_00ns_ha_XY = new TH2F("Muon_Endcap_00ns_ha_XY", "Hadron Simhits with tof < 50ns in XY :: Muon Endcap",  n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_50ns_el_XY = new TH2F("Muon_Endcap_50ns_el_XY", "Electron Simhits with 50 < tof < 250ns in XY :: Muon Endcap",n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_50ns_mu_XY = new TH2F("Muon_Endcap_50ns_mu_XY", "Muon_Endcap Simhits with 50 < tof < 250ns in XY :: Muon Endcap",    n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  Muon_Endcap_50ns_ha_XY = new TH2F("Muon_Endcap_50ns_ha_XY", "Hadron Simhits with 50 < tof < 250ns in XY :: Muon Endcap",  n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);


  RPCb_Muons_SHPT     = new TH1F("RPCb_Muons_SHPT", "RPCb :: All Muon Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCf_Muons_SHPT     = new TH1F("RPCf_Muons_SHPT", "RPCf :: All Muon Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  CSC_Muons_SHPT      = new TH1F("CSC_Muons_SHPT",  "CSC :: All Muon Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  DT_Muons_SHPT       = new TH1F("DT_Muons_SHPT",   "DT :: All Muon Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  GEM_Muons_SHPT      = new TH1F("GEM_Muons_SHPT",  "GEM :: All Muon Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  ME0_Muons_SHPT      = new TH1F("ME0_Muons_SHPT",  "ME0 :: All Muon Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);

  RPCb_Hadrons_SHPT   = new TH1F("RPCb_Hadrons_SHPT", "RPCb :: All Hadron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCf_Hadrons_SHPT   = new TH1F("RPCf_Hadrons_SHPT", "RPCf :: All Hadron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  CSC_Hadrons_SHPT    = new TH1F("CSC_Hadrons_SHPT",  "CSC :: All Hadron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  DT_Hadrons_SHPT     = new TH1F("DT_Hadrons_SHPT",   "DT :: All Hadron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  GEM_Hadrons_SHPT    = new TH1F("GEM_Hadrons_SHPT",  "GEM :: All Hadron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  ME0_Hadrons_SHPT    = new TH1F("ME0_Hadrons_SHPT",  "ME0 :: All Hadron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);

  RPCb_Electrons_SHPT = new TH1F("RPCb_Electrons_SHPT", "RPCb :: All Electron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCf_Electrons_SHPT = new TH1F("RPCf_Electrons_SHPT", "RPCf :: All Electron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  CSC_Electrons_SHPT  = new TH1F("CSC_Electrons_SHPT",  "CSC :: All Electron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  DT_Electrons_SHPT   = new TH1F("DT_Electrons_SHPT",   "DT :: All Electron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  GEM_Electrons_SHPT  = new TH1F("GEM_Electrons_SHPT",  "GEM :: All Electron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  ME0_Electrons_SHPT  = new TH1F("ME0_Electrons_SHPT",  "ME0 :: All Electron Hits :: SimHit Process Type", n_pro, n1_pro, n2_pro);

  RPCb_Electrons_000ns_SHPT = new TH1F("RPCb_Electrons_000ns_SHPT", "RPCb :: Electron Hits with tof < 50 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCf_Electrons_000ns_SHPT = new TH1F("RPCf_Electrons_000ns_SHPT", "RPCf :: Electron Hits with tof < 50 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  CSC_Electrons_000ns_SHPT  = new TH1F("CSC_Electrons_000ns_SHPT",  "CSC ::  Electron Hits with tof < 50 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  DT_Electrons_000ns_SHPT   = new TH1F("DT_Electrons_000ns_SHPT",   "DT ::   Electron Hits with tof < 50 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  GEM_Electrons_000ns_SHPT  = new TH1F("GEM_Electrons_000ns_SHPT",  "GEM ::  Electron Hits with tof < 50 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  ME0_Electrons_000ns_SHPT  = new TH1F("ME0_Electrons_000ns_SHPT",  "ME0 ::  Electron Hits with tof < 50 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCb_Electrons_050ns_SHPT = new TH1F("RPCb_Electrons_050ns_SHPT", "RPCb :: Electron Hits with 50 < tof < 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCf_Electrons_050ns_SHPT = new TH1F("RPCf_Electrons_050ns_SHPT", "RPCf :: Electron Hits with 50 < tof < 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  CSC_Electrons_050ns_SHPT  = new TH1F("CSC_Electrons_050ns_SHPT",  "CSC ::  Electron Hits with 50 < tof < 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  DT_Electrons_050ns_SHPT   = new TH1F("DT_Electrons_050ns_SHPT",   "DT ::   Electron Hits with 50 < tof < 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  GEM_Electrons_050ns_SHPT  = new TH1F("GEM_Electrons_050ns_SHPT",  "GEM ::  Electron Hits with 50 < tof < 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  ME0_Electrons_050ns_SHPT  = new TH1F("ME0_Electrons_050ns_SHPT",  "ME0 ::  Electron Hits with 50 < tof < 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCb_Electrons_250ns_SHPT = new TH1F("RPCb_Electrons_250ns_SHPT", "RPCb :: Electron Hits with tof > 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  RPCf_Electrons_250ns_SHPT = new TH1F("RPCf_Electrons_250ns_SHPT", "RPCf :: Electron Hits with tof > 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  CSC_Electrons_250ns_SHPT  = new TH1F("CSC_Electrons_250ns_SHPT",  "CSC ::  Electron Hits with tof > 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  DT_Electrons_250ns_SHPT   = new TH1F("DT_Electrons_250ns_SHPT",   "DT ::   Electron Hits with tof > 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  GEM_Electrons_250ns_SHPT  = new TH1F("GEM_Electrons_250ns_SHPT",  "GEM ::  Electron Hits with tof > 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);
  ME0_Electrons_250ns_SHPT  = new TH1F("ME0_Electrons_250ns_SHPT",  "ME0 ::  Electron Hits with tof > 250 ns :: SimHit Process Type", n_pro, n1_pro, n2_pro);

  RPCb_EntryExit_All_Glob_dr = new TH1F("RPCb_EntryExit_All_Glob_dr", "RPCb :: #Delta R (Global Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  DT_EntryExit_All_Glob_dr   = new TH1F("DT_EntryExit_All_Glob_dr",   "DT   :: #Delta R (Global Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  RPCf_EntryExit_All_Glob_dz = new TH1F("RPCf_EntryExit_All_Glob_dz", "RPCf :: #Delta z (Global Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  CSC_EntryExit_All_Glob_dz  = new TH1F("CSC_EntryExit_All_Glob_dz",  "CSC  :: #Delta z (Global Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  GEM_EntryExit_All_Glob_dz  = new TH1F("GEM_EntryExit_All_Glob_dz",  "GEM  :: #Delta z (Global Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  ME0_EntryExit_All_Glob_dz  = new TH1F("ME0_EntryExit_All_Glob_dz",  "ME0  :: #Delta z (Global Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);

  CSC_EntryExit_Electrons_Glob_dz  = new TH1F("CSC_EntryExit_Electrons_Glob_dz",  "CSC  :: #Delta z (Global Coords) of Electron SimHits Entry and Exit point", n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Muons_Glob_dz      = new TH1F("CSC_EntryExit_Muons_Glob_dz",      "CSC  :: #Delta z (Global Coords) of Muon SimHits Entry and Exit point", n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Hadrons_Glob_dz    = new TH1F("CSC_EntryExit_Hadrons_Glob_dz",    "CSC  :: #Delta z (Global Coords) of Hadron SimHits Entry and Exit point", n_dz, n1_dz, n2_dz);
  /*
  CSC_EntryExit_All_Deposit_dz       = new TH2F("CSC_EntryExit_All_Deposit_dz",       "CSC :: #Delta z (Global Coords) vs EE_{deposit}",              n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Electrons_Deposit_dz = new TH2F("CSC_EntryExit_Electrons_Deposit_dz", "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Electrons", n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Muons_Deposit_dz     = new TH2F("CSC_EntryExit_Muons_Deposit_dz",     "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Muons",     n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Hadrons_Deposit_dz   = new TH2F("CSC_EntryExit_Hadrons_Deposit_dz",   "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Hadrons",   n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  */
  CSC_EntryExit_All_Glob_dGap  = new TH1F("CSC_EntryExit_All_Glob_dGap",              "CSC  :: #Delta Gap (Global Coords) of SimHit Entry and Exit point", n_dG, n1_dG, n2_dG);
  CSC_EntryExit_Electrons_Glob_dGap  = new TH1F("CSC_EntryExit_Electrons_Glob_dGap",  "CSC  :: #Delta Gap (Global Coords) of Electron SimHits Entry and Exit point", n_dG, n1_dG, n2_dG);
  CSC_EntryExit_Muons_Glob_dGap      = new TH1F("CSC_EntryExit_Muons_Glob_dGap",      "CSC  :: #Delta Gap (Global Coords) of Muon SimHits Entry and Exit point", n_dG, n1_dG, n2_dG);
  CSC_EntryExit_Hadrons_Glob_dGap    = new TH1F("CSC_EntryExit_Hadrons_Glob_dGap",    "CSC  :: #Delta Gap (Global Coords) of Hadron SimHits Entry and Exit point", n_dG, n1_dG, n2_dG);

  // --------- dZ for different deposits ---------
  CSC_EntryExit_el_Deposit_dz = new TH2F("CSC_EntryExit_el_Deposit_dz", "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Electrons", n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_mu_Deposit_dz = new TH2F("CSC_EntryExit_mu_Deposit_dz", "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Muons",     n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_pi_Deposit_dz = new TH2F("CSC_EntryExit_pi_Deposit_dz", "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Pions",     n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_ka_Deposit_dz = new TH2F("CSC_EntryExit_ka_Deposit_dz", "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Kaons",     n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_p_Deposit_dz  = new TH2F("CSC_EntryExit_p_Deposit_dz",  "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Protons",   n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_n_Deposit_dz  = new TH2F("CSC_EntryExit_n_Deposit_dz",  "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Neugrons",  n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_g_Deposit_dz  = new TH2F("CSC_EntryExit_g_Deposit_dz",  "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Photons",   n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_N_Deposit_dz  = new TH2F("CSC_EntryExit_N_Deposit_dz",  "CSC :: #Delta z (Global Coords) vs E_{deposit} :: Nuclei",    n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);

  CSC_EntryExit_el_KinEn_dz = new TH2F("CSC_EntryExit_el_KinEn_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Electrons", n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_mu_KinEn_dz = new TH2F("CSC_EntryExit_mu_KinEn_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Muons",     n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_pi_KinEn_dz = new TH2F("CSC_EntryExit_pi_KinEn_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Pions",     n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_ka_KinEn_dz = new TH2F("CSC_EntryExit_ka_KinEn_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Kaons",     n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_p_KinEn_dz  = new TH2F("CSC_EntryExit_p_KinEn_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Protons",   n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_n_KinEn_dz  = new TH2F("CSC_EntryExit_n_KinEn_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Neugrons",  n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_g_KinEn_dz  = new TH2F("CSC_EntryExit_g_KinEn_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Photons",   n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_N_KinEn_dz  = new TH2F("CSC_EntryExit_N_KinEn_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Nuclei",    n_E, n1_E, n2_E, n_dz, n1_dz, n2_dz);

  CSC_EntryExit_el_Time_dz = new TH2F("CSC_EntryExit_el_Time_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Electrons", n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_mu_Time_dz = new TH2F("CSC_EntryExit_mu_Time_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Muons",     n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_pi_Time_dz = new TH2F("CSC_EntryExit_pi_Time_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Pions",     n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_ka_Time_dz = new TH2F("CSC_EntryExit_ka_Time_dz", "CSC :: #Delta z (Global Coords) vs E_{kin} :: Kaons",     n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_p_Time_dz  = new TH2F("CSC_EntryExit_p_Time_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Protons",   n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_n_Time_dz  = new TH2F("CSC_EntryExit_n_Time_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Neugrons",  n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_g_Time_dz  = new TH2F("CSC_EntryExit_g_Time_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Photons",   n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);
  CSC_EntryExit_N_Time_dz  = new TH2F("CSC_EntryExit_N_Time_dz",  "CSC :: #Delta z (Global Coords) vs E_{kin} :: Nuclei",    n_dz, n1_dz, n2_dz, n_tof, n1_tof, n2_tof);


  // --------- dR or eventually pi*dR^2 ---------
  CSC_EntryExit_el_Deposit_dR = new TH2F("CSC_EntryExit_el_Deposit_dR", "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Electrons", n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_mu_Deposit_dR = new TH2F("CSC_EntryExit_mu_Deposit_dR", "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Muons",     n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_pi_Deposit_dR = new TH2F("CSC_EntryExit_pi_Deposit_dR", "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Pions",     n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_ka_Deposit_dR = new TH2F("CSC_EntryExit_ka_Deposit_dR", "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Kaons",     n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_p_Deposit_dR  = new TH2F("CSC_EntryExit_p_Deposit_dR",  "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Protons",   n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_n_Deposit_dR  = new TH2F("CSC_EntryExit_n_Deposit_dR",  "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Neugrons",  n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_g_Deposit_dR  = new TH2F("CSC_EntryExit_g_Deposit_dR",  "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Photons",   n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_N_Deposit_dR  = new TH2F("CSC_EntryExit_N_Deposit_dR",  "CSC :: #Delta R (Global Coords) vs E_{deposit} :: Nuclei",    n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR);

  CSC_EntryExit_el_KinEn_dR = new TH2F("CSC_EntryExit_el_KinEn_dR", "CSC :: #Delta R (Global Coords) vs E_{kin} :: Electrons", n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_mu_KinEn_dR = new TH2F("CSC_EntryExit_mu_KinEn_dR", "CSC :: #Delta R (Global Coords) vs E_{kin} :: Muons",     n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_pi_KinEn_dR = new TH2F("CSC_EntryExit_pi_KinEn_dR", "CSC :: #Delta R (Global Coords) vs E_{kin} :: Pions",     n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_ka_KinEn_dR = new TH2F("CSC_EntryExit_ka_KinEn_dR", "CSC :: #Delta R (Global Coords) vs E_{kin} :: Kaons",     n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_p_KinEn_dR  = new TH2F("CSC_EntryExit_p_KinEn_dR",  "CSC :: #Delta R (Global Coords) vs E_{kin} :: Protons",   n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_n_KinEn_dR  = new TH2F("CSC_EntryExit_n_KinEn_dR",  "CSC :: #Delta R (Global Coords) vs E_{kin} :: Neugrons",  n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_g_KinEn_dR  = new TH2F("CSC_EntryExit_g_KinEn_dR",  "CSC :: #Delta R (Global Coords) vs E_{kin} :: Photons",   n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_N_KinEn_dR  = new TH2F("CSC_EntryExit_N_KinEn_dR",  "CSC :: #Delta R (Global Coords) vs E_{kin} :: Nuclei",    n_E, n1_E, n2_E, n_dR, n1_dR, n2_dR);

  RPCb_EntryExit_All_Loc_dz = new TH1F("RPCb_EntryExit_All_Loc_dz", "RPCb :: #Delta z (Local Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  RPCf_EntryExit_All_Loc_dz = new TH1F("RPCf_EntryExit_All_Loc_dz", "RPCf :: #Delta z (Local Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  CSC_EntryExit_All_Loc_dz  = new TH1F("CSC_EntryExit_All_Loc_dz",  "CSC  :: #Delta z (Local Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  DT_EntryExit_All_Loc_dz   = new TH1F("DT_EntryExit_All_Loc_dz",   "DT   :: #Delta z (Local Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  GEM_EntryExit_All_Loc_dz  = new TH1F("GEM_EntryExit_All_Loc_dz",  "GEM  :: #Delta z (Local Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);
  ME0_EntryExit_All_Loc_dz  = new TH1F("ME0_EntryExit_All_Loc_dz",  "ME0  :: #Delta z (Local Coords) of SimHit Entry and Exit point", n_dz, n1_dz, n2_dz);

  CSC_EntryExit_Electrons_Loc_dz  = new TH1F("CSC_EntryExit_Electrons_Loc_dz",  "CSC  :: #Delta z (Local Coords) of Electron SimHits Entry and Exit point", n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Muons_Loc_dz      = new TH1F("CSC_EntryExit_Muons_Loc_dz",      "CSC  :: #Delta z (Local Coords) of Muon SimHits Entry and Exit point", n_dz, n1_dz, n2_dz);
  CSC_EntryExit_Hadrons_Loc_dz    = new TH1F("CSC_EntryExit_Hadrons_Loc_dz",    "CSC  :: #Delta z (Local Coords) of Hadron SimHits Entry and Exit point", n_dz, n1_dz, n2_dz);

  // --------- dz vs dR ---------
  CSC_EntryExit_el_dz_dR = new TH2F("CSC_EntryExit_el_dz_dR", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Electrons",  n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_mu_dz_dR = new TH2F("CSC_EntryExit_mu_dz_dR", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Muons",      n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_pi_dz_dR = new TH2F("CSC_EntryExit_pi_dz_dR", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Pions",      n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_ka_dz_dR = new TH2F("CSC_EntryExit_ka_dz_dR", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Kaons",      n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_p_dz_dR  = new TH2F("CSC_EntryExit_p_dz_dR",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Protons",    n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_n_dz_dR  = new TH2F("CSC_EntryExit_n_dz_dR",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Neugrons",   n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_g_dz_dR  = new TH2F("CSC_EntryExit_g_dz_dR",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Photons",    n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);
  CSC_EntryExit_N_dz_dR  = new TH2F("CSC_EntryExit_N_dz_dR",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Nuclei",     n_dz, n1_dz, n2_dz, n_dR, n1_dR, n2_dR);

  CSC_EntryExit_el_dz_dR_detail = new TH2F("CSC_EntryExit_el_dz_dR_detail", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Electrons",  n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_mu_dz_dR_detail = new TH2F("CSC_EntryExit_mu_dz_dR_detail", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Muons",      n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_pi_dz_dR_detail = new TH2F("CSC_EntryExit_pi_dz_dR_detail", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Pions",      n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_ka_dz_dR_detail = new TH2F("CSC_EntryExit_ka_dz_dR_detail", "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Kaons",      n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_p_dz_dR_detail  = new TH2F("CSC_EntryExit_p_dz_dR_detail",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Protons",    n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_n_dz_dR_detail  = new TH2F("CSC_EntryExit_n_dz_dR_detail",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Neugrons",   n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_g_dz_dR_detail  = new TH2F("CSC_EntryExit_g_dz_dR_detail",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Photons",    n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_N_dz_dR_detail  = new TH2F("CSC_EntryExit_N_dz_dR_detail",  "CSC :: #Delta z (Global Coords) vs #Delta R (Global Coords) :: Nuclei",     n_dz, n1_dz, n2_dz2, n_dR, n1_dR, n2_dR2);

  CSC_EntryExit_el_GapLength_Deposit = new TH2F("CSC_EntryExit_el_GapLength_Deposit", "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Electrons",  n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_mu_GapLength_Deposit = new TH2F("CSC_EntryExit_mu_GapLength_Deposit", "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Muons",      n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_pi_GapLength_Deposit = new TH2F("CSC_EntryExit_pi_GapLength_Deposit", "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Pions",      n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_ka_GapLength_Deposit = new TH2F("CSC_EntryExit_ka_GapLength_Deposit", "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Kaons",      n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_p_GapLength_Deposit  = new TH2F("CSC_EntryExit_p_GapLength_Deposit",  "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Protons",    n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_n_GapLength_Deposit  = new TH2F("CSC_EntryExit_n_GapLength_Deposit",  "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Neugrons",   n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_g_GapLength_Deposit  = new TH2F("CSC_EntryExit_g_GapLength_Deposit",  "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Photons",    n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);
  CSC_EntryExit_N_GapLength_Deposit  = new TH2F("CSC_EntryExit_N_GapLength_Deposit",  "CSC :: #Delta Gap (Global Coords) vs E_{deposit} :: Nuclei",     n_D, n1_D, n2_D, n_dG, n1_dG, n2_dG);

  CSC_EntryExit_el_Deposit_pidR2 = new TH2F("CSC_EntryExit_el_Deposit_pidR2", "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Electrons",  n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_mu_Deposit_pidR2 = new TH2F("CSC_EntryExit_mu_Deposit_pidR2", "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Muons",      n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_pi_Deposit_pidR2 = new TH2F("CSC_EntryExit_pi_Deposit_pidR2", "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Pions",      n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_ka_Deposit_pidR2 = new TH2F("CSC_EntryExit_ka_Deposit_pidR2", "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Kaons",      n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_p_Deposit_pidR2  = new TH2F("CSC_EntryExit_p_Deposit_pidR2",  "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Protons",    n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_n_Deposit_pidR2  = new TH2F("CSC_EntryExit_n_Deposit_pidR2",  "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Neugrons",   n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_g_Deposit_pidR2  = new TH2F("CSC_EntryExit_g_Deposit_pidR2",  "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Photons",    n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);
  CSC_EntryExit_N_Deposit_pidR2  = new TH2F("CSC_EntryExit_N_Deposit_pidR2",  "CSC :: #pi (#Delta R / 2)^{2} (Global Coords) vs E_{deposit} :: Nuclei",     n_D, n1_D, n2_D, n_dR, n1_dR, n2_dR2);

  CSC_EntryExit_el_Deposit_dRdz = new TH2F("CSC_EntryExit_el_Deposit_dRdz", "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Electrons",  n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_mu_Deposit_dRdz = new TH2F("CSC_EntryExit_mu_Deposit_dRdz", "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Muons",      n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_pi_Deposit_dRdz = new TH2F("CSC_EntryExit_pi_Deposit_dRdz", "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Pions",      n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_ka_Deposit_dRdz = new TH2F("CSC_EntryExit_ka_Deposit_dRdz", "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Kaons",      n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_p_Deposit_dRdz  = new TH2F("CSC_EntryExit_p_Deposit_dRdz",  "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Protons",    n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_n_Deposit_dRdz  = new TH2F("CSC_EntryExit_n_Deposit_dRdz",  "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Neugrons",   n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_g_Deposit_dRdz  = new TH2F("CSC_EntryExit_g_Deposit_dRdz",  "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Photons",    n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);
  CSC_EntryExit_N_Deposit_dRdz  = new TH2F("CSC_EntryExit_N_Deposit_dRdz",  "CSC :: #sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: Nuclei",     n_D, n1_D, n2_D, n_dz, n1_dz, n2_dz);



  SimVertices_RZ      = new TH2F("SimVertices_RZ",      "Vertices in RZ",                           n_zrc_z, n_zrc_z1, n_zrc_z2, n_zrc_r, n_zrc_r1, n_zrc_r2);
  SimVertices_Muon_RZ = new TH2F("SimVertices_Muon_RZ", "Vertices in RZ in Muon System and Cavern", n_zrc_z, n_zrc_z1, n_zrc_z2, n_zrc_r, n_zrc_r1, n_zrc_r2);
  PrimVertices_Z = new TH1F("PrimVertices_Z", "Primary Vertices distribution along Z", n_pv_z, n1_pv_z, n2_pv_z);
  PrimVertices_R = new TH1F("PrimVertices_R", "Primary Vertices distribution along R", n_pv_r, n1_pv_r, n2_pv_r);

  RPC_hits  = new TH1F("RPC_hits",  "Hits in different parts of the RPC system", n_cat, n1_cat, n2_cat);
  RPC_area  = new TH1F("RPC_area",  "Surface of different parts of the RPC system", n_cat, n1_cat, n2_cat);
  RPC_rates = new TH1F("RPC_rates", "Rates in different parts of the RPC system", n_cat, n1_cat, n2_cat);

  for(int i=0; i<2; ++i) { 
    for(int j=0; j<4; ++j) {
      for(int k=0; k<3; ++k) {
	RPC_hits_array[i][j][k] = 0; 
	RPC_area_array[i][j][k] = 0.0;
	RPC_rates_array[i][j][k]= 0.0;
      }
    }
  }

  RPCb_hits_tof = new TH1F("RPCb_hits_tof", "Simhit time :: RPCb", m_tof, m1_tof, m2_tof);
  RPCb_hits_eta = new TH1F("RPCb_hits_eta", "Simhit time :: RPCb", m_eta, m1_eta, m2_eta);
  RPCb_hits_phi = new TH1F("RPCb_hits_phi", "Simhit time :: RPCb", m_phi, m1_phi, m2_phi);
  RPCb_hits_lin = new TH1F("RPCb_hits_lin", "Simhit time :: RPCb", m_lin, m1_lin, m2_lin);

  RPCf_hits_tof = new TH1F("RPCf_hits_tof", "Simhit time :: RPCf", m_tof, m1_tof, m2_tof);
  RPCf_hits_eta = new TH1F("RPCf_hits_eta", "Simhit time :: RPCf", m_eta, m1_eta, m2_eta);
  RPCf_hits_phi = new TH1F("RPCf_hits_phi", "Simhit time :: RPCf", m_phi, m1_phi, m2_phi);
  RPCf_hits_lin = new TH1F("RPCf_hits_lin", "Simhit time :: RPCf", m_lin, m1_lin, m2_lin);

  CSC_hits_tof  = new TH1F("CSC_hits_tof",  "Simhit time :: CSC",  m_tof, m1_tof, m2_tof);
  CSC_hits_eta  = new TH1F("CSC_hits_eta",  "Simhit time :: CSC",  m_eta, m1_eta, m2_eta);
  CSC_hits_phi  = new TH1F("CSC_hits_phi",  "Simhit time :: CSC",  m_phi, m1_phi, m2_phi);
  CSC_hits_lin  = new TH1F("CSC_hits_lin",  "Simhit time :: CSC",  m_lin, m1_lin, m2_lin);

  DT_hits_tof   = new TH1F("DT_hits_tof",   "Simhit time :: DT",   m_tof, m1_tof, m2_tof);
  DT_hits_eta   = new TH1F("DT_hits_eta",   "Simhit time :: DT",   m_eta, m1_eta, m2_eta);
  DT_hits_phi   = new TH1F("DT_hits_phi",   "Simhit time :: DT",   m_phi, m1_phi, m2_phi);
  DT_hits_lin   = new TH1F("DT_hits_lin",   "Simhit time :: DT",   m_lin, m1_lin, m2_lin);

  GEM_hits_tof  = new TH1F("GEM_hits_tof",  "Simhit time :: GEM",  m_tof, m1_tof, m2_tof);
  GEM_hits_eta  = new TH1F("GEM_hits_eta",  "Simhit time :: GEM",  m_eta, m1_eta, m2_eta);
  GEM_hits_phi  = new TH1F("GEM_hits_phi",  "Simhit time :: GEM",  m_phi, m1_phi, m2_phi);
  GEM_hits_lin  = new TH1F("GEM_hits_lin",  "Simhit time :: GEM",  m_lin, m1_lin, m2_lin);

  ME0_hits_tof  = new TH1F("ME0_hits_tof",  "Simhit time :: ME0",  m_tof, m1_tof, m2_tof);
  ME0_hits_eta  = new TH1F("ME0_hits_eta",  "Simhit time :: ME0",  m_eta, m1_eta, m2_eta);
  ME0_hits_phi  = new TH1F("ME0_hits_phi",  "Simhit time :: ME0",  m_phi, m1_phi, m2_phi);
  ME0_hits_lin  = new TH1F("ME0_hits_lin",  "Simhit time :: ME0",  m_lin, m1_lin, m2_lin);

  /*
  RB4_hits_tof = new TH1F("RB4_hits_tof", "Simhit time :: RB4", m_tof, m1_tof, m2_tof);
  RE4_hits_tof = new TH1F("RE4_hits_tof", "Simhit time :: RE4", m_tof, m1_tof, m2_tof);
  MB4_hits_tof = new TH1F("MB4_hits_tof", "Simhit time :: MB4", m_tof, m1_tof, m2_tof);
  ME4_hits_tof = new TH1F("ME4_hits_tof", "Simhit time :: ME4", m_tof, m1_tof, m2_tof);
  */
  /*
  RB4_hits_phi = new TH1F("RB4_hits_phi", "Simhit time :: RB4", m_phi, m1_phi, m2_phi);
  RE4_hits_phi = new TH1F("RE4_hits_phi", "Simhit time :: RE4", m_phi, m1_phi, m2_phi);
  MB4_hits_phi = new TH1F("MB4_hits_phi", "Simhit time :: MB4", m_phi, m1_phi, m2_phi);
  ME4_hits_phi = new TH1F("ME4_hits_phi", "Simhit time :: ME4", m_phi, m1_phi, m2_phi);
  */
  MB1_hits_phi = new TH1F("MB1_hits_phi", "Simhit phi :: MB1", m_phi, m1_phi, m2_phi);
  MB2_hits_phi = new TH1F("MB2_hits_phi", "Simhit phi :: MB2", m_phi, m1_phi, m2_phi);
  MB3_hits_phi = new TH1F("MB3_hits_phi", "Simhit phi :: MB3", m_phi, m1_phi, m2_phi);
  MB4_hits_phi = new TH1F("MB4_hits_phi", "Simhit phi :: MB4", m_phi, m1_phi, m2_phi);
  RB1_hits_phi = new TH1F("RB1_hits_phi", "Simhit phi :: RB1", m_phi, m1_phi, m2_phi);
  RB2_hits_phi = new TH1F("RB2_hits_phi", "Simhit phi :: RB2", m_phi, m1_phi, m2_phi);
  RB3_hits_phi = new TH1F("RB3_hits_phi", "Simhit phi :: RB3", m_phi, m1_phi, m2_phi);
  RB4_hits_phi = new TH1F("RB4_hits_phi", "Simhit phi :: RB4", m_phi, m1_phi, m2_phi);

  ME11_hits_phi = new TH1F("ME11_hits_phi", "Simhit phi :: ME11", m_phi, m1_phi, m2_phi);
  ME12_hits_phi = new TH1F("ME12_hits_phi", "Simhit phi :: ME12", m_phi, m1_phi, m2_phi);
  ME13_hits_phi = new TH1F("ME13_hits_phi", "Simhit phi :: ME13", m_phi, m1_phi, m2_phi);
  ME21_hits_phi = new TH1F("ME21_hits_phi", "Simhit phi :: ME21", m_phi, m1_phi, m2_phi);
  ME22_hits_phi = new TH1F("ME22_hits_phi", "Simhit phi :: ME22", m_phi, m1_phi, m2_phi);
  ME31_hits_phi = new TH1F("ME31_hits_phi", "Simhit phi :: ME31", m_phi, m1_phi, m2_phi);
  ME32_hits_phi = new TH1F("ME32_hits_phi", "Simhit phi :: ME32", m_phi, m1_phi, m2_phi);
  ME41_hits_phi = new TH1F("ME41_hits_phi", "Simhit phi :: ME41", m_phi, m1_phi, m2_phi);
  ME42_hits_phi = new TH1F("ME42_hits_phi", "Simhit phi :: ME42", m_phi, m1_phi, m2_phi);
  RE12_hits_phi = new TH1F("RE12_hits_phi", "Simhit phi :: RE12", m_phi, m1_phi, m2_phi);
  RE13_hits_phi = new TH1F("RE13_hits_phi", "Simhit phi :: RE13", m_phi, m1_phi, m2_phi);
  RE22_hits_phi = new TH1F("RE22_hits_phi", "Simhit phi :: RE21", m_phi, m1_phi, m2_phi);
  RE23_hits_phi = new TH1F("RE23_hits_phi", "Simhit phi :: RE22", m_phi, m1_phi, m2_phi);
  RE32_hits_phi = new TH1F("RE32_hits_phi", "Simhit phi :: RE32", m_phi, m1_phi, m2_phi);
  RE33_hits_phi = new TH1F("RE33_hits_phi", "Simhit phi :: RE33", m_phi, m1_phi, m2_phi);
  RE42_hits_phi = new TH1F("RE42_hits_phi", "Simhit phi :: RE42", m_phi, m1_phi, m2_phi);
  RE43_hits_phi = new TH1F("RE43_hits_phi", "Simhit phi :: RE43", m_phi, m1_phi, m2_phi);

  RPCb_HPC      = new TH1F("RPCb_el_HPC", "Simhits per Chamber :: RPCb", n_hits, n1_hits, n2_hits);
  RPCf_HPC      = new TH1F("RPCf_el_HPC", "Simhits per Chamber :: RPCf", n_hits, n1_hits, n2_hits);
  CSC_HPC       = new TH1F("CSC_el_HPC",  "Simhits per Chamber :: CSC",  n_hits, n1_hits, n2_hits);
  DT_HPC        = new TH1F("DT_el_HPC",   "Simhits per Chamber :: DT",   n_hits, n1_hits, n2_hits);
  RPCb_el_HPC   = new TH1F("RPCb_el_HPC", "Simhits per Chamber :: RPCb", n_hits, n1_hits, n2_hits);
  RPCf_el_HPC   = new TH1F("RPCf_el_HPC", "Simhits per Chamber :: RPCf", n_hits, n1_hits, n2_hits);
  CSC_el_HPC    = new TH1F("CSC_el_HPC",  "Simhits per Chamber :: CSC",  n_hits, n1_hits, n2_hits);
  DT_el_HPC     = new TH1F("DT_el_HPC",   "Simhits per Chamber :: DT",   n_hits, n1_hits, n2_hits);

  RPCb_el_HPL   = new TH1F("RPCb_el_HPL", "Layers hit by Electron  :: RPCb", n_lay, n1_lay, n2_lay);
  RPCf_el_HPL   = new TH1F("RPCf_el_HPL", "Layers hit by Electron  :: RPCf", n_lay, n1_lay, n2_lay);
  CSC_el_HPL    = new TH1F("CSC_el_HPL",  "Layers hit by Electron  :: CSC",  n_lay, n1_lay, n2_lay);
  DT_el_HPL     = new TH1F("DT_el_HPL",   "Layers hit by Electron  :: DT",   n_lay, n1_lay, n2_lay);
  RPCb_mu_HPL   = new TH1F("RPCb_mu_HPL", "Layers hit by Muon  :: RPCb", n_lay, n1_lay, n2_lay);
  RPCf_mu_HPL   = new TH1F("RPCf_mu_HPL", "Layers hit by Muon  :: RPCf", n_lay, n1_lay, n2_lay);
  CSC_mu_HPL    = new TH1F("CSC_mu_HPL",  "Layers hit by Muon  :: CSC",  n_lay, n1_lay, n2_lay);
  DT_mu_HPL     = new TH1F("DT_mu_HPL",   "Layers hit by Muons :: DT",   n_lay, n1_lay, n2_lay);
  RPCb_ha_HPL   = new TH1F("RPCb_ha_HPL", "Layers hit by Hadrons  :: RPCb", n_lay, n1_lay, n2_lay);
  RPCf_ha_HPL   = new TH1F("RPCf_ha_HPL", "Layers hit by Hadrons  :: RPCf", n_lay, n1_lay, n2_lay);
  CSC_ha_HPL    = new TH1F("CSC_ha_HPL",  "Layers hit by Hadrons  :: CSC",  n_lay, n1_lay, n2_lay);
  DT_ha_HPL     = new TH1F("DT_ha_HPL",   "Layers hit by Hadrons  :: DT",   n_lay, n1_lay, n2_lay);

  MB1_el_HPL = new TH1F("MB1_el_HPL",   "Layers hit by Electrons  :: DT MB1",   n_lay, n1_lay, n2_lay);
  MB2_el_HPL = new TH1F("MB2_el_HPL",   "Layers hit by Electrons  :: DT MB2",   n_lay, n1_lay, n2_lay);
  MB3_el_HPL = new TH1F("MB3_el_HPL",   "Layers hit by Electrons  :: DT MB3",   n_lay, n1_lay, n2_lay);
  MB4_el_HPL = new TH1F("MB4_el_HPL",   "Layers hit by Electrons  :: DT MB4",   n_lay, n1_lay, n2_lay);
  MB1_mu_HPL = new TH1F("MB1_mu_HPL",   "Layers hit by Muons  :: DT MB1",   n_lay, n1_lay, n2_lay);
  MB2_mu_HPL = new TH1F("MB2_mu_HPL",   "Layers hit by Muons  :: DT MB2",   n_lay, n1_lay, n2_lay);
  MB3_mu_HPL = new TH1F("MB3_mu_HPL",   "Layers hit by Muons  :: DT MB3",   n_lay, n1_lay, n2_lay);
  MB4_mu_HPL = new TH1F("MB4_mu_HPL",   "Layers hit by Muons  :: DT MB4",   n_lay, n1_lay, n2_lay);
  MB1_ha_HPL = new TH1F("MB1_ha_HPL",   "Layers hit by Hadrons  :: DT MB1",   n_lay, n1_lay, n2_lay);
  MB2_ha_HPL = new TH1F("MB2_ha_HPL",   "Layers hit by Hadrons  :: DT MB2",   n_lay, n1_lay, n2_lay);
  MB3_ha_HPL = new TH1F("MB3_ha_HPL",   "Layers hit by Hadrons  :: DT MB3",   n_lay, n1_lay, n2_lay);
  MB4_ha_HPL = new TH1F("MB4_ha_HPL",   "Layers hit by Hadrons  :: DT MB4",   n_lay, n1_lay, n2_lay);

  ME1_el_HPL = new TH1F("ME1_el_HPL",   "Layers hit by Electrons  :: CSC ME1",   n_lay, n1_lay, n2_lay);
  ME2_el_HPL = new TH1F("ME2_el_HPL",   "Layers hit by Electrons  :: CSC ME2",   n_lay, n1_lay, n2_lay);
  ME3_el_HPL = new TH1F("ME3_el_HPL",   "Layers hit by Electrons  :: CSC ME3",   n_lay, n1_lay, n2_lay);
  ME4_el_HPL = new TH1F("ME4_el_HPL",   "Layers hit by Electrons  :: CSC ME4",   n_lay, n1_lay, n2_lay);
  ME1_mu_HPL = new TH1F("ME1_mu_HPL",   "Layers hit by Muons  :: CSC ME1",   n_lay, n1_lay, n2_lay);
  ME2_mu_HPL = new TH1F("ME2_mu_HPL",   "Layers hit by Muons  :: CSC ME2",   n_lay, n1_lay, n2_lay);
  ME3_mu_HPL = new TH1F("ME3_mu_HPL",   "Layers hit by Muons  :: CSC ME3",   n_lay, n1_lay, n2_lay);
  ME4_mu_HPL = new TH1F("ME4_mu_HPL",   "Layers hit by Muons  :: CSC ME4",   n_lay, n1_lay, n2_lay);
  ME1_ha_HPL = new TH1F("ME1_ha_HPL",   "Layers hit by Hadrons  :: CSC ME1",   n_lay, n1_lay, n2_lay);
  ME2_ha_HPL = new TH1F("ME2_ha_HPL",   "Layers hit by Hadrons  :: CSC ME2",   n_lay, n1_lay, n2_lay);
  ME3_ha_HPL = new TH1F("ME3_ha_HPL",   "Layers hit by Hadrons  :: CSC ME3",   n_lay, n1_lay, n2_lay);
  ME4_ha_HPL = new TH1F("ME4_ha_HPL",   "Layers hit by Hadrons  :: CSC ME4",   n_lay, n1_lay, n2_lay);

  CSC_ME1_area  = new TH1F("CSC_ME1_area",  "Area of CSC ME1 as function of R",  r_CSC, r1_CSC, r2_CSC);
  CSC_ME2_area  = new TH1F("CSC_ME2_area",  "Area of CSC ME2 as function of R",  r_CSC, r1_CSC, r2_CSC);
  CSC_ME3_area  = new TH1F("CSC_ME3_area",  "Area of CSC ME3 as function of R",  r_CSC, r1_CSC, r2_CSC);
  CSC_ME4_area  = new TH1F("CSC_ME4_area",  "Area of CSC ME4 as function of R",  r_CSC, r1_CSC, r2_CSC);

  CSC_ME1_all_hits    = new TH1F("CSC_ME1_all_hits",    "Hits in CSC ME1 as function of R [All]",       r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_000ns_hits  = new TH1F("CSC_ME1_000ns_hits",  "Hits in CSC ME1 as function of R [<250ns]",    r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_250ns_hits  = new TH1F("CSC_ME1_250ns_hits",  "Hits in CSC ME1 as function of R [>250ns]",    r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_00ns_hits   = new TH1F("CSC_ME1_00ns_hits",   "Hits in CSC ME1 as function of R [0-50ns]",    r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_50ns_hits   = new TH1F("CSC_ME1_50ns_hits",   "Hits in CSC ME1 as function of R [50-250ns]",  r_CSC, r1_CSC, r2_CSC);

  CSC_ME1_all_rates   = new TH1F("CSC_ME1_all_rates",   "Rates of CSC ME1 as function of R [All]",      r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_000ns_rates = new TH1F("CSC_ME1_000ns_rates", "Rates of CSC ME1 as function of R [<250ns]",   r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_250ns_rates = new TH1F("CSC_ME1_250ns_rates", "Rates of CSC ME1 as function of R [>250ns]",   r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_00ns_rates  = new TH1F("CSC_ME1_00ns_rates",  "Rates of CSC ME1 as function of R [0-50ns]",   r_CSC, r1_CSC, r2_CSC);
  CSC_ME1_50ns_rates  = new TH1F("CSC_ME1_50ns_rates",  "Rates of CSC ME1 as function of R [50-250ns]", r_CSC, r1_CSC, r2_CSC);

  CSC_ME2_all_hits    = new TH1F("CSC_ME2_all_hits",    "Hits in CSC ME2 as function of R [All]",       r_CSC, r1_CSC, r2_CSC);
  CSC_ME2_all_rates   = new TH1F("CSC_ME2_all_rates",   "Rates of CSC ME2 as function of R [All]",      r_CSC, r1_CSC, r2_CSC);

  CSC_ME3_all_hits    = new TH1F("CSC_ME3_all_hits",    "Hits in CSC ME3 as function of R [All]",       r_CSC, r1_CSC, r2_CSC);
  CSC_ME3_all_rates   = new TH1F("CSC_ME3_all_rates",   "Rates of CSC ME3 as function of R [All]",      r_CSC, r1_CSC, r2_CSC);

  CSC_ME4_all_hits    = new TH1F("CSC_ME4_all_hits",    "Hits in CSC ME4 as function of R [All]",       r_CSC, r1_CSC, r2_CSC);
  CSC_ME4_all_rates   = new TH1F("CSC_ME4_all_rates",   "Rates of CSC ME4 as function of R [All]",      r_CSC, r1_CSC, r2_CSC);

  CSC_Geom_ME1_area  = new TH1F("CSC_Geom_ME1_area",  "Area of CSC ME1 as function of R",  51, r_ME1_geom);
  CSC_Geom_ME2_area  = new TH1F("CSC_Geom_ME2_area",  "Area of CSC ME2 as function of R",  53, r_MEX_geom);
  CSC_Geom_ME3_area  = new TH1F("CSC_Geom_ME3_area",  "Area of CSC ME3 as function of R",  53, r_MEX_geom);
  CSC_Geom_ME4_area  = new TH1F("CSC_Geom_ME4_area",  "Area of CSC ME4 as function of R",  53, r_MEX_geom);

  CSC_Geom_ME1_all_hits    = new TH1F("CSC_Geom_ME1_all_hits",    "Hits in CSC ME1 as function of R [All]", 51, r_ME1_geom);
  CSC_Geom_ME2_all_hits    = new TH1F("CSC_Geom_ME2_all_hits",    "Hits in CSC ME2 as function of R [All]", 53, r_MEX_geom);
  CSC_Geom_ME3_all_hits    = new TH1F("CSC_Geom_ME3_all_hits",    "Hits in CSC ME3 as function of R [All]", 53, r_MEX_geom);
  CSC_Geom_ME4_all_hits    = new TH1F("CSC_Geom_ME4_all_hits",    "Hits in CSC ME4 as function of R [All]", 53, r_MEX_geom);

  CSC_Geom_ME1_all_rates   = new TH1F("CSC_Geom_ME1_all_rates",   "Rates of CSC ME1 as function of R [All]", 51, r_ME1_geom);
  CSC_Geom_ME2_all_rates   = new TH1F("CSC_Geom_ME2_all_rates",   "Rates of CSC ME2 as function of R [All]", 53, r_MEX_geom);
  CSC_Geom_ME3_all_rates   = new TH1F("CSC_Geom_ME3_all_rates",   "Rates of CSC ME3 as function of R [All]", 53, r_MEX_geom);
  CSC_Geom_ME4_all_rates   = new TH1F("CSC_Geom_ME4_all_rates",   "Rates of CSC ME4 as function of R [All]", 53, r_MEX_geom);
}


MyNeutronSimHitAnalyzer::~MyNeutronSimHitAnalyzer(){

  if(tech_debug) std::cout<<"[MyNeutronSimHitAnalyzer :: Destructor]"<<std::endl; 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  TStyle *plain  = new TStyle("Plain","Plain Style (no colours/fill areas)");
  plain->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  TCanvas * Dummy = new TCanvas("dummy", "dummy", 600, 600);

  if(pdf_output) {
    // First Plot for PDF File :: print empty Dummy 
    pdfFileName = pdfFileNameBase + ".pdf[";
    Dummy->Print(pdfFileName.c_str());
    // Name for next plot for PDF File
    pdfFileName = pdfFileNameBase + ".pdf";
  }

  // Center of Mass Energy Label
  std::stringstream comlabelss; comlabelss<<"CMS Simulation #sqrt{s} = "<<comenergy<<" TeV"; std::string comlabel = comlabelss.str();

  // Legends
  double l1_x1, l1_y1, l1_x2, l1_y2;
  double l2_x1, l2_y1, l2_x2, l2_y2;
  double l3_x1, l3_y1, l3_x2, l3_y2;
  double l4_x1, l4_y1, l4_x2, l4_y2;
  // First version
  // l1_x1 = 0.60; l1_x2 = 0.85; l1_y1 = 0.60; l1_y2 = 0.85;
  // Second version
  l1_x1 = 0.65; l1_x2 = 0.85; l1_y1 = 0.65; l1_y2 = 0.85;
  TLegend *l1 = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  l1->SetLineColor(0);    l1->SetLineStyle(0);  l1->SetLineWidth(0);
  l1->SetFillColor(4000); l1->SetBorderSize(0); l1->SetNColumns(2);
  // First version, with neutrons and gammas included
  /*
  l1->AddEntry(RPCb_el_hits, "e","p");
  l1->AddEntry(RPCb_p_hits,  "p","p");
  l1->AddEntry(RPCb_mu_hits, "#mu","p");
  l1->AddEntry(RPCb_n_hits,  "n","p");
  l1->AddEntry(RPCb_pi_hits, "#pi","p");
  l1->AddEntry(RPCb_g_hits,  "#gamma","p");
  l1->AddEntry(RPCb_ka_hits, "K","p");
  l1->AddEntry(RPCb_N_hits,  "Nuclei","p");
  */
  // second version 
  l1->AddEntry(RPCb_el_hits, "e","p");
  l1->AddEntry(RPCb_pi_hits, "#pi","p");
  l1->AddEntry(RPCb_mu_hits, "#mu","p");
  l1->AddEntry(RPCb_ka_hits, "K","p");
  l1->AddEntry(RPCb_p_hits,  "p","p");
  l1->AddEntry(RPCb_N_hits,  "Nuclei","p");

  TLegend *l2 = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  l2->SetLineColor(0);    l2->SetLineStyle(0);  l2->SetLineWidth(0);
  l2->SetFillColor(4000); l2->SetBorderSize(0); l2->SetNColumns(1);
  l2->AddEntry(RPCb_el_deps, "e","l");
  l2->AddEntry(RPCb_mu_deps, "#mu","l");
  l2->AddEntry(RPCb_ha_deps, "had","l");


  // Some visual help and comments
  double log250ns =  log10(250);
  TH1F * line_250ns_hits = new TH1F("line_250ns_hits", "", n_E, n1_E, n2_E); for(int i=0; i<n_E; ++i) { line_250ns_hits->SetBinContent(i+1,  log250ns); }
  TH1F * line_250ns_deps = new TH1F("line_250ns_deps", "", n_D, n1_D, n2_D); for(int i=0; i<n_D; ++i) { line_250ns_deps->SetBinContent(i+1,  log250ns); }
  line_250ns_hits->SetLineColor(kBlack); line_250ns_hits->SetLineStyle(2); line_250ns_hits->SetLineWidth(1);
  line_250ns_deps->SetLineColor(kBlack); line_250ns_deps->SetLineStyle(2); line_250ns_deps->SetLineWidth(1);
  double logMax =  log10(maxsimtime);
  TH1F * line_max_hits = new TH1F("line_max_hits", "", n_E, n1_E, n2_E); for(int i=0; i<n_E; ++i) { line_max_hits->SetBinContent(i+1,  logMax); }
  TH1F * line_max_deps = new TH1F("line_max_deps", "", n_D, n1_D, n2_D); for(int i=0; i<n_D; ++i) { line_max_deps->SetBinContent(i+1,  logMax); }
  line_max_hits->SetLineColor(kRed); line_max_hits->SetLineStyle(2); line_max_hits->SetLineWidth(1);
  line_max_deps->SetLineColor(kRed); line_max_deps->SetLineStyle(2); line_max_deps->SetLineWidth(1);
  TLatex latex_left;   latex_left.SetNDC();   latex_left.SetTextSize(0.03);    latex_left.SetTextAlign(11);    latex_left.SetTextColor(2);
  TLatex latex_right;  latex_right.SetNDC();  latex_right.SetTextSize(0.03);   latex_right.SetTextAlign(31);
  TLatex latex_cmslab; latex_cmslab.SetNDC(); latex_cmslab.SetTextSize(0.035); latex_cmslab.SetTextAlign(11);
  TLatex latex_legend; latex_legend.SetNDC(); latex_legend.SetTextSize(0.04);  latex_legend.SetTextAlign(31);

  // Old Settings ... ranging log(TOF) from 1 to 8 (100ms)
  // double lab_250ns_x  = 0.09, lab_250ns_y  = 0.250;
  // double lab_1ms_x    = 0.06, lab_1ms_y    = 0.665;
  // double lab_1MeV_x   = 0.50, lab_1MeV_y   = 0.025;
  // double lab_prompt_x = 0.89, lab_prompt_y = 0.225; 
  // double lab_neutr_x  = 0.89, lab_neutr_y  = 0.275; 
  // New Settings ... ranging log(TOF) from 1 to 11 (100s)
  // double lab_250ns_x  = 0.09, lab_250ns_y  = 0.205;
  // double lab_1ms_x    = 0.06, lab_1ms_y    = 0.495;
  // double lab_1MeV_x   = 0.50, lab_1MeV_y   = 0.025;
  // double lab_prompt_x = 0.89, lab_prompt_y = 0.185; 
  // double lab_neutr_x  = 0.89, lab_neutr_y  = 0.225; 
  // Newest Settings ... ranging log(TOF) from 1 to 12 (1000s)
  double lab_250ns_x  = 0.09, lab_250ns_y  = 0.1925;
  double lab_1ms_x    = 0.07, lab_1ms_y    = 0.4585;
  double lab_1MeV_x   = 0.50, lab_1MeV_y   = 0.0250;
  double lab_prompt_x = 0.89, lab_prompt_y = 0.1750; 
  double lab_neutr_x  = 0.89, lab_neutr_y  = 0.2100; 
  // 
  double lab_max_x = 0.0, lab_max_y = 0.0; std::string lab_time = "";
  // New Settings ... ranging log(TOF) from 1 to 11 (100s)
  // if(maxsimtime==100000000000)  {lab_max_x = 0.15; lab_max_y = 0.825; lab_time = "#font[12]{max simulation time = 100s}";}  // case 100s 
  // if(maxsimtime==10000000000)   {lab_max_x = 0.15; lab_max_y = 0.825; lab_time = "#font[12]{max simulation time = 10s}";}   // case 10s 
  // if(maxsimtime==1000000000)    {lab_max_x = 0.15; lab_max_y = 0.750; lab_time = "#font[12]{max simulation time = 1s}";}    // case 1s 
  // if(maxsimtime==100000000)     {lab_max_x = 0.15; lab_max_y = 0.675; lab_time = "#font[12]{max simulation time = 100ms}";} // case 100ms 
  // Newest Settings ... ranging log(TOF) from 1 to 12 (1000s)
  if(maxsimtime==1000000000000) {lab_max_x = 0.125; lab_max_y = 0.830; lab_time = "#font[12]{max simulation time = 1000s}";}  // case 1000s 
  if(maxsimtime==100000000000)  {lab_max_x = 0.125; lab_max_y = 0.800; lab_time = "#font[12]{max simulation time = 100s}";}   // case 100s 
  if(maxsimtime==10000000000)   {lab_max_x = 0.125; lab_max_y = 0.730; lab_time = "#font[12]{max simulation time = 10s}";}    // case 10s 
  if(maxsimtime==1000000000)    {lab_max_x = 0.125; lab_max_y = 0.685; lab_time = "#font[12]{max simulation time = 1s}";}     // case 1s 
  if(maxsimtime==100000000)     {lab_max_x = 0.125; lab_max_y = 0.620; lab_time = "#font[12]{max simulation time = 100ms}";}  // case 100ms 

  Canvas_RPCb_hits = new TCanvas("Canvas_RPCb_hits", "Simhit time vs E_{kin} :: RPCb", 600, 600);
  Canvas_RPCf_hits = new TCanvas("Canvas_RPCf_hits", "Simhit time vs E_{kin} :: RPCf", 600, 600);
  Canvas_CSC_hits  = new TCanvas("Canvas_CSC_hits",  "Simhit time vs E_{kin} :: CSC",  600, 600);
  Canvas_DT_hits   = new TCanvas("Canvas_DT_hits",   "Simhit time vs E_{kin} :: DT",   600, 600);
  Canvas_GEM_hits  = new TCanvas("Canvas_GEM_hits",  "Simhit time vs E_{kin} :: GEM",   600, 600);
  Canvas_ME0_hits  = new TCanvas("Canvas_ME0_hits",  "Simhit time vs E_{kin} :: ME0",   600, 600);

  // Combine histograms in single plot
  Canvas_RPCb_hits->cd();
  RPCb_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  RPCb_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  RPCb_el_hits->GetXaxis()->SetTitleOffset(1.2);
  RPCb_el_hits->SetTitle("SimHit time vs E_{kin} :: RPCb");
  RPCb_el_hits->SetMarkerStyle(7);  RPCb_el_hits->SetMarkerColor(kBlack);   RPCb_el_hits->SetMarkerSize(1);  RPCb_el_hits->Draw("P");
  RPCb_mu_hits->SetMarkerStyle(24); RPCb_mu_hits->SetMarkerColor(kBlue);    RPCb_mu_hits->SetMarkerSize(1);  RPCb_mu_hits->Draw("PSame");
  RPCb_pi_hits->SetMarkerStyle(33); RPCb_pi_hits->SetMarkerColor(kGreen);   RPCb_pi_hits->SetMarkerSize(1);  RPCb_pi_hits->Draw("PSame");
  RPCb_ka_hits->SetMarkerStyle(5);  RPCb_ka_hits->SetMarkerColor(kOrange);  RPCb_ka_hits->SetMarkerSize(1);  RPCb_ka_hits->Draw("PSame");
  RPCb_p_hits->SetMarkerStyle(26);  RPCb_p_hits->SetMarkerColor(kMagenta);  RPCb_p_hits->SetMarkerSize(1);   RPCb_p_hits->Draw("PSame");
  RPCb_n_hits->SetMarkerStyle(32);  RPCb_n_hits->SetMarkerColor(kViolet);   RPCb_n_hits->SetMarkerSize(1);   RPCb_n_hits->Draw("PSame");
  RPCb_g_hits->SetMarkerStyle(30);  RPCb_g_hits->SetMarkerColor(kCyan);     RPCb_g_hits->SetMarkerSize(1);   RPCb_g_hits->Draw("PSame");
  RPCb_N_hits->SetMarkerStyle(2);   RPCb_N_hits->SetMarkerColor(kRed);      RPCb_N_hits->SetMarkerSize(1);   RPCb_N_hits->Draw("PSame");
  line_250ns_hits->Draw("][same");  line_max_hits->Draw("][same");
  l1->Draw();

  // latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCb");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");                  // these values will have to be rechecked
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // indicate also 1 second
  latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());

  Canvas_RPCb_hits->SetTicks(1,1);
  Canvas_RPCb_hits->Write();
  if(pdf_output) {Canvas_RPCb_hits->Print(pdfFileName.c_str());}

  Canvas_RPCf_hits->cd();
  RPCf_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  RPCf_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  RPCf_el_hits->GetXaxis()->SetTitleOffset(1.2);
  RPCf_el_hits->SetTitle("SimHit time vs E_{kin} :: RPCf");
  RPCf_el_hits->SetMarkerStyle(7);  RPCf_el_hits->SetMarkerColor(kBlack);   RPCf_el_hits->SetMarkerSize(1);  RPCf_el_hits->Draw("P");
  RPCf_mu_hits->SetMarkerStyle(24); RPCf_mu_hits->SetMarkerColor(kBlue);    RPCf_mu_hits->SetMarkerSize(1);  RPCf_mu_hits->Draw("PSame");
  RPCf_pi_hits->SetMarkerStyle(33); RPCf_pi_hits->SetMarkerColor(kGreen);   RPCf_pi_hits->SetMarkerSize(1);  RPCf_pi_hits->Draw("PSame");
  RPCf_ka_hits->SetMarkerStyle(5);  RPCf_ka_hits->SetMarkerColor(kOrange);  RPCf_ka_hits->SetMarkerSize(1);  RPCf_ka_hits->Draw("PSame");
  RPCf_p_hits->SetMarkerStyle(26);  RPCf_p_hits->SetMarkerColor(kMagenta);  RPCf_p_hits->SetMarkerSize(1);   RPCf_p_hits->Draw("PSame");
  RPCf_n_hits->SetMarkerStyle(32);  RPCf_n_hits->SetMarkerColor(kViolet);   RPCf_n_hits->SetMarkerSize(1);   RPCf_n_hits->Draw("PSame");
  RPCf_g_hits->SetMarkerStyle(30);  RPCf_g_hits->SetMarkerColor(kCyan);     RPCf_g_hits->SetMarkerSize(1);   RPCf_g_hits->Draw("PSame");
  RPCf_N_hits->SetMarkerStyle(2);   RPCf_N_hits->SetMarkerColor(kRed);      RPCf_N_hits->SetMarkerSize(1);   RPCf_N_hits->Draw("PSame");
  line_250ns_hits->Draw("][same");  line_max_hits->Draw("][same");
  l1->Draw();
  // latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCf");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_RPCf_hits->SetTicks(1,1);
  Canvas_RPCf_hits->Write();
  if(pdf_output) {Canvas_RPCf_hits->Print(pdfFileName.c_str());}

  Canvas_CSC_hits->cd();
  CSC_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  CSC_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  CSC_el_hits->GetXaxis()->SetTitleOffset(1.2);
  CSC_el_hits->SetTitle("SimHit time vs E_{kin} :: CSC");
  CSC_el_hits->SetMarkerStyle(7);  CSC_el_hits->SetMarkerColor(kBlack);   CSC_el_hits->SetMarkerSize(1);  CSC_el_hits->Draw("P");
  CSC_mu_hits->SetMarkerStyle(24); CSC_mu_hits->SetMarkerColor(kBlue);    CSC_mu_hits->SetMarkerSize(1);  CSC_mu_hits->Draw("PSame");
  CSC_pi_hits->SetMarkerStyle(33); CSC_pi_hits->SetMarkerColor(kGreen);   CSC_pi_hits->SetMarkerSize(1);  CSC_pi_hits->Draw("PSame");
  CSC_ka_hits->SetMarkerStyle(5);  CSC_ka_hits->SetMarkerColor(kOrange);  CSC_ka_hits->SetMarkerSize(1);  CSC_ka_hits->Draw("PSame");
  CSC_p_hits->SetMarkerStyle(26);  CSC_p_hits->SetMarkerColor(kMagenta);  CSC_p_hits->SetMarkerSize(1);   CSC_p_hits->Draw("PSame");
  CSC_n_hits->SetMarkerStyle(32);  CSC_n_hits->SetMarkerColor(kViolet);   CSC_n_hits->SetMarkerSize(1);   CSC_n_hits->Draw("PSame");
  CSC_g_hits->SetMarkerStyle(30);  CSC_g_hits->SetMarkerColor(kCyan);     CSC_g_hits->SetMarkerSize(1);   CSC_g_hits->Draw("PSame");
  CSC_N_hits->SetMarkerStyle(2);   CSC_N_hits->SetMarkerColor(kRed);      CSC_N_hits->SetMarkerSize(1);   CSC_N_hits->Draw("PSame");
  line_250ns_hits->Draw("][same"); line_max_hits->Draw("][same");
  l1->Draw();
  // latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_hits->SetTicks(1,1);
  Canvas_CSC_hits->Write();
  if(pdf_output) {Canvas_CSC_hits->Print(pdfFileName.c_str());}


  Canvas_DT_hits->cd();
  DT_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  DT_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  DT_el_hits->GetXaxis()->SetTitleOffset(1.2);
  DT_el_hits->SetTitle("SimHit time vs E_{kin} :: DT");
  DT_el_hits->SetMarkerStyle(7);  DT_el_hits->SetMarkerColor(kBlack);   DT_el_hits->SetMarkerSize(1);  DT_el_hits->Draw("P");
  DT_mu_hits->SetMarkerStyle(24); DT_mu_hits->SetMarkerColor(kBlue);    DT_mu_hits->SetMarkerSize(1);  DT_mu_hits->Draw("PSame");
  DT_pi_hits->SetMarkerStyle(33); DT_pi_hits->SetMarkerColor(kGreen);   DT_pi_hits->SetMarkerSize(1);  DT_pi_hits->Draw("PSame");
  DT_ka_hits->SetMarkerStyle(5);  DT_ka_hits->SetMarkerColor(kOrange);  DT_ka_hits->SetMarkerSize(1);  DT_ka_hits->Draw("PSame");
  DT_p_hits->SetMarkerStyle(26);  DT_p_hits->SetMarkerColor(kMagenta);  DT_p_hits->SetMarkerSize(1);   DT_p_hits->Draw("PSame");
  DT_n_hits->SetMarkerStyle(32);  DT_n_hits->SetMarkerColor(kViolet);   DT_n_hits->SetMarkerSize(1);   DT_n_hits->Draw("PSame");
  DT_g_hits->SetMarkerStyle(30);  DT_g_hits->SetMarkerColor(kCyan);     DT_g_hits->SetMarkerSize(1);   DT_g_hits->Draw("PSame");
  DT_N_hits->SetMarkerStyle(2);   DT_N_hits->SetMarkerColor(kRed);      DT_N_hits->SetMarkerSize(1);   DT_N_hits->Draw("PSame");
  line_250ns_hits->Draw("][same"); line_max_hits->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"DT");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_DT_hits->SetTicks(1,1);
  Canvas_DT_hits->Write();
  if(pdf_output) {Canvas_DT_hits->Print(pdfFileName.c_str());}

  Canvas_GEM_hits->cd();
  GEM_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  GEM_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  GEM_el_hits->GetXaxis()->SetTitleOffset(1.2);
  GEM_el_hits->SetTitle("SimHit time vs E_{kin} :: GEM");
  GEM_el_hits->SetMarkerStyle(7);  GEM_el_hits->SetMarkerColor(kBlack);   GEM_el_hits->SetMarkerSize(1);  GEM_el_hits->Draw("P");
  GEM_mu_hits->SetMarkerStyle(24); GEM_mu_hits->SetMarkerColor(kBlue);    GEM_mu_hits->SetMarkerSize(1);  GEM_mu_hits->Draw("PSame");
  GEM_pi_hits->SetMarkerStyle(33); GEM_pi_hits->SetMarkerColor(kGreen);   GEM_pi_hits->SetMarkerSize(1);  GEM_pi_hits->Draw("PSame");
  GEM_ka_hits->SetMarkerStyle(5);  GEM_ka_hits->SetMarkerColor(kOrange);  GEM_ka_hits->SetMarkerSize(1);  GEM_ka_hits->Draw("PSame");
  GEM_p_hits->SetMarkerStyle(26);  GEM_p_hits->SetMarkerColor(kMagenta);  GEM_p_hits->SetMarkerSize(1);   GEM_p_hits->Draw("PSame");
  GEM_n_hits->SetMarkerStyle(32);  GEM_n_hits->SetMarkerColor(kViolet);   GEM_n_hits->SetMarkerSize(1);   GEM_n_hits->Draw("PSame");
  GEM_g_hits->SetMarkerStyle(30);  GEM_g_hits->SetMarkerColor(kCyan);     GEM_g_hits->SetMarkerSize(1);   GEM_g_hits->Draw("PSame");
  GEM_N_hits->SetMarkerStyle(2);   GEM_N_hits->SetMarkerColor(kRed);      GEM_N_hits->SetMarkerSize(1);   GEM_N_hits->Draw("PSame");
  line_250ns_hits->Draw("][same"); line_max_hits->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"GEM");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_GEM_hits->SetTicks(1,1);
  Canvas_GEM_hits->Write();
  if(pdf_output) {Canvas_GEM_hits->Print(pdfFileName.c_str());}

  Canvas_ME0_hits->cd();
  ME0_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  ME0_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  ME0_el_hits->GetXaxis()->SetTitleOffset(1.2);
  ME0_el_hits->SetTitle("SimHit time vs E_{kin} :: ME0");
  ME0_el_hits->SetMarkerStyle(7);  ME0_el_hits->SetMarkerColor(kBlack);   ME0_el_hits->SetMarkerSize(1);  ME0_el_hits->Draw("P");
  ME0_mu_hits->SetMarkerStyle(24); ME0_mu_hits->SetMarkerColor(kBlue);    ME0_mu_hits->SetMarkerSize(1);  ME0_mu_hits->Draw("PSame");
  ME0_pi_hits->SetMarkerStyle(33); ME0_pi_hits->SetMarkerColor(kGreen);   ME0_pi_hits->SetMarkerSize(1);  ME0_pi_hits->Draw("PSame");
  ME0_ka_hits->SetMarkerStyle(5);  ME0_ka_hits->SetMarkerColor(kOrange);  ME0_ka_hits->SetMarkerSize(1);  ME0_ka_hits->Draw("PSame");
  ME0_p_hits->SetMarkerStyle(26);  ME0_p_hits->SetMarkerColor(kMagenta);  ME0_p_hits->SetMarkerSize(1);   ME0_p_hits->Draw("PSame");
  ME0_n_hits->SetMarkerStyle(32);  ME0_n_hits->SetMarkerColor(kViolet);   ME0_n_hits->SetMarkerSize(1);   ME0_n_hits->Draw("PSame");
  ME0_g_hits->SetMarkerStyle(30);  ME0_g_hits->SetMarkerColor(kCyan);     ME0_g_hits->SetMarkerSize(1);   ME0_g_hits->Draw("PSame");
  ME0_N_hits->SetMarkerStyle(2);   ME0_N_hits->SetMarkerColor(kRed);      ME0_N_hits->SetMarkerSize(1);   ME0_N_hits->Draw("PSame");
  line_250ns_hits->Draw("][same"); line_max_hits->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"ME0");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_ME0_hits->SetTicks(1,1);
  Canvas_ME0_hits->Write();
  if(pdf_output) {Canvas_ME0_hits->Print(pdfFileName.c_str());}

  Canvas_RPCb_hits_fancy = new TCanvas("Canvas_RPCb_hits_fancy", "Simhit time vs E_{kin} :: RPCb", 600, 600);
  Canvas_RPCf_hits_fancy = new TCanvas("Canvas_RPCf_hits_fancy", "Simhit time vs E_{kin} :: RPCf", 600, 600);
  Canvas_CSC_hits_fancy  = new TCanvas("Canvas_CSC_hits_fancy",  "Simhit time vs E_{kin} :: CSC",  600, 600);
  Canvas_DT_hits_fancy   = new TCanvas("Canvas_DT_hits_fancy",   "Simhit time vs E_{kin} :: DT",   600, 600);
  Canvas_GEM_hits_fancy  = new TCanvas("Canvas_GEM_hits_fancy",  "Simhit time vs E_{kin} :: GEM",   600, 600);
  Canvas_ME0_hits_fancy  = new TCanvas("Canvas_ME0_hits_fancy",  "Simhit time vs E_{kin} :: ME0",   600, 600);

  Canvas_CSC_hits_fancy->cd();
  CSC_el_hits->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  CSC_el_hits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  CSC_el_hits->GetXaxis()->SetTitleOffset(1.2);
  CSC_el_hits->SetTitle("SimHit time vs E_{kin} :: CSC");
  CSC_el_hits->SetMarkerStyle(7);  CSC_el_hits->SetMarkerColor(kBlack);   CSC_el_hits->SetMarkerSize(1);  CSC_el_hits->Draw("P");
  CSC_mu_hits->SetMarkerStyle(24); CSC_mu_hits->SetMarkerColor(kBlue);    CSC_mu_hits->SetMarkerSize(1);  CSC_mu_hits->Draw("PSame");
  CSC_pi_hits->SetMarkerStyle(33); CSC_pi_hits->SetMarkerColor(kGreen);   CSC_pi_hits->SetMarkerSize(1);  CSC_pi_hits->Draw("PSame");
  CSC_ka_hits->SetMarkerStyle(5);  CSC_ka_hits->SetMarkerColor(kOrange);  CSC_ka_hits->SetMarkerSize(1);  CSC_ka_hits->Draw("PSame");
  CSC_p_hits->SetMarkerStyle(26);  CSC_p_hits->SetMarkerColor(kMagenta);  CSC_p_hits->SetMarkerSize(1);   CSC_p_hits->Draw("PSame");
  CSC_n_hits->SetMarkerStyle(32);  CSC_n_hits->SetMarkerColor(kViolet);   CSC_n_hits->SetMarkerSize(1);   CSC_n_hits->Draw("PSame");
  CSC_g_hits->SetMarkerStyle(30);  CSC_g_hits->SetMarkerColor(kCyan);     CSC_g_hits->SetMarkerSize(1);   CSC_g_hits->Draw("PSame");
  CSC_N_hits->SetMarkerStyle(2);   CSC_N_hits->SetMarkerColor(kRed);      CSC_N_hits->SetMarkerSize(1);   CSC_N_hits->Draw("PSame");
  // line_250ns_hits->Draw("][same"); line_max_hits->Draw("][same");
  // l1->Draw();
  // latex_cmslab.DrawLatex(0.10, 0.925,"CMS Simulation #sqrt{s} = 8 TeV");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  // latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_1MeV_x,lab_1MeV_y,    "#font[12]{1 MeV}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_hits->SetTicks(1,1);
  // kill X-axis label on big plot
  CSC_el_hits->GetXaxis()->SetTitle("");
  CSC_el_hits->GetXaxis()->SetLabelSize(0.001);
  double canvasratio = 0.25;
  // bottom canvas
  Canvas_CSC_hits_fancy->SetBottomMargin(canvasratio + (1-canvasratio)*Canvas_CSC_hits_fancy->GetBottomMargin()-canvasratio*Canvas_CSC_hits_fancy->GetTopMargin());
  TPad *Canvas_CSC_hits_fancy_1 = new TPad("BottomPad","",0,0,1,1);
  Canvas_CSC_hits_fancy_1->SetTopMargin( (1-canvasratio) - (1-canvasratio)*Canvas_CSC_hits_fancy_1->GetBottomMargin() + canvasratio*Canvas_CSC_hits_fancy_1->GetTopMargin());
  Canvas_CSC_hits_fancy_1->SetRightMargin( canvasratio   + (1-canvasratio)*Canvas_CSC_hits_fancy_1->GetRightMargin()  - canvasratio*Canvas_CSC_hits_fancy_1->GetLeftMargin());
  // original and new values
  // double originalbottom = Canvas_CSC_hits_fancy->GetBottomMargin();
  // double originalright = Canvas_CSC_hits_fancy->GetRightMargin();
  // double newbottom = canvasratio + (1-canvasratio)*originalbottom-canvasratio*Canvas_CSC_hits_fancy_1->GetTopMargin();
  // double newright = canvasratio + (1-canvasratio)*originalright-canvasratio*Canvas_CSC_hits_fancy_1->GetLeftMargin();
  // right canvas
  Canvas_CSC_hits_fancy->SetRightMargin(canvasratio + (1-canvasratio)*Canvas_CSC_hits_fancy->GetRightMargin()-canvasratio*Canvas_CSC_hits_fancy->GetLeftMargin());
  TPad *Canvas_CSC_hits_fancy_2 = new TPad("RightPad","",0,0,1,1);
  Canvas_CSC_hits_fancy_2->SetBottomMargin(Canvas_CSC_hits_fancy->GetBottomMargin());
  Canvas_CSC_hits_fancy_2->SetLeftMargin((1-canvasratio) - (1-canvasratio)*Canvas_CSC_hits_fancy_2->GetRightMargin() + canvasratio*Canvas_CSC_hits_fancy_2->GetLeftMargin());
  Canvas_CSC_hits_fancy_1->SetFillStyle(4000);         Canvas_CSC_hits_fancy_2->SetFillStyle(4000);
  Canvas_CSC_hits_fancy_1->SetFillColor(4000);         Canvas_CSC_hits_fancy_2->SetFillColor(4000);
  Canvas_CSC_hits_fancy_1->SetFrameFillColor(4000);    Canvas_CSC_hits_fancy_2->SetFrameFillColor(4000);
  Canvas_CSC_hits_fancy_1->SetFrameFillStyle(4000);    Canvas_CSC_hits_fancy_2->SetFrameFillStyle(4000);
  Canvas_CSC_hits_fancy_1->SetFrameBorderMode(0);      Canvas_CSC_hits_fancy_2->SetFrameBorderMode(0);
  Canvas_CSC_hits_fancy_1->SetTicks(1,1);              Canvas_CSC_hits_fancy_2->SetTicks(1,1);
  Canvas_CSC_hits_fancy_1->Draw();                     Canvas_CSC_hits_fancy_2->Draw();
  // bottom canvas :: E_{kin}
  Canvas_CSC_hits_fancy_1->cd(); 
  // Canvas_CSC_hits_fancy_1->SetLogy();
  TH1F * CSC_kins = (TH1F*) CSC_el_kins->Clone(); CSC_kins->Add(CSC_mu_kins); CSC_kins->Add(CSC_ha_kins); // check whether everything is already saved ... 
  CSC_kins->Rebin(10);                                                                                   // check whether originals are already written to ROOT file ...
  CSC_kins->SetFillColor(kBlue); /*CSC_kins->SetFillStyle(3001);*/ CSC_kins->SetLineColor(kBlue); CSC_kins->SetLineStyle(1); CSC_kins->SetLineWidth(0); 
  CSC_kins->Draw();
  // CSC_kins->GetXaxis()->SetLabelSize(0.001);
  CSC_kins->GetYaxis()->SetLabelSize(0.75);    
  CSC_kins->GetYaxis()->SetTickLength(CSC_el_hits->GetYaxis()->GetTickLength());
  CSC_kins->GetYaxis()->SetNdivisions(-504);
  CSC_kins->GetYaxis()->SetRangeUser(0,2500);
  Canvas_CSC_hits_fancy_1->RedrawAxis();
  Canvas_CSC_hits_fancy_1->Update();
  // right canvas :: TOF
  Canvas_CSC_hits_fancy_2->cd();
  // Canvas_CSC_hits_fancy_2->SetLogy();
  // CSC_hits_tof->Rebin(10);
  CSC_hits_tof->SetFillColor(kBlue); /*CSC_hits_tof->SetFillStyle(3001);*/ CSC_hits_tof->SetLineColor(kBlue); CSC_hits_tof->SetLineStyle(1); CSC_hits_tof->SetLineWidth(1);
  CSC_hits_tof->Draw("hbar");
  CSC_hits_tof->GetXaxis()->SetLabelSize(0.001);
  CSC_hits_tof->GetYaxis()->SetLabelSize(0.75); // drawing a histogram with horizontal bar option (hbar) inverts x- and y-axis
  CSC_hits_tof->GetYaxis()->SetNdivisions(-504);
  CSC_hits_tof->GetYaxis()->SetRangeUser(0,2500);
  // Canvas_CSC_hits_fancy_2->RedrawAxis(); 
  Canvas_CSC_hits_fancy_2->Update();
  Canvas_CSC_hits_fancy->cd();
  Canvas_CSC_hits_fancy->SetTicks(1,1);
  Canvas_CSC_hits->RedrawAxis();
  Canvas_CSC_hits_fancy->Write();
  if(pdf_output) {Canvas_CSC_hits_fancy->Print(pdfFileName.c_str());}
  /*
  TCanvas * test = new TCanvas("test", "test", 600, 600);
  test->cd();
  test->SetTicks(1,1);
  CSC_hits_tof->Draw("hbar");
  test->Update();
  test->RedrawAxis();
  test->Write();
  if(pdf_output) {test->Print(pdfFileName.c_str());}
  */

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

  Canvas_RPCb_deposits = new TCanvas("Canvas_RPCb_deposits", "Simhit time vs E_{deposit} :: RPCb", 600, 600);
  Canvas_RPCf_deposits = new TCanvas("Canvas_RPCf_deposits", "Simhit time vs E_{deposit} :: RPCf", 600, 600);
  Canvas_CSC_deposits  = new TCanvas("Canvas_CSC_deposits",  "Simhit time vs E_{deposit} :: CSC",  600, 600);
  Canvas_DT_deposits   = new TCanvas("Canvas_DT_deposits",   "Simhit time vs E_{deposit} :: DT",   600, 600);
  Canvas_GEM_deposits  = new TCanvas("Canvas_GEM_deposits",  "Simhit time vs E_{deposit} :: GEM",  600, 600);
  Canvas_ME0_deposits  = new TCanvas("Canvas_ME0_deposits",  "Simhit time vs E_{deposit} :: ME0",  600, 600);

  // SimHit time vs Energy Deposit
  // Combine histograms in single plot
  Canvas_RPCb_deposits->cd();
  RPCb_el_deposits->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  RPCb_el_deposits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  RPCb_el_deposits->GetXaxis()->SetTitleOffset(1.2);
  RPCb_el_deposits->SetTitle("SimHit time vs E_{deposit} :: RPCb");
  RPCb_el_deposits->SetMarkerStyle(7);  RPCb_el_deposits->SetMarkerColor(kBlack);   RPCb_el_deposits->SetMarkerSize(1);  RPCb_el_deposits->Draw("P");
  RPCb_mu_deposits->SetMarkerStyle(24); RPCb_mu_deposits->SetMarkerColor(kBlue);    RPCb_mu_deposits->SetMarkerSize(1);  RPCb_mu_deposits->Draw("PSame");
  RPCb_pi_deposits->SetMarkerStyle(33); RPCb_pi_deposits->SetMarkerColor(kGreen);   RPCb_pi_deposits->SetMarkerSize(1);  RPCb_pi_deposits->Draw("PSame");
  RPCb_ka_deposits->SetMarkerStyle(5);  RPCb_ka_deposits->SetMarkerColor(kOrange);  RPCb_ka_deposits->SetMarkerSize(1);  RPCb_ka_deposits->Draw("PSame");
  RPCb_p_deposits->SetMarkerStyle(26);  RPCb_p_deposits->SetMarkerColor(kMagenta);  RPCb_p_deposits->SetMarkerSize(1);   RPCb_p_deposits->Draw("PSame");
  RPCb_n_deposits->SetMarkerStyle(32);  RPCb_n_deposits->SetMarkerColor(kViolet);   RPCb_n_deposits->SetMarkerSize(1);   RPCb_n_deposits->Draw("PSame");
  RPCb_g_deposits->SetMarkerStyle(30);  RPCb_g_deposits->SetMarkerColor(kCyan);   RPCb_g_deposits->SetMarkerSize(1);   RPCb_g_deposits->Draw("PSame");
  RPCb_N_deposits->SetMarkerStyle(2);   RPCb_N_deposits->SetMarkerColor(kRed);      RPCb_N_deposits->SetMarkerSize(1);   RPCb_N_deposits->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");       
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCb");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(0.40, 0.025,              "#font[12]{1 keV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_RPCb_deposits->SetTicks(1,1);
  Canvas_RPCb_deposits->Write();
  if(pdf_output) {Canvas_RPCb_deposits->Print(pdfFileName.c_str());}

  Canvas_RPCf_deposits->cd();
  RPCf_el_deposits->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  RPCf_el_deposits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  RPCf_el_deposits->GetXaxis()->SetTitleOffset(1.2);
  RPCf_el_deposits->SetTitle("SimHit time vs E_{deposit} :: RPCf");
  RPCf_el_deposits->SetMarkerStyle(7);  RPCf_el_deposits->SetMarkerColor(kBlack);   RPCf_el_deposits->SetMarkerSize(1);  RPCf_el_deposits->Draw("P");
  RPCf_mu_deposits->SetMarkerStyle(24); RPCf_mu_deposits->SetMarkerColor(kBlue);    RPCf_mu_deposits->SetMarkerSize(1);  RPCf_mu_deposits->Draw("PSame");
  RPCf_pi_deposits->SetMarkerStyle(33); RPCf_pi_deposits->SetMarkerColor(kGreen);   RPCf_pi_deposits->SetMarkerSize(1);  RPCf_pi_deposits->Draw("PSame");
  RPCf_ka_deposits->SetMarkerStyle(5);  RPCf_ka_deposits->SetMarkerColor(kOrange);  RPCf_ka_deposits->SetMarkerSize(1);  RPCf_ka_deposits->Draw("PSame");
  RPCf_p_deposits->SetMarkerStyle(26);  RPCf_p_deposits->SetMarkerColor(kMagenta);  RPCf_p_deposits->SetMarkerSize(1);   RPCf_p_deposits->Draw("PSame");
  RPCf_n_deposits->SetMarkerStyle(32);  RPCf_n_deposits->SetMarkerColor(kViolet);   RPCf_n_deposits->SetMarkerSize(1);   RPCf_n_deposits->Draw("PSame");
  RPCf_g_deposits->SetMarkerStyle(30);  RPCf_g_deposits->SetMarkerColor(kCyan);     RPCf_g_deposits->SetMarkerSize(1);   RPCf_g_deposits->Draw("PSame");
  RPCf_N_deposits->SetMarkerStyle(2);   RPCf_N_deposits->SetMarkerColor(kRed);      RPCf_N_deposits->SetMarkerSize(1);   RPCf_N_deposits->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCf");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(0.40, 0.025,              "#font[12]{1 keV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_RPCf_deposits->SetTicks(1,1);
  Canvas_RPCf_deposits->Write();
  if(pdf_output) {Canvas_RPCf_deposits->Print(pdfFileName.c_str());}

  Canvas_CSC_deposits->cd();
  CSC_el_deposits->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  CSC_el_deposits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  CSC_el_deposits->GetXaxis()->SetTitleOffset(1.2);
  CSC_el_deposits->SetTitle("SimHit time vs E_{deposit} :: CSC");
  CSC_el_deposits->SetMarkerStyle(7);  CSC_el_deposits->SetMarkerColor(kBlack);   CSC_el_deposits->SetMarkerSize(1);  CSC_el_deposits->Draw("P");
  CSC_mu_deposits->SetMarkerStyle(24); CSC_mu_deposits->SetMarkerColor(kBlue);    CSC_mu_deposits->SetMarkerSize(1);  CSC_mu_deposits->Draw("PSame");
  CSC_pi_deposits->SetMarkerStyle(33); CSC_pi_deposits->SetMarkerColor(kGreen);   CSC_pi_deposits->SetMarkerSize(1);  CSC_pi_deposits->Draw("PSame");
  CSC_ka_deposits->SetMarkerStyle(5);  CSC_ka_deposits->SetMarkerColor(kOrange);  CSC_ka_deposits->SetMarkerSize(1);  CSC_ka_deposits->Draw("PSame");
  CSC_p_deposits->SetMarkerStyle(26);  CSC_p_deposits->SetMarkerColor(kMagenta);  CSC_p_deposits->SetMarkerSize(1);   CSC_p_deposits->Draw("PSame");
  CSC_n_deposits->SetMarkerStyle(32);  CSC_n_deposits->SetMarkerColor(kViolet);   CSC_n_deposits->SetMarkerSize(1);   CSC_n_deposits->Draw("PSame");
  CSC_g_deposits->SetMarkerStyle(30);  CSC_g_deposits->SetMarkerColor(kCyan);     CSC_g_deposits->SetMarkerSize(1);   CSC_g_deposits->Draw("PSame");
  CSC_N_deposits->SetMarkerStyle(2);   CSC_N_deposits->SetMarkerColor(kRed);      CSC_N_deposits->SetMarkerSize(1);   CSC_N_deposits->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_deposits->SetTicks(1,1);
  Canvas_CSC_deposits->Write();
  if(pdf_output) {Canvas_CSC_deposits->Print(pdfFileName.c_str());}

  Canvas_DT_deposits->cd();
  DT_el_deposits->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  DT_el_deposits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  DT_el_deposits->GetXaxis()->SetTitleOffset(1.2);
  DT_el_deposits->SetTitle("SimHit time vs E_{deposit} :: DT");
  DT_el_deposits->SetMarkerStyle(7);  DT_el_deposits->SetMarkerColor(kBlack);   DT_el_deposits->SetMarkerSize(1);  DT_el_deposits->Draw("P");
  DT_mu_deposits->SetMarkerStyle(24); DT_mu_deposits->SetMarkerColor(kBlue);    DT_mu_deposits->SetMarkerSize(1);  DT_mu_deposits->Draw("PSame");
  DT_pi_deposits->SetMarkerStyle(33); DT_pi_deposits->SetMarkerColor(kGreen);   DT_pi_deposits->SetMarkerSize(1);  DT_pi_deposits->Draw("PSame");
  DT_ka_deposits->SetMarkerStyle(5);  DT_ka_deposits->SetMarkerColor(kOrange);  DT_ka_deposits->SetMarkerSize(1);  DT_ka_deposits->Draw("PSame");
  DT_p_deposits->SetMarkerStyle(26);  DT_p_deposits->SetMarkerColor(kMagenta);  DT_p_deposits->SetMarkerSize(1);   DT_p_deposits->Draw("PSame");
  DT_n_deposits->SetMarkerStyle(32);  DT_n_deposits->SetMarkerColor(kViolet);   DT_n_deposits->SetMarkerSize(1);   DT_n_deposits->Draw("PSame");
  DT_g_deposits->SetMarkerStyle(30);  DT_g_deposits->SetMarkerColor(kCyan);     DT_g_deposits->SetMarkerSize(1);   DT_g_deposits->Draw("PSame");
  DT_N_deposits->SetMarkerStyle(2);   DT_N_deposits->SetMarkerColor(kRed);      DT_N_deposits->SetMarkerSize(1);   DT_N_deposits->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"DT");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(0.40, 0.025,              "#font[12]{1 keV}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_DT_deposits->SetTicks(1,1);
  Canvas_DT_deposits->Write();
  if(pdf_output) {Canvas_DT_deposits->Print(pdfFileName.c_str());}

  Canvas_GEM_deposits->cd();
  GEM_el_deposits->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  GEM_el_deposits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  GEM_el_deposits->GetXaxis()->SetTitleOffset(1.2);
  GEM_el_deposits->SetTitle("SimHit time vs E_{deposit} :: GEM");
  GEM_el_deposits->SetMarkerStyle(7);  GEM_el_deposits->SetMarkerColor(kBlack);   GEM_el_deposits->SetMarkerSize(1);  GEM_el_deposits->Draw("P");
  GEM_mu_deposits->SetMarkerStyle(24); GEM_mu_deposits->SetMarkerColor(kBlue);    GEM_mu_deposits->SetMarkerSize(1);  GEM_mu_deposits->Draw("PSame");
  GEM_pi_deposits->SetMarkerStyle(33); GEM_pi_deposits->SetMarkerColor(kGreen);   GEM_pi_deposits->SetMarkerSize(1);  GEM_pi_deposits->Draw("PSame");
  GEM_ka_deposits->SetMarkerStyle(5);  GEM_ka_deposits->SetMarkerColor(kOrange);  GEM_ka_deposits->SetMarkerSize(1);  GEM_ka_deposits->Draw("PSame");
  GEM_p_deposits->SetMarkerStyle(26);  GEM_p_deposits->SetMarkerColor(kMagenta);  GEM_p_deposits->SetMarkerSize(1);   GEM_p_deposits->Draw("PSame");
  GEM_n_deposits->SetMarkerStyle(32);  GEM_n_deposits->SetMarkerColor(kViolet);   GEM_n_deposits->SetMarkerSize(1);   GEM_n_deposits->Draw("PSame");
  GEM_g_deposits->SetMarkerStyle(30);  GEM_g_deposits->SetMarkerColor(kCyan);     GEM_g_deposits->SetMarkerSize(1);   GEM_g_deposits->Draw("PSame");
  GEM_N_deposits->SetMarkerStyle(2);   GEM_N_deposits->SetMarkerColor(kRed);      GEM_N_deposits->SetMarkerSize(1);   GEM_N_deposits->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"GEM");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_GEM_deposits->SetTicks(1,1);
  Canvas_GEM_deposits->Write();
  if(pdf_output) {Canvas_GEM_deposits->Print(pdfFileName.c_str());}

  Canvas_ME0_deposits->cd();
  ME0_el_deposits->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  ME0_el_deposits->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  ME0_el_deposits->GetXaxis()->SetTitleOffset(1.2);
  ME0_el_deposits->SetTitle("SimHit time vs E_{deposit} :: ME0");
  ME0_el_deposits->SetMarkerStyle(7);  ME0_el_deposits->SetMarkerColor(kBlack);   ME0_el_deposits->SetMarkerSize(1);  ME0_el_deposits->Draw("P");
  ME0_mu_deposits->SetMarkerStyle(24); ME0_mu_deposits->SetMarkerColor(kBlue);    ME0_mu_deposits->SetMarkerSize(1);  ME0_mu_deposits->Draw("PSame");
  ME0_pi_deposits->SetMarkerStyle(33); ME0_pi_deposits->SetMarkerColor(kGreen);   ME0_pi_deposits->SetMarkerSize(1);  ME0_pi_deposits->Draw("PSame");
  ME0_ka_deposits->SetMarkerStyle(5);  ME0_ka_deposits->SetMarkerColor(kOrange);  ME0_ka_deposits->SetMarkerSize(1);  ME0_ka_deposits->Draw("PSame");
  ME0_p_deposits->SetMarkerStyle(26);  ME0_p_deposits->SetMarkerColor(kMagenta);  ME0_p_deposits->SetMarkerSize(1);   ME0_p_deposits->Draw("PSame");
  ME0_n_deposits->SetMarkerStyle(32);  ME0_n_deposits->SetMarkerColor(kViolet);   ME0_n_deposits->SetMarkerSize(1);   ME0_n_deposits->Draw("PSame");
  ME0_g_deposits->SetMarkerStyle(30);  ME0_g_deposits->SetMarkerColor(kCyan);     ME0_g_deposits->SetMarkerSize(1);   ME0_g_deposits->Draw("PSame");
  ME0_N_deposits->SetMarkerStyle(2);   ME0_N_deposits->SetMarkerColor(kRed);      ME0_N_deposits->SetMarkerSize(1);   ME0_N_deposits->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"ME0");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,"#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,"#font[12]{1 ms}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,"#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_ME0_deposits->SetTicks(1,1);
  Canvas_ME0_deposits->Write();
  if(pdf_output) {Canvas_ME0_deposits->Print(pdfFileName.c_str());}

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

  Canvas_RPCb_1D_deps = new TCanvas("Canvas_RPCb_1D_deps", "E_{deposit} :: RPCb", 600, 600);
  Canvas_RPCf_1D_deps = new TCanvas("Canvas_RPCf_1D_deps", "E_{deposit} :: RPCf", 600, 600);
  Canvas_CSC_1D_deps  = new TCanvas("Canvas_CSC_1D_deps",  "E_{deposit} :: CSC",  600, 600);
  Canvas_DT_1D_deps   = new TCanvas("Canvas_DT_1D_deps",   "E_{deposit} :: DT",   600, 600);
  Canvas_GEM_1D_deps  = new TCanvas("Canvas_GEM_1D_deps",  "E_{deposit} :: GEM",  600, 600);
  Canvas_ME0_1D_deps  = new TCanvas("Canvas_ME0_1D_deps",  "E_{deposit} :: ME0",  600, 600);

  Canvas_RPCb_1D_deps->cd();
  RPCb_el_deps->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  RPCb_el_deps->GetYaxis()->SetTitle("Hits");
  RPCb_el_deps->SetTitle("E_{deposit} :: RPCb");
  RPCb_el_deps->SetLineStyle(1);  RPCb_el_deps->SetLineColor(kBlack);   RPCb_el_deps->SetLineWidth(1);  RPCb_el_deps->Draw("H");
  RPCb_mu_deps->SetLineStyle(1); RPCb_mu_deps->SetLineColor(kBlue);    RPCb_mu_deps->SetLineWidth(1);  RPCb_mu_deps->Draw("HSame");
  RPCb_ha_deps->SetLineStyle(1); RPCb_ha_deps->SetLineColor(kGreen);   RPCb_ha_deps->SetLineWidth(1);  RPCb_ha_deps->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCb");
  latex_right.DrawLatex(0.40, 0.025,"#font[12]{1 keV}");
  Canvas_RPCb_1D_deps->SetTicks(1,1);
  Canvas_RPCb_1D_deps->Write();
  if(pdf_output) {Canvas_RPCb_1D_deps->Print(pdfFileName.c_str());}

  Canvas_RPCf_1D_deps->cd();
  RPCf_el_deps->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  RPCf_el_deps->GetYaxis()->SetTitle("Hits");
  RPCf_el_deps->SetTitle("E_{deposit} :: RPCf");
  RPCf_el_deps->SetLineStyle(1);  RPCf_el_deps->SetLineColor(kBlack);   RPCf_el_deps->SetLineWidth(1);  RPCf_el_deps->Draw("H");
  RPCf_mu_deps->SetLineStyle(1); RPCf_mu_deps->SetLineColor(kBlue);    RPCf_mu_deps->SetLineWidth(1);  RPCf_mu_deps->Draw("HSame");
  RPCf_ha_deps->SetLineStyle(1); RPCf_ha_deps->SetLineColor(kGreen);   RPCf_ha_deps->SetLineWidth(1);  RPCf_ha_deps->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCf");
  latex_right.DrawLatex(0.40, 0.025,"#font[12]{1 keV}");
  Canvas_RPCf_1D_deps->SetTicks(1,1);
  Canvas_RPCf_1D_deps->Write();
  if(pdf_output) {Canvas_RPCf_1D_deps->Print(pdfFileName.c_str());}

  Canvas_CSC_1D_deps->cd();
  CSC_el_deps->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  CSC_el_deps->GetYaxis()->SetTitle("Hits");
  CSC_el_deps->SetTitle("E_{deposit} :: CSC");
  CSC_el_deps->SetLineStyle(1);  CSC_el_deps->SetLineColor(kBlack);   CSC_el_deps->SetLineWidth(1);  CSC_el_deps->Draw("H");
  CSC_mu_deps->SetLineStyle(1); CSC_mu_deps->SetLineColor(kBlue);    CSC_mu_deps->SetLineWidth(1);  CSC_mu_deps->Draw("HSame");
  CSC_ha_deps->SetLineStyle(1); CSC_ha_deps->SetLineColor(kGreen);   CSC_ha_deps->SetLineWidth(1);  CSC_ha_deps->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_right.DrawLatex(0.40, 0.025,"#font[12]{1 keV}");
  Canvas_CSC_1D_deps->SetTicks(1,1);
  Canvas_CSC_1D_deps->Write();
  if(pdf_output) {Canvas_CSC_1D_deps->Print(pdfFileName.c_str());}

  Canvas_DT_1D_deps->cd();
  DT_el_deps->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  DT_el_deps->GetYaxis()->SetTitle("Hits");
  DT_el_deps->SetTitle("E_{deposit} :: DT");
  DT_el_deps->SetLineStyle(1);  DT_el_deps->SetLineColor(kBlack);   DT_el_deps->SetLineWidth(1);  DT_el_deps->Draw("H");
  DT_mu_deps->SetLineStyle(1); DT_mu_deps->SetLineColor(kBlue);    DT_mu_deps->SetLineWidth(1);  DT_mu_deps->Draw("HSame");
  DT_ha_deps->SetLineStyle(1); DT_ha_deps->SetLineColor(kGreen);   DT_ha_deps->SetLineWidth(1);  DT_ha_deps->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"DT");
  latex_right.DrawLatex(0.40, 0.025,"#font[12]{1 keV}");
  Canvas_DT_1D_deps->SetTicks(1,1);
  Canvas_DT_1D_deps->Write();
  if(pdf_output) {Canvas_DT_1D_deps->Print(pdfFileName.c_str());}

  Canvas_GEM_1D_deps->cd();
  GEM_el_deps->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  GEM_el_deps->GetYaxis()->SetTitle("Hits");
  GEM_el_deps->SetTitle("E_{deposit} :: GEM");
  GEM_el_deps->SetLineStyle(1);  GEM_el_deps->SetLineColor(kBlack);   GEM_el_deps->SetLineWidth(1);  GEM_el_deps->Draw("H");
  GEM_mu_deps->SetLineStyle(1); GEM_mu_deps->SetLineColor(kBlue);    GEM_mu_deps->SetLineWidth(1);  GEM_mu_deps->Draw("HSame");
  GEM_ha_deps->SetLineStyle(1); GEM_ha_deps->SetLineColor(kGreen);   GEM_ha_deps->SetLineWidth(1);  GEM_ha_deps->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"GEM");
  latex_right.DrawLatex(0.40, 0.025,"#font[12]{1 keV}");
  Canvas_GEM_1D_deps->SetTicks(1,1);
  Canvas_GEM_1D_deps->Write();
  if(pdf_output) {Canvas_GEM_1D_deps->Print(pdfFileName.c_str());}

  Canvas_ME0_1D_deps->cd();
  ME0_el_deps->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  ME0_el_deps->GetYaxis()->SetTitle("Hits");
  ME0_el_deps->SetTitle("E_{deposit} :: ME0");
  ME0_el_deps->SetLineStyle(1);  ME0_el_deps->SetLineColor(kBlack);   ME0_el_deps->SetLineWidth(1);  ME0_el_deps->Draw("H");
  ME0_mu_deps->SetLineStyle(1); ME0_mu_deps->SetLineColor(kBlue);    ME0_mu_deps->SetLineWidth(1);  ME0_mu_deps->Draw("HSame");
  ME0_ha_deps->SetLineStyle(1); ME0_ha_deps->SetLineColor(kGreen);   ME0_ha_deps->SetLineWidth(1);  ME0_ha_deps->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"ME0");
  latex_right.DrawLatex(0.40, 0.025,"#font[12]{1 keV}");
  Canvas_ME0_1D_deps->SetTicks(1,1);
  Canvas_ME0_1D_deps->Write();
  if(pdf_output) {Canvas_ME0_1D_deps->Print(pdfFileName.c_str());}

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

  Canvas_RPCb_1D_tof = new TCanvas("Canvas_RPCb_1D_tof", "Time Of Flight :: RPCb", 600, 600);
  Canvas_RPCf_1D_tof = new TCanvas("Canvas_RPCf_1D_tof", "Time Of Flight :: RPCf", 600, 600);
  Canvas_CSC_1D_tof  = new TCanvas("Canvas_CSC_1D_tof",  "Time Of Flight :: CSC",  600, 600);
  Canvas_DT_1D_tof   = new TCanvas("Canvas_DT_1D_tof",   "Time Of Flight :: DT",   600, 600);
  Canvas_GEM_1D_tof  = new TCanvas("Canvas_GEM_1D_tof",  "Time Of Flight :: GEM",  600, 600);
  Canvas_ME0_1D_tof  = new TCanvas("Canvas_ME0_1D_tof",  "Time Of Flight :: ME0",  600, 600);

  Canvas_RPCb_1D_tof->cd();
  RPCb_el_tof->GetXaxis()->SetTitle("^{10}log TOF (ns)");
  RPCb_el_tof->GetYaxis()->SetTitle("Hits");
  RPCb_el_tof->SetTitle("Time Of Flight :: RPCb");
  RPCb_el_tof->SetLineStyle(1);  RPCb_el_tof->SetLineColor(kBlack);   RPCb_el_tof->SetLineWidth(1);  RPCb_el_tof->Draw("H");
  RPCb_mu_tof->SetLineStyle(1); RPCb_mu_tof->SetLineColor(kBlue);    RPCb_mu_tof->SetLineWidth(1);  RPCb_mu_tof->Draw("HSame");
  RPCb_ha_tof->SetLineStyle(1); RPCb_ha_tof->SetLineColor(kGreen);   RPCb_ha_tof->SetLineWidth(1);  RPCb_ha_tof->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCb");
  Canvas_RPCb_1D_tof->SetTicks(1,1);
  Canvas_RPCb_1D_tof->Write();
  if(pdf_output) {Canvas_RPCb_1D_tof->Print(pdfFileName.c_str());}

  Canvas_RPCf_1D_tof->cd();
  RPCf_el_tof->GetXaxis()->SetTitle("^{10}log TOF (ns)");
  RPCf_el_tof->GetYaxis()->SetTitle("Hits");
  RPCf_el_tof->SetTitle("Time Of Flight :: RPCf");
  RPCf_el_tof->SetLineStyle(1);  RPCf_el_tof->SetLineColor(kBlack);   RPCf_el_tof->SetLineWidth(1);  RPCf_el_tof->Draw("H");
  RPCf_mu_tof->SetLineStyle(1); RPCf_mu_tof->SetLineColor(kBlue);    RPCf_mu_tof->SetLineWidth(1);  RPCf_mu_tof->Draw("HSame");
  RPCf_ha_tof->SetLineStyle(1); RPCf_ha_tof->SetLineColor(kGreen);   RPCf_ha_tof->SetLineWidth(1);  RPCf_ha_tof->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCf");
  Canvas_RPCf_1D_tof->SetTicks(1,1);
  Canvas_RPCf_1D_tof->Write();
  if(pdf_output) {Canvas_RPCf_1D_tof->Print(pdfFileName.c_str());}

  Canvas_CSC_1D_tof->cd();
  CSC_el_tof->GetXaxis()->SetTitle("^{10}log TOF (ns)");
  CSC_el_tof->GetYaxis()->SetTitle("Hits");
  CSC_el_tof->SetTitle("Time Of Flight :: CSC");
  CSC_el_tof->SetLineStyle(1);  CSC_el_tof->SetLineColor(kBlack);   CSC_el_tof->SetLineWidth(1);  CSC_el_tof->Draw("H");
  CSC_mu_tof->SetLineStyle(1); CSC_mu_tof->SetLineColor(kBlue);    CSC_mu_tof->SetLineWidth(1);  CSC_mu_tof->Draw("HSame");
  CSC_ha_tof->SetLineStyle(1); CSC_ha_tof->SetLineColor(kGreen);   CSC_ha_tof->SetLineWidth(1);  CSC_ha_tof->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  Canvas_CSC_1D_tof->SetTicks(1,1);
  Canvas_CSC_1D_tof->Write();
  if(pdf_output) {Canvas_CSC_1D_tof->Print(pdfFileName.c_str());}

  Canvas_DT_1D_tof->cd();
  DT_el_tof->GetXaxis()->SetTitle("^{10}log TOF (ns)");
  DT_el_tof->GetYaxis()->SetTitle("Hits");
  DT_el_tof->SetTitle("Time Of Flight :: DT");
  DT_el_tof->SetLineStyle(1);  DT_el_tof->SetLineColor(kBlack);   DT_el_tof->SetLineWidth(1);  DT_el_tof->Draw("H");
  DT_mu_tof->SetLineStyle(1); DT_mu_tof->SetLineColor(kBlue);    DT_mu_tof->SetLineWidth(1);  DT_mu_tof->Draw("HSame");
  DT_ha_tof->SetLineStyle(1); DT_ha_tof->SetLineColor(kGreen);   DT_ha_tof->SetLineWidth(1);  DT_ha_tof->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"DT");
  Canvas_DT_1D_tof->SetTicks(1,1);
  Canvas_DT_1D_tof->Write();
  if(pdf_output) {Canvas_DT_1D_tof->Print(pdfFileName.c_str());}

  Canvas_GEM_1D_tof->cd();
  GEM_el_tof->GetXaxis()->SetTitle("^{10}log TOF (ns)");
  GEM_el_tof->GetYaxis()->SetTitle("Hits");
  GEM_el_tof->SetTitle("Time Of Flight :: GEM");
  GEM_el_tof->SetLineStyle(1);  GEM_el_tof->SetLineColor(kBlack);   GEM_el_tof->SetLineWidth(1);  GEM_el_tof->Draw("H");
  GEM_mu_tof->SetLineStyle(1); GEM_mu_tof->SetLineColor(kBlue);    GEM_mu_tof->SetLineWidth(1);  GEM_mu_tof->Draw("HSame");
  GEM_ha_tof->SetLineStyle(1); GEM_ha_tof->SetLineColor(kGreen);   GEM_ha_tof->SetLineWidth(1);  GEM_ha_tof->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"GEM");
  Canvas_GEM_1D_tof->SetTicks(1,1);
  Canvas_GEM_1D_tof->Write();
  if(pdf_output) {Canvas_GEM_1D_tof->Print(pdfFileName.c_str());}

  Canvas_ME0_1D_tof->cd();
  ME0_el_tof->GetXaxis()->SetTitle("^{10}log TOF (ns)");
  ME0_el_tof->GetYaxis()->SetTitle("Hits");
  ME0_el_tof->SetTitle("Time Of Flight :: ME0");
  ME0_el_tof->SetLineStyle(1);  ME0_el_tof->SetLineColor(kBlack);   ME0_el_tof->SetLineWidth(1);  ME0_el_tof->Draw("H");
  ME0_mu_tof->SetLineStyle(1); ME0_mu_tof->SetLineColor(kBlue);    ME0_mu_tof->SetLineWidth(1);  ME0_mu_tof->Draw("HSame");
  ME0_ha_tof->SetLineStyle(1); ME0_ha_tof->SetLineColor(kGreen);   ME0_ha_tof->SetLineWidth(1);  ME0_ha_tof->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"ME0");
  Canvas_ME0_1D_tof->SetTicks(1,1);
  Canvas_ME0_1D_tof->Write();
  if(pdf_output) {Canvas_ME0_1D_tof->Print(pdfFileName.c_str());}

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

  DT_MB1_all_hits_tof->Write();
  DT_MB2_all_hits_tof->Write();
  DT_MB3_all_hits_tof->Write();
  DT_MB4_all_hits_tof->Write();
  DT_MB1_W2_hits_tof->Write();
  DT_MB1_W1_hits_tof->Write();
  DT_MB1_W0_hits_tof->Write();
  DT_MBX_all_hits_tof->Write();
  CSC_ME11_hits_tof->Write();
  CSC_ME12_hits_tof->Write();
  CSC_ME13_hits_tof->Write();
  CSC_ME21_hits_tof->Write();
  CSC_ME22_hits_tof->Write();
  CSC_ME31_hits_tof->Write();
  CSC_ME32_hits_tof->Write();
  CSC_ME41_hits_tof->Write();
  CSC_ME42_hits_tof->Write();
  CSC_MEX1_hits_tof->Write();
  CSC_MEX2_hits_tof->Write();  
  CSC_MEXX_hits_tof->Write();




  // --------------------------- 
  outputfile->cd();

  Canvas_RPCb_hits_phi = new TCanvas("Canvas_RPCb_hits_phi", "Hits vs phi :: RPCb", 600, 600);
  Canvas_RPCf_hits_phi = new TCanvas("Canvas_RPCf_hits_phi", "Hits vs phi :: RPCf", 600, 600);
  Canvas_CSC_hits_phi  = new TCanvas("Canvas_CSC_hits_phi",  "Hits vs phi :: CSC", 600, 600);
  Canvas_DT_hits_phi   = new TCanvas("Canvas_DT_hits_phi",   "Hits vs phi :: DT", 600, 600);

  Canvas_RPCb_hits_phi->cd();
  RB4_hits_phi->GetXaxis()->SetTitle("phi (-)");
  RB4_hits_phi->GetYaxis()->SetTitle("Hits");
  RB4_hits_phi->SetTitle("Hits vs phi :: RPCb");
  RB4_hits_phi->SetLineStyle(1); RB4_hits_phi->SetLineColor(kRed);     RB4_hits_phi->SetLineWidth(1);  RB4_hits_phi->SetMarkerStyle(20); RB4_hits_phi->SetMarkerColor(kRed);     RB4_hits_phi->Draw("HP");
  RB3_hits_phi->SetLineStyle(1); RB3_hits_phi->SetLineColor(kMagenta); RB3_hits_phi->SetLineWidth(1);  RB3_hits_phi->SetMarkerStyle(24); RB3_hits_phi->SetMarkerColor(kMagenta); RB3_hits_phi->Draw("HPSame");
  RB2_hits_phi->SetLineStyle(1); RB2_hits_phi->SetLineColor(kBlue);    RB2_hits_phi->SetLineWidth(1);  RB2_hits_phi->SetMarkerStyle(27); RB2_hits_phi->SetMarkerColor(kBlue);    RB2_hits_phi->Draw("HPSame");
  RB1_hits_phi->SetLineStyle(1); RB1_hits_phi->SetLineColor(kOrange);  RB1_hits_phi->SetLineWidth(1);  RB1_hits_phi->SetMarkerStyle(29); RB1_hits_phi->SetMarkerColor(kOrange);  RB1_hits_phi->Draw("HPSame");
  l1_x1 = 0.725; l1_y1 = 0.725; l1_x2 = 0.875; l1_y2 = 0.875;
  TLegend *lRPCb = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  lRPCb->SetLineColor(0);    lRPCb->SetLineStyle(0);  lRPCb->SetLineWidth(0);
  lRPCb->SetFillColor(4000); lRPCb->SetBorderSize(0); lRPCb->SetNColumns(2);
  lRPCb->AddEntry(RB1_hits_phi, "RB1","lp");
  lRPCb->AddEntry(RB2_hits_phi, "RB2","lp");
  lRPCb->AddEntry(RB3_hits_phi, "RB3","lp");
  lRPCb->AddEntry(RB4_hits_phi, "RB4","lp");
  lRPCb->Draw();
  RB4_hits_phi->GetYaxis()->SetRangeUser(0,1.25*RB4_hits_phi->GetBinContent(RB4_hits_phi->GetMaximumBin()));
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCb");
  Canvas_RPCb_hits_phi->SetTicks(1,1);
  Canvas_RPCb_hits_phi->Write();
  if(pdf_output) {Canvas_RPCb_hits_phi->Print(pdfFileName.c_str());}

  Canvas_DT_hits_phi->cd();
  MB4_hits_phi->GetXaxis()->SetTitle("phi (-)");
  MB4_hits_phi->GetYaxis()->SetTitle("Hits");
  MB4_hits_phi->SetTitle("Hits vs phi :: DT");
  MB4_hits_phi->SetLineStyle(1); MB4_hits_phi->SetLineColor(kRed);     MB4_hits_phi->SetLineWidth(1);  MB4_hits_phi->SetMarkerStyle(20); MB4_hits_phi->SetMarkerColor(kRed);     MB4_hits_phi->Draw("HP");
  MB3_hits_phi->SetLineStyle(1); MB3_hits_phi->SetLineColor(kMagenta); MB3_hits_phi->SetLineWidth(1);  MB3_hits_phi->SetMarkerStyle(24); MB3_hits_phi->SetMarkerColor(kMagenta); MB3_hits_phi->Draw("HPSame");
  MB2_hits_phi->SetLineStyle(1); MB2_hits_phi->SetLineColor(kBlue);    MB2_hits_phi->SetLineWidth(1);  MB2_hits_phi->SetMarkerStyle(27); MB2_hits_phi->SetMarkerColor(kBlue);    MB2_hits_phi->Draw("HPSame");
  MB1_hits_phi->SetLineStyle(1); MB1_hits_phi->SetLineColor(kOrange);  MB1_hits_phi->SetLineWidth(1);  MB1_hits_phi->SetMarkerStyle(29); MB1_hits_phi->SetMarkerColor(kOrange);  MB1_hits_phi->Draw("HPSame");
  l1_x1 = 0.725; l1_y1 = 0.725; l1_x2 = 0.875; l1_y2 = 0.875;
  TLegend *lDT = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  lDT->SetLineColor(0);    lDT->SetLineStyle(0);  lDT->SetLineWidth(0);
  lDT->SetFillColor(4000); lDT->SetBorderSize(0); lDT->SetNColumns(2);
  lDT->AddEntry(MB1_hits_phi, "MB1","lp");
  lDT->AddEntry(MB2_hits_phi, "MB2","lp");
  lDT->AddEntry(MB3_hits_phi, "MB3","lp");
  lDT->AddEntry(MB4_hits_phi, "MB4","lp");
  lDT->Draw();
  MB4_hits_phi->GetYaxis()->SetRangeUser(0,1.25*MB4_hits_phi->GetBinContent(MB4_hits_phi->GetMaximumBin()));
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"DT");
  Canvas_DT_hits_phi->SetTicks(1,1);
  Canvas_DT_hits_phi->Write();
  if(pdf_output) {Canvas_DT_hits_phi->Print(pdfFileName.c_str());}

  Canvas_RPCf_hits_phi->cd();
  RE12_hits_phi->GetXaxis()->SetTitle("phi (-)");
  RE12_hits_phi->GetYaxis()->SetTitle("Hits");
  RE12_hits_phi->SetTitle("Hits vs phi :: RPCf");
  RE12_hits_phi->SetLineStyle(1); RE12_hits_phi->SetLineColor(kRed);     RE12_hits_phi->SetLineWidth(1);  RE12_hits_phi->SetMarkerStyle(26); RE12_hits_phi->SetMarkerSize(1);    RE12_hits_phi->SetMarkerColor(kRed);     RE12_hits_phi->Draw("HP");
  RE13_hits_phi->SetLineStyle(2); RE13_hits_phi->SetLineColor(kRed);     RE13_hits_phi->SetLineWidth(1);  RE13_hits_phi->SetMarkerStyle(20); RE13_hits_phi->SetMarkerSize(0.75); RE13_hits_phi->SetMarkerColor(kRed);     RE13_hits_phi->Draw("HPSame");
  RE22_hits_phi->SetLineStyle(1); RE22_hits_phi->SetLineColor(kMagenta); RE22_hits_phi->SetLineWidth(1);  RE22_hits_phi->SetMarkerStyle(27); RE22_hits_phi->SetMarkerSize(1);    RE22_hits_phi->SetMarkerColor(kMagenta); RE22_hits_phi->Draw("HPSame");
  RE23_hits_phi->SetLineStyle(2); RE23_hits_phi->SetLineColor(kMagenta); RE23_hits_phi->SetLineWidth(1);  RE23_hits_phi->SetMarkerStyle(25); RE23_hits_phi->SetMarkerSize(1);    RE23_hits_phi->SetMarkerColor(kMagenta); RE23_hits_phi->Draw("HPSame");
  RE32_hits_phi->SetLineStyle(1); RE32_hits_phi->SetLineColor(kBlue);    RE32_hits_phi->SetLineWidth(1);  RE32_hits_phi->SetMarkerStyle(24); RE32_hits_phi->SetMarkerSize(1);    RE32_hits_phi->SetMarkerColor(kBlue);    RE32_hits_phi->Draw("HPSame");
  RE33_hits_phi->SetLineStyle(2); RE33_hits_phi->SetLineColor(kBlue);    RE33_hits_phi->SetLineWidth(1);  RE33_hits_phi->SetMarkerStyle(23); RE33_hits_phi->SetMarkerSize(1);    RE33_hits_phi->SetMarkerColor(kBlue);    RE33_hits_phi->Draw("HPSame");
  RE42_hits_phi->SetLineStyle(1); RE42_hits_phi->SetLineColor(kOrange);  RE42_hits_phi->SetLineWidth(1);  RE42_hits_phi->SetMarkerStyle(20); RE42_hits_phi->SetMarkerSize(1.25); RE42_hits_phi->SetMarkerColor(kOrange);  RE42_hits_phi->Draw("HPSame");
  RE43_hits_phi->SetLineStyle(2); RE43_hits_phi->SetLineColor(kOrange);  RE43_hits_phi->SetLineWidth(1);  RE43_hits_phi->SetMarkerStyle(28); RE43_hits_phi->SetMarkerSize(1);    RE43_hits_phi->SetMarkerColor(kOrange);  RE43_hits_phi->Draw("HPSame");
  l1_x1 = 0.575; l1_y1 = 0.575; l1_x2 = 0.875; l1_y2 = 0.875;
  TLegend *lRPCf = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  lRPCf->SetLineColor(0);    lRPCf->SetLineStyle(0);  lRPCf->SetLineWidth(0);
  lRPCf->SetFillColor(4000); lRPCf->SetBorderSize(0); lRPCf->SetNColumns(2);
  lRPCf->AddEntry(RE12_hits_phi, "RE12","lp");
  lRPCf->AddEntry(RE13_hits_phi, "RE13","lp");
  lRPCf->AddEntry(RE22_hits_phi, "RE22","lp");
  lRPCf->AddEntry(RE23_hits_phi, "RE23","lp");
  lRPCf->AddEntry(RE32_hits_phi, "RE32","lp");
  lRPCf->AddEntry(RE33_hits_phi, "RE33","lp");
  lRPCf->AddEntry(RE42_hits_phi, "RE42","lp");
  lRPCf->AddEntry(RE43_hits_phi, "RE43","lp");
  lRPCf->Draw();
  RE12_hits_phi->GetYaxis()->SetRangeUser(0,1.25*RE43_hits_phi->GetBinContent(RE43_hits_phi->GetMaximumBin()));
  std::cout<<"RE43 Max Bin = "<<RE43_hits_phi->GetMaximumBin()<<" has bincontent = "<<RE43_hits_phi->GetBinContent(RE43_hits_phi->GetMaximumBin())<<" ==> max height = "<<1.25*RE43_hits_phi->GetBinContent(RE43_hits_phi->GetMaximumBin())<<std::endl;
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"RPCf");
  Canvas_RPCf_hits_phi->SetTicks(1,1);
  Canvas_RPCf_hits_phi->Write();
  if(pdf_output) {Canvas_RPCf_hits_phi->Print(pdfFileName.c_str());}

  Canvas_CSC_hits_phi->cd();
  ME11_hits_phi->GetXaxis()->SetTitle("phi (-)");
  ME11_hits_phi->GetYaxis()->SetTitle("Hits");
  ME11_hits_phi->SetTitle("Hits vs phi :: CSC");
  ME11_hits_phi->SetLineStyle(1); ME11_hits_phi->SetLineColor(kRed);     ME11_hits_phi->SetLineWidth(1);  ME11_hits_phi->SetMarkerStyle(29); ME11_hits_phi->SetMarkerSize(1);    ME11_hits_phi->SetMarkerColor(kRed);     ME11_hits_phi->Draw("HP");
  ME12_hits_phi->SetLineStyle(2); ME12_hits_phi->SetLineColor(kRed);     ME12_hits_phi->SetLineWidth(1);  ME12_hits_phi->SetMarkerStyle(26); ME12_hits_phi->SetMarkerSize(1);    ME12_hits_phi->SetMarkerColor(kRed);     ME12_hits_phi->Draw("HPSame");
  ME13_hits_phi->SetLineStyle(3); ME13_hits_phi->SetLineColor(kRed);     ME13_hits_phi->SetLineWidth(1);  ME13_hits_phi->SetMarkerStyle(20); ME13_hits_phi->SetMarkerSize(0.75); ME13_hits_phi->SetMarkerColor(kRed);     ME13_hits_phi->Draw("HPSame");
  ME21_hits_phi->SetLineStyle(1); ME21_hits_phi->SetLineColor(kMagenta); ME21_hits_phi->SetLineWidth(1);  ME21_hits_phi->SetMarkerStyle(27); ME21_hits_phi->SetMarkerSize(1);    ME21_hits_phi->SetMarkerColor(kMagenta); ME21_hits_phi->Draw("HPSame");
  ME22_hits_phi->SetLineStyle(2); ME22_hits_phi->SetLineColor(kMagenta); ME22_hits_phi->SetLineWidth(1);  ME22_hits_phi->SetMarkerStyle(25); ME22_hits_phi->SetMarkerSize(1);    ME22_hits_phi->SetMarkerColor(kMagenta); ME22_hits_phi->Draw("HPSame");
  ME31_hits_phi->SetLineStyle(1); ME31_hits_phi->SetLineColor(kBlue);    ME31_hits_phi->SetLineWidth(1);  ME31_hits_phi->SetMarkerStyle(24); ME31_hits_phi->SetMarkerSize(1);    ME31_hits_phi->SetMarkerColor(kBlue);    ME31_hits_phi->Draw("HPSame");
  ME32_hits_phi->SetLineStyle(2); ME32_hits_phi->SetLineColor(kBlue);    ME32_hits_phi->SetLineWidth(1);  ME32_hits_phi->SetMarkerStyle(23); ME32_hits_phi->SetMarkerSize(1);    ME32_hits_phi->SetMarkerColor(kBlue);    ME32_hits_phi->Draw("HPSame");
  ME41_hits_phi->SetLineStyle(1); ME41_hits_phi->SetLineColor(kOrange);  ME41_hits_phi->SetLineWidth(1);  ME41_hits_phi->SetMarkerStyle(20); ME41_hits_phi->SetMarkerSize(1.25); ME41_hits_phi->SetMarkerColor(kOrange);  ME41_hits_phi->Draw("HPSame");
  ME42_hits_phi->SetLineStyle(2); ME42_hits_phi->SetLineColor(kOrange);  ME42_hits_phi->SetLineWidth(1);  ME42_hits_phi->SetMarkerStyle(28); ME42_hits_phi->SetMarkerSize(1);    ME42_hits_phi->SetMarkerColor(kOrange);  ME42_hits_phi->Draw("HPSame");
  l1_x1 = 0.575; l1_y1 = 0.575; l1_x2 = 0.875; l1_y2 = 0.875;
  TLegend *lCSC = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  lCSC->SetLineColor(0);    lCSC->SetLineStyle(0);  lCSC->SetLineWidth(0);
  lCSC->SetFillColor(4000); lCSC->SetBorderSize(0); lCSC->SetNColumns(3);
  lCSC->AddEntry(ME11_hits_phi, "ME11","lp");
  lCSC->AddEntry(ME12_hits_phi, "ME12","lp");
  lCSC->AddEntry(ME13_hits_phi, "ME13","lp");
  lCSC->AddEntry(ME21_hits_phi, "ME21","lp");
  lCSC->AddEntry(ME22_hits_phi, "ME22","lp");
  lCSC->AddEntry(ME31_hits_phi, "ME31","lp");
  lCSC->AddEntry(ME32_hits_phi, "ME32","lp");
  lCSC->AddEntry(ME41_hits_phi, "ME41","lp");
  lCSC->AddEntry(ME42_hits_phi, "ME42","lp");
  lCSC->Draw();
  ME11_hits_phi->GetYaxis()->SetRangeUser(0,1.25*std::max(ME11_hits_phi->GetBinContent(ME11_hits_phi->GetMaximumBin()),ME31_hits_phi->GetBinContent(ME31_hits_phi->GetMaximumBin())));
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  Canvas_CSC_hits_phi->SetTicks(1,1);
  Canvas_CSC_hits_phi->Write();
  if(pdf_output) {Canvas_CSC_hits_phi->Print(pdfFileName.c_str());}




  // -----    Muon XZ plots    ----- //
  Canvas_Muon_RZ = new TCanvas("Canvas_Muon_RZ", "RZ-view of SimHits :: Muon", 600, 600);
  Muon_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_RZ->SetTitle("RZ-view of SimHits :: Muon");
  Muon_RZ->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"All SimHits");
  Canvas_Muon_RZ->SetTicks(1,1);
  Canvas_Muon_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_RZ->Print(pdfFileName.c_str());}

  Canvas_Muon_000ns_RZ = new TCanvas("Canvas_Muon_000ns_RZ", "RZ-view of SimHits with tof < 250 ns :: Muon", 600, 600);
  Muon_000ns_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_000ns_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_000ns_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_000ns_RZ->SetTitle("RZ-view of SimHits with tof < 250 ns :: Muon");
  Muon_000ns_RZ->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 250 ns");
  Canvas_Muon_000ns_RZ->SetTicks(1,1);
  Canvas_Muon_000ns_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_000ns_RZ->Print(pdfFileName.c_str());}

  Canvas_Muon_250ns_RZ = new TCanvas("Canvas_Muon_250ns_RZ", "RZ-view of SimHits with tof > 250 ns :: Muon", 600, 600);
  Muon_250ns_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_250ns_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_250ns_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_250ns_RZ->SetTitle("RZ-view of SimHits with tof > 250 ns :: Muon");
  Muon_250ns_RZ->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF > 250 ns");
  Canvas_Muon_250ns_RZ->SetTicks(1,1);
  Canvas_Muon_250ns_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_250ns_RZ->Print(pdfFileName.c_str());}

  Muon_000ns_el_RZ->SetMarkerColor(kBlack);  Muon_000ns_mu_RZ->SetMarkerColor(kBlue);  Muon_000ns_ha_RZ->SetMarkerColor(kRed);
  l2_x1 = 0.15; l2_x2 = 0.35; l2_y1 = 0.20; l2_y2 = 0.40;
  TLegend *l2a = new TLegend(l2_x1, l2_y1,l2_x2,l2_y2,NULL,"brNDC");
  l2a->SetLineColor(0);    l2a->SetLineStyle(0);  l2a->SetLineWidth(0);
  l2a->SetFillColor(4000); l2a->SetBorderSize(0); l2a->SetNColumns(1);
  l2a->AddEntry(Muon_000ns_el_RZ, "e","p");
  l2a->AddEntry(Muon_000ns_mu_RZ, "#mu","p");
  l2a->AddEntry(Muon_000ns_ha_RZ, "had","p");

  Canvas_Muon_000ns_Cont_RZ = new TCanvas("Canvas_Muon_000ns_Cont_RZ", "RZ-view of SimHits with tof < 250 ns :: Muon", 600, 600);
  Muon_000ns_el_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_000ns_el_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_000ns_el_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_000ns_el_RZ->SetTitle("RZ-view of SimHits with tof < 250 ns :: Muon");
  Muon_000ns_el_RZ->SetMarkerStyle(7); Muon_000ns_el_RZ->SetMarkerColor(kBlack);   Muon_000ns_el_RZ->SetMarkerSize(1);  Muon_000ns_el_RZ->Draw("P");
  Muon_000ns_mu_RZ->SetMarkerStyle(7); Muon_000ns_mu_RZ->SetMarkerColor(kBlue);    Muon_000ns_mu_RZ->SetMarkerSize(1);  Muon_000ns_mu_RZ->Draw("PSame");
  Muon_000ns_ha_RZ->SetMarkerStyle(7); Muon_000ns_ha_RZ->SetMarkerColor(kRed);     Muon_000ns_ha_RZ->SetMarkerSize(1);  Muon_000ns_ha_RZ->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 250 ns");
  l2a->Draw();
  Canvas_Muon_000ns_Cont_RZ->SetTicks(1,1);
  Canvas_Muon_000ns_Cont_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_000ns_Cont_RZ->Print(pdfFileName.c_str());}

  Canvas_Muon_250ns_Cont_RZ = new TCanvas("Canvas_Muon_250ns_Cont_RZ", "RZ-view of SimHits with tof > 250 ns :: Muon", 600, 600);
  Muon_250ns_el_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_250ns_el_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_250ns_el_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_250ns_el_RZ->SetTitle("RZ-view of SimHits with tof > 250 ns :: Muon");
  Muon_250ns_el_RZ->SetMarkerStyle(7); Muon_250ns_el_RZ->SetMarkerColor(kBlack);   Muon_250ns_el_RZ->SetMarkerSize(1);  Muon_250ns_el_RZ->Draw("P");
  Muon_250ns_mu_RZ->SetMarkerStyle(7); Muon_250ns_mu_RZ->SetMarkerColor(kBlue);    Muon_250ns_mu_RZ->SetMarkerSize(1);  Muon_250ns_mu_RZ->Draw("PSame");
  Muon_250ns_ha_RZ->SetMarkerStyle(7); Muon_250ns_ha_RZ->SetMarkerColor(kRed);     Muon_250ns_ha_RZ->SetMarkerSize(1);  Muon_250ns_ha_RZ->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF > 250 ns");
  l2a->Draw();
  Canvas_Muon_250ns_Cont_RZ->SetTicks(1,1);
  Canvas_Muon_250ns_Cont_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_250ns_Cont_RZ->Print(pdfFileName.c_str());}

  Canvas_Muon_00ns_Cont_RZ = new TCanvas("Canvas_Muon_00ns_Cont_RZ", "RZ-view of SimHits with tof < 50 ns :: Muon", 600, 600);
  Muon_00ns_el_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_00ns_el_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_00ns_el_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_00ns_el_RZ->SetTitle("RZ-view of SimHits with tof < 50 ns :: Muon");
  Muon_00ns_el_RZ->SetMarkerStyle(7); Muon_00ns_el_RZ->SetMarkerColor(kBlack);   Muon_00ns_el_RZ->SetMarkerSize(1);  Muon_00ns_el_RZ->Draw("P");
  Muon_00ns_mu_RZ->SetMarkerStyle(7); Muon_00ns_mu_RZ->SetMarkerColor(kBlue);    Muon_00ns_mu_RZ->SetMarkerSize(1);  Muon_00ns_mu_RZ->Draw("PSame");
  Muon_00ns_ha_RZ->SetMarkerStyle(7); Muon_00ns_ha_RZ->SetMarkerColor(kRed);     Muon_00ns_ha_RZ->SetMarkerSize(1);  Muon_00ns_ha_RZ->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 50 ns");
  l2a->Draw();
  Canvas_Muon_00ns_Cont_RZ->SetTicks(1,1);
  Canvas_Muon_00ns_Cont_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_00ns_Cont_RZ->Print(pdfFileName.c_str());}

  Canvas_Muon_50ns_Cont_RZ = new TCanvas("Canvas_Muon_50ns_Cont_RZ", "RZ-view of SimHits with 50 < tof < 250 ns :: Muon", 600, 600);
  Muon_50ns_el_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_50ns_el_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_50ns_el_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_50ns_el_RZ->SetTitle("RZ-view of SimHits with 50 < tof < 250 ns :: Muon");
  Muon_50ns_el_RZ->SetMarkerStyle(7); Muon_50ns_el_RZ->SetMarkerColor(kBlack);   Muon_50ns_el_RZ->SetMarkerSize(1);  Muon_50ns_el_RZ->Draw("P");
  Muon_50ns_mu_RZ->SetMarkerStyle(7); Muon_50ns_mu_RZ->SetMarkerColor(kBlue);    Muon_50ns_mu_RZ->SetMarkerSize(1);  Muon_50ns_mu_RZ->Draw("PSame");
  Muon_50ns_ha_RZ->SetMarkerStyle(7); Muon_50ns_ha_RZ->SetMarkerColor(kRed);     Muon_50ns_ha_RZ->SetMarkerSize(1);  Muon_50ns_ha_RZ->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits 50 < TOF < 250 ns");
  l2a->Draw();
  Canvas_Muon_50ns_Cont_RZ->SetTicks(1,1);
  Canvas_Muon_50ns_Cont_RZ->Write(); 
  if(pdf_output) {Canvas_Muon_50ns_Cont_RZ->Print(pdfFileName.c_str());}
  // -------------------------------- //

  // ----- Muon Barrel XY plots ----- //
  Canvas_Muon_Barrel_XY = new TCanvas("Canvas_Muon_Barrel_XY", "XY-view of SimHits :: Muon Barrel", 600, 600);
  Muon_Barrel_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_XY->GetYaxis()->SetTitle("Y (cm)");
  Muon_Barrel_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_XY->SetTitle("XY-view of SimHits :: Muon Barrel");
  Muon_Barrel_XY->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"All SimHits");
  Canvas_Muon_Barrel_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Barrel_000ns_XY = new TCanvas("Canvas_Muon_Barrel_000ns_XY", "XY-view of SimHits with tof < 250 ns :: Muon Barrel", 600, 600);
  Muon_Barrel_000ns_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_000ns_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Barrel_000ns_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_000ns_XY->SetTitle("XY-view of SimHits with tof < 250 ns :: Muon Barrel");
  Muon_Barrel_000ns_XY->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 250 ns");
  Canvas_Muon_Barrel_000ns_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_000ns_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_000ns_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Barrel_250ns_XY = new TCanvas("Canvas_Muon_Barrel_250ns_XY", "XY-view of SimHits with tof > 250 ns :: Muon Barrel", 600, 600);
  Muon_Barrel_250ns_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_250ns_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Barrel_250ns_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_250ns_XY->SetTitle("XY-view of SimHits with tof > 250 ns :: Muon Barrel");
  Muon_Barrel_250ns_XY->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF > 250 ns");
  Canvas_Muon_Barrel_250ns_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_250ns_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_250ns_XY->Print(pdfFileName.c_str());}

  Muon_Barrel_000ns_el_XY->SetMarkerColor(kBlack);  Muon_Barrel_000ns_mu_XY->SetMarkerColor(kBlue);  Muon_Barrel_000ns_ha_XY->SetMarkerColor(kRed);
  l2_x1 = 0.15; l2_x2 = 0.35; l2_y1 = 0.825; l2_y2 = 0.875;
  TLegend *l2b = new TLegend(l2_x1, l2_y1,l2_x2,l2_y2,NULL,"brNDC");
  l2b->SetLineColor(0);    l2b->SetLineStyle(0);  l2b->SetLineWidth(0);
  l2b->SetFillColor(4000); l2b->SetBorderSize(0); l2b->SetNColumns(3);
  l2b->AddEntry(Muon_000ns_el_RZ, "e","p");
  l2b->AddEntry(Muon_000ns_mu_RZ, "#mu","p");
  l2b->AddEntry(Muon_000ns_ha_RZ, "had","p");

  Canvas_Muon_Barrel_000ns_Cont_XY = new TCanvas("Canvas_Muon_Barrel_000ns_Cont_XY", "XY-view of SimHits with tof < 250 ns :: Muon Barrel", 600, 600);
  Muon_Barrel_000ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_000ns_el_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Barrel_000ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_000ns_el_XY->SetTitle("XY-view of SimHits with tof < 250 ns :: Muon Barrel");
  Muon_Barrel_000ns_el_XY->SetMarkerStyle(7); Muon_Barrel_000ns_el_XY->SetMarkerColor(kBlack);   Muon_Barrel_000ns_el_XY->SetMarkerSize(1);  Muon_Barrel_000ns_el_XY->Draw("P");
  Muon_Barrel_000ns_mu_XY->SetMarkerStyle(7); Muon_Barrel_000ns_mu_XY->SetMarkerColor(kBlue);    Muon_Barrel_000ns_mu_XY->SetMarkerSize(1);  Muon_Barrel_000ns_mu_XY->Draw("PSame");
  Muon_Barrel_000ns_ha_XY->SetMarkerStyle(7); Muon_Barrel_000ns_ha_XY->SetMarkerColor(kRed);     Muon_Barrel_000ns_ha_XY->SetMarkerSize(1);  Muon_Barrel_000ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 250 ns");
  l2b->Draw();
  Canvas_Muon_Barrel_000ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_000ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_000ns_Cont_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Barrel_250ns_Cont_XY = new TCanvas("Canvas_Muon_Barrel_250ns_Cont_XY", "XY-view of SimHits with tof > 250 ns :: Muon Barrel", 600, 600);
  Muon_Barrel_250ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_250ns_el_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Barrel_250ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_250ns_el_XY->SetTitle("XY-view of SimHits with tof > 250 ns :: Muon Barrel");
  Muon_Barrel_250ns_el_XY->SetMarkerStyle(7); Muon_Barrel_250ns_el_XY->SetMarkerColor(kBlack);   Muon_Barrel_250ns_el_XY->SetMarkerSize(1);  Muon_Barrel_250ns_el_XY->Draw("P");
  Muon_Barrel_250ns_mu_XY->SetMarkerStyle(7); Muon_Barrel_250ns_mu_XY->SetMarkerColor(kBlue);    Muon_Barrel_250ns_mu_XY->SetMarkerSize(1);  Muon_Barrel_250ns_mu_XY->Draw("PSame");
  Muon_Barrel_250ns_ha_XY->SetMarkerStyle(7); Muon_Barrel_250ns_ha_XY->SetMarkerColor(kRed);     Muon_Barrel_250ns_ha_XY->SetMarkerSize(1);  Muon_Barrel_250ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF > 250 ns");
  l2b->Draw();
  Canvas_Muon_Barrel_250ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_250ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_250ns_Cont_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Barrel_00ns_Cont_XY = new TCanvas("Canvas_Muon_Barrel_00ns_Cont_XY", "XY-view of SimHits with tof < 50 ns :: Muon Barrel", 600, 600);
  Muon_Barrel_00ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_00ns_el_XY->GetYaxis()->SetTitle("Y (cm)");
  Muon_Barrel_00ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_00ns_el_XY->SetTitle("XY-view of SimHits with tof < 50 ns :: Muon Barrel");
  Muon_Barrel_00ns_el_XY->SetMarkerStyle(7); Muon_Barrel_00ns_el_XY->SetMarkerColor(kBlack);   Muon_Barrel_00ns_el_XY->SetMarkerSize(1);  Muon_Barrel_00ns_el_XY->Draw("P");
  Muon_Barrel_00ns_mu_XY->SetMarkerStyle(7); Muon_Barrel_00ns_mu_XY->SetMarkerColor(kBlue);    Muon_Barrel_00ns_mu_XY->SetMarkerSize(1);  Muon_Barrel_00ns_mu_XY->Draw("PSame");
  Muon_Barrel_00ns_ha_XY->SetMarkerStyle(7); Muon_Barrel_00ns_ha_XY->SetMarkerColor(kRed);     Muon_Barrel_00ns_ha_XY->SetMarkerSize(1);  Muon_Barrel_00ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 50 ns");
  l2b->Draw();
  Canvas_Muon_Barrel_00ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_00ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_00ns_Cont_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Barrel_50ns_Cont_XY = new TCanvas("Canvas_Muon_Barrel_50ns_Cont_XY", "XY-view of SimHits with 50 < tof < 250 ns :: Muon Barrel", 600, 600);
  Muon_Barrel_50ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Barrel_50ns_el_XY->GetYaxis()->SetTitle("Y (cm)");
  Muon_Barrel_50ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Barrel_50ns_el_XY->SetTitle("XY-view of SimHits with 50 < tof < 250 ns :: Muon Barrel");
  Muon_Barrel_50ns_el_XY->SetMarkerStyle(7); Muon_Barrel_50ns_el_XY->SetMarkerColor(kBlack);   Muon_Barrel_50ns_el_XY->SetMarkerSize(1);  Muon_Barrel_50ns_el_XY->Draw("P");
  Muon_Barrel_50ns_mu_XY->SetMarkerStyle(7); Muon_Barrel_50ns_mu_XY->SetMarkerColor(kBlue);    Muon_Barrel_50ns_mu_XY->SetMarkerSize(1);  Muon_Barrel_50ns_mu_XY->Draw("PSame");
  Muon_Barrel_50ns_ha_XY->SetMarkerStyle(7); Muon_Barrel_50ns_ha_XY->SetMarkerColor(kRed);     Muon_Barrel_50ns_ha_XY->SetMarkerSize(1);  Muon_Barrel_50ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits 50 < TOF < 250 ns");
  l2b->Draw();
  Canvas_Muon_Barrel_50ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Barrel_50ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Barrel_50ns_Cont_XY->Print(pdfFileName.c_str());}
  // -------------------------------- //

  // ----- Muon Endcap XY plots ----- //
  Canvas_Muon_Endcap_XY = new TCanvas("Canvas_Muon_Endcap_XY", "XY-view of SimHits :: Muon Endcap", 600, 600);
  Muon_Endcap_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_XY->GetYaxis()->SetTitle("Y (cm)");
  Muon_Endcap_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_XY->SetTitle("XY-view of SimHits :: Muon Endcap");
  Muon_Endcap_XY->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"All SimHits");
  Canvas_Muon_Endcap_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Endcap_000ns_XY = new TCanvas("Canvas_Muon_Endcap_000ns_XY", "XY-view of SimHits with tof < 250 ns :: Muon Endcap", 600, 600);
  Muon_Endcap_000ns_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_000ns_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Endcap_000ns_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_000ns_XY->SetTitle("XY-view of SimHits with tof < 250 ns :: Muon Endcap");
  Muon_Endcap_000ns_XY->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 250 ns");
  Canvas_Muon_Endcap_000ns_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_000ns_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_000ns_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Endcap_250ns_XY = new TCanvas("Canvas_Muon_Endcap_250ns_XY", "XY-view of SimHits with tof > 250 ns :: Muon Endcap", 600, 600);
  Muon_Endcap_250ns_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_250ns_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Endcap_250ns_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_250ns_XY->SetTitle("XY-view of SimHits with tof > 250 ns :: Muon Endcap");
  Muon_Endcap_250ns_XY->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF > 250 ns");
  Canvas_Muon_Endcap_250ns_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_250ns_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_250ns_XY->Print(pdfFileName.c_str());}

  Muon_Endcap_000ns_el_XY->SetMarkerColor(kBlack);  Muon_Endcap_000ns_mu_XY->SetMarkerColor(kBlue);  Muon_Endcap_000ns_ha_XY->SetMarkerColor(kRed);

  Canvas_Muon_Endcap_000ns_Cont_XY = new TCanvas("Canvas_Muon_Endcap_000ns_Cont_XY", "XY-view of SimHits with tof < 250 ns :: Muon Endcap", 600, 600);
  Muon_Endcap_000ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_000ns_el_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Endcap_000ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_000ns_el_XY->SetTitle("XY-view of SimHits with tof < 250 ns :: Muon Endcap");
  Muon_Endcap_000ns_el_XY->SetMarkerStyle(7); Muon_Endcap_000ns_el_XY->SetMarkerColor(kBlack);   Muon_Endcap_000ns_el_XY->SetMarkerSize(1);  Muon_Endcap_000ns_el_XY->Draw("P");
  Muon_Endcap_000ns_mu_XY->SetMarkerStyle(7); Muon_Endcap_000ns_mu_XY->SetMarkerColor(kBlue);    Muon_Endcap_000ns_mu_XY->SetMarkerSize(1);  Muon_Endcap_000ns_mu_XY->Draw("PSame");
  Muon_Endcap_000ns_ha_XY->SetMarkerStyle(7); Muon_Endcap_000ns_ha_XY->SetMarkerColor(kRed);     Muon_Endcap_000ns_ha_XY->SetMarkerSize(1);  Muon_Endcap_000ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 250 ns");
  l2b->Draw();
  Canvas_Muon_Endcap_000ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_000ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_000ns_Cont_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Endcap_250ns_Cont_XY = new TCanvas("Canvas_Muon_Endcap_250ns_Cont_XY", "XY-view of SimHits with tof > 250 ns :: Muon Endcap", 600, 600);
  Muon_Endcap_250ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_250ns_el_XY->GetYaxis()->SetTitle("R (cm)");
  Muon_Endcap_250ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_250ns_el_XY->SetTitle("XY-view of SimHits with tof > 250 ns :: Muon Endcap");
  Muon_Endcap_250ns_el_XY->SetMarkerStyle(7); Muon_Endcap_250ns_el_XY->SetMarkerColor(kBlack);   Muon_Endcap_250ns_el_XY->SetMarkerSize(1);  Muon_Endcap_250ns_el_XY->Draw("P");
  Muon_Endcap_250ns_mu_XY->SetMarkerStyle(7); Muon_Endcap_250ns_mu_XY->SetMarkerColor(kBlue);    Muon_Endcap_250ns_mu_XY->SetMarkerSize(1);  Muon_Endcap_250ns_mu_XY->Draw("PSame");
  Muon_Endcap_250ns_ha_XY->SetMarkerStyle(7); Muon_Endcap_250ns_ha_XY->SetMarkerColor(kRed);     Muon_Endcap_250ns_ha_XY->SetMarkerSize(1);  Muon_Endcap_250ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF > 250 ns");
  l2b->Draw();
  Canvas_Muon_Endcap_250ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_250ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_250ns_Cont_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Endcap_00ns_Cont_XY = new TCanvas("Canvas_Muon_Endcap_00ns_Cont_XY", "XY-view of SimHits with tof < 50 ns :: Muon Endcap", 600, 600);
  Muon_Endcap_00ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_00ns_el_XY->GetYaxis()->SetTitle("Y (cm)");
  Muon_Endcap_00ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_00ns_el_XY->SetTitle("XY-view of SimHits with tof < 50 ns :: Muon Endcap");
  Muon_Endcap_00ns_el_XY->SetMarkerStyle(7); Muon_Endcap_00ns_el_XY->SetMarkerColor(kBlack);   Muon_Endcap_00ns_el_XY->SetMarkerSize(1);  Muon_Endcap_00ns_el_XY->Draw("P");
  Muon_Endcap_00ns_mu_XY->SetMarkerStyle(7); Muon_Endcap_00ns_mu_XY->SetMarkerColor(kBlue);    Muon_Endcap_00ns_mu_XY->SetMarkerSize(1);  Muon_Endcap_00ns_mu_XY->Draw("PSame");
  Muon_Endcap_00ns_ha_XY->SetMarkerStyle(7); Muon_Endcap_00ns_ha_XY->SetMarkerColor(kRed);     Muon_Endcap_00ns_ha_XY->SetMarkerSize(1);  Muon_Endcap_00ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits TOF < 50 ns");
  l2b->Draw();
  Canvas_Muon_Endcap_00ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_00ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_00ns_Cont_XY->Print(pdfFileName.c_str());}

  Canvas_Muon_Endcap_50ns_Cont_XY = new TCanvas("Canvas_Muon_Endcap_50ns_Cont_XY", "XY-view of SimHits with 50 < tof < 250 ns :: Muon Endcap", 600, 600);
  Muon_Endcap_50ns_el_XY->GetXaxis()->SetTitle("X (cm)");
  Muon_Endcap_50ns_el_XY->GetYaxis()->SetTitle("Y (cm)");
  Muon_Endcap_50ns_el_XY->GetYaxis()->SetTitleOffset(1.30);
  Muon_Endcap_50ns_el_XY->SetTitle("XY-view of SimHits with 50 < tof < 250 ns :: Muon Endcap");
  Muon_Endcap_50ns_el_XY->SetMarkerStyle(7); Muon_Endcap_50ns_el_XY->SetMarkerColor(kBlack);   Muon_Endcap_50ns_el_XY->SetMarkerSize(1);  Muon_Endcap_50ns_el_XY->Draw("P");
  Muon_Endcap_50ns_mu_XY->SetMarkerStyle(7); Muon_Endcap_50ns_mu_XY->SetMarkerColor(kBlue);    Muon_Endcap_50ns_mu_XY->SetMarkerSize(1);  Muon_Endcap_50ns_mu_XY->Draw("PSame");
  Muon_Endcap_50ns_ha_XY->SetMarkerStyle(7); Muon_Endcap_50ns_ha_XY->SetMarkerColor(kRed);     Muon_Endcap_50ns_ha_XY->SetMarkerSize(1);  Muon_Endcap_50ns_ha_XY->Draw("PSame");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.15, 0.15,"SimHits 50 < TOF < 250 ns");
  l2b->Draw();
  Canvas_Muon_Endcap_50ns_Cont_XY->SetTicks(1,1);
  Canvas_Muon_Endcap_50ns_Cont_XY->Write(); 
  if(pdf_output) {Canvas_Muon_Endcap_50ns_Cont_XY->Print(pdfFileName.c_str());}
  // -------------------------------- //


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

  TDir_Muon_XY_RZ_views->cd();
  // --------------------------- 
  Muon_00ns_el_RZ->Write();
  Muon_00ns_mu_RZ->Write();
  Muon_00ns_ha_RZ->Write();
  Muon_50ns_el_RZ->Write();
  Muon_50ns_mu_RZ->Write();
  Muon_50ns_ha_RZ->Write();
  // --------------------------- 

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

  /*
  Canvas_SimVertices_RZ = new TCanvas("Canvas_SimVertices_RZ", "RZ-view of Sim Vertices", 600, 600);
  SimVertices_RZ->GetXaxis()->SetTitle("Z (cm)");
  SimVertices_RZ->GetYaxis()->SetTitle("R (cm)");
  SimVertices_RZ->GetYaxis()->SetTitleOffset(1.30);
  SimVertices_RZ->SetTitle("RZ-view of Sim Vertices");
  SimVertices_RZ->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_SimVertices_RZ->SetTicks(1,1);
  Canvas_SimVertices_RZ->Write(); 
  if(pdf_output) {Canvas_SimVertices_RZ->Print(pdfFileName.c_str());}

  Canvas_SimVertices_Muon_RZ = new TCanvas("Canvas_SimVertices_Muon_RZ", "RZ-view of Sim Vertices in Muon System and Cavern", 600, 600);
  SimVertices_Muon_RZ->GetXaxis()->SetTitle("Z (cm)");
  SimVertices_Muon_RZ->GetYaxis()->SetTitle("R (cm)");
  SimVertices_Muon_RZ->GetYaxis()->SetTitleOffset(1.30);
  SimVertices_Muon_RZ->SetTitle("RZ-view of Sim Vertices in Muon System and Cavern");
  // SimVertices_Muon_RZ->SetMarkerStyle(5);
  // SimVertices_Muon_RZ->SetMarkerColor(kBlack);
  // SimVertices_Muon_RZ->SetMarkerSize(1.25);
  // SimVertices_Muon_RZ->Draw("P");
  SimVertices_Muon_RZ->Draw("colz");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_SimVertices_Muon_RZ->SetTicks(1,1);
  Canvas_SimVertices_Muon_RZ->Write(); 
  if(pdf_output) {Canvas_SimVertices_Muon_RZ->Print(pdfFileName.c_str());}
  PrimVertices_Z->Write();
  PrimVertices_R->Write();
  */

  TDir_Muon_XY_RZ_views->cd();
  // --------------------------- 
  RPCb_XY->Write();
  RPCb_RZ->Write();
  RPCf_XY->Write();
  RPCf_RZ->Write();
  CSC_XY->Write();
  CSC_RZ->Write();
  DT_XY->Write();
  DT_RZ->Write();
  GEM_XY->Write();
  GEM_RZ->Write();
  ME0_XY->Write();
  ME0_RZ->Write();
  Muon_Barrel_XY->Write();
  Muon_Endcap_XY->Write();
  Muon_RZ->Write();

  RPCb_000ns_XY->Write();
  RPCb_000ns_RZ->Write();
  RPCf_000ns_XY->Write();
  RPCf_000ns_RZ->Write();
  CSC_000ns_XY->Write();
  CSC_000ns_RZ->Write();
  DT_000ns_XY->Write();
  DT_000ns_RZ->Write();
  GEM_000ns_XY->Write();
  GEM_000ns_RZ->Write();
  ME0_000ns_XY->Write();
  ME0_000ns_RZ->Write();
  Muon_Barrel_000ns_XY->Write();
  Muon_Endcap_000ns_XY->Write();
  Muon_000ns_RZ->Write();

  RPCb_250ns_XY->Write();
  RPCb_250ns_RZ->Write();
  RPCf_250ns_XY->Write();
  RPCf_250ns_RZ->Write();
  CSC_250ns_XY->Write();
  CSC_250ns_RZ->Write();
  DT_250ns_XY->Write();
  DT_250ns_RZ->Write();
  GEM_250ns_XY->Write();
  GEM_250ns_RZ->Write();
  ME0_250ns_XY->Write();
  ME0_250ns_RZ->Write();
  Muon_Barrel_250ns_XY->Write();
  Muon_Endcap_250ns_XY->Write();
  Muon_250ns_RZ->Write();
  // --------------------------- 
  outputfile->cd();

  for(int m=0; m<15; ++m) {

    RPCb_Muons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    RPCf_Muons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    CSC_Muons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    DT_Muons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    GEM_Muons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    ME0_Muons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());

    RPCb_Hadrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    RPCf_Hadrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    CSC_Hadrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    DT_Hadrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    GEM_Hadrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    ME0_Hadrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());

    RPCb_Electrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    RPCf_Electrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    CSC_Electrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    DT_Electrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    GEM_Electrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    ME0_Electrons_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());

    RPCb_Electrons_250ns_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    RPCf_Electrons_250ns_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    CSC_Electrons_250ns_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    DT_Electrons_250ns_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    GEM_Electrons_250ns_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
    ME0_Electrons_250ns_SHPT->GetXaxis()->SetBinLabel( m+1, proc[m].c_str());
  }

  TDir_Muon_hit_info->cd();
  // --------------------------- 
  RPCb_Muons_SHPT->Write();
  RPCf_Muons_SHPT->Write();
  CSC_Muons_SHPT->Write();
  DT_Muons_SHPT->Write();
  GEM_Muons_SHPT->Write();
  ME0_Muons_SHPT->Write();

  RPCb_Hadrons_SHPT->Write();
  RPCf_Hadrons_SHPT->Write();
  CSC_Hadrons_SHPT->Write();
  DT_Hadrons_SHPT->Write();
  GEM_Hadrons_SHPT->Write();
  ME0_Hadrons_SHPT->Write();

  RPCb_Electrons_SHPT->Write();
  RPCf_Electrons_SHPT->Write();
  CSC_Electrons_SHPT->Write();
  DT_Electrons_SHPT->Write();
  GEM_Electrons_SHPT->Write();
  ME0_Electrons_SHPT->Write();

  RPCb_Electrons_000ns_SHPT->Write();
  RPCf_Electrons_000ns_SHPT->Write();
  CSC_Electrons_000ns_SHPT->Write();
  DT_Electrons_000ns_SHPT->Write();
  GEM_Electrons_000ns_SHPT->Write();
  ME0_Electrons_000ns_SHPT->Write();

  RPCb_Electrons_050ns_SHPT->Write();
  RPCf_Electrons_050ns_SHPT->Write();
  CSC_Electrons_050ns_SHPT->Write();
  DT_Electrons_050ns_SHPT->Write();
  GEM_Electrons_050ns_SHPT->Write();
  ME0_Electrons_050ns_SHPT->Write();

  RPCb_Electrons_250ns_SHPT->Write();
  RPCf_Electrons_250ns_SHPT->Write();
  CSC_Electrons_250ns_SHPT->Write();
  DT_Electrons_250ns_SHPT->Write();
  GEM_Electrons_250ns_SHPT->Write();
  ME0_Electrons_250ns_SHPT->Write();
  // ---------------------------
  outputfile->cd();


  TDir_Muon_entry_info->cd();
  RPCb_EntryExit_All_Glob_dr->Write();
  DT_EntryExit_All_Glob_dr->Write();
  RPCf_EntryExit_All_Glob_dz->Write();
  CSC_EntryExit_All_Glob_dz->Write();
  GEM_EntryExit_All_Glob_dz->Write();
  ME0_EntryExit_All_Glob_dz->Write();

  CSC_EntryExit_Electrons_Glob_dz->Write();
  CSC_EntryExit_Muons_Glob_dz->Write();
  CSC_EntryExit_Hadrons_Glob_dz->Write();

  CSC_EntryExit_All_Glob_dGap->Write();
  CSC_EntryExit_Electrons_Glob_dGap->Write();
  CSC_EntryExit_Muons_Glob_dGap->Write();
  CSC_EntryExit_Hadrons_Glob_dGap->Write();

  RPCb_EntryExit_All_Loc_dz->Write();
  RPCf_EntryExit_All_Loc_dz->Write();
  CSC_EntryExit_All_Loc_dz->Write();
  DT_EntryExit_All_Loc_dz->Write();
  GEM_EntryExit_All_Loc_dz->Write();
  ME0_EntryExit_All_Loc_dz->Write();

  CSC_EntryExit_Electrons_Loc_dz->Write();
  CSC_EntryExit_Muons_Loc_dz->Write();
  CSC_EntryExit_Hadrons_Loc_dz->Write();

  /*
  CSC_EntryExit_All_Deposit_dz->Write();
  CSC_EntryExit_Electrons_Deposit_dz->Write();
  CSC_EntryExit_Muons_Deposit_dz->Write();
  CSC_EntryExit_Hadrons_Deposit_dz->Write();
  */
  CSC_EntryExit_el_Deposit_dz->Write();
  CSC_EntryExit_mu_Deposit_dz->Write();
  CSC_EntryExit_pi_Deposit_dz->Write();
  CSC_EntryExit_ka_Deposit_dz->Write();
  CSC_EntryExit_p_Deposit_dz->Write();
  CSC_EntryExit_n_Deposit_dz->Write();
  CSC_EntryExit_g_Deposit_dz->Write();
  CSC_EntryExit_N_Deposit_dz->Write();

  CSC_EntryExit_el_KinEn_dz->Write();
  CSC_EntryExit_mu_KinEn_dz->Write();
  CSC_EntryExit_pi_KinEn_dz->Write();
  CSC_EntryExit_ka_KinEn_dz->Write();
  CSC_EntryExit_p_KinEn_dz->Write();
  CSC_EntryExit_n_KinEn_dz->Write();
  CSC_EntryExit_g_KinEn_dz->Write();
  CSC_EntryExit_N_KinEn_dz->Write();

  CSC_EntryExit_el_Time_dz->Write();
  CSC_EntryExit_mu_Time_dz->Write();
  CSC_EntryExit_pi_Time_dz->Write();
  CSC_EntryExit_ka_Time_dz->Write();
  CSC_EntryExit_p_Time_dz->Write();
  CSC_EntryExit_n_Time_dz->Write();
  CSC_EntryExit_g_Time_dz->Write();
  CSC_EntryExit_N_Time_dz->Write();

  CSC_EntryExit_el_Deposit_dR->Write();
  CSC_EntryExit_mu_Deposit_dR->Write();
  CSC_EntryExit_pi_Deposit_dR->Write();
  CSC_EntryExit_ka_Deposit_dR->Write();
  CSC_EntryExit_p_Deposit_dR->Write();
  CSC_EntryExit_n_Deposit_dR->Write();
  CSC_EntryExit_g_Deposit_dR->Write();
  CSC_EntryExit_N_Deposit_dR->Write();

  CSC_EntryExit_el_KinEn_dR->Write();
  CSC_EntryExit_mu_KinEn_dR->Write();
  CSC_EntryExit_pi_KinEn_dR->Write();
  CSC_EntryExit_ka_KinEn_dR->Write();
  CSC_EntryExit_p_KinEn_dR->Write();
  CSC_EntryExit_n_KinEn_dR->Write();
  CSC_EntryExit_g_KinEn_dR->Write();
  CSC_EntryExit_N_KinEn_dR->Write();

  CSC_EntryExit_el_dz_dR->Write();
  CSC_EntryExit_mu_dz_dR->Write();
  CSC_EntryExit_pi_dz_dR->Write();
  CSC_EntryExit_ka_dz_dR->Write();
  CSC_EntryExit_p_dz_dR->Write();
  CSC_EntryExit_n_dz_dR->Write();
  CSC_EntryExit_g_dz_dR->Write();
  CSC_EntryExit_N_dz_dR->Write();

  CSC_EntryExit_el_dz_dR_detail->Write();
  CSC_EntryExit_mu_dz_dR_detail->Write();
  CSC_EntryExit_pi_dz_dR_detail->Write();
  CSC_EntryExit_ka_dz_dR_detail->Write();
  CSC_EntryExit_p_dz_dR_detail->Write();
  CSC_EntryExit_n_dz_dR_detail->Write();
  CSC_EntryExit_g_dz_dR_detail->Write();
  CSC_EntryExit_N_dz_dR_detail->Write();

  CSC_EntryExit_el_GapLength_Deposit->Write();
  CSC_EntryExit_mu_GapLength_Deposit->Write();
  CSC_EntryExit_pi_GapLength_Deposit->Write();
  CSC_EntryExit_ka_GapLength_Deposit->Write();
  CSC_EntryExit_p_GapLength_Deposit->Write();
  CSC_EntryExit_n_GapLength_Deposit->Write();
  CSC_EntryExit_g_GapLength_Deposit->Write();
  CSC_EntryExit_N_GapLength_Deposit->Write();

  CSC_EntryExit_el_Deposit_pidR2->Write();
  CSC_EntryExit_mu_Deposit_pidR2->Write();
  CSC_EntryExit_pi_Deposit_pidR2->Write();
  CSC_EntryExit_ka_Deposit_pidR2->Write();
  CSC_EntryExit_p_Deposit_pidR2->Write();
  CSC_EntryExit_n_Deposit_pidR2->Write();
  CSC_EntryExit_g_Deposit_pidR2->Write();
  CSC_EntryExit_N_Deposit_pidR2->Write();

  CSC_EntryExit_el_Deposit_dRdz->Write();
  CSC_EntryExit_mu_Deposit_dRdz->Write();
  CSC_EntryExit_pi_Deposit_dRdz->Write();
  CSC_EntryExit_ka_Deposit_dRdz->Write();
  CSC_EntryExit_p_Deposit_dRdz->Write();
  CSC_EntryExit_n_Deposit_dRdz->Write();
  CSC_EntryExit_g_Deposit_dRdz->Write();
  CSC_EntryExit_N_Deposit_dRdz->Write();


  // --------------------------- 
  outputfile->cd();

  Canvas_CSC_EntryExit_vs_deposits  = new TCanvas("Canvas_CSC_EntryExit_vs_deposits",  "#Delta z (Global Coords) vs E_{deposit} :: CSC", 600, 600);
  Canvas_CSC_EntryExit_vs_kinenergy = new TCanvas("Canvas_CSC_EntryExit_vs_kinenergy", "#Delta z (Global Coords) vs E_{kin} :: CSC", 600, 600);
  Canvas_CSC_EntryExit_vs_time      = new TCanvas("Canvas_CSC_EntryExit_vs_time",      "#Delta z (Global Coords) vs TOF :: CSC", 600, 600);

  Canvas_CSC_EntryExit_vs_deposits->cd();
  CSC_EntryExit_el_Deposit_dz->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  CSC_EntryExit_el_Deposit_dz->GetYaxis()->SetTitle("#Delta z (Global Coords) (cm)");
  CSC_EntryExit_el_Deposit_dz->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_Deposit_dz->SetTitleOffset(1.2);
  CSC_EntryExit_el_Deposit_dz->SetTitle("SimHit #Delta z vs E_{deposit} :: CSC");
  CSC_EntryExit_el_Deposit_dz->SetMarkerStyle(7);  CSC_EntryExit_el_Deposit_dz->SetMarkerColor(kBlack);   CSC_EntryExit_el_Deposit_dz->SetMarkerSize(1);  CSC_EntryExit_el_Deposit_dz->Draw("P");
  CSC_EntryExit_mu_Deposit_dz->SetMarkerStyle(24); CSC_EntryExit_mu_Deposit_dz->SetMarkerColor(kBlue);    CSC_EntryExit_mu_Deposit_dz->SetMarkerSize(1);  CSC_EntryExit_mu_Deposit_dz->Draw("PSame");
  CSC_EntryExit_pi_Deposit_dz->SetMarkerStyle(33); CSC_EntryExit_pi_Deposit_dz->SetMarkerColor(kGreen);   CSC_EntryExit_pi_Deposit_dz->SetMarkerSize(1);  CSC_EntryExit_pi_Deposit_dz->Draw("PSame");
  CSC_EntryExit_ka_Deposit_dz->SetMarkerStyle(5);  CSC_EntryExit_ka_Deposit_dz->SetMarkerColor(kOrange);  CSC_EntryExit_ka_Deposit_dz->SetMarkerSize(1);  CSC_EntryExit_ka_Deposit_dz->Draw("PSame");
  CSC_EntryExit_p_Deposit_dz->SetMarkerStyle(26);  CSC_EntryExit_p_Deposit_dz->SetMarkerColor(kMagenta);  CSC_EntryExit_p_Deposit_dz->SetMarkerSize(1);   CSC_EntryExit_p_Deposit_dz->Draw("PSame");
  CSC_EntryExit_n_Deposit_dz->SetMarkerStyle(32);  CSC_EntryExit_n_Deposit_dz->SetMarkerColor(kViolet);   CSC_EntryExit_n_Deposit_dz->SetMarkerSize(1);   CSC_EntryExit_n_Deposit_dz->Draw("PSame");
  CSC_EntryExit_g_Deposit_dz->SetMarkerStyle(30);  CSC_EntryExit_g_Deposit_dz->SetMarkerColor(kCyan);     CSC_EntryExit_g_Deposit_dz->SetMarkerSize(1);   CSC_EntryExit_g_Deposit_dz->Draw("PSame");
  CSC_EntryExit_N_Deposit_dz->SetMarkerStyle(2);   CSC_EntryExit_N_Deposit_dz->SetMarkerColor(kRed);      CSC_EntryExit_N_Deposit_dz->SetMarkerSize(1);   CSC_EntryExit_N_Deposit_dz->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_EntryExit_vs_deposits->SetTicks(1,1);
  Canvas_CSC_EntryExit_vs_deposits->Write();
  if(pdf_output) {Canvas_CSC_EntryExit_vs_deposits->Print(pdfFileName.c_str());}

  Canvas_CSC_EntryExit_vs_kinenergy->cd();
  CSC_EntryExit_el_KinEn_dz->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  CSC_EntryExit_el_KinEn_dz->GetYaxis()->SetTitle("#Delta z (Global Coords) (cm)");
  CSC_EntryExit_el_KinEn_dz->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_KinEn_dz->SetTitleOffset(1.2);
  CSC_EntryExit_el_KinEn_dz->SetTitle("SimHit #Delta z vs E_{kin} :: CSC");
  CSC_EntryExit_el_KinEn_dz->SetMarkerStyle(7);  CSC_EntryExit_el_KinEn_dz->SetMarkerColor(kBlack);   CSC_EntryExit_el_KinEn_dz->SetMarkerSize(1);  CSC_EntryExit_el_KinEn_dz->Draw("P");
  CSC_EntryExit_mu_KinEn_dz->SetMarkerStyle(24); CSC_EntryExit_mu_KinEn_dz->SetMarkerColor(kBlue);    CSC_EntryExit_mu_KinEn_dz->SetMarkerSize(1);  CSC_EntryExit_mu_KinEn_dz->Draw("PSame");
  CSC_EntryExit_pi_KinEn_dz->SetMarkerStyle(33); CSC_EntryExit_pi_KinEn_dz->SetMarkerColor(kGreen);   CSC_EntryExit_pi_KinEn_dz->SetMarkerSize(1);  CSC_EntryExit_pi_KinEn_dz->Draw("PSame");
  CSC_EntryExit_ka_KinEn_dz->SetMarkerStyle(5);  CSC_EntryExit_ka_KinEn_dz->SetMarkerColor(kOrange);  CSC_EntryExit_ka_KinEn_dz->SetMarkerSize(1);  CSC_EntryExit_ka_KinEn_dz->Draw("PSame");
  CSC_EntryExit_p_KinEn_dz->SetMarkerStyle(26);  CSC_EntryExit_p_KinEn_dz->SetMarkerColor(kMagenta);  CSC_EntryExit_p_KinEn_dz->SetMarkerSize(1);   CSC_EntryExit_p_KinEn_dz->Draw("PSame");
  CSC_EntryExit_n_KinEn_dz->SetMarkerStyle(32);  CSC_EntryExit_n_KinEn_dz->SetMarkerColor(kViolet);   CSC_EntryExit_n_KinEn_dz->SetMarkerSize(1);   CSC_EntryExit_n_KinEn_dz->Draw("PSame");
  CSC_EntryExit_g_KinEn_dz->SetMarkerStyle(30);  CSC_EntryExit_g_KinEn_dz->SetMarkerColor(kCyan);     CSC_EntryExit_g_KinEn_dz->SetMarkerSize(1);   CSC_EntryExit_g_KinEn_dz->Draw("PSame");
  CSC_EntryExit_N_KinEn_dz->SetMarkerStyle(2);   CSC_EntryExit_N_KinEn_dz->SetMarkerColor(kRed);      CSC_EntryExit_N_KinEn_dz->SetMarkerSize(1);   CSC_EntryExit_N_KinEn_dz->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_EntryExit_vs_kinenergy->SetTicks(1,1);
  Canvas_CSC_EntryExit_vs_kinenergy->Write();
  if(pdf_output) {Canvas_CSC_EntryExit_vs_kinenergy->Print(pdfFileName.c_str());}

  Canvas_CSC_EntryExit_vs_time->cd();
  CSC_EntryExit_el_Time_dz->GetXaxis()->SetTitle("#Delta z (Global Coords) (cm)");
  CSC_EntryExit_el_Time_dz->GetYaxis()->SetTitle("^{10}log TOF (ns)");
  CSC_EntryExit_el_Time_dz->GetXaxis()->SetTitleOffset(1.2);
  CSC_EntryExit_el_Time_dz->SetTitle("SimHit time vs E_{deposit} :: CSC");
  CSC_EntryExit_el_Time_dz->SetMarkerStyle(7);  CSC_EntryExit_el_Time_dz->SetMarkerColor(kBlack);   CSC_EntryExit_el_Time_dz->SetMarkerSize(1);  CSC_EntryExit_el_Time_dz->Draw("P");
  CSC_EntryExit_mu_Time_dz->SetMarkerStyle(24); CSC_EntryExit_mu_Time_dz->SetMarkerColor(kBlue);    CSC_EntryExit_mu_Time_dz->SetMarkerSize(1);  CSC_EntryExit_mu_Time_dz->Draw("PSame");
  CSC_EntryExit_pi_Time_dz->SetMarkerStyle(33); CSC_EntryExit_pi_Time_dz->SetMarkerColor(kGreen);   CSC_EntryExit_pi_Time_dz->SetMarkerSize(1);  CSC_EntryExit_pi_Time_dz->Draw("PSame");
  CSC_EntryExit_ka_Time_dz->SetMarkerStyle(5);  CSC_EntryExit_ka_Time_dz->SetMarkerColor(kOrange);  CSC_EntryExit_ka_Time_dz->SetMarkerSize(1);  CSC_EntryExit_ka_Time_dz->Draw("PSame");
  CSC_EntryExit_p_Time_dz->SetMarkerStyle(26);  CSC_EntryExit_p_Time_dz->SetMarkerColor(kMagenta);  CSC_EntryExit_p_Time_dz->SetMarkerSize(1);   CSC_EntryExit_p_Time_dz->Draw("PSame");
  CSC_EntryExit_n_Time_dz->SetMarkerStyle(32);  CSC_EntryExit_n_Time_dz->SetMarkerColor(kViolet);   CSC_EntryExit_n_Time_dz->SetMarkerSize(1);   CSC_EntryExit_n_Time_dz->Draw("PSame");
  CSC_EntryExit_g_Time_dz->SetMarkerStyle(30);  CSC_EntryExit_g_Time_dz->SetMarkerColor(kCyan);     CSC_EntryExit_g_Time_dz->SetMarkerSize(1);   CSC_EntryExit_g_Time_dz->Draw("PSame");
  CSC_EntryExit_N_Time_dz->SetMarkerStyle(2);   CSC_EntryExit_N_Time_dz->SetMarkerColor(kRed);      CSC_EntryExit_N_Time_dz->SetMarkerSize(1);   CSC_EntryExit_N_Time_dz->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,"#font[12]{250 ns}");
  latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,"#font[12]{1 ms}");
  latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,"#font[12]{neutron background}");
  latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_EntryExit_vs_time->SetTicks(1,1);
  Canvas_CSC_EntryExit_vs_time->Write();
  if(pdf_output) {Canvas_CSC_EntryExit_vs_time->Print(pdfFileName.c_str());}


  Canvas_CSC_PathLength_vs_deposits  = new TCanvas("Canvas_CSC_PathLength_vs_deposits",  "#Delta z (Global Coords) vs E_{deposit} :: CSC", 600, 600);
  Canvas_CSC_PathLength_vs_kinenergy = new TCanvas("Canvas_CSC_PathLength_vs_kinenergy", "#Delta z (Global Coords) vs E_{kin} :: CSC", 600, 600);

  Canvas_CSC_PathLength_vs_deposits->cd();
  CSC_EntryExit_el_Deposit_dR->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  CSC_EntryExit_el_Deposit_dR->GetYaxis()->SetTitle("#Delta R (Global Coords) (cm)");
  CSC_EntryExit_el_Deposit_dR->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_Deposit_dR->SetTitleOffset(1.2);
  CSC_EntryExit_el_Deposit_dR->SetTitle("SimHit #Delta z vs E_{deposit} :: CSC");
  CSC_EntryExit_el_Deposit_dR->SetMarkerStyle(7);  CSC_EntryExit_el_Deposit_dR->SetMarkerColor(kBlack);   CSC_EntryExit_el_Deposit_dR->SetMarkerSize(1);  CSC_EntryExit_el_Deposit_dR->Draw("P");
  CSC_EntryExit_mu_Deposit_dR->SetMarkerStyle(24); CSC_EntryExit_mu_Deposit_dR->SetMarkerColor(kBlue);    CSC_EntryExit_mu_Deposit_dR->SetMarkerSize(1);  CSC_EntryExit_mu_Deposit_dR->Draw("PSame");
  CSC_EntryExit_pi_Deposit_dR->SetMarkerStyle(33); CSC_EntryExit_pi_Deposit_dR->SetMarkerColor(kGreen);   CSC_EntryExit_pi_Deposit_dR->SetMarkerSize(1);  CSC_EntryExit_pi_Deposit_dR->Draw("PSame");
  CSC_EntryExit_ka_Deposit_dR->SetMarkerStyle(5);  CSC_EntryExit_ka_Deposit_dR->SetMarkerColor(kOrange);  CSC_EntryExit_ka_Deposit_dR->SetMarkerSize(1);  CSC_EntryExit_ka_Deposit_dR->Draw("PSame");
  CSC_EntryExit_p_Deposit_dR->SetMarkerStyle(26);  CSC_EntryExit_p_Deposit_dR->SetMarkerColor(kMagenta);  CSC_EntryExit_p_Deposit_dR->SetMarkerSize(1);   CSC_EntryExit_p_Deposit_dR->Draw("PSame");
  CSC_EntryExit_n_Deposit_dR->SetMarkerStyle(32);  CSC_EntryExit_n_Deposit_dR->SetMarkerColor(kViolet);   CSC_EntryExit_n_Deposit_dR->SetMarkerSize(1);   CSC_EntryExit_n_Deposit_dR->Draw("PSame");
  CSC_EntryExit_g_Deposit_dR->SetMarkerStyle(30);  CSC_EntryExit_g_Deposit_dR->SetMarkerColor(kCyan);     CSC_EntryExit_g_Deposit_dR->SetMarkerSize(1);   CSC_EntryExit_g_Deposit_dR->Draw("PSame");
  CSC_EntryExit_N_Deposit_dR->SetMarkerStyle(2);   CSC_EntryExit_N_Deposit_dR->SetMarkerColor(kRed);      CSC_EntryExit_N_Deposit_dR->SetMarkerSize(1);   CSC_EntryExit_N_Deposit_dR->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_PathLength_vs_deposits->SetTicks(1,1);
  Canvas_CSC_PathLength_vs_deposits->Write();
  if(pdf_output) {Canvas_CSC_PathLength_vs_deposits->Print(pdfFileName.c_str());}


  Canvas_CSC_PathLength_vs_kinenergy->cd();
  CSC_EntryExit_el_KinEn_dR->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  CSC_EntryExit_el_KinEn_dR->GetYaxis()->SetTitle("#Delta R (Global Coords) (cm)");
  CSC_EntryExit_el_KinEn_dR->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_KinEn_dR->SetTitleOffset(1.2);
  CSC_EntryExit_el_KinEn_dR->SetTitle("SimHit #Delta z vs E_{kin} :: CSC");
  CSC_EntryExit_el_KinEn_dR->SetMarkerStyle(7);  CSC_EntryExit_el_KinEn_dR->SetMarkerColor(kBlack);   CSC_EntryExit_el_KinEn_dR->SetMarkerSize(1);  CSC_EntryExit_el_KinEn_dR->Draw("P");
  CSC_EntryExit_mu_KinEn_dR->SetMarkerStyle(24); CSC_EntryExit_mu_KinEn_dR->SetMarkerColor(kBlue);    CSC_EntryExit_mu_KinEn_dR->SetMarkerSize(1);  CSC_EntryExit_mu_KinEn_dR->Draw("PSame");
  CSC_EntryExit_pi_KinEn_dR->SetMarkerStyle(33); CSC_EntryExit_pi_KinEn_dR->SetMarkerColor(kGreen);   CSC_EntryExit_pi_KinEn_dR->SetMarkerSize(1);  CSC_EntryExit_pi_KinEn_dR->Draw("PSame");
  CSC_EntryExit_ka_KinEn_dR->SetMarkerStyle(5);  CSC_EntryExit_ka_KinEn_dR->SetMarkerColor(kOrange);  CSC_EntryExit_ka_KinEn_dR->SetMarkerSize(1);  CSC_EntryExit_ka_KinEn_dR->Draw("PSame");
  CSC_EntryExit_p_KinEn_dR->SetMarkerStyle(26);  CSC_EntryExit_p_KinEn_dR->SetMarkerColor(kMagenta);  CSC_EntryExit_p_KinEn_dR->SetMarkerSize(1);   CSC_EntryExit_p_KinEn_dR->Draw("PSame");
  CSC_EntryExit_n_KinEn_dR->SetMarkerStyle(32);  CSC_EntryExit_n_KinEn_dR->SetMarkerColor(kViolet);   CSC_EntryExit_n_KinEn_dR->SetMarkerSize(1);   CSC_EntryExit_n_KinEn_dR->Draw("PSame");
  CSC_EntryExit_g_KinEn_dR->SetMarkerStyle(30);  CSC_EntryExit_g_KinEn_dR->SetMarkerColor(kCyan);     CSC_EntryExit_g_KinEn_dR->SetMarkerSize(1);   CSC_EntryExit_g_KinEn_dR->Draw("PSame");
  CSC_EntryExit_N_KinEn_dR->SetMarkerStyle(2);   CSC_EntryExit_N_KinEn_dR->SetMarkerColor(kRed);      CSC_EntryExit_N_KinEn_dR->SetMarkerSize(1);   CSC_EntryExit_N_KinEn_dR->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_PathLength_vs_kinenergy->SetTicks(1,1);
  Canvas_CSC_PathLength_vs_kinenergy->Write();
  if(pdf_output) {Canvas_CSC_PathLength_vs_kinenergy->Print(pdfFileName.c_str());}


  Canvas_CSC_dz_vs_dR              = new TCanvas("Canvas_CSC_dz_vs_dR",              "#Delta R (Global Coords) vs dz (Global Coords) :: CSC", 600, 600);
  Canvas_CSC_dz_vs_dR_detail       = new TCanvas("Canvas_CSC_dz_vs_dR_detail",       "#Delta R (Global Coords) vs dz (Global Coords) :: CSC", 600, 600);
  Canvas_CSC_deposits_vs_GapLength = new TCanvas("Canvas_CSC_deposits_vs_GapLength", "#Delta Gap (Global Coords) vs E_{deposit} :: CSC", 600, 600);
  Canvas_CSC_deposits_vs_pidR2     = new TCanvas("Canvas_CSC_deposits_vs_pidR2",     "#pi (#frac{#Delta R}{2})^{2} (Global Coords) vs E_{deposit} :: CSC (dz==0)", 600, 600);
  Canvas_CSC_deposits_vs_dRdz      = new TCanvas("Canvas_CSC_deposits_vs_dRdz",      "#sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) vs E_{deposit} :: CSC (dz!=0)", 600, 600);

  Canvas_CSC_dz_vs_dR->cd();
  CSC_EntryExit_el_dz_dR->GetXaxis()->SetTitle("#Delta z (Global Coords) (cm)");
  CSC_EntryExit_el_dz_dR->GetYaxis()->SetTitle("#Delta R (Global Coords) (cm)");
  CSC_EntryExit_el_dz_dR->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_dz_dR->SetTitleOffset(1.2);
  CSC_EntryExit_el_dz_dR->SetTitle("SimHit #Delta R vs #Delta z :: CSC");
  CSC_EntryExit_el_dz_dR->SetMarkerStyle(7);  CSC_EntryExit_el_dz_dR->SetMarkerColor(kBlack);   CSC_EntryExit_el_dz_dR->SetMarkerSize(1);  CSC_EntryExit_el_dz_dR->Draw("P");
  CSC_EntryExit_mu_dz_dR->SetMarkerStyle(24); CSC_EntryExit_mu_dz_dR->SetMarkerColor(kBlue);    CSC_EntryExit_mu_dz_dR->SetMarkerSize(1);  CSC_EntryExit_mu_dz_dR->Draw("PSame");
  CSC_EntryExit_pi_dz_dR->SetMarkerStyle(33); CSC_EntryExit_pi_dz_dR->SetMarkerColor(kGreen);   CSC_EntryExit_pi_dz_dR->SetMarkerSize(1);  CSC_EntryExit_pi_dz_dR->Draw("PSame");
  CSC_EntryExit_ka_dz_dR->SetMarkerStyle(5);  CSC_EntryExit_ka_dz_dR->SetMarkerColor(kOrange);  CSC_EntryExit_ka_dz_dR->SetMarkerSize(1);  CSC_EntryExit_ka_dz_dR->Draw("PSame");
  CSC_EntryExit_p_dz_dR->SetMarkerStyle(26);  CSC_EntryExit_p_dz_dR->SetMarkerColor(kMagenta);  CSC_EntryExit_p_dz_dR->SetMarkerSize(1);   CSC_EntryExit_p_dz_dR->Draw("PSame");
  CSC_EntryExit_n_dz_dR->SetMarkerStyle(32);  CSC_EntryExit_n_dz_dR->SetMarkerColor(kViolet);   CSC_EntryExit_n_dz_dR->SetMarkerSize(1);   CSC_EntryExit_n_dz_dR->Draw("PSame");
  CSC_EntryExit_g_dz_dR->SetMarkerStyle(30);  CSC_EntryExit_g_dz_dR->SetMarkerColor(kCyan);     CSC_EntryExit_g_dz_dR->SetMarkerSize(1);   CSC_EntryExit_g_dz_dR->Draw("PSame");
  CSC_EntryExit_N_dz_dR->SetMarkerStyle(2);   CSC_EntryExit_N_dz_dR->SetMarkerColor(kRed);      CSC_EntryExit_N_dz_dR->SetMarkerSize(1);   CSC_EntryExit_N_dz_dR->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_dz_vs_dR->SetTicks(1,1);
  Canvas_CSC_dz_vs_dR->Write();
  if(pdf_output) {Canvas_CSC_dz_vs_dR->Print(pdfFileName.c_str());}

  Canvas_CSC_dz_vs_dR_detail->cd();
  CSC_EntryExit_el_dz_dR_detail->GetXaxis()->SetTitle("#Delta z (Global Coords) (cm)");
  CSC_EntryExit_el_dz_dR_detail->GetYaxis()->SetTitle("#Delta R (Global Coords) (cm)");
  CSC_EntryExit_el_dz_dR_detail->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_dz_dR_detail->SetTitleOffset(1.2);
  CSC_EntryExit_el_dz_dR_detail->SetTitle("SimHit #Delta R vs #Delta z :: CSC");
  CSC_EntryExit_el_dz_dR_detail->SetMarkerStyle(7);  CSC_EntryExit_el_dz_dR_detail->SetMarkerColor(kBlack);   CSC_EntryExit_el_dz_dR_detail->SetMarkerSize(1);  CSC_EntryExit_el_dz_dR_detail->Draw("P");
  CSC_EntryExit_mu_dz_dR_detail->SetMarkerStyle(24); CSC_EntryExit_mu_dz_dR_detail->SetMarkerColor(kBlue);    CSC_EntryExit_mu_dz_dR_detail->SetMarkerSize(1);  CSC_EntryExit_mu_dz_dR_detail->Draw("PSame");
  CSC_EntryExit_pi_dz_dR_detail->SetMarkerStyle(33); CSC_EntryExit_pi_dz_dR_detail->SetMarkerColor(kGreen);   CSC_EntryExit_pi_dz_dR_detail->SetMarkerSize(1);  CSC_EntryExit_pi_dz_dR_detail->Draw("PSame");
  CSC_EntryExit_ka_dz_dR_detail->SetMarkerStyle(5);  CSC_EntryExit_ka_dz_dR_detail->SetMarkerColor(kOrange);  CSC_EntryExit_ka_dz_dR_detail->SetMarkerSize(1);  CSC_EntryExit_ka_dz_dR_detail->Draw("PSame");
  CSC_EntryExit_p_dz_dR_detail->SetMarkerStyle(26);  CSC_EntryExit_p_dz_dR_detail->SetMarkerColor(kMagenta);  CSC_EntryExit_p_dz_dR_detail->SetMarkerSize(1);   CSC_EntryExit_p_dz_dR_detail->Draw("PSame");
  CSC_EntryExit_n_dz_dR_detail->SetMarkerStyle(32);  CSC_EntryExit_n_dz_dR_detail->SetMarkerColor(kViolet);   CSC_EntryExit_n_dz_dR_detail->SetMarkerSize(1);   CSC_EntryExit_n_dz_dR_detail->Draw("PSame");
  CSC_EntryExit_g_dz_dR_detail->SetMarkerStyle(30);  CSC_EntryExit_g_dz_dR_detail->SetMarkerColor(kCyan);     CSC_EntryExit_g_dz_dR_detail->SetMarkerSize(1);   CSC_EntryExit_g_dz_dR_detail->Draw("PSame");
  CSC_EntryExit_N_dz_dR_detail->SetMarkerStyle(2);   CSC_EntryExit_N_dz_dR_detail->SetMarkerColor(kRed);      CSC_EntryExit_N_dz_dR_detail->SetMarkerSize(1);   CSC_EntryExit_N_dz_dR_detail->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_dz_vs_dR_detail->SetTicks(1,1);
  Canvas_CSC_dz_vs_dR_detail->Write();
  if(pdf_output) {Canvas_CSC_dz_vs_dR_detail->Print(pdfFileName.c_str());}

  Canvas_CSC_deposits_vs_GapLength->cd();
  CSC_EntryExit_el_GapLength_Deposit->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV) (cm)");
  CSC_EntryExit_el_GapLength_Deposit->GetYaxis()->SetTitle("#Delta Gap (Global Coords) (cm)");
  CSC_EntryExit_el_GapLength_Deposit->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_GapLength_Deposit->SetTitleOffset(1.2);
  CSC_EntryExit_el_GapLength_Deposit->SetTitle("SimHit #Delta Gap vs E_{deposit} :: CSC");
  CSC_EntryExit_el_GapLength_Deposit->SetMarkerStyle(7);  CSC_EntryExit_el_GapLength_Deposit->SetMarkerColor(kBlack);   CSC_EntryExit_el_GapLength_Deposit->SetMarkerSize(1);  CSC_EntryExit_el_GapLength_Deposit->Draw("P");
  CSC_EntryExit_mu_GapLength_Deposit->SetMarkerStyle(24); CSC_EntryExit_mu_GapLength_Deposit->SetMarkerColor(kBlue);    CSC_EntryExit_mu_GapLength_Deposit->SetMarkerSize(1);  CSC_EntryExit_mu_GapLength_Deposit->Draw("PSame");
  CSC_EntryExit_pi_GapLength_Deposit->SetMarkerStyle(33); CSC_EntryExit_pi_GapLength_Deposit->SetMarkerColor(kGreen);   CSC_EntryExit_pi_GapLength_Deposit->SetMarkerSize(1);  CSC_EntryExit_pi_GapLength_Deposit->Draw("PSame");
  CSC_EntryExit_ka_GapLength_Deposit->SetMarkerStyle(5);  CSC_EntryExit_ka_GapLength_Deposit->SetMarkerColor(kOrange);  CSC_EntryExit_ka_GapLength_Deposit->SetMarkerSize(1);  CSC_EntryExit_ka_GapLength_Deposit->Draw("PSame");
  CSC_EntryExit_p_GapLength_Deposit->SetMarkerStyle(26);  CSC_EntryExit_p_GapLength_Deposit->SetMarkerColor(kMagenta);  CSC_EntryExit_p_GapLength_Deposit->SetMarkerSize(1);   CSC_EntryExit_p_GapLength_Deposit->Draw("PSame");
  CSC_EntryExit_n_GapLength_Deposit->SetMarkerStyle(32);  CSC_EntryExit_n_GapLength_Deposit->SetMarkerColor(kViolet);   CSC_EntryExit_n_GapLength_Deposit->SetMarkerSize(1);   CSC_EntryExit_n_GapLength_Deposit->Draw("PSame");
  CSC_EntryExit_g_GapLength_Deposit->SetMarkerStyle(30);  CSC_EntryExit_g_GapLength_Deposit->SetMarkerColor(kCyan);     CSC_EntryExit_g_GapLength_Deposit->SetMarkerSize(1);   CSC_EntryExit_g_GapLength_Deposit->Draw("PSame");
  CSC_EntryExit_N_GapLength_Deposit->SetMarkerStyle(2);   CSC_EntryExit_N_GapLength_Deposit->SetMarkerColor(kRed);      CSC_EntryExit_N_GapLength_Deposit->SetMarkerSize(1);   CSC_EntryExit_N_GapLength_Deposit->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
 Canvas_CSC_deposits_vs_GapLength->SetTicks(1,1);
 Canvas_CSC_deposits_vs_GapLength->Write();
 if(pdf_output) {Canvas_CSC_deposits_vs_GapLength->Print(pdfFileName.c_str());}

  Canvas_CSC_deposits_vs_pidR2->cd();
  CSC_EntryExit_el_Deposit_pidR2->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  CSC_EntryExit_el_Deposit_pidR2->GetYaxis()->SetTitle("#pi {(#frac{#Delta R}{2})}^{2} (Global Coords) (cm)");
  CSC_EntryExit_el_Deposit_pidR2->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_Deposit_pidR2->SetTitleOffset(1.2);
  CSC_EntryExit_el_Deposit_pidR2->SetTitle("SimHit #pi #Delta R^{2} vs E_{deposit} :: CSC");
  CSC_EntryExit_el_Deposit_pidR2->SetMarkerStyle(7);  CSC_EntryExit_el_Deposit_pidR2->SetMarkerColor(kBlack);   CSC_EntryExit_el_Deposit_pidR2->SetMarkerSize(1);  CSC_EntryExit_el_Deposit_pidR2->Draw("P");
  CSC_EntryExit_mu_Deposit_pidR2->SetMarkerStyle(24); CSC_EntryExit_mu_Deposit_pidR2->SetMarkerColor(kBlue);    CSC_EntryExit_mu_Deposit_pidR2->SetMarkerSize(1);  CSC_EntryExit_mu_Deposit_pidR2->Draw("PSame");
  CSC_EntryExit_pi_Deposit_pidR2->SetMarkerStyle(33); CSC_EntryExit_pi_Deposit_pidR2->SetMarkerColor(kGreen);   CSC_EntryExit_pi_Deposit_pidR2->SetMarkerSize(1);  CSC_EntryExit_pi_Deposit_pidR2->Draw("PSame");
  CSC_EntryExit_ka_Deposit_pidR2->SetMarkerStyle(5);  CSC_EntryExit_ka_Deposit_pidR2->SetMarkerColor(kOrange);  CSC_EntryExit_ka_Deposit_pidR2->SetMarkerSize(1);  CSC_EntryExit_ka_Deposit_pidR2->Draw("PSame");
  CSC_EntryExit_p_Deposit_pidR2->SetMarkerStyle(26);  CSC_EntryExit_p_Deposit_pidR2->SetMarkerColor(kMagenta);  CSC_EntryExit_p_Deposit_pidR2->SetMarkerSize(1);   CSC_EntryExit_p_Deposit_pidR2->Draw("PSame");
  CSC_EntryExit_n_Deposit_pidR2->SetMarkerStyle(32);  CSC_EntryExit_n_Deposit_pidR2->SetMarkerColor(kViolet);   CSC_EntryExit_n_Deposit_pidR2->SetMarkerSize(1);   CSC_EntryExit_n_Deposit_pidR2->Draw("PSame");
  CSC_EntryExit_g_Deposit_pidR2->SetMarkerStyle(30);  CSC_EntryExit_g_Deposit_pidR2->SetMarkerColor(kCyan);     CSC_EntryExit_g_Deposit_pidR2->SetMarkerSize(1);   CSC_EntryExit_g_Deposit_pidR2->Draw("PSame");
  CSC_EntryExit_N_Deposit_pidR2->SetMarkerStyle(2);   CSC_EntryExit_N_Deposit_pidR2->SetMarkerColor(kRed);      CSC_EntryExit_N_Deposit_pidR2->SetMarkerSize(1);   CSC_EntryExit_N_Deposit_pidR2->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  // latex_left.DrawLatex(lab_max_x,lab_max_y,       lab_time.c_str());
  Canvas_CSC_deposits_vs_pidR2->SetTicks(1,1);
  Canvas_CSC_deposits_vs_pidR2->Write();
  if(pdf_output) {Canvas_CSC_deposits_vs_pidR2->Print(pdfFileName.c_str());}

  Canvas_CSC_deposits_vs_dRdz->cd();
  CSC_EntryExit_el_Deposit_dRdz->GetXaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  CSC_EntryExit_el_Deposit_dRdz->GetYaxis()->SetTitle("#sqrt{#Delta R^{2} + #Delta z^{2}} (Global Coords) (cm)");
  CSC_EntryExit_el_Deposit_dRdz->GetYaxis()->SetTitleOffset(1.15);
  CSC_EntryExit_el_Deposit_dRdz->SetTitleOffset(1.2);
  CSC_EntryExit_el_Deposit_dRdz->SetTitle("SimHit #sqrt{#Delta R^{2} + #Delta z^{2}} vs E_{deposit} :: CSC");
  CSC_EntryExit_el_Deposit_dRdz->SetMarkerStyle(7);  CSC_EntryExit_el_Deposit_dRdz->SetMarkerColor(kBlack);   CSC_EntryExit_el_Deposit_dRdz->SetMarkerSize(1);  CSC_EntryExit_el_Deposit_dRdz->Draw("P");
  CSC_EntryExit_mu_Deposit_dRdz->SetMarkerStyle(24); CSC_EntryExit_mu_Deposit_dRdz->SetMarkerColor(kBlue);    CSC_EntryExit_mu_Deposit_dRdz->SetMarkerSize(1);  CSC_EntryExit_mu_Deposit_dRdz->Draw("PSame");
  CSC_EntryExit_pi_Deposit_dRdz->SetMarkerStyle(33); CSC_EntryExit_pi_Deposit_dRdz->SetMarkerColor(kGreen);   CSC_EntryExit_pi_Deposit_dRdz->SetMarkerSize(1);  CSC_EntryExit_pi_Deposit_dRdz->Draw("PSame");
  CSC_EntryExit_ka_Deposit_dRdz->SetMarkerStyle(5);  CSC_EntryExit_ka_Deposit_dRdz->SetMarkerColor(kOrange);  CSC_EntryExit_ka_Deposit_dRdz->SetMarkerSize(1);  CSC_EntryExit_ka_Deposit_dRdz->Draw("PSame");
  CSC_EntryExit_p_Deposit_dRdz->SetMarkerStyle(26);  CSC_EntryExit_p_Deposit_dRdz->SetMarkerColor(kMagenta);  CSC_EntryExit_p_Deposit_dRdz->SetMarkerSize(1);   CSC_EntryExit_p_Deposit_dRdz->Draw("PSame");
  CSC_EntryExit_n_Deposit_dRdz->SetMarkerStyle(32);  CSC_EntryExit_n_Deposit_dRdz->SetMarkerColor(kViolet);   CSC_EntryExit_n_Deposit_dRdz->SetMarkerSize(1);   CSC_EntryExit_n_Deposit_dRdz->Draw("PSame");
  CSC_EntryExit_g_Deposit_dRdz->SetMarkerStyle(30);  CSC_EntryExit_g_Deposit_dRdz->SetMarkerColor(kCyan);     CSC_EntryExit_g_Deposit_dRdz->SetMarkerSize(1);   CSC_EntryExit_g_Deposit_dRdz->Draw("PSame");
  CSC_EntryExit_N_Deposit_dRdz->SetMarkerStyle(2);   CSC_EntryExit_N_Deposit_dRdz->SetMarkerColor(kRed);      CSC_EntryExit_N_Deposit_dRdz->SetMarkerSize(1);   CSC_EntryExit_N_Deposit_dRdz->Draw("PSame");
  // line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"CSC");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  // latex_right.DrawLatex(lab_250ns_x,lab_250ns_y,  "#font[12]{250 ns}");
  // latex_right.DrawLatex(lab_1ms_x,lab_1ms_y,      "#font[12]{1 ms}");
  // latex_right.DrawLatex(lab_prompt_x,lab_prompt_y,"#font[12]{prompt and decay}");
  // latex_right.DrawLatex(lab_neutr_x,lab_neutr_y,  "#font[12]{neutron background}");
  //   latex_right.DrawLatex(lab_max_x,lab_max_y,    lab_time.c_str());
  Canvas_CSC_deposits_vs_dRdz->SetTicks(1,1);
  Canvas_CSC_deposits_vs_dRdz->Write();
  if(pdf_output) {Canvas_CSC_deposits_vs_dRdz->Print(pdfFileName.c_str());}




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

  
  TDir_Muon_hit_info->cd();
  // ------------------- 
  RPCb_HPC->Write();
  RPCf_HPC->Write();
  CSC_HPC->Write();
  DT_HPC->Write();
  RPCb_el_HPC->Write();
  RPCf_el_HPC->Write();
  CSC_el_HPC->Write();
  DT_el_HPC->Write();
  RPCb_el_HPL->Write();
  RPCf_el_HPL->Write();
  CSC_el_HPL->Write();
  DT_el_HPL->Write();
  RPCb_mu_HPL->Write();
  RPCf_mu_HPL->Write();
  CSC_mu_HPL->Write();
  DT_mu_HPL->Write();
  RPCb_ha_HPL->Write();
  RPCf_ha_HPL->Write();
  CSC_ha_HPL->Write();
  DT_ha_HPL->Write();

  MB1_el_HPL->Write();
  MB2_el_HPL->Write();
  MB3_el_HPL->Write();
  MB4_el_HPL->Write();
  MB1_mu_HPL->Write();
  MB2_mu_HPL->Write();
  MB3_mu_HPL->Write();
  MB4_mu_HPL->Write();
  MB1_ha_HPL->Write();
  MB2_ha_HPL->Write();
  MB3_ha_HPL->Write();
  MB4_ha_HPL->Write();
  ME1_el_HPL->Write();
  ME2_el_HPL->Write();
  ME3_el_HPL->Write();
  ME4_el_HPL->Write();
  ME1_mu_HPL->Write();
  ME2_mu_HPL->Write();
  ME3_mu_HPL->Write();
  ME4_mu_HPL->Write();
  ME1_ha_HPL->Write();
  ME2_ha_HPL->Write();
  ME3_ha_HPL->Write();
  ME4_ha_HPL->Write();
  // ------------------- 
  outputfile->cd();
  
  
  Canvas_RPCb_Layers = new TCanvas("Canvas_RPCb_Layers", "Number of Layers hit :: RPCb", 600, 600);
  Canvas_RPCf_Layers = new TCanvas("Canvas_RPCf_Layers", "Number of Layers hit :: RPCf", 600, 600);
  Canvas_CSC_Layers  = new TCanvas("Canvas_CSC_Layers",  "Number of Layers hit :: CSC",  600, 600);
  Canvas_DT_Layers   = new TCanvas("Canvas_DT_Layers",   "Number of Layers hit :: DT",   600, 600);

  Canvas_RPCb_Layers->cd();
  // RPCb_el_HPL->GetXaxis()->SetTitle("Number of Layers hit in RPCb system"); // original intention
  RPCb_el_HPL->GetXaxis()->SetTitle("Hit Distribution in Layers in RPCb system");
  RPCb_el_HPL->GetYaxis()->SetTitle("Events");
  RPCb_el_HPL->SetTitle("Layers :: RPCb");
  RPCb_el_HPL->SetLineStyle(1);  RPCb_el_HPL->SetLineColor(kBlack);   RPCb_el_HPL->SetLineWidth(1);  RPCb_el_HPL->Draw("H");
  RPCb_mu_HPL->SetLineStyle(1);  RPCb_mu_HPL->SetLineColor(kBlue);    RPCb_mu_HPL->SetLineWidth(1);  RPCb_mu_HPL->Draw("HSame");
  RPCb_ha_HPL->SetLineStyle(1);  RPCb_ha_HPL->SetLineColor(kRed);     RPCb_ha_HPL->SetLineWidth(1);  RPCb_ha_HPL->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPCb_Layers->SetTicks(1,1);
  Canvas_RPCb_Layers->Write();
  if(pdf_output) {Canvas_RPCb_Layers->Print(pdfFileName.c_str());}

  Canvas_RPCf_Layers->cd();
  // RPCf_el_HPL->GetXaxis()->SetTitle("Number of Layers hit in RPCf system"); // original intention
  RPCf_el_HPL->GetXaxis()->SetTitle("Hit Distribution in Layers in RPCf system");
  RPCf_el_HPL->GetYaxis()->SetTitle("Events");
  RPCf_el_HPL->SetTitle("Layers :: RPCf");
  RPCf_el_HPL->SetLineStyle(1);  RPCf_el_HPL->SetLineColor(kBlack);   RPCf_el_HPL->SetLineWidth(1);  RPCf_el_HPL->Draw("H");
  RPCf_mu_HPL->SetLineStyle(1);  RPCf_mu_HPL->SetLineColor(kBlue);    RPCf_mu_HPL->SetLineWidth(1);  RPCf_mu_HPL->Draw("HSame");
  RPCf_ha_HPL->SetLineStyle(1);  RPCf_ha_HPL->SetLineColor(kRed);     RPCf_ha_HPL->SetLineWidth(1);  RPCf_ha_HPL->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPCf_Layers->SetTicks(1,1);
  Canvas_RPCf_Layers->Write();
  if(pdf_output) {Canvas_RPCf_Layers->Print(pdfFileName.c_str());}
  

  Canvas_CSC_Layers->cd();
  // CSC_el_HPL->GetXaxis()->SetTitle("Number of Layers hit in CSC system"); // original intention
  CSC_el_HPL->GetXaxis()->SetTitle("Hit Distribution in Layers in CSC system");
  CSC_el_HPL->GetYaxis()->SetTitle("Events");
  CSC_el_HPL->SetTitle("Layers :: CSC");
  CSC_el_HPL->SetLineStyle(1);  CSC_el_HPL->SetLineColor(kBlack);   CSC_el_HPL->SetLineWidth(1);  CSC_el_HPL->Draw("H");
  CSC_mu_HPL->SetLineStyle(1);  CSC_mu_HPL->SetLineColor(kBlue);    CSC_mu_HPL->SetLineWidth(1);  CSC_mu_HPL->Draw("HSame");
  CSC_ha_HPL->SetLineStyle(1);  CSC_ha_HPL->SetLineColor(kRed);     CSC_ha_HPL->SetLineWidth(1);  CSC_ha_HPL->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_CSC_Layers->SetTicks(1,1);
  Canvas_CSC_Layers->Write();
  if(pdf_output) {Canvas_CSC_Layers->Print(pdfFileName.c_str());}

  Canvas_DT_Layers->cd();
  // DT_el_HPL->GetXaxis()->SetTitle("Number of Layers hit in DT system"); // original intention
  DT_el_HPL->GetXaxis()->SetTitle("Hit Distribution in Layers DT system");
  DT_el_HPL->GetYaxis()->SetTitle("Events");
  DT_el_HPL->SetTitle("Layers :: DT");
  DT_el_HPL->SetLineStyle(1);  DT_el_HPL->SetLineColor(kBlack);   DT_el_HPL->SetLineWidth(1);  DT_el_HPL->Draw("H");
  DT_mu_HPL->SetLineStyle(1);  DT_mu_HPL->SetLineColor(kBlue);    DT_mu_HPL->SetLineWidth(1);  DT_mu_HPL->Draw("HSame");
  DT_ha_HPL->SetLineStyle(1);  DT_ha_HPL->SetLineColor(kRed);     DT_ha_HPL->SetLineWidth(1);  DT_ha_HPL->Draw("HSame");
  l2->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_DT_Layers->SetTicks(1,1);
  Canvas_DT_Layers->Write();
  if(pdf_output) {Canvas_DT_Layers->Print(pdfFileName.c_str());}


  int rpc_barrel_hits = 0;
  int rpc_endcap_hits = 0;
  for(int j=0; j<4; ++j) {
    for(int k=0; k<3; ++k) {
      rpc_barrel_hits += RPC_hits_array[0][j][k];
      rpc_endcap_hits += RPC_hits_array[1][j][k];
      // std::cout<<"RPC_hits_array[region=0][station="<<j+1<<"][ring="<<k<<"] = "  <<RPC_hits_array[0][j][k]<<" ==> rpc_barrel_hits = "<<rpc_barrel_hits<<std::endl; 
      // std::cout<<"RPC_hits_array[region=1][station="<<j+1<<"][ring="<<k+1<<"] = "<<RPC_hits_array[1][j][k]<<" ==> rpc_endcap_hits = "<<rpc_endcap_hits<<std::endl; 
    }
  }

  // Relation Instantaneous Luminosity and PU interactions

  // calculation ::
  // --------------
  // # collisions / second = L x sigma 
  // # collisions / 25 ns = L x sigma / 40*10^6
  // # collisions / 50 ns = L x sigma / 20*10^6
  // ------------------------------------------
  // Luminosity is in units of cm^-2 s^-1
  // Cross section is in units of barn (10^-28 m^-2
  // 1 mb = 10^3 x 10^-28 m^-2 = 10^3 x 10^-32 cm^-2 = 10^-27 cm^-2
  // ------------------------------------------
  // not all bunches in LHC are colliding
  // 3564 bunches
  // design: 2808 (2760?) filled, 2662 colliding ==> fraction = 0.75
  // 2012:   1380 filled, 1380 colliding         ==> fraction = 0.39

  // 2012
  // 8 TeV
  // sigma_inel = 73 mb ( 8 TeV)
  // L = 10^33  (peak :: 8 * 10^33)
  // L x sigma_inel = 73 x 10^6 s^-1 / 20*10^6 = 3.65 
  // L x sigma_inel = 73 x 10^6 s^-1 / 40*10^6 = 7.30
  // 7.30 / 0.39 = 18.7 (collisions happen only in 40% of all possible 25ns timestamps)

  // int pu = 21; // measured ...
  int pu = 24; // normalized to 5 * 10^33 cm^-2 s^-1

  // 2015
  // com = 14 TeV
  // sigma_inel = 80 mb (14 TeV)
  // L = 10^34 1/cm^2 1/s (design)
  // int pu = 25;

  // old calculation
  int entries = 2500;
  // int bx_space = 25;
  // double corr_fact = pow(10,9) * 1.0/bx_space * 1.0/(entries * 1.0/pu);
  double sg_corr_fact = 2;  // simulation has only single gas volumes implemented instead of double gas volumes
  // std::cout<<"Rate Correction factor =  pow(10,9) * 1.0/(entries * 1.0/pu * bx_space) = "<<pow(10,9) * 1.0/(entries * 1.0/pu * bx_space)<<std::endl;
  // std::cout<<"Rate Correction factor because of simulation of Single Gas Layers = "<<sg_corr_fact<<std::endl;
  // std::cout<<"RPC Hits in Barrel = "<<rpc_barrel_hits<<" || RPC Hits in Endcap = "<<rpc_endcap_hits<<std::endl;
  // std::cout<<"RPC Barrel Active area = "<<rpc_barrel_area<<" cm2 || RPC Endcap Active area = "<<rpc_endcap_area<<" cm2"<<std::endl;

  // new calculation
  // ---------------
  // average amount of hits per event = (hits / entries)
  // ( average amount of hits per event ) / area
  // amount of collisions in 1s = 11245 (revolution frequency) * 1380 (colliding bunches) * pu (no of collisions per crossing)
  // => Rate [Hz/cm^2] = ( average amount of hits per event ) / area * amount of collisions in 1s
  // => Rate [Hz/cm^2] = N /entries / area * 11245 * 1380 * 21

  // getLumi(pu, space, com)
  for(int i=0; i<50; ++i) {
    // std::cout<<"Instantaneous Luminosity for "<<i<<" PU interactions at "<<bunchspacing<<"ns bunch spacing and at "<<comenergy<<"TeV is "<<getLumi(i,bunchspacing,comenergy)<<" * 10^34 cm^-2 s^-1 "<<std::endl;
    InstLumi[i] = getLumi(i,bunchspacing,comenergy);
  }
  // getPU(lumi, space, com)
  for(int i=0; i<60; ++i) {
    // std::cout<<"Pile Up for Instantaneous Luminosity "<<i*1.0/100<<" * 10^34 cm^-2 s^-1 at "<<bunchspacing<<"ns bunch spacing and at "<<comenergy<<"TeV is "<<getPU(i*1.0/100,bunchspacing,comenergy)<<" PU"<<std::endl;
  }

  double bunches;
  if(bunchspacing==50) bunches = 1380;
  if(bunchspacing==25) bunches = 2808;
  double corr_fact = 1.0/entries * 11245 * bunches * pu;

  // Rates as function of pu / instantaneous luminosity
  const int max_pu = 33;
  for(int i=0; i<max_pu; ++i) {
    RPC_rates_Summary[0][i] =  (rpc_barrel_hits+rpc_endcap_hits) * 1.0/(rpc_barrel_area+rpc_endcap_area) * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
    RPC_rates_Summary[1][i] =  (rpc_barrel_hits)                 * 1.0/(rpc_barrel_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
    RPC_rates_Summary[2][i] =  (rpc_endcap_hits)                 * 1.0/(rpc_endcap_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
    RPC_uncer_Rate[0][i] =  sqrt(rpc_barrel_hits+rpc_endcap_hits) * 1.0/(rpc_barrel_area+rpc_endcap_area) * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
    RPC_uncer_Rate[1][i] =  sqrt(rpc_barrel_hits)                 * 1.0/(rpc_barrel_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
    RPC_uncer_Rate[2][i] =  sqrt(rpc_endcap_hits)                 * 1.0/(rpc_endcap_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
    RPC_uncer_Lumi[0][i] =  0;
    RPC_uncer_Lumi[1][i] =  0;
    RPC_uncer_Lumi[2][i] =  0;
  }

  gr_RPC_Rates_All    = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[0], RPC_uncer_Lumi[0], RPC_uncer_Rate[0]); 
  gr_RPC_Rates_Barrel = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[1], RPC_uncer_Lumi[1], RPC_uncer_Rate[1]); 
  gr_RPC_Rates_Endcap = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[2], RPC_uncer_Lumi[2], RPC_uncer_Rate[2]); 

  gr_RPC_Rates_All->SetMarkerStyle(21);    gr_RPC_Rates_All->SetMarkerColor(kCyan);     gr_RPC_Rates_All->SetLineColor(kCyan);
  gr_RPC_Rates_Barrel->SetMarkerStyle(21); gr_RPC_Rates_Barrel->SetMarkerColor(kRed);   gr_RPC_Rates_Barrel->SetLineColor(kRed);
  gr_RPC_Rates_Endcap->SetMarkerStyle(21); gr_RPC_Rates_Endcap->SetMarkerColor(kBlack); gr_RPC_Rates_Endcap->SetLineColor(kBlack);

  l3_x1 = 0.20; l3_x2 = 0.45; l3_y1 = 0.65; l3_y2 = 0.85;
  l4_x1 = 0.20; l4_x2 = 0.35; l4_y1 = 0.75; l4_y2 = 0.85;
  TLegend *l3 = new TLegend(l3_x1,l3_y1,l3_x2,l3_y2,NULL,"brNDC");
  l3->SetLineColor(0);    l3->SetLineStyle(0);  l3->SetLineWidth(0);
  l3->SetFillColor(4000); l3->SetBorderSize(0); l3->SetNColumns(1);
  l3->AddEntry(gr_RPC_Rates_Endcap, "Endcap","pl");
  l3->AddEntry(gr_RPC_Rates_All,    "Barrel + Endcap","pl");
  l3->AddEntry(gr_RPC_Rates_Barrel, "Barrel","pl");

  Canvas_RPC_Rates = new TCanvas("Canvas_RPC_Rates", "Rates in RPC System", 600, 600);
  gr_RPC_Rates_Endcap->GetXaxis()->SetTitle("Instantaneous Luminosity #times 10^{34} (cm^{-2}s^{-1})");
  gr_RPC_Rates_Endcap->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  gr_RPC_Rates_Endcap->GetYaxis()->SetTitleOffset(1.30);
  gr_RPC_Rates_Endcap->GetXaxis()->SetRangeUser(0.00,0.68);
  gr_RPC_Rates_Endcap->GetYaxis()->SetRangeUser(0.00,6.00);
  gr_RPC_Rates_Endcap->SetTitle("Rates in RPC System");
  gr_RPC_Rates_Endcap->Draw("PA");
  gr_RPC_Rates_Barrel->Draw("Psame");
  gr_RPC_Rates_All->Draw("Psame");
  l3->Draw("same");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPC_Rates->SetTicks(1,1);
  Canvas_RPC_Rates->Write(); 
  if(pdf_output) {Canvas_RPC_Rates->Print(pdfFileName.c_str());}


  // Barrel + Endcap @ 24 pu
  RPC_hits->SetBinContent(1,  rpc_barrel_hits+rpc_endcap_hits);
  RPC_area->SetBinContent(1,  rpc_barrel_area+rpc_endcap_area);
  if((rpc_barrel_hits > 0 || rpc_endcap_hits > 0) && (rpc_barrel_area > 0.0 || rpc_endcap_area > 0.0)) {
    RPC_rates->SetBinContent(1,  (rpc_barrel_hits+rpc_endcap_hits)*1.0/(rpc_barrel_area+rpc_endcap_area)*corr_fact*sg_corr_fact);
    RPC_rates->SetBinError(  1,  sqrt((rpc_barrel_hits+rpc_endcap_hits))*1.0/(rpc_barrel_area+rpc_endcap_area)*corr_fact*sg_corr_fact);
  }
  // Barrel @ 24 pu 
  RPC_hits->SetBinContent(2,  rpc_barrel_hits);
  RPC_area->SetBinContent(2,  rpc_barrel_area);
  if(rpc_barrel_hits > 0 && rpc_barrel_area > 0.0) {
    RPC_rates->SetBinContent(2,  rpc_barrel_hits*1.0/rpc_barrel_area*corr_fact*sg_corr_fact);
    RPC_rates->SetBinError(  2,  sqrt(rpc_barrel_hits)*1.0/rpc_barrel_area*corr_fact*sg_corr_fact);
  }
  for(int k=0; k<3; ++k) {  // Barrel :: ring == wheel
   for(int j=0; j<4; ++j) {
      RPC_hits->SetBinContent(2+4*(k)+(j+1),RPC_hits_array[0][j][k]);
      RPC_area->SetBinContent(2+4*(k)+(j+1),RPC_area_array[0][j][k]);
      // std::cout<<"Bin "<<1+4*(k)+(j+1)<<" filled with RPC_area_array[region="<<0<<"][station"<<j+1<<"]["<<k<<"] = "<<RPC_area_array[0][j][k]<<" cm2"<<std::endl;
      if(RPC_hits_array[0][j][k] >  0 && RPC_area_array[0][j][k] > 0.0) {
	RPC_rates->SetBinContent(2+4*(k)+(j+1), RPC_hits_array[0][j][k]*1.0/RPC_area_array[0][j][k]*corr_fact*sg_corr_fact);
	RPC_rates->SetBinError(  2+4*(k)+(j+1), sqrt(RPC_hits_array[0][j][k])*1.0/RPC_area_array[0][j][k]*corr_fact*sg_corr_fact);
      }
    }
  }
  // Endcap @ 24 pu
  RPC_hits->SetBinContent(15, rpc_endcap_hits);
  RPC_area->SetBinContent(15, rpc_endcap_area);
  if(rpc_endcap_hits > 0 && rpc_endcap_area > 0.0) {
    RPC_rates->SetBinContent(15,  rpc_endcap_hits*1.0/rpc_endcap_area*corr_fact*sg_corr_fact);
    RPC_rates->SetBinError(  15,  sqrt(rpc_endcap_hits)*1.0/rpc_endcap_area*corr_fact*sg_corr_fact);
  }
  for(int j=0; j<4; ++j) {  // Endcap :: station == disk
    for(int k=0; k<3; ++k) {
      RPC_hits->SetBinContent(15+3*j+(k+1),RPC_hits_array[1][j][k]);
      RPC_area->SetBinContent(15+3*j+(k+1),RPC_area_array[1][j][k]);
      // std::cout<<"Bin "<<14+3*j+(k+1)<<" filled with RPC_area_array[region="<<1<<"][station"<<j+1<<"]["<<k+1<<"] = "<<RPC_area_array[1][j][k]<<" cm2"<<std::endl;
      if(RPC_hits_array[1][j][k] >  0 && RPC_area_array[1][j][k] > 0.0) {
	RPC_rates->SetBinContent(15+3*j+(k+1), RPC_hits_array[1][j][k]*1.0/RPC_area_array[1][j][k]*corr_fact*sg_corr_fact);
	RPC_rates->SetBinError(  15+3*j+(k+1), sqrt(RPC_hits_array[1][j][k])*1.0/RPC_area_array[1][j][k]*corr_fact*sg_corr_fact);
      }
    }
  }

  for(int m=0; m<27; ++m) {
    RPC_hits->GetXaxis()->SetBinLabel( m+1, allcat[m].c_str());
    RPC_area->GetXaxis()->SetBinLabel( m+1, allcat[m].c_str());
    RPC_rates->GetXaxis()->SetBinLabel(m+1, allcat[m].c_str());
  }
  /*
  for(int i=0; i<2;++i) {
    for(int j=0; j<4; ++j) {
      for(int k=0; k<3; ++k) {
	int region = i;
        int station = j+1;
        int ring; if(region==0) ring = k; else ring = k+1;
	std::cout<<"RPC_area_array[region="<<region<<"][station="<<station<<"][ring="<<ring<<"] = "<<RPC_area_array[i][j][k]<<" cm2"<<std::endl;
      }
    }
  }
  for(int i=0; i<2;++i) {
    for(int j=0; j<4; ++j) {
      for(int k=0; k<3; ++k) {
	int region = i;
	int station = j+1;
	int ring; if(region==0) ring = k; else ring = k+1;
	std::cout<<"RPC_hits_array[region="<<region<<"][station="<<station<<"][ring="<<ring<<"] = "<<RPC_hits_array[i][j][k]<<" hits"<<std::endl;
      }
    }
  }
  */
  TDir_Muon_rates->cd();
  // ------------------- 
  RPC_hits->Write();
  RPC_area->Write();
  RPC_rates->Write();

  CSC_ME1_area->Write();
  CSC_ME2_area->Write();
  CSC_ME3_area->Write();
  CSC_ME4_area->Write();

  CSC_ME1_all_hits->Write();
  CSC_ME1_000ns_hits->Write();
  CSC_ME1_250ns_hits->Write();
  CSC_ME1_00ns_hits->Write();
  CSC_ME1_50ns_hits->Write();

  CSC_ME2_all_hits->Write();
  CSC_ME3_all_hits->Write();
  CSC_ME4_all_hits->Write();

  for(int i=0; i<r_CSC; ++i) {
    // all hits
    if(CSC_ME1_all_hits->GetBinContent(i+1) > 0.0 && CSC_ME1_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME1_all_hits->GetBinContent(i+1)*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME1_all_hits->GetBinContent(i+1))*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME1_all_rates->SetBinContent(i+1, rate);
      CSC_ME1_all_rates->SetBinError(i+1, uncert);
    }
    // 0 ns < hits < 250 ns
    if(CSC_ME1_000ns_hits->GetBinContent(i+1) > 0.0 && CSC_ME1_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME1_000ns_hits->GetBinContent(i+1)*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME1_000ns_hits->GetBinContent(i+1))*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME1_000ns_rates->SetBinContent(i+1, rate);
      CSC_ME1_000ns_rates->SetBinError(i+1, uncert);
    }
    // hits > 250 ns        :: slow neutrons (neutron capture)
    if(CSC_ME1_250ns_hits->GetBinContent(i+1) > 0.0 && CSC_ME1_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME1_250ns_hits->GetBinContent(i+1)*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME1_250ns_hits->GetBinContent(i+1))*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME1_250ns_rates->SetBinContent(i+1, rate);
      CSC_ME1_250ns_rates->SetBinError(i+1, uncert);
    }
    // 0 ns < hits < 50 ns   :: prompt & punch through
    if(CSC_ME1_00ns_hits->GetBinContent(i+1) > 0.0 && CSC_ME1_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME1_00ns_hits->GetBinContent(i+1)*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME1_00ns_hits->GetBinContent(i+1))*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME1_00ns_rates->SetBinContent(i+1, rate);
      CSC_ME1_00ns_rates->SetBinError(i+1, uncert);
    }
    // 50 ns < hits < 250 ns :: fast neutrons (proton knock out)
    if(CSC_ME1_50ns_hits->GetBinContent(i+1) > 0.0 && CSC_ME1_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME1_50ns_hits->GetBinContent(i+1)*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME1_50ns_hits->GetBinContent(i+1))*1.0/CSC_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME1_50ns_rates->SetBinContent(i+1, rate);
      CSC_ME1_50ns_rates->SetBinError(i+1, uncert);
    }
  }
  CSC_ME1_all_rates->Write();
  CSC_ME1_000ns_rates->Write();
  CSC_ME1_250ns_rates->Write();
  CSC_ME1_00ns_rates->Write();
  CSC_ME1_50ns_rates->Write();

  for(int i=0; i<r_CSC; ++i) {
    // all hits
    if(CSC_ME2_all_hits->GetBinContent(i+1) > 0.0 && CSC_ME2_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME2_all_hits->GetBinContent(i+1)*1.0/CSC_ME2_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME2_all_hits->GetBinContent(i+1))*1.0/CSC_ME2_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME2_all_rates->SetBinContent(i+1, rate);
      CSC_ME2_all_rates->SetBinError(i+1, uncert);
    }
    if(CSC_ME3_all_hits->GetBinContent(i+1) > 0.0 && CSC_ME3_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME3_all_hits->GetBinContent(i+1)*1.0/CSC_ME3_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME3_all_hits->GetBinContent(i+1))*1.0/CSC_ME3_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME3_all_rates->SetBinContent(i+1, rate);
      CSC_ME3_all_rates->SetBinError(i+1, uncert);
    }
    if(CSC_ME4_all_hits->GetBinContent(i+1) > 0.0 && CSC_ME4_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_ME4_all_hits->GetBinContent(i+1)*1.0/CSC_ME4_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_ME4_all_hits->GetBinContent(i+1))*1.0/CSC_ME4_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_ME4_all_rates->SetBinContent(i+1, rate);
      CSC_ME4_all_rates->SetBinError(i+1, uncert);
    }
  }
  CSC_ME2_all_rates->Write();
  CSC_ME3_all_rates->Write();
  CSC_ME4_all_rates->Write();

  CSC_Geom_ME1_area->Write();
  CSC_Geom_ME2_area->Write();
  CSC_Geom_ME3_area->Write();
  CSC_Geom_ME4_area->Write();

  CSC_Geom_ME1_all_hits->Write();
  CSC_Geom_ME2_all_hits->Write();
  CSC_Geom_ME3_all_hits->Write();
  CSC_Geom_ME4_all_hits->Write();

  for(int i=0; i<51; ++i) {
    if(CSC_Geom_ME1_all_hits->GetBinContent(i+1) > 0.0 && CSC_Geom_ME1_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_Geom_ME1_all_hits->GetBinContent(i+1)*1.0/CSC_Geom_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_Geom_ME1_all_hits->GetBinContent(i+1))*1.0/CSC_Geom_ME1_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_Geom_ME1_all_rates->SetBinContent(i+1, rate);
      CSC_Geom_ME1_all_rates->SetBinError(i+1, uncert);
    }
  }
  for(int i=0; i<53; ++i) {
    if(CSC_Geom_ME2_all_hits->GetBinContent(i+1) > 0.0 && CSC_Geom_ME2_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_Geom_ME2_all_hits->GetBinContent(i+1)*1.0/CSC_Geom_ME2_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_Geom_ME2_all_hits->GetBinContent(i+1))*1.0/CSC_Geom_ME2_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_Geom_ME2_all_rates->SetBinContent(i+1, rate);
      CSC_Geom_ME2_all_rates->SetBinError(i+1, uncert);
    }
    if(CSC_Geom_ME3_all_hits->GetBinContent(i+1) > 0.0 && CSC_Geom_ME3_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_Geom_ME3_all_hits->GetBinContent(i+1)*1.0/CSC_Geom_ME3_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_Geom_ME3_all_hits->GetBinContent(i+1))*1.0/CSC_Geom_ME3_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_Geom_ME3_all_rates->SetBinContent(i+1, rate);
      CSC_Geom_ME3_all_rates->SetBinError(i+1, uncert);
    }
    if(CSC_Geom_ME4_all_hits->GetBinContent(i+1) > 0.0 && CSC_Geom_ME4_area->GetBinContent(i+1) > 0.0) {
      double rate = CSC_Geom_ME4_all_hits->GetBinContent(i+1)*1.0/CSC_Geom_ME4_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      double uncert = sqrt(CSC_Geom_ME4_all_hits->GetBinContent(i+1))*1.0/CSC_Geom_ME4_area->GetBinContent(i+1)*corr_fact/n_layers_csc;
      CSC_Geom_ME4_all_rates->SetBinContent(i+1, rate);
      CSC_Geom_ME4_all_rates->SetBinError(i+1, uncert);
    }
  }

  CSC_Geom_ME1_all_rates->Write();
  CSC_Geom_ME2_all_rates->Write();
  CSC_Geom_ME3_all_rates->Write();
  CSC_Geom_ME4_all_rates->Write();

  // ------------------- 
  outputfile->cd();

  for(int i=0; i<max_pu; ++i) {
    for(int j=0; j<4; ++j) {  
      for(int k=0; k<3; ++k) {                 /* RPC_hits_array[1][j][k] */                /* RPC_area_array[1][j][k] */
	if(RPC_hits->GetBinContent(15+3*j+(k+1)) > 0.0 || RPC_area->GetBinContent(15+3*j+(k+1))) {
	  RPC_rates_Summary[2+3*j+(k+1)][i] =       RPC_hits->GetBinContent(15+3*j+(k+1))  * 1.0/RPC_area->GetBinContent(15+3*j+(k+1)) * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
	  RPC_uncer_Rate[2+3*j+(k+1)][i]    =  sqrt(RPC_hits->GetBinContent(15+3*j+(k+1))) * 1.0/RPC_area->GetBinContent(15+3*j+(k+1)) * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
	  RPC_uncer_Lumi[2+3*j+(k+1)][i]    =  0;
	}
	else {
	  RPC_rates_Summary[2+3*j+(k+1)][i] = 0;
	  RPC_uncer_Rate[2+3*j+(k+1)][i]    = 0;
	  RPC_uncer_Lumi[2+3*j+(k+1)][i]    = 0;
	}
      }
    }
  }
  // RE11 --> RPC_rates_Summary[3]
  gr_RPC_Rates_RE12 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[4], RPC_uncer_Lumi[4], RPC_uncer_Rate[4]); 
  gr_RPC_Rates_RE13 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[5], RPC_uncer_Lumi[5], RPC_uncer_Rate[5]); 
  gr_RPC_Rates_RE12->SetMarkerStyle(21);    gr_RPC_Rates_RE12->SetMarkerColor(kRed);      gr_RPC_Rates_RE12->SetLineColor(kRed);
  gr_RPC_Rates_RE13->SetMarkerStyle(21);    gr_RPC_Rates_RE13->SetMarkerColor(kCyan);     gr_RPC_Rates_RE13->SetLineColor(kCyan);
  // RE21 --> RPC_rates_Summary[6]
  gr_RPC_Rates_RE22 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[7], RPC_uncer_Lumi[7], RPC_uncer_Rate[7]); 
  gr_RPC_Rates_RE23 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[8], RPC_uncer_Lumi[8], RPC_uncer_Rate[8]); 
  gr_RPC_Rates_RE22->SetMarkerStyle(21);    gr_RPC_Rates_RE22->SetMarkerColor(kRed);      gr_RPC_Rates_RE22->SetLineColor(kRed);
  gr_RPC_Rates_RE23->SetMarkerStyle(21);    gr_RPC_Rates_RE23->SetMarkerColor(kCyan);     gr_RPC_Rates_RE23->SetLineColor(kCyan);
  // RE31 --> RPC_rates_Summary[9]
  gr_RPC_Rates_RE32 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[10], RPC_uncer_Lumi[10], RPC_uncer_Rate[10]); 
  gr_RPC_Rates_RE33 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[11], RPC_uncer_Lumi[11], RPC_uncer_Rate[11]); 
  gr_RPC_Rates_RE32->SetMarkerStyle(21);    gr_RPC_Rates_RE32->SetMarkerColor(kRed);      gr_RPC_Rates_RE32->SetLineColor(kRed);
  gr_RPC_Rates_RE33->SetMarkerStyle(21);    gr_RPC_Rates_RE33->SetMarkerColor(kCyan);     gr_RPC_Rates_RE33->SetLineColor(kCyan);
  // RE41 --> RPC_rates_Summary[12]
  gr_RPC_Rates_RE42 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[13], RPC_uncer_Lumi[13], RPC_uncer_Rate[13]); 
  gr_RPC_Rates_RE43 = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[14], RPC_uncer_Lumi[14], RPC_uncer_Rate[14]); 
  gr_RPC_Rates_RE42->SetMarkerStyle(21);    gr_RPC_Rates_RE42->SetMarkerColor(kRed);      gr_RPC_Rates_RE42->SetLineColor(kRed);
  gr_RPC_Rates_RE43->SetMarkerStyle(21);    gr_RPC_Rates_RE43->SetMarkerColor(kCyan);     gr_RPC_Rates_RE43->SetLineColor(kCyan);

  TLegend *l3a = new TLegend(l4_x1,l4_y1,l4_x2,l4_y2,NULL,"brNDC");
  l3a->SetLineColor(0);    l3a->SetLineStyle(0);  l3a->SetLineWidth(0);
  l3a->SetFillColor(4000); l3a->SetBorderSize(0); l3a->SetNColumns(1);
  l3a->AddEntry(gr_RPC_Rates_RE12, "RE1/2","pl");
  l3a->AddEntry(gr_RPC_Rates_RE13, "RE1/3","pl");
  TLegend *l3b = new TLegend(l4_x1,l4_y1,l4_x2,l4_y2,NULL,"brNDC");
  l3b->SetLineColor(0);    l3b->SetLineStyle(0);  l3b->SetLineWidth(0);
  l3b->SetFillColor(4000); l3b->SetBorderSize(0); l3b->SetNColumns(1);
  l3b->AddEntry(gr_RPC_Rates_RE12, "RE2/2","pl");
  l3b->AddEntry(gr_RPC_Rates_RE13, "RE2/3","pl");
  TLegend *l3c = new TLegend(l4_x1,l4_y1,l4_x2,l4_y2,NULL,"brNDC");
  l3c->SetLineColor(0);    l3c->SetLineStyle(0);  l3c->SetLineWidth(0);
  l3c->SetFillColor(4000); l3c->SetBorderSize(0); l3c->SetNColumns(1);
  l3c->AddEntry(gr_RPC_Rates_RE12, "RE3/2","pl");
  l3c->AddEntry(gr_RPC_Rates_RE13, "RE3/3","pl");
  TLegend *l3d = new TLegend(l4_x1,l4_y1,l4_x2,l4_y2,NULL,"brNDC");
  l3d->SetLineColor(0);    l3d->SetLineStyle(0);  l3d->SetLineWidth(0);
  l3d->SetFillColor(4000); l3d->SetBorderSize(0); l3d->SetNColumns(1);
  l3d->AddEntry(gr_RPC_Rates_RE12, "RE4/2","pl");
  l3d->AddEntry(gr_RPC_Rates_RE13, "RE4/3","pl");

  Canvas_RPC_Rates_RE1 = new TCanvas("Canvas_RPC_Rates_RE1", "Rates in RPC System :: First Endcap Station", 600, 600);
  gr_RPC_Rates_RE12->GetXaxis()->SetTitle("Instantaneous Luminosity #times 10^{34} (cm^{-2}s^{-1})");
  gr_RPC_Rates_RE12->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  gr_RPC_Rates_RE12->GetYaxis()->SetTitleOffset(1.30);
  gr_RPC_Rates_RE12->GetXaxis()->SetRangeUser(0.00,0.68);
  gr_RPC_Rates_RE12->GetYaxis()->SetRangeUser(0.00,4.00);
  gr_RPC_Rates_RE12->SetTitle("Rates in RPC System");
  gr_RPC_Rates_RE12->Draw("PA");
  gr_RPC_Rates_RE13->Draw("Psame");
  l3a->Draw("same");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPC_Rates_RE1->SetTicks(1,1);
  Canvas_RPC_Rates_RE1->Write(); 
  if(pdf_output) {Canvas_RPC_Rates_RE1->Print(pdfFileName.c_str());}

  Canvas_RPC_Rates_RE2 = new TCanvas("Canvas_RPC_Rates_RE2", "Rates in RPC System :: First Endcap Station", 600, 600);
  gr_RPC_Rates_RE22->GetXaxis()->SetTitle("Instantaneous Luminosity #times 10^{34} (cm^{-2}s^{-1})");
  gr_RPC_Rates_RE22->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  gr_RPC_Rates_RE22->GetYaxis()->SetTitleOffset(1.30);
  gr_RPC_Rates_RE22->GetXaxis()->SetRangeUser(0.00,0.68);
  gr_RPC_Rates_RE22->GetYaxis()->SetRangeUser(0.00,15.0);
  gr_RPC_Rates_RE22->SetTitle("Rates in RPC System");
  gr_RPC_Rates_RE22->Draw("PA");
  gr_RPC_Rates_RE23->Draw("Psame");
  l3b->Draw("same");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPC_Rates_RE2->SetTicks(1,1);
  Canvas_RPC_Rates_RE2->Write(); 
  if(pdf_output) {Canvas_RPC_Rates_RE2->Print(pdfFileName.c_str());}

  Canvas_RPC_Rates_RE3 = new TCanvas("Canvas_RPC_Rates_RE3", "Rates in RPC System :: First Endcap Station", 600, 600);
  gr_RPC_Rates_RE32->GetXaxis()->SetTitle("Instantaneous Luminosity #times 10^{34} (cm^{-2}s^{-1})");
  gr_RPC_Rates_RE32->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  gr_RPC_Rates_RE32->GetYaxis()->SetTitleOffset(1.30);
  gr_RPC_Rates_RE32->GetXaxis()->SetRangeUser(0.00,0.68);
  gr_RPC_Rates_RE32->GetYaxis()->SetRangeUser(0.00,11.0);
  gr_RPC_Rates_RE32->SetTitle("Rates in RPC System");
  gr_RPC_Rates_RE32->Draw("PA");
  gr_RPC_Rates_RE33->Draw("Psame");
  l3c->Draw("same");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPC_Rates_RE3->SetTicks(1,1);
  Canvas_RPC_Rates_RE3->Write(); 
  if(pdf_output) {Canvas_RPC_Rates_RE3->Print(pdfFileName.c_str());}

  Canvas_RPC_Rates_RE4 = new TCanvas("Canvas_RPC_Rates_RE4", "Rates in RPC System :: First Endcap Station", 600, 600);
  gr_RPC_Rates_RE42->GetXaxis()->SetTitle("Instantaneous Luminosity #times 10^{34} (cm^{-2}s^{-1})");
  gr_RPC_Rates_RE42->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  gr_RPC_Rates_RE42->GetYaxis()->SetTitleOffset(1.30);
  gr_RPC_Rates_RE42->GetXaxis()->SetRangeUser(0.00,0.68);
  gr_RPC_Rates_RE42->GetYaxis()->SetRangeUser(0.00,15.0);
  gr_RPC_Rates_RE42->SetTitle("Rates in RPC System");
  gr_RPC_Rates_RE42->Draw("PA");
  gr_RPC_Rates_RE43->Draw("Psame");
  l3d->Draw("same");
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  Canvas_RPC_Rates_RE4->SetTicks(1,1);
  Canvas_RPC_Rates_RE4->Write(); 
  if(pdf_output) {Canvas_RPC_Rates_RE4->Print(pdfFileName.c_str());}

  /*
    TDir_Muon_rates->cd();
  // ------------------------ 
  gr_RPC_Rates_RE12->Write();
  gr_RPC_Rates_RE13->Write();
  gr_RPC_Rates_RE22->Write();
  gr_RPC_Rates_RE23->Write();
  gr_RPC_Rates_RE32->Write();
  gr_RPC_Rates_RE33->Write();
  gr_RPC_Rates_RE42->Write();
  gr_RPC_Rates_RE43->Write();
  // ------------------------
  */
  if(pdf_output) {
  // last plot for PDF File :: print empty Dummy
    pdfFileName = pdfFileNameBase + ".pdf]";
    Dummy->Print(pdfFileName.c_str());
  }
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyNeutronSimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(tech_debug) std::cout<<"[MyNeutronSimHitAnalyzer :: Analyze]"<<std::endl;

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
    double ms_time     = (*iHit).timeOfFlight()*1.0/1000000; // time ns->ms

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
      int station = rollId.station();
      int ring = rollId.ring();

      bool isME11 = false, isME11odd = false, isME11even = false;
      if(rollId.station() == 1 && (rollId.ring() == 1 || rollId.ring() == 4) ) {
	isME11 = true;
	if(rollId.chamber()%2==1) isME11odd = true; else isME11even = true;
      }

      CSC_MEXX_hits_tof->Fill(ms_time);
      if(station==1) {
	if(ring==1)      CSC_ME11_hits_tof->Fill(ms_time); 
	else if(ring==2) CSC_ME12_hits_tof->Fill(ms_time);
	else if(ring==3) CSC_ME13_hits_tof->Fill(ms_time);
	else if(ring==4) CSC_ME11_hits_tof->Fill(ms_time);
	else {}
      }
      if(station==2) {
	if(ring==1) CSC_ME21_hits_tof->Fill(ms_time);
	else        CSC_ME22_hits_tof->Fill(ms_time);
      }
      if(station==3) {
	if(ring==1) CSC_ME31_hits_tof->Fill(ms_time);
	else        CSC_ME32_hits_tof->Fill(ms_time);
      }
      if(station==4) {
	if(ring==1) CSC_ME41_hits_tof->Fill(ms_time);
	else        CSC_ME42_hits_tof->Fill(ms_time);
      }
      if(ring==1 || ring==4) { CSC_MEX1_hits_tof->Fill(ms_time); }
      if(ring==2 || ring==3) { CSC_MEX2_hits_tof->Fill(ms_time); }

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
      int station = wireId.station();
      int wheel = wireId.wheel();

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

      DT_MBX_all_hits_tof->Fill(ms_time);
      if(station==1) {
	DT_MB1_all_hits_tof->Fill(ms_time);
	if(abs(wheel)==2)      DT_MB1_W2_hits_tof->Fill(ms_time);
	else if(abs(wheel)==1) DT_MB1_W1_hits_tof->Fill(ms_time);
	else                   DT_MB1_W0_hits_tof->Fill(ms_time);
      }
      else if(station==2) { DT_MB2_all_hits_tof->Fill(ms_time); }
      else if(station==3) { DT_MB3_all_hits_tof->Fill(ms_time); }
      else if(station==4) { DT_MB4_all_hits_tof->Fill(ms_time); }
      else {}
      
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
MyNeutronSimHitAnalyzer::getLumi(int pu, int space, int com)
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

double MyNeutronSimHitAnalyzer::getPU(double lumi, int space, int com)
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
MyNeutronSimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyNeutronSimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyNeutronSimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  if(tech_debug) std::cout<<"[MyNeutronSimHitAnalyzer :: BeginRun]"<<std::endl;

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
MyNeutronSimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{

}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyNeutronSimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyNeutronSimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyNeutronSimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyNeutronSimHitAnalyzer);

//  LocalWords:  str
