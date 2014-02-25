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
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
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
#include "FastSimulation/Tracking/test/FastTrackAnalyzer.h"
// DetIds
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include "DataFormats/MuonDetId/interface/DTWireId.h"
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

class MyNeutronSimHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyNeutronSimHitAnalyzer(const edm::ParameterSet&);
      ~MyNeutronSimHitAnalyzer();
  edm::ESHandle <RPCGeometry> rpcGeom;
  edm::ESHandle <CSCGeometry> cscGeom;
  edm::ESHandle <DTGeometry>  dtGeom;

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
  std::string pdfFileNameBase, pdfFileName, rootFileName;
  bool debug;
  TFile * outputfile;

  TH2F * RPCb_el_hits, * RPCb_mu_hits, * RPCb_pi_hits, * RPCb_ka_hits, * RPCb_p_hits, * RPCb_n_hits, * RPCb_g_hits, * RPCb_N_hits;
  TH2F * RPCf_el_hits, * RPCf_mu_hits, * RPCf_pi_hits, * RPCf_ka_hits, * RPCf_p_hits, * RPCf_n_hits, * RPCf_g_hits, * RPCf_N_hits;
  TH2F * CSC_el_hits,  * CSC_mu_hits,  * CSC_pi_hits,  * CSC_ka_hits,  * CSC_p_hits,  * CSC_n_hits,  * CSC_g_hits,  * CSC_N_hits;
  TH2F * DT_el_hits,   * DT_mu_hits,   * DT_pi_hits,   * DT_ka_hits,   * DT_p_hits,   * DT_n_hits,   * DT_g_hits,   * DT_N_hits;

  TH2F * RPCb_el_deposits, * RPCb_mu_deposits, * RPCb_pi_deposits, * RPCb_ka_deposits, * RPCb_p_deposits, * RPCb_n_deposits, * RPCb_g_deposits, * RPCb_N_deposits;
  TH2F * RPCf_el_deposits, * RPCf_mu_deposits, * RPCf_pi_deposits, * RPCf_ka_deposits, * RPCf_p_deposits, * RPCf_n_deposits, * RPCf_g_deposits, * RPCf_N_deposits;
  TH2F * CSC_el_deposits,  * CSC_mu_deposits,  * CSC_pi_deposits,  * CSC_ka_deposits,  * CSC_p_deposits,  * CSC_n_deposits,  * CSC_g_deposits,  * CSC_N_deposits;
  TH2F * DT_el_deposits,   * DT_mu_deposits,   * DT_pi_deposits,   * DT_ka_deposits,   * DT_p_deposits,   * DT_n_deposits,   * DT_g_deposits,   * DT_N_deposits;

  TH1F * RPCb_hits_tof, * RPCf_hits_tof, * CSC_hits_tof, * DT_hits_tof;
  TH1F * RPCb_hits_eta, * RPCf_hits_eta, * CSC_hits_eta, * DT_hits_eta;
  TH1F * RPCb_hits_phi, * RPCf_hits_phi, * CSC_hits_phi, * DT_hits_phi;

  TCanvas * Canvas_RPCb_hits, * Canvas_RPCf_hits, * Canvas_CSC_hits, * Canvas_DT_hits;
  TCanvas * Canvas_RPCb_deposits, * Canvas_RPCf_deposits, * Canvas_CSC_deposits, * Canvas_DT_deposits;

  TH2F    * RPCb_XY, * RPCb_RZ, * RPCf_XY, * RPCf_RZ, *CSC_XY, * CSC_RZ, * DT_XY, * DT_RZ, * Muon_RZ;
  TCanvas * Canvas_RPC_XY, * Canvas_RPC_RZ, * Canvas_CSC_XY, * Canvas_CSC_RZ, * Canvas_DT_XY, * Canvas_DT_RZ, * Canvas_Muon_RZ;

  TH2F    * RPCb_250ns_XY, * RPCb_250ns_RZ, * RPCf_250ns_XY, * RPCf_250ns_RZ, *CSC_250ns_XY, * CSC_250ns_RZ, * DT_250ns_XY, * DT_250ns_RZ, * Muon_250ns_RZ;
  TCanvas * Canvas_RPC_250ns_XY, * Canvas_RPC_250ns_RZ, * Canvas_CSC_250ns_XY, * Canvas_CSC_250ns_RZ, * Canvas_DT_250ns_XY, * Canvas_DT_250ns_RZ, * Canvas_Muon_250ns_RZ;

  TH2F * SimVertices_RZ, * SimVertices_Muon_RZ;
  TCanvas * Canvas_SimVertices_RZ, * Canvas_SimVertices_Muon_RZ;

  TH1F * PrimVertices_Z, * PrimVertices_R;
  TH1F * RPC_hits, * RPC_area, * RPC_rates;
  int RPC_hits_array[2][4][3];
  double RPC_area_array[2][4][3], RPC_rates_array[2][4][3]; 
  double rpc_barrel_area, rpc_endcap_area = 0.0;
};

//
// constants, enums and typedefs
//

int n_tof = 700,  n1_tof = 1,  n2_tof = 8;

int m_tof = 70,   m1_tof = 1,  m2_tof = 8;
int m_eta = 50;   double m1_eta =  0.0,  m2_eta = 2.5;
int m_phi = 63;   double m1_phi = -3.15, m2_phi = 3.15;

int n_D   = 900,  n1_D   = -3, n2_D   = 6;
int n_E   = 900,  n1_E   = -4, n2_E   = 5;


int n_xy_x =  800; double n_xy_x1 = -800;  double n_xy_x2 = +800;
int n_xy_y =  800; double n_xy_y1 = -800;  double n_xy_y2 = +800;
int n_zr_z =  550; double n_zr_z1 =    0;  double n_zr_z2 = 1100;
int n_zr_r =  400; double n_zr_r1 =    0;  double n_zr_r2 =  800;

int n_zrc_z = 1350; double n_zrc_z1 = 0; double n_zrc_z2 = 2700;
int n_zrc_r =  650; double n_zrc_r1 = 0; double n_zrc_r2 = 1300;

int n_pv_z = 50, n1_pv_z = -25, n2_pv_z = 25;
int n_pv_r = 10, n1_pv_r = 0,   n2_pv_r = 2;
int n_cat  = 26; double n1_cat  = 0.5, n2_cat  = 26.5;

std::string cat[26] = {"Barrel", "W0_RB1", "W0_RB2", "W0_RB3", "W0_RB4", "W1_RB1", "W1_RB2", "W1_RB3", "W1_RB4", "W2_RB1", "W2_RB2", "W2_RB3", "W2_RB4", 
                  "Endcap", "RE11", "RE12", "RE13", "RE21", "RE22", "RE23", "RE31", "RE32", "RE33", "RE41", "RE42", "RE43"};
//
// static data member definitions
//

//
// constructors and destructor
//
MyNeutronSimHitAnalyzer::MyNeutronSimHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  pdfFileNameBase = iConfig.getUntrackedParameter<std::string>("PdfFileNameBase");
  rootFileName    = iConfig.getUntrackedParameter<std::string>("RootFileName");
  debug           = iConfig.getUntrackedParameter<bool>("Debug");

  outputfile      = new TFile(rootFileName.c_str(), "RECREATE" );

  // Simhit Time vs Ekin
  RPCb_el_hits = new TH2F("RPCb_el_hits", "Simhit time vs E_{kin} :: RPCb :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_mu_hits = new TH2F("RPCb_mu_hits", "Simhit time vs E_{kin} :: RPCb :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_pi_hits = new TH2F("RPCb_pi_hits", "Simhit time vs E_{kin} :: RPCb :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_ka_hits = new TH2F("RPCb_ka_hits", "Simhit time vs E_{kin} :: RPCb :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_p_hits  = new TH2F("RPCb_p_hits",  "Simhit time vs E_{kin} :: RPCb :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_n_hits  = new TH2F("RPCb_n_hits",  "Simhit time vs E_{kin} :: RPCb :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_g_hits  = new TH2F("RPCb_g_hits",  "Simhit time vs E_{kin} :: RPCb :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCb_N_hits  = new TH2F("RPCb_N_hits",  "Simhit time vs E_{kin} :: RPCb :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  RPCf_el_hits = new TH2F("RPCf_el_hits", "Simhit time vs E_{kin} :: RPCf :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_mu_hits = new TH2F("RPCf_mu_hits", "Simhit time vs E_{kin} :: RPCf :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_pi_hits = new TH2F("RPCf_pi_hits", "Simhit time vs E_{kin} :: RPCf :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_ka_hits = new TH2F("RPCf_ka_hits", "Simhit time vs E_{kin} :: RPCf :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_p_hits  = new TH2F("RPCf_p_hits",  "Simhit time vs E_{kin} :: RPCf :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_n_hits  = new TH2F("RPCf_n_hits",  "Simhit time vs E_{kin} :: RPCf :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_g_hits  = new TH2F("RPCf_g_hits",  "Simhit time vs E_{kin} :: RPCf :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  RPCf_N_hits  = new TH2F("RPCf_N_hits",  "Simhit time vs E_{kin} :: RPCf :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  CSC_el_hits = new TH2F("CSC_el_hits", "Simhit time vs E_{kin} :: CSC :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_mu_hits = new TH2F("CSC_mu_hits", "Simhit time vs E_{kin} :: CSC :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_pi_hits = new TH2F("CSC_pi_hits", "Simhit time vs E_{kin} :: CSC :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_ka_hits = new TH2F("CSC_ka_hits", "Simhit time vs E_{kin} :: CSC :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_p_hits  = new TH2F("CSC_p_hits",  "Simhit time vs E_{kin} :: CSC :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_n_hits  = new TH2F("CSC_n_hits",  "Simhit time vs E_{kin} :: CSC :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_g_hits  = new TH2F("CSC_g_hits",  "Simhit time vs E_{kin} :: CSC :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  CSC_N_hits  = new TH2F("CSC_N_hits",  "Simhit time vs E_{kin} :: CSC :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  DT_el_hits = new TH2F("DT_el_hits", "Simhit time vs E_{kin} :: DT :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_mu_hits = new TH2F("DT_mu_hits", "Simhit time vs E_{kin} :: DT :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_pi_hits = new TH2F("DT_pi_hits", "Simhit time vs E_{kin} :: DT :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_ka_hits = new TH2F("DT_ka_hits", "Simhit time vs E_{kin} :: DT :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_p_hits  = new TH2F("DT_p_hits",  "Simhit time vs E_{kin} :: DT :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);\
  DT_n_hits  = new TH2F("DT_n_hits",  "Simhit time vs E_{kin} :: DT :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_g_hits  = new TH2F("DT_g_hits",  "Simhit time vs E_{kin} :: DT :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  DT_N_hits  = new TH2F("DT_N_hits",  "Simhit time vs E_{kin} :: DT :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  // Simhit Time vs E deposit
  RPCb_el_deposits = new TH2F("RPCb_el_deposits", "Simhit time vs E_{deposit} :: RPCb :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_mu_deposits = new TH2F("RPCb_mu_deposits", "Simhit time vs E_{deposit} :: RPCb :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_pi_deposits = new TH2F("RPCb_pi_deposits", "Simhit time vs E_{deposit} :: RPCb :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_ka_deposits = new TH2F("RPCb_ka_deposits", "Simhit time vs E_{deposit} :: RPCb :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_p_deposits  = new TH2F("RPCb_p_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_n_deposits  = new TH2F("RPCb_n_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_g_deposits  = new TH2F("RPCb_g_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCb_N_deposits  = new TH2F("RPCb_N_deposits",  "Simhit time vs E_{deposit} :: RPCb :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  RPCf_el_deposits = new TH2F("RPCf_el_deposits", "Simhit time vs E_{deposit} :: RPCf :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_mu_deposits = new TH2F("RPCf_mu_deposits", "Simhit time vs E_{deposit} :: RPCf :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_pi_deposits = new TH2F("RPCf_pi_deposits", "Simhit time vs E_{deposit} :: RPCf :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_ka_deposits = new TH2F("RPCf_ka_deposits", "Simhit time vs E_{deposit} :: RPCf :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_p_deposits  = new TH2F("RPCf_p_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_n_deposits  = new TH2F("RPCf_n_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_g_deposits  = new TH2F("RPCf_g_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  RPCf_N_deposits  = new TH2F("RPCf_N_deposits",  "Simhit time vs E_{deposit} :: RPCf :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  CSC_el_deposits = new TH2F("CSC_el_deposits", "Simhit time vs E_{deposit} :: CSC :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_mu_deposits = new TH2F("CSC_mu_deposits", "Simhit time vs E_{deposit} :: CSC :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_pi_deposits = new TH2F("CSC_pi_deposits", "Simhit time vs E_{deposit} :: CSC :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_ka_deposits = new TH2F("CSC_ka_deposits", "Simhit time vs E_{deposit} :: CSC :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_p_deposits  = new TH2F("CSC_p_deposits",  "Simhit time vs E_{deposit} :: CSC :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_n_deposits  = new TH2F("CSC_n_deposits",  "Simhit time vs E_{deposit} :: CSC :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_g_deposits  = new TH2F("CSC_g_deposits",  "Simhit time vs E_{deposit} :: CSC :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  CSC_N_deposits  = new TH2F("CSC_N_deposits",  "Simhit time vs E_{deposit} :: CSC :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  DT_el_deposits = new TH2F("DT_el_deposits", "Simhit time vs E_{deposit} :: DT :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_mu_deposits = new TH2F("DT_mu_deposits", "Simhit time vs E_{deposit} :: DT :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_pi_deposits = new TH2F("DT_pi_deposits", "Simhit time vs E_{deposit} :: DT :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_ka_deposits = new TH2F("DT_ka_deposits", "Simhit time vs E_{deposit} :: DT :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_p_deposits  = new TH2F("DT_p_deposits",  "Simhit time vs E_{deposit} :: DT :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_n_deposits  = new TH2F("DT_n_deposits",  "Simhit time vs E_{deposit} :: DT :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_g_deposits  = new TH2F("DT_g_deposits",  "Simhit time vs E_{deposit} :: DT :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  DT_N_deposits  = new TH2F("DT_N_deposits",  "Simhit time vs E_{deposit} :: DT :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  RPCb_XY = new TH2F("RPCb_XY", "Simhits in XY :: RPCb", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCb_RZ = new TH2F("RPCb_RZ", "Simhits in RZ :: RPCb", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  RPCf_XY = new TH2F("RPCf_XY", "Simhits in XY :: RPCf", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCf_RZ = new TH2F("RPCf_RZ", "Simhits in RZ :: RPCf", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  CSC_XY  = new TH2F("CSC_XY",  "Simhits in XY :: CSC", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  CSC_RZ  = new TH2F("CSC_RZ",  "Simhits in RZ :: CSC", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  DT_XY   = new TH2F("DT_XY",   "Simhits in XY :: DT", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  DT_RZ   = new TH2F("DT_RZ",   "Simhits in RZ :: DT", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_RZ = new TH2F("Muon_RZ", "Simhits in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);

  RPCb_250ns_XY = new TH2F("RPCb_250ns_XY", "Simhits with tof > 250ns in XY :: RPCb", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCb_250ns_RZ = new TH2F("RPCb_250ns_RZ", "Simhits with tof > 250ns in RZ :: RPCb", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  RPCf_250ns_XY = new TH2F("RPCf_250ns_XY", "Simhits with tof > 250ns in XY :: RPCf", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  RPCf_250ns_RZ = new TH2F("RPCf_250ns_RZ", "Simhits with tof > 250ns in RZ :: RPCf", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  CSC_250ns_XY  = new TH2F("CSC_250ns_XY",  "Simhits with tof > 250ns in XY :: CSC", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  CSC_250ns_RZ  = new TH2F("CSC_250ns_RZ",  "Simhits with tof > 250ns in RZ :: CSC", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  DT_250ns_XY   = new TH2F("DT_250ns_XY",   "Simhits with tof > 250ns in XY :: DT", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  DT_250ns_RZ   = new TH2F("DT_250ns_RZ",   "Simhits with tof > 250ns in RZ :: DT", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  Muon_250ns_RZ = new TH2F("Muon_250ns_RZ", "Simhits with tof > 250ns in RZ :: Muon", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);

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

  RPCf_hits_tof = new TH1F("RPCf_hits_tof", "Simhit time :: RPCf", m_tof, m1_tof, m2_tof);
  RPCf_hits_eta = new TH1F("RPCf_hits_eta", "Simhit time :: RPCf", m_eta, m1_eta, m2_eta);
  RPCf_hits_phi = new TH1F("RPCf_hits_phi", "Simhit time :: RPCf", m_phi, m1_phi, m2_phi);

  CSC_hits_tof  = new TH1F("CSC_hits_tof",  "Simhit time :: CSC",  m_tof, m1_tof, m2_tof);
  CSC_hits_eta  = new TH1F("CSC_hits_eta",  "Simhit time :: CSC",  m_eta, m1_eta, m2_eta);
  CSC_hits_phi  = new TH1F("CSC_hits_phi",  "Simhit time :: CSC",  m_phi, m1_phi, m2_phi);

  DT_hits_tof   = new TH1F("DT_hits_tof",   "Simhit time :: DT",   m_tof, m1_tof, m2_tof);
  DT_hits_eta   = new TH1F("DT_hits_eta",   "Simhit time :: DT",   m_eta, m1_eta, m2_eta);
  DT_hits_phi   = new TH1F("DT_hits_phi",   "Simhit time :: DT",   m_phi, m1_phi, m2_phi);

}


MyNeutronSimHitAnalyzer::~MyNeutronSimHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  TStyle *plain  = new TStyle("Plain","Plain Style (no colours/fill areas)");
  plain->cd();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  TCanvas * Dummy = new TCanvas("dummy", "dummy", 600, 600);
  // First Plot for PDF File :: print empty Dummy 
  pdfFileName = pdfFileNameBase + ".pdf[";
  Dummy->Print(pdfFileName.c_str());
  // Name for next plot for PDF File
  pdfFileName = pdfFileNameBase + ".pdf";


  // Legend
  double l1_x1, l1_y1, l1_x2, l1_y2;
  // First version
  // l1_x1 = 0.60; l1_x2 = 0.85; l1_y1 = 0.60; l1_y2 = 0.85;
  // Second version
  l1_x1 = 0.60; l1_x2 = 0.85; l1_y1 = 0.65; l1_y2 = 0.85;
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

  Canvas_RPCb_hits = new TCanvas("Canvas_RPCb_hits", "Simhit time vs E_{kin} :: RPCb", 600, 600);
  Canvas_RPCf_hits = new TCanvas("Canvas_RPCf_hits", "Simhit time vs E_{kin} :: RPCf", 600, 600);
  Canvas_CSC_hits  = new TCanvas("Canvas_CSC_hits",  "Simhit time vs E_{kin} :: CSC",  600, 600);
  Canvas_DT_hits   = new TCanvas("Canvas_DT_hits",   "Simhit time vs E_{kin} :: DT",   600, 600);

  // Combine histograms in single plot
  Canvas_RPCb_hits->cd();
  RPCb_el_hits->GetXaxis()->SetTitle("log_{10} E_{kin} (MeV)");
  RPCb_el_hits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  RPCb_el_hits->SetTitle("SimHit time vs E_{kin} :: RPCb");
  RPCb_el_hits->SetMarkerStyle(7);  RPCb_el_hits->SetMarkerColor(kBlack);   RPCb_el_hits->SetMarkerSize(1);  RPCb_el_hits->Draw("P");
  RPCb_mu_hits->SetMarkerStyle(24); RPCb_mu_hits->SetMarkerColor(kBlue);    RPCb_mu_hits->SetMarkerSize(1);  RPCb_mu_hits->Draw("PSame");
  RPCb_pi_hits->SetMarkerStyle(33); RPCb_pi_hits->SetMarkerColor(kGreen);   RPCb_pi_hits->SetMarkerSize(1);  RPCb_pi_hits->Draw("PSame");
  RPCb_ka_hits->SetMarkerStyle(5);  RPCb_ka_hits->SetMarkerColor(kOrange);  RPCb_ka_hits->SetMarkerSize(1);  RPCb_ka_hits->Draw("PSame");
  RPCb_p_hits->SetMarkerStyle(26);  RPCb_p_hits->SetMarkerColor(kMagenta);  RPCb_p_hits->SetMarkerSize(1);   RPCb_p_hits->Draw("PSame");
  RPCb_n_hits->SetMarkerStyle(32);  RPCb_n_hits->SetMarkerColor(kViolet);   RPCb_n_hits->SetMarkerSize(1);   RPCb_n_hits->Draw("PSame");
  RPCb_g_hits->SetMarkerStyle(30);  RPCb_g_hits->SetMarkerColor(kCyan);     RPCb_g_hits->SetMarkerSize(1);   RPCb_g_hits->Draw("PSame");
  RPCb_N_hits->SetMarkerStyle(2);   RPCb_N_hits->SetMarkerColor(kRed);      RPCb_N_hits->SetMarkerSize(1);   RPCb_N_hits->Draw("PSame");
  l1->Draw();
  Canvas_RPCb_hits->Write();
  Canvas_RPCb_hits->Print(pdfFileName.c_str());

  Canvas_RPCf_hits->cd();
  RPCf_el_hits->GetXaxis()->SetTitle("log_{10} E_{kin} (MeV)");
  RPCf_el_hits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  RPCf_el_hits->SetTitle("SimHit time vs E_{kin} :: RPCf");
  RPCf_el_hits->SetMarkerStyle(7);  RPCf_el_hits->SetMarkerColor(kBlack);   RPCf_el_hits->SetMarkerSize(1);  RPCf_el_hits->Draw("P");
  RPCf_mu_hits->SetMarkerStyle(24); RPCf_mu_hits->SetMarkerColor(kBlue);    RPCf_mu_hits->SetMarkerSize(1);  RPCf_mu_hits->Draw("PSame");
  RPCf_pi_hits->SetMarkerStyle(33); RPCf_pi_hits->SetMarkerColor(kGreen);   RPCf_pi_hits->SetMarkerSize(1);  RPCf_pi_hits->Draw("PSame");
  RPCf_ka_hits->SetMarkerStyle(5);  RPCf_ka_hits->SetMarkerColor(kOrange);  RPCf_ka_hits->SetMarkerSize(1);  RPCf_ka_hits->Draw("PSame");
  RPCf_p_hits->SetMarkerStyle(26);  RPCf_p_hits->SetMarkerColor(kMagenta);  RPCf_p_hits->SetMarkerSize(1);   RPCf_p_hits->Draw("PSame");
  RPCf_n_hits->SetMarkerStyle(32);  RPCf_n_hits->SetMarkerColor(kViolet);   RPCf_n_hits->SetMarkerSize(1);   RPCf_n_hits->Draw("PSame");
  RPCf_g_hits->SetMarkerStyle(30);  RPCf_g_hits->SetMarkerColor(kCyan);     RPCf_g_hits->SetMarkerSize(1);   RPCf_g_hits->Draw("PSame");
  RPCf_N_hits->SetMarkerStyle(2);   RPCf_N_hits->SetMarkerColor(kRed);      RPCf_N_hits->SetMarkerSize(1);   RPCf_N_hits->Draw("PSame");
  l1->Draw();
  Canvas_RPCf_hits->Write();
  Canvas_RPCf_hits->Print(pdfFileName.c_str());


  Canvas_CSC_hits->cd();
  CSC_el_hits->GetXaxis()->SetTitle("log_{10} E_{kin} (MeV)");
  CSC_el_hits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  CSC_el_hits->SetTitle("SimHit time vs E_{kin} :: CSC");
  CSC_el_hits->SetMarkerStyle(7);  CSC_el_hits->SetMarkerColor(kBlack);   CSC_el_hits->SetMarkerSize(1);  CSC_el_hits->Draw("P");
  CSC_mu_hits->SetMarkerStyle(24); CSC_mu_hits->SetMarkerColor(kBlue);    CSC_mu_hits->SetMarkerSize(1);  CSC_mu_hits->Draw("PSame");
  CSC_pi_hits->SetMarkerStyle(33); CSC_pi_hits->SetMarkerColor(kGreen);   CSC_pi_hits->SetMarkerSize(1);  CSC_pi_hits->Draw("PSame");
  CSC_ka_hits->SetMarkerStyle(5);  CSC_ka_hits->SetMarkerColor(kOrange);  CSC_ka_hits->SetMarkerSize(1);  CSC_ka_hits->Draw("PSame");
  CSC_p_hits->SetMarkerStyle(26);  CSC_p_hits->SetMarkerColor(kMagenta);  CSC_p_hits->SetMarkerSize(1);   CSC_p_hits->Draw("PSame");
  CSC_n_hits->SetMarkerStyle(32);  CSC_n_hits->SetMarkerColor(kViolet);   CSC_n_hits->SetMarkerSize(1);   CSC_n_hits->Draw("PSame");
  CSC_g_hits->SetMarkerStyle(30);  CSC_g_hits->SetMarkerColor(kCyan);     CSC_g_hits->SetMarkerSize(1);   CSC_g_hits->Draw("PSame");
  CSC_N_hits->SetMarkerStyle(2);   CSC_N_hits->SetMarkerColor(kRed);      CSC_N_hits->SetMarkerSize(1);   CSC_N_hits->Draw("PSame");
  l1->Draw();
  Canvas_CSC_hits->Write();
  Canvas_CSC_hits->Print(pdfFileName.c_str());


  Canvas_DT_hits->cd();
  DT_el_hits->GetXaxis()->SetTitle("log_{10} E_{kin} (MeV)");
  DT_el_hits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  DT_el_hits->SetTitle("SimHit time vs E_{kin} :: DT");
  DT_el_hits->SetMarkerStyle(7);  DT_el_hits->SetMarkerColor(kBlack);   DT_el_hits->SetMarkerSize(1);  DT_el_hits->Draw("P");
  DT_mu_hits->SetMarkerStyle(24); DT_mu_hits->SetMarkerColor(kBlue);    DT_mu_hits->SetMarkerSize(1);  DT_mu_hits->Draw("PSame");
  DT_pi_hits->SetMarkerStyle(33); DT_pi_hits->SetMarkerColor(kGreen);   DT_pi_hits->SetMarkerSize(1);  DT_pi_hits->Draw("PSame");
  DT_ka_hits->SetMarkerStyle(5);  DT_ka_hits->SetMarkerColor(kOrange);  DT_ka_hits->SetMarkerSize(1);  DT_ka_hits->Draw("PSame");
  DT_p_hits->SetMarkerStyle(26);  DT_p_hits->SetMarkerColor(kMagenta);  DT_p_hits->SetMarkerSize(1);   DT_p_hits->Draw("PSame");
  DT_n_hits->SetMarkerStyle(32);  DT_n_hits->SetMarkerColor(kViolet);   DT_n_hits->SetMarkerSize(1);   DT_n_hits->Draw("PSame");
  DT_g_hits->SetMarkerStyle(30);  DT_g_hits->SetMarkerColor(kCyan);     DT_g_hits->SetMarkerSize(1);   DT_g_hits->Draw("PSame");
  DT_N_hits->SetMarkerStyle(2);   DT_N_hits->SetMarkerColor(kRed);      DT_N_hits->SetMarkerSize(1);   DT_N_hits->Draw("PSame");
  l1->Draw();
  Canvas_DT_hits->Write();
  Canvas_DT_hits->Print(pdfFileName.c_str());


  RPCf_el_hits->Write(); 
  RPCf_mu_hits->Write(); 
  RPCf_pi_hits->Write(); 
  RPCf_ka_hits->Write(); 
  RPCf_p_hits ->Write(); 
  RPCf_n_hits ->Write(); 
  RPCf_g_hits ->Write(); 
  RPCf_N_hits ->Write(); 

  RPCb_el_hits->Write(); 
  RPCb_mu_hits->Write(); 
  RPCb_pi_hits->Write(); 
  RPCb_ka_hits->Write(); 
  RPCb_p_hits ->Write(); 
  RPCb_n_hits ->Write(); 
  RPCb_g_hits ->Write(); 
  RPCb_N_hits ->Write(); 

  CSC_el_hits->Write(); 
  CSC_mu_hits->Write(); 
  CSC_pi_hits->Write(); 
  CSC_ka_hits->Write(); 
  CSC_p_hits ->Write(); 
  CSC_n_hits ->Write(); 
  CSC_g_hits ->Write(); 
  CSC_N_hits ->Write(); 

  DT_el_hits->Write(); 
  DT_mu_hits->Write(); 
  DT_pi_hits->Write(); 
  DT_ka_hits->Write(); 
  DT_p_hits ->Write(); 
  DT_n_hits ->Write(); 
  DT_g_hits ->Write(); 
  DT_N_hits ->Write(); 

  Canvas_RPCb_deposits = new TCanvas("Canvas_RPCb_deposits", "Simhit time vs E_{kin} :: RPCb", 600, 600);
  Canvas_RPCf_deposits = new TCanvas("Canvas_RPCf_deposits", "Simhit time vs E_{kin} :: RPCf", 600, 600);
  Canvas_CSC_deposits  = new TCanvas("Canvas_CSC_deposits",  "Simhit time vs E_{kin} :: CSC",  600, 600);
  Canvas_DT_deposits   = new TCanvas("Canvas_DT_deposits",   "Simhit time vs E_{kin} :: DT",   600, 600);

  // SimHit time vs Energy Deposit
  // Combine histograms in single plot
  Canvas_RPCb_deposits->cd();
  RPCb_el_deposits->GetXaxis()->SetTitle("log_{10} E_{deposit} (keV)");
  RPCb_el_deposits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  RPCb_el_deposits->SetTitle("SimHit time vs E_{deposit} :: RPCb");
  RPCb_el_deposits->SetMarkerStyle(7);  RPCb_el_deposits->SetMarkerColor(kBlack);   RPCb_el_deposits->SetMarkerSize(1);  RPCb_el_deposits->Draw("P");
  RPCb_mu_deposits->SetMarkerStyle(24); RPCb_mu_deposits->SetMarkerColor(kBlue);    RPCb_mu_deposits->SetMarkerSize(1);  RPCb_mu_deposits->Draw("PSame");
  RPCb_pi_deposits->SetMarkerStyle(33); RPCb_pi_deposits->SetMarkerColor(kGreen);   RPCb_pi_deposits->SetMarkerSize(1);  RPCb_pi_deposits->Draw("PSame");
  RPCb_ka_deposits->SetMarkerStyle(5);  RPCb_ka_deposits->SetMarkerColor(kOrange);  RPCb_ka_deposits->SetMarkerSize(1);  RPCb_ka_deposits->Draw("PSame");
  RPCb_p_deposits->SetMarkerStyle(26);  RPCb_p_deposits->SetMarkerColor(kMagenta);  RPCb_p_deposits->SetMarkerSize(1);   RPCb_p_deposits->Draw("PSame");
  RPCb_n_deposits->SetMarkerStyle(32);  RPCb_n_deposits->SetMarkerColor(kViolet);   RPCb_n_deposits->SetMarkerSize(1);   RPCb_n_deposits->Draw("PSame");
  RPCb_g_deposits->SetMarkerStyle(30);  RPCb_g_deposits->SetMarkerColor(kCyan);   RPCb_g_deposits->SetMarkerSize(1);   RPCb_g_deposits->Draw("PSame");
  RPCb_N_deposits->SetMarkerStyle(2);   RPCb_N_deposits->SetMarkerColor(kRed);      RPCb_N_deposits->SetMarkerSize(1);   RPCb_N_deposits->Draw("PSame");
  l1->Draw();
  Canvas_RPCb_deposits->Write();
  Canvas_RPCb_deposits->Print(pdfFileName.c_str());


  Canvas_RPCf_deposits->cd();
  RPCf_el_deposits->GetXaxis()->SetTitle("log_{10} E_{deposit} (keV)");
  RPCf_el_deposits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  RPCf_el_deposits->SetTitle("SimHit time vs E_{deposit} :: RPCf");
  RPCf_el_deposits->SetMarkerStyle(7);  RPCf_el_deposits->SetMarkerColor(kBlack);   RPCf_el_deposits->SetMarkerSize(1);  RPCf_el_deposits->Draw("P");
  RPCf_mu_deposits->SetMarkerStyle(24); RPCf_mu_deposits->SetMarkerColor(kBlue);    RPCf_mu_deposits->SetMarkerSize(1);  RPCf_mu_deposits->Draw("PSame");
  RPCf_pi_deposits->SetMarkerStyle(33); RPCf_pi_deposits->SetMarkerColor(kGreen);   RPCf_pi_deposits->SetMarkerSize(1);  RPCf_pi_deposits->Draw("PSame");
  RPCf_ka_deposits->SetMarkerStyle(5);  RPCf_ka_deposits->SetMarkerColor(kOrange);  RPCf_ka_deposits->SetMarkerSize(1);  RPCf_ka_deposits->Draw("PSame");
  RPCf_p_deposits->SetMarkerStyle(26);  RPCf_p_deposits->SetMarkerColor(kMagenta);  RPCf_p_deposits->SetMarkerSize(1);   RPCf_p_deposits->Draw("PSame");
  RPCf_n_deposits->SetMarkerStyle(32);  RPCf_n_deposits->SetMarkerColor(kViolet);   RPCf_n_deposits->SetMarkerSize(1);   RPCf_n_deposits->Draw("PSame");
  RPCf_g_deposits->SetMarkerStyle(30);  RPCf_g_deposits->SetMarkerColor(kCyan);     RPCf_g_deposits->SetMarkerSize(1);   RPCf_g_deposits->Draw("PSame");
  RPCf_N_deposits->SetMarkerStyle(2);   RPCf_N_deposits->SetMarkerColor(kRed);      RPCf_N_deposits->SetMarkerSize(1);   RPCf_N_deposits->Draw("PSame");
  l1->Draw();
  Canvas_RPCf_deposits->Write();
  Canvas_RPCf_deposits->Print(pdfFileName.c_str());


  Canvas_CSC_deposits->cd();
  CSC_el_deposits->GetXaxis()->SetTitle("log_{10} E_{deposit} (keV)");
  CSC_el_deposits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  CSC_el_deposits->SetTitle("SimHit time vs E_{deposit} :: CSC");
  CSC_el_deposits->SetMarkerStyle(7);  CSC_el_deposits->SetMarkerColor(kBlack);   CSC_el_deposits->SetMarkerSize(1);  CSC_el_deposits->Draw("P");
  CSC_mu_deposits->SetMarkerStyle(24); CSC_mu_deposits->SetMarkerColor(kBlue);    CSC_mu_deposits->SetMarkerSize(1);  CSC_mu_deposits->Draw("PSame");
  CSC_pi_deposits->SetMarkerStyle(33); CSC_pi_deposits->SetMarkerColor(kGreen);   CSC_pi_deposits->SetMarkerSize(1);  CSC_pi_deposits->Draw("PSame");
  CSC_ka_deposits->SetMarkerStyle(5);  CSC_ka_deposits->SetMarkerColor(kOrange);  CSC_ka_deposits->SetMarkerSize(1);  CSC_ka_deposits->Draw("PSame");
  CSC_p_deposits->SetMarkerStyle(26);  CSC_p_deposits->SetMarkerColor(kMagenta);  CSC_p_deposits->SetMarkerSize(1);   CSC_p_deposits->Draw("PSame");
  CSC_n_deposits->SetMarkerStyle(32);  CSC_n_deposits->SetMarkerColor(kViolet);   CSC_n_deposits->SetMarkerSize(1);   CSC_n_deposits->Draw("PSame");
  CSC_g_deposits->SetMarkerStyle(30);  CSC_g_deposits->SetMarkerColor(kCyan);     CSC_g_deposits->SetMarkerSize(1);   CSC_g_deposits->Draw("PSame");
  CSC_N_deposits->SetMarkerStyle(2);   CSC_N_deposits->SetMarkerColor(kRed);      CSC_N_deposits->SetMarkerSize(1);   CSC_N_deposits->Draw("PSame");
  l1->Draw();
  Canvas_CSC_deposits->Write();
  Canvas_CSC_deposits->Print(pdfFileName.c_str());


  Canvas_DT_deposits->cd();
  DT_el_deposits->GetXaxis()->SetTitle("log_{10} E_{deposit} (keV)");
  DT_el_deposits->GetYaxis()->SetTitle("log_{10} TOF (ns)");
  DT_el_deposits->SetTitle("SimHit time vs E_{deposit} :: DT");
  DT_el_deposits->SetMarkerStyle(7);  DT_el_deposits->SetMarkerColor(kBlack);   DT_el_deposits->SetMarkerSize(1);  DT_el_deposits->Draw("P");
  DT_mu_deposits->SetMarkerStyle(24); DT_mu_deposits->SetMarkerColor(kBlue);    DT_mu_deposits->SetMarkerSize(1);  DT_mu_deposits->Draw("PSame");
  DT_pi_deposits->SetMarkerStyle(33); DT_pi_deposits->SetMarkerColor(kGreen);   DT_pi_deposits->SetMarkerSize(1);  DT_pi_deposits->Draw("PSame");
  DT_ka_deposits->SetMarkerStyle(5);  DT_ka_deposits->SetMarkerColor(kOrange);  DT_ka_deposits->SetMarkerSize(1);  DT_ka_deposits->Draw("PSame");
  DT_p_deposits->SetMarkerStyle(26);  DT_p_deposits->SetMarkerColor(kMagenta);  DT_p_deposits->SetMarkerSize(1);   DT_p_deposits->Draw("PSame");
  DT_n_deposits->SetMarkerStyle(32);  DT_n_deposits->SetMarkerColor(kViolet);   DT_n_deposits->SetMarkerSize(1);   DT_n_deposits->Draw("PSame");
  DT_g_deposits->SetMarkerStyle(30);  DT_g_deposits->SetMarkerColor(kCyan);     DT_g_deposits->SetMarkerSize(1);   DT_g_deposits->Draw("PSame");
  DT_N_deposits->SetMarkerStyle(2);   DT_N_deposits->SetMarkerColor(kRed);      DT_N_deposits->SetMarkerSize(1);   DT_N_deposits->Draw("PSame");
  l1->Draw();
  Canvas_DT_deposits->Write();
  Canvas_DT_deposits->Print(pdfFileName.c_str());


  RPCf_el_deposits->Write(); 
  RPCf_mu_deposits->Write(); 
  RPCf_pi_deposits->Write(); 
  RPCf_ka_deposits->Write(); 
  RPCf_p_deposits ->Write(); 
  RPCf_n_deposits ->Write(); 
  RPCf_g_deposits ->Write(); 
  RPCf_N_deposits ->Write(); 

  RPCb_el_deposits->Write(); 
  RPCb_mu_deposits->Write(); 
  RPCb_pi_deposits->Write(); 
  RPCb_ka_deposits->Write(); 
  RPCb_p_deposits ->Write(); 
  RPCb_n_deposits ->Write(); 
  RPCb_g_deposits ->Write(); 
  RPCb_N_deposits ->Write(); 

  CSC_el_deposits->Write(); 
  CSC_mu_deposits->Write(); 
  CSC_pi_deposits->Write(); 
  CSC_ka_deposits->Write(); 
  CSC_p_deposits ->Write(); 
  CSC_n_deposits ->Write(); 
  CSC_g_deposits ->Write(); 
  CSC_N_deposits ->Write(); 

  DT_el_deposits->Write(); 
  DT_mu_deposits->Write(); 
  DT_pi_deposits->Write(); 
  DT_ka_deposits->Write(); 
  DT_p_deposits ->Write(); 
  DT_n_deposits ->Write(); 
  DT_g_deposits ->Write(); 
  DT_N_deposits ->Write(); 


  Canvas_Muon_RZ = new TCanvas("Canvas_Muon_RZ", "RZ-view of SimHits :: Muon", 600, 600);
  Muon_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_RZ->SetTitle("RZ-view of SimHits :: Muon");
  Muon_RZ->Draw("colz");
  Canvas_Muon_RZ->Write(); 
  Canvas_Muon_RZ->Print(pdfFileName.c_str());

  Canvas_Muon_250ns_RZ = new TCanvas("Canvas_Muon_250ns_RZ", "RZ-view of SimHits with tof > 250 ns :: Muon", 600, 600);
  Muon_250ns_RZ->GetXaxis()->SetTitle("Z (cm)");
  Muon_250ns_RZ->GetYaxis()->SetTitle("R (cm)");
  Muon_250ns_RZ->GetYaxis()->SetTitleOffset(1.30);
  Muon_250ns_RZ->SetTitle("RZ-view of SimHits with tof > 250 ns :: Muon");
  Muon_250ns_RZ->Draw("colz");
  Canvas_Muon_250ns_RZ->Write(); 
  Canvas_Muon_250ns_RZ->Print(pdfFileName.c_str());

  Canvas_SimVertices_RZ = new TCanvas("Canvas_SimVertices_RZ", "RZ-view of Sim Vertices", 600, 600);
  SimVertices_RZ->GetXaxis()->SetTitle("Z (cm)");
  SimVertices_RZ->GetYaxis()->SetTitle("R (cm)");
  SimVertices_RZ->GetYaxis()->SetTitleOffset(1.30);
  SimVertices_RZ->SetTitle("RZ-view of Sim Vertices");
  SimVertices_RZ->Draw("colz");
  Canvas_SimVertices_RZ->Write(); 
  Canvas_SimVertices_RZ->Print(pdfFileName.c_str());

  Canvas_SimVertices_Muon_RZ = new TCanvas("Canvas_SimVertices_Muon_RZ", "RZ-view of Sim Vertices in Muon System and Cavern", 600, 600);
  SimVertices_Muon_RZ->GetXaxis()->SetTitle("Z (cm)");
  SimVertices_Muon_RZ->GetYaxis()->SetTitle("R (cm)");
  SimVertices_Muon_RZ->GetYaxis()->SetTitleOffset(1.30);
  SimVertices_Muon_RZ->SetTitle("RZ-view of Sim Vertices in Muon System and Cavern");
  SimVertices_Muon_RZ->SetMarkerStyle(5);
  SimVertices_Muon_RZ->SetMarkerColor(kBlack);
  SimVertices_Muon_RZ->SetMarkerSize(1.25);
  SimVertices_Muon_RZ->Draw("P");
  Canvas_SimVertices_Muon_RZ->Write(); 
  Canvas_SimVertices_Muon_RZ->Print(pdfFileName.c_str());

  // last plot for PDF File :: print empty Dummy
  pdfFileName = pdfFileNameBase + ".pdf]";
  Dummy->Print(pdfFileName.c_str());


  PrimVertices_Z->Write();
  PrimVertices_R->Write();

  RPCb_XY->Write();
  RPCb_RZ->Write();
  RPCf_XY->Write();
  RPCf_RZ->Write();
  CSC_XY->Write();
  CSC_RZ->Write();
  DT_XY->Write();
  DT_RZ->Write();
  Muon_RZ->Write();

  RPCb_250ns_XY->Write();
  RPCb_250ns_RZ->Write();
  RPCf_250ns_XY->Write();
  RPCf_250ns_RZ->Write();
  CSC_250ns_XY->Write();
  CSC_250ns_RZ->Write();
  DT_250ns_XY->Write();
  DT_250ns_RZ->Write();
  Muon_250ns_RZ->Write();

  RPCb_hits_tof->Write();
  RPCb_hits_eta->Write();
  RPCb_hits_phi->Write();

  RPCf_hits_tof->Write();
  RPCf_hits_eta->Write();
  RPCf_hits_phi->Write();

  CSC_hits_tof->Write();
  CSC_hits_eta->Write();
  CSC_hits_phi->Write();

  DT_hits_tof->Write();
  DT_hits_eta->Write();
  DT_hits_phi->Write();

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
  int pu = 21; // measured ...

  // 2015
  // com = 14 TeV
  // sigma_inel = 80 mb (14 TeV)
  // L = 10^34 1/cm^2 1/s (design)
  // int pu = 25;

  int entries = 2500;
  int bx_space = 50;
  double corr_fact = pow(10,9) * 1.0/(entries * 1.0/pu * bx_space);
  double sg_corr_fact = 2;  // simulation has only single gas volumes implemented instead of double gas volumes
  std::cout<<"Rate Correction factor =  pow(10,9) * 1.0/(entries * 1.0/pu * bx_space) = "<<pow(10,9) * 1.0/(entries * 1.0/pu * bx_space)<<std::endl;
  std::cout<<"Rate Correction factor because of simulation of Single Gas Layers = "<<sg_corr_fact<<std::endl;
  std::cout<<"RPC Hits in Barrel = "<<rpc_barrel_hits<<" || RPC Hits in Endcap = "<<rpc_endcap_hits<<std::endl;
  std::cout<<"RPC Barrel Active area = "<<rpc_barrel_area<<" cm2 || RPC Endcap Active area = "<<rpc_endcap_area<<" cm2"<<std::endl;

  RPC_hits->SetBinContent(1,  rpc_barrel_hits);
  RPC_area->SetBinContent(1,  rpc_barrel_area);
  if(rpc_barrel_hits > 0 && rpc_barrel_area > 0.0) {
    RPC_rates->SetBinContent(1,  rpc_barrel_hits*1.0/rpc_barrel_area*corr_fact*sg_corr_fact);
    RPC_rates->SetBinError(  1,  sqrt(rpc_barrel_hits)*1.0/rpc_barrel_area*corr_fact*sg_corr_fact);
  }
  for(int k=0; k<3; ++k) {  // Barrel :: ring == wheel
   for(int j=0; j<4; ++j) {
      RPC_hits->SetBinContent(1+4*(k)+(j+1),RPC_hits_array[0][j][k]);
      RPC_area->SetBinContent(1+4*(k)+(j+1),RPC_area_array[0][j][k]);
      // std::cout<<"Bin "<<1+4*(k)+(j+1)<<" filled with RPC_area_array[region="<<0<<"][station"<<j+1<<"]["<<k<<"] = "<<RPC_area_array[0][j][k]<<" cm2"<<std::endl;
      if(RPC_hits_array[0][j][k] >  0 && RPC_area_array[0][j][k] > 0.0) {
	RPC_rates->SetBinContent(1+4*(k)+(j+1), RPC_hits_array[0][j][k]*1.0/RPC_area_array[0][j][k]*corr_fact*sg_corr_fact);
	RPC_rates->SetBinError(  1+4*(k)+(j+1), sqrt(RPC_hits_array[0][j][k])*1.0/RPC_area_array[0][j][k]*corr_fact*sg_corr_fact);
      }
    }
  }
  RPC_hits->SetBinContent(14, rpc_endcap_hits);
  RPC_area->SetBinContent(14, rpc_endcap_area);
  if(rpc_endcap_hits > 0 && rpc_endcap_area > 0.0) {
    RPC_rates->SetBinContent(14,  rpc_endcap_hits*1.0/rpc_endcap_area*corr_fact*sg_corr_fact);
    RPC_rates->SetBinError(  14,  sqrt(rpc_endcap_hits)*1.0/rpc_endcap_area*corr_fact*sg_corr_fact);
  }
  for(int j=0; j<4; ++j) {  // Endcap :: station == disk
    for(int k=0; k<3; ++k) {
      RPC_hits->SetBinContent(14+3*j+(k+1),RPC_hits_array[1][j][k]);
      RPC_area->SetBinContent(14+3*j+(k+1),RPC_area_array[1][j][k]);
      // std::cout<<"Bin "<<14+3*j+(k+1)<<" filled with RPC_area_array[region="<<1<<"][station"<<j+1<<"]["<<k+1<<"] = "<<RPC_area_array[1][j][k]<<" cm2"<<std::endl;
      if(RPC_hits_array[1][j][k] >  0 && RPC_area_array[1][j][k] > 0.0) {
	RPC_rates->SetBinContent(14+3*j+(k+1), RPC_hits_array[1][j][k]*1.0/RPC_area_array[1][j][k]*corr_fact*sg_corr_fact);
	RPC_rates->SetBinError(  14+3*j+(k+1), sqrt(RPC_hits_array[1][j][k])*1.0/RPC_area_array[1][j][k]*corr_fact*sg_corr_fact);
      }
    }
  }
  for(int m=0; m<26; ++m) {
    RPC_hits->GetXaxis()->SetBinLabel( m+1, cat[m].c_str());
    RPC_area->GetXaxis()->SetBinLabel( m+1, cat[m].c_str());
    RPC_rates->GetXaxis()->SetBinLabel(m+1, cat[m].c_str());
  }

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

  RPC_hits->Write();
  RPC_area->Write();
  RPC_rates->Write();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyNeutronSimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // ===============================
  //      General Collections
  // ===============================
  
  // SimHits
  // std::cout << " Getting the SimHits " <<std::endl;
  std::vector<edm::Handle<edm::PSimHitContainer> > theSimHitContainers;
  iEvent.getManyByType(theSimHitContainers);
  if(debug) std::cout << " The number of SimHit Containers is  " << theSimHitContainers.size() <<std::endl;
  std::vector<PSimHit> theSimHits;
  for (int i = 0; i < int(theSimHitContainers.size()); ++i) {
    theSimHits.insert(theSimHits.end(),theSimHitContainers.at(i)->begin(),theSimHitContainers.at(i)->end());
  }
  // SimTracks
  std::vector<SimTrack> theSimTracks;
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel("g4SimHits",SimTk);
  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  if(debug) std::cout << "This Event has " <<  theSimTracks.size() << " sim tracks " << std::endl;
  // SimVertices
  std::vector<SimVertex> theSimVertices; 
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel("g4SimHits",SimVtx);
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  if(debug) std::cout << "This Event has " <<  theSimVertices.size() << " sim vertices " << std::endl;

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
    if(vtx_r > 400 || vtx_z > 500 ) {
      SimVertices_Muon_RZ->Fill(fabs(vtx_z),vtx_r);
    }
  }

  
  // create a map associating geant particle id and position in the 
  // event SimTrack vector
  std::map<unsigned, unsigned> geantToIndex;
  for( unsigned it=0; it<theSimTracks.size(); ++it ) {
    geantToIndex[ theSimTracks[it].trackId() ] = it;
  }  

  // ===============================
  //      Loop over the SimHits
  // ===============================

  for (std::vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); ++iHit) {
    DetId theDetUnitId((*iHit).detUnitId());
    DetId simdetid= DetId((*iHit).detUnitId());


    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){ // Only RPCs

      // RPC Geometry
      // ============
      RPCDetId rollId(theDetUnitId);
      RPCGeomServ rpcsrv(rollId);
      const RPCRoll* rollasociated = rpcGeom->roll(rollId);
      const BoundPlane & RPCSurface = rollasociated->surface();
      GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal((*iHit).localPosition());
      if(debug) {
	std::cout<<"RPC SimHit in "<<std::setw(12)<<(int)rollId<< /*" a.k.a. "<<std::setw(24)<<rpcsrv.name()<<" details: "<<*/ std::setw(24)<<rollId;
	std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
	std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta();
	std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
      }
      double RPC_GlobalPoint_R = sqrt(pow(RPCGlobalPoint.x(),2)+pow(RPCGlobalPoint.y(),2));

      // SIMHIT Investigations
      // =====================
      int pid           = (*iHit).particleType();

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

      // time and energy
      // ---------------
      double log_time   = log10((*iHit).timeOfFlight());
      double log_energy = log10((*iHit).momentumAtEntry().perp()*1000);
      double log_deposit = log10((*iHit).energyLoss()*1000000);

      // count hits
      // ----------
      int region = abs(rollId.region());
      int station = abs(rollId.station())-1;
      int ring = (region==0)? abs(rollId.ring()) : abs(rollId.ring())-1;
      ++RPC_hits_array[region][station][ring];

      // RPCb
      // ----
      if(rollId.region() == 0) {
	RPCb_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());             if((*iHit).timeOfFlight()>250) { RPCb_250ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y()); }
	RPCb_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));  if((*iHit).timeOfFlight()>250) { RPCb_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); }
	Muon_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));  if((*iHit).timeOfFlight()>250) { Muon_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); }
	// std::cout<<"RPCb :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and deposit deposit "<<(*iHit).energyLoss()<< " [GeV]";
	// std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
	if(abs(pid)==2212)      {RPCb_p_hits->Fill(log_energy,log_time); RPCb_p_deposits->Fill(log_deposit,log_time);}
	else if(abs(pid)==2112) {RPCb_n_hits->Fill(log_energy,log_time); RPCb_n_deposits->Fill(log_deposit,log_time);}
	else { 
	  switch (abs(pid)%1000) {
	    // leptons
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
	  default:   {
	    RPCb_N_hits->Fill(log_energy,log_time); RPCb_N_deposits->Fill(log_deposit,log_time);
	    std::cout<<"RPCb :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and deposit deposit "<<(*iHit).energyLoss()<< " [GeV]";
	    std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV]"<<std::endl;
	    break;
	  }
	  }
	}
	RPCb_hits_tof->Fill(log_time);
	RPCb_hits_eta->Fill(fabs(RPCGlobalPoint.eta()));
	RPCb_hits_phi->Fill(RPCGlobalPoint.phi());
      }
      if(rollId.region() != 0) {
	RPCf_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());            if((*iHit).timeOfFlight()>250) {RPCf_250ns_XY->Fill(RPCGlobalPoint.x(), RPCGlobalPoint.y());}
	RPCf_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); if((*iHit).timeOfFlight()>250) {RPCf_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));}
	Muon_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R)); if((*iHit).timeOfFlight()>250) {Muon_250ns_RZ->Fill(fabs(RPCGlobalPoint.z()), fabs(RPC_GlobalPoint_R));}
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
	    RPCf_N_hits->Fill(log_energy,log_time); RPCf_N_deposits->Fill(log_deposit,log_time); 
	    std::cout<<"RPCf :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<<" [GeV]";
	    std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV]"<<std::endl;
	    break;
	  }
	  }
	}
	RPCf_hits_tof->Fill(log_time);
	RPCf_hits_eta->Fill(fabs(RPCGlobalPoint.eta()));
	RPCf_hits_phi->Fill(RPCGlobalPoint.phi());
      }
    }

    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::CSC){ // Only CSCs
      CSCDetId rollId(theDetUnitId);
      GlobalPoint CSCGlobalPoint = cscGeom->idToDet(rollId)->toGlobal((*iHit).localPosition());
      if (debug) {
       std::cout<<"CSC SimHit in "<<std::setw(24)<<rollId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
       std::cout<<" | r = "<<std::setw(12)<<CSCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<CSCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<CSCGlobalPoint.eta();
       std::cout<<" | global position = "<<CSCGlobalPoint<<std::endl;
      }
      double CSC_GlobalPoint_R = sqrt(pow(CSCGlobalPoint.x(),2)+pow(CSCGlobalPoint.y(),2));
      CSC_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());             if((*iHit).timeOfFlight()>250) {CSC_250ns_XY->Fill(CSCGlobalPoint.x(), CSCGlobalPoint.y());}
      CSC_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));  if((*iHit).timeOfFlight()>250) {CSC_250ns_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));}
      Muon_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R)); if((*iHit).timeOfFlight()>250) {Muon_250ns_RZ->Fill(fabs(CSCGlobalPoint.z()), fabs(CSC_GlobalPoint_R));}

      int pid           = (*iHit).particleType();
      // get track
      int trackId = (*iHit).trackId();
      const SimTrack& track = theSimTracks[trackId];
      std::cout<<"TrackId = "<<trackId<<" track = "<<&track<<std::endl;
      // int trackType = track.type();
      // std::cout<<"TrackId = "<<trackId<<" type = "<<trackType<<std::endl;

      /*
      // get vertex
      int vertexId = track.vertIndex();
      const SimVertex& vertex = theSimVertices[vertexId];
      std::cout<<"TrackId = "<<trackId<<" VertexId = "<<vertexId<<std::endl;
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
      double log_time   = log10((*iHit).timeOfFlight());                 
      double log_energy = log10((*iHit).momentumAtEntry().perp()*1000);  
      double log_deposit = log10((*iHit).energyLoss()*1000000);          
      // std::cout<<"CSC :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
      // std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [MeV]"<<std::endl;
      if(abs(pid)==2212)      {CSC_p_hits->Fill(log_energy,log_time); CSC_p_deposits->Fill(log_deposit,log_time);}
      else if(abs(pid)==2112) {CSC_n_hits->Fill(log_energy,log_time); CSC_n_deposits->Fill(log_deposit,log_time);}
      else {
	switch (abs(pid)%1000) {
	  // leptons
	case 11:  CSC_el_hits->Fill(log_energy,log_time); CSC_el_deposits->Fill(log_deposit,log_time); break;
	case 13:  CSC_mu_hits->Fill(log_energy,log_time); CSC_mu_deposits->Fill(log_deposit,log_time); break;
	  // Pions
	case 111: CSC_pi_hits->Fill(log_energy,log_time); CSC_pi_deposits->Fill(log_deposit,log_time); break;
	case 211: CSC_pi_hits->Fill(log_energy,log_time); CSC_pi_deposits->Fill(log_deposit,log_time); break;
	case 130: CSC_pi_hits->Fill(log_energy,log_time); CSC_pi_deposits->Fill(log_deposit,log_time); break;
	  // Kaons
	case 310: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 311: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 321: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 313: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 323: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 315: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 325: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 317: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 327: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 319: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	case 329: CSC_ka_hits->Fill(log_energy,log_time); CSC_ka_deposits->Fill(log_deposit,log_time); break;
	  // Protons
	  // case 2212: CSC_p_hits->Fill(log_energy,log_time); CSC_p_deposits->Fill(log_deposit,log_time); break;
	  // Neutrons
	  // case 2112: CSC_n_hits->Fill(log_energy,log_time); CSC_n_deposits->Fill(log_deposit,log_time); break;
	  // Photons
	case 22:   CSC_g_hits->Fill(log_energy,log_time); CSC_g_deposits->Fill(log_deposit,log_time); break;
	  // Nucleons
	default:   {
	  CSC_N_hits->Fill(log_energy,log_time); CSC_N_deposits->Fill(log_deposit,log_time); 
	  std::cout<<"CSC :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV]"<<std::endl;
	  break;
	}
	}
      }
      CSC_hits_tof->Fill(log_time);
      CSC_hits_eta->Fill(fabs(CSCGlobalPoint.eta()));
      CSC_hits_phi->Fill(CSCGlobalPoint.phi());
    }
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::DT){
      DTWireId wireId(theDetUnitId);
      GlobalPoint DTGlobalPoint = dtGeom->idToDet(wireId)->toGlobal((*iHit).localPosition());
      // dtGeometry->idToDet(id)->surface().toGlobal(LocalPoint(point->displacement.x(), point->displacement.y(), point->displacement.z()));
      if (debug) {
       std::cout<<"DT SimHit in "<<std::setw(24)<<wireId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<DTGlobalPoint.z();
       std::cout<<" | r = "<<std::setw(12)<<DTGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<DTGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<DTGlobalPoint.eta();
       std::cout<<" | global position = "<<DTGlobalPoint<<std::endl;
      }
      double DT_GlobalPoint_R = sqrt(pow(DTGlobalPoint.x(),2)+pow(DTGlobalPoint.y(),2));
      DT_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());               if((*iHit).timeOfFlight()>250) {DT_250ns_XY->Fill(DTGlobalPoint.x(), DTGlobalPoint.y());}
      DT_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));    if((*iHit).timeOfFlight()>250) {DT_250ns_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));}
      Muon_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));  if((*iHit).timeOfFlight()>250) {Muon_250ns_RZ->Fill(fabs(DTGlobalPoint.z()), fabs(DT_GlobalPoint_R));}

      int pid           = (*iHit).particleType();
      int trackid = (*iHit).trackId();
      std::cout<<"TrackId = "<<trackid<<std::endl;
      double log_time   = log10((*iHit).timeOfFlight());
      double log_energy = log10((*iHit).momentumAtEntry().perp()*1000);
      double log_deposit = log10((*iHit).energyLoss()*1000000);
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
	  DT_N_hits->Fill(log_energy,log_time); DT_N_deposits->Fill(log_deposit,log_time); 
	  std::cout<<"DT :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	  std::cout<<" 10 log (tof) = "<<log_time<<" [ns] and 10 log (E) = "<<log_deposit<<" [keV]"<<std::endl;
	  break;
	}
	}
      }
      DT_hits_tof->Fill(log_time);
      DT_hits_eta->Fill(fabs(DTGlobalPoint.eta()));
      DT_hits_phi->Fill(DTGlobalPoint.phi());
    }
  }
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
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
  iSetup.get<MuonGeometryRecord>().get(dtGeom);

  // Loop over RPC Geometry
  // std::cout <<"Analyze RPC Geometry :: Loop over RPC Chambers"<<std::endl;
  for (TrackingGeometry::DetContainer::const_iterator it=rpcGeom->dets().begin(); it<rpcGeom->dets().end(); ++it) {
    if( dynamic_cast< RPCChamber* >( *it ) != 0 ){
      RPCChamber* ch = dynamic_cast< RPCChamber* >( *it );
      std::vector< const RPCRoll*> rolls = (ch->rolls());
      for(std::vector<const RPCRoll*>::const_iterator r = rolls.begin();r != rolls.end(); ++r) {
	RPCDetId rpcId = (*r)->id();
	int n_strips=(*r)->nstrips();
	RPCGeomServ rpcsrv(rpcId);

	int region = abs(rpcId.region());
	int station = abs(rpcId.station())-1;
	int ring = (region==0)? abs(rpcId.ring()) : abs(rpcId.ring())-1; // region == 0 ? Barrel : Endcap

	if (region == 0) {
	  const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&((*r)->topology()));
	  float stripl = top_->stripLength();
	  float stripw = top_->pitch();
	  rpc_barrel_area += stripw*stripl*n_strips;
	  RPC_area_array[region][station][ring] += stripw*stripl*n_strips;
	  // std::cout<<(int) rpcId<<" aka "<<rpcId<<" = "<< stripw*stripl*n_strips <<" cm2";
	  // std::cout<<"==> RPC_area_array[region="<<region<<"][station="<<station+1<<"][ring="<<ring<<"] = "<<RPC_area_array[region][station][ring]<<" cm2"<<std::endl;
	}
	if(region == 1) {
	  const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&((*r)->topology()));
	  float stripl = top_->stripLength();
	  float stripw = top_->pitch();
	  rpc_endcap_area += stripw*stripl*n_strips;
	  RPC_area_array[region][station][ring] += stripw*stripl*n_strips;
	  // std::cout<<(int) rpcId<<" aka "<<rpcId<<" = "<< stripw*stripl*n_strips <<" cm2";
	  // std::cout<<" ==> RPC_area_array[region="<<region<<"][station="<<station+1<<"][ring"<<ring+1<<"] = "<<RPC_area_array[region][station][ring]<<" cm2"<<std::endl;
	}
      }
    }
  }
  // std::cout<<"RPC Barrel Active area = "<<rpc_barrel_area<<" cm2 || RPC Endcap Active area = "<<rpc_endcap_area<<" cm2"<<std::endl;
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
