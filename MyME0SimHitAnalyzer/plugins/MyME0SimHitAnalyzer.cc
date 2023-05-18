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
// #include "FWCore/Framework/interface/EDAnalyzer.h"   // pre-13X
#include "FWCore/Framework/interface/one/EDAnalyzer.h"  // 13X
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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
// DetIds
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include "DataFormats/MuonDetId/interface/DTWireId.h"
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

// OLD class MyME0SimHitAnalyzer : public edm::EDAnalyzer {
class MyME0SimHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit MyME0SimHitAnalyzer(const edm::ParameterSet&);
      ~MyME0SimHitAnalyzer();
  // edm::ESHandle <RPCGeometry> rpcGeom;
  edm::ESHandle <CSCGeometry> cscGeom;
  // edm::ESHandle <DTGeometry>  dtGeom;
  edm::ESHandle <GEMGeometry> gemGeom;
  // edm::ESHandle <ME0Geometry> me0Geom;

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
  double edepMin, hipMin;
  bool phys_debug, tech_debug, edep_30eV; // gem_only_ge11, 
  TFile * outputfile;


  // required for 7XY: registration of the data access
  // see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent
  edm::EDGetTokenT<edm::SimTrackContainer>  SIMTrack_Token;
  edm::EDGetTokenT<edm::SimVertexContainer> SIMVertex_Token;

  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeom_Token;
  edm::ESGetToken<DTGeometry,  MuonGeometryRecord> dtGeom_Token;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeom_Token;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeom_Token;

  TH1F * EventCounter;

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_radius, * TDir_Muon_250ns_radius, * TDir_Muon_25ns_radius;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * ME0_el_hits_R,  * ME0_mu_hits_R,  * ME0_pi_hits_R,  * ME0_ka_hits_R,  * ME0_p_hits_R,  * ME0_n_hits_R,  * ME0_g_hits_R,  * ME0_N_hits_R, * ME0_OH_hits_R;
  TH1F * G11_el_hits_R,  * G11_mu_hits_R,  * G11_pi_hits_R,  * G11_ka_hits_R,  * G11_p_hits_R,  * G11_n_hits_R,  * G11_g_hits_R,  * G11_N_hits_R, * G11_OH_hits_R;
  TH1F * G21_el_hits_R,  * G21_mu_hits_R,  * G21_pi_hits_R,  * G21_ka_hits_R,  * G21_p_hits_R,  * G21_n_hits_R,  * G21_g_hits_R,  * G21_N_hits_R, * G21_OH_hits_R;

  TH1F * ME0_All_hits_R, * ME0_All_hits_000_R, * ME0_All_hits_250_R, * ME0_All_hits_00_R, * ME0_All_hits_25_R, * ME0_HIP_hits_R;
  TH1F * G11_All_hits_R, * G11_All_hits_000_R, * G11_All_hits_250_R, * G11_All_hits_00_R, * G11_All_hits_25_R, * G11_HIP_hits_R; 
  TH1F * G21_All_hits_R, * G21_All_hits_000_R, * G21_All_hits_250_R, * G21_All_hits_00_R, * G21_All_hits_25_R, * G21_HIP_hits_R; 

  // future :: need to split GE11 in short & long
  TH1F * G11_L1_el_hits_R,  * G11_L1_mu_hits_R,  * G11_L1_pi_hits_R,  * G11_L1_ka_hits_R,  * G11_L1_p_hits_R,  * G11_L1_n_hits_R,  * G11_L1_g_hits_R,  * G11_L1_N_hits_R, * G11_L1_OH_hits_R,  * G11_L1_All_hits_R,  * G11_L1_HIP_hits_R;
  TH1F * G11_L2_el_hits_R,  * G11_L2_mu_hits_R,  * G11_L2_pi_hits_R,  * G11_L2_ka_hits_R,  * G11_L2_p_hits_R,  * G11_L2_n_hits_R,  * G11_L2_g_hits_R,  * G11_L2_N_hits_R, * G11_L2_OH_hits_R,  * G11_L2_All_hits_R,  * G11_L2_HIP_hits_R;
  TH1F * G11_Od_el_hits_R,  * G11_Od_mu_hits_R,  * G11_Od_pi_hits_R,  * G11_Od_ka_hits_R,  * G11_Od_p_hits_R,  * G11_Od_n_hits_R,  * G11_Od_g_hits_R,  * G11_Od_N_hits_R, * G11_Od_OH_hits_R,  * G11_Od_All_hits_R,  * G11_Od_HIP_hits_R;
  TH1F * G11_Ev_el_hits_R,  * G11_Ev_mu_hits_R,  * G11_Ev_pi_hits_R,  * G11_Ev_ka_hits_R,  * G11_Ev_p_hits_R,  * G11_Ev_n_hits_R,  * G11_Ev_g_hits_R,  * G11_Ev_N_hits_R, * G11_Ev_OH_hits_R,  * G11_Ev_All_hits_R,  * G11_Ev_HIP_hits_R;
  // TH1F * M11_Od_el_hits_R,  * M11_Od_mu_hits_R,  * M11_Od_pi_hits_R,  * M11_Od_ka_hits_R,  * M11_Od_p_hits_R,  * M11_Od_n_hits_R,  * M11_Od_g_hits_R,  * M11_Od_N_hits_R, * M11_Od_OH_hits_R,  * M11_Od_All_hits_R,  * M11_Od_HIP_hits_R;
  // TH1F * M11_Ev_el_hits_R,  * M11_Ev_mu_hits_R,  * M11_Ev_pi_hits_R,  * M11_Ev_ka_hits_R,  * M11_Ev_p_hits_R,  * M11_Ev_n_hits_R,  * M11_Ev_g_hits_R,  * M11_Ev_N_hits_R, * M11_Ev_OH_hits_R,  * M11_Ev_All_hits_R,  * M11_Ev_HIP_hits_R;
  // TH1F * M11_Od_L1_All_hits_R, * M11_Od_L2_All_hits_R, * M11_Od_L3_All_hits_R, * M11_Od_L4_All_hits_R, * M11_Od_L5_All_hits_R, * M11_Od_L6_All_hits_R;
  // TH1F * M11_Ev_L1_All_hits_R, * M11_Ev_L2_All_hits_R, * M11_Ev_L3_All_hits_R, * M11_Ev_L4_All_hits_R, * M11_Ev_L5_All_hits_R, * M11_Ev_L6_All_hits_R;


  // Separation [0 - 250ns] and [250ns - beyond]
  TH1F * ME0_el_hits_000_R,  * ME0_mu_hits_000_R,  * ME0_pi_hits_000_R,  * ME0_ka_hits_000_R,  * ME0_p_hits_000_R,  * ME0_n_hits_000_R,  * ME0_g_hits_000_R,  * ME0_N_hits_000_R, * ME0_OH_hits_000_R;
  TH1F * ME0_el_hits_250_R,  * ME0_mu_hits_250_R,  * ME0_pi_hits_250_R,  * ME0_ka_hits_250_R,  * ME0_p_hits_250_R,  * ME0_n_hits_250_R,  * ME0_g_hits_250_R,  * ME0_N_hits_250_R, * ME0_OH_hits_250_R;
  TH1F * G11_el_hits_000_R,  * G11_mu_hits_000_R,  * G11_pi_hits_000_R,  * G11_ka_hits_000_R,  * G11_p_hits_000_R,  * G11_n_hits_000_R,  * G11_g_hits_000_R,  * G11_N_hits_000_R, * G11_OH_hits_000_R;
  TH1F * G11_el_hits_250_R,  * G11_mu_hits_250_R,  * G11_pi_hits_250_R,  * G11_ka_hits_250_R,  * G11_p_hits_250_R,  * G11_n_hits_250_R,  * G11_g_hits_250_R,  * G11_N_hits_250_R, * G11_OH_hits_250_R;
  TH1F * G21_el_hits_000_R,  * G21_mu_hits_000_R,  * G21_pi_hits_000_R,  * G21_ka_hits_000_R,  * G21_p_hits_000_R,  * G21_n_hits_000_R,  * G21_g_hits_000_R,  * G21_N_hits_000_R, * G21_OH_hits_000_R;
  TH1F * G21_el_hits_250_R,  * G21_mu_hits_250_R,  * G21_pi_hits_250_R,  * G21_ka_hits_250_R,  * G21_p_hits_250_R,  * G21_n_hits_250_R,  * G21_g_hits_250_R,  * G21_N_hits_250_R, * G21_OH_hits_250_R;

  // Separation [0 - 25ns] and [25ns - beyond]
  TH1F * ME0_el_hits_00_R,  * ME0_mu_hits_00_R,  * ME0_pi_hits_00_R,  * ME0_ka_hits_00_R,  * ME0_p_hits_00_R,  * ME0_n_hits_00_R,  * ME0_g_hits_00_R,  * ME0_N_hits_00_R, * ME0_OH_hits_00_R;
  TH1F * ME0_el_hits_25_R,  * ME0_mu_hits_25_R,  * ME0_pi_hits_25_R,  * ME0_ka_hits_25_R,  * ME0_p_hits_25_R,  * ME0_n_hits_25_R,  * ME0_g_hits_25_R,  * ME0_N_hits_25_R, * ME0_OH_hits_25_R;
  TH1F * G11_el_hits_00_R,  * G11_mu_hits_00_R,  * G11_pi_hits_00_R,  * G11_ka_hits_00_R,  * G11_p_hits_00_R,  * G11_n_hits_00_R,  * G11_g_hits_00_R,  * G11_N_hits_00_R, * G11_OH_hits_00_R;
  TH1F * G11_el_hits_25_R,  * G11_mu_hits_25_R,  * G11_pi_hits_25_R,  * G11_ka_hits_25_R,  * G11_p_hits_25_R,  * G11_n_hits_25_R,  * G11_g_hits_25_R,  * G11_N_hits_25_R, * G11_OH_hits_25_R;
  TH1F * G21_el_hits_00_R,  * G21_mu_hits_00_R,  * G21_pi_hits_00_R,  * G21_ka_hits_00_R,  * G21_p_hits_00_R,  * G21_n_hits_00_R,  * G21_g_hits_00_R,  * G21_N_hits_00_R, * G21_OH_hits_00_R;
  TH1F * G21_el_hits_25_R,  * G21_mu_hits_25_R,  * G21_pi_hits_25_R,  * G21_ka_hits_25_R,  * G21_p_hits_25_R,  * G21_n_hits_25_R,  * G21_g_hits_25_R,  * G21_N_hits_25_R, * G21_OH_hits_25_R;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 


  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_etapart; // * TDir_Muon_250ns_etapart, * TDir_Muon_25ns_hits_etapart;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * ME0_el_hits_E,  * ME0_mu_hits_E,  * ME0_pi_hits_E,  * ME0_ka_hits_E,  * ME0_p_hits_E,  * ME0_n_hits_E,  * ME0_g_hits_E,  * ME0_N_hits_E, * ME0_OH_hits_E;
  TH1F * G11_el_hits_E,  * G11_mu_hits_E,  * G11_pi_hits_E,  * G11_ka_hits_E,  * G11_p_hits_E,  * G11_n_hits_E,  * G11_g_hits_E,  * G11_N_hits_E, * G11_OH_hits_E;
  TH1F * G21_el_hits_E,  * G21_mu_hits_E,  * G21_pi_hits_E,  * G21_ka_hits_E,  * G21_p_hits_E,  * G21_n_hits_E,  * G21_g_hits_E,  * G21_N_hits_E, * G21_OH_hits_E;

  TH1F * ME0_All_hits_E, * ME0_All_hits_000_E, * ME0_All_hits_250_E, * ME0_All_hits_00_E, * ME0_All_hits_25_E, * ME0_HIP_hits_E; 
  TH1F * G11_All_hits_E, * G11_All_hits_000_E, * G11_All_hits_250_E, * G11_All_hits_00_E, * G11_All_hits_25_E, * G11_HIP_hits_E; 
  TH1F * G21_All_hits_E, * G21_All_hits_000_E, * G21_All_hits_250_E, * G21_All_hits_00_E, * G21_All_hits_25_E, * G21_HIP_hits_E; 

  TH1F * G11_L1Odd_el_hits_E,  * G11_L1Odd_mu_hits_E,   * G11_L1Odd_pi_hits_E,   * G11_L1Odd_ka_hits_E,   * G11_L1Odd_p_hits_E,    * G11_L1Odd_n_hits_E;
  TH1F * G11_L1Odd_g_hits_E,   * G11_L1Odd_N_hits_E,    * G11_L1Odd_OH_hits_E,   * G11_L1Odd_All_hits_E,  * G11_L1Odd_HIP_hits_E;
  TH1F * G11_L1Even_el_hits_E, * G11_L1Even_mu_hits_E,  * G11_L1Even_pi_hits_E,  * G11_L1Even_ka_hits_E,  * G11_L1Even_p_hits_E,   * G11_L1Even_n_hits_E;
  TH1F * G11_L1Even_g_hits_E,  * G11_L1Even_N_hits_E,   * G11_L1Even_OH_hits_E,  * G11_L1Even_All_hits_E, * G11_L1Even_HIP_hits_E;
  TH1F * G11_L2Odd_el_hits_E,  * G11_L2Odd_mu_hits_E,   * G11_L2Odd_pi_hits_E,   * G11_L2Odd_ka_hits_E,   * G11_L2Odd_p_hits_E,    * G11_L2Odd_n_hits_E;
  TH1F * G11_L2Odd_g_hits_E,   * G11_L2Odd_N_hits_E,    * G11_L2Odd_OH_hits_E,   * G11_L2Odd_All_hits_E,  * G11_L2Odd_HIP_hits_E;  
  TH1F * G11_L2Even_el_hits_E, * G11_L2Even_mu_hits_E,  * G11_L2Even_pi_hits_E,  * G11_L2Even_ka_hits_E,  * G11_L2Even_p_hits_E,   * G11_L2Even_n_hits_E;
  TH1F * G11_L2Even_g_hits_E,  * G11_L2Even_N_hits_E,   * G11_L2Even_OH_hits_E,  * G11_L2Even_All_hits_E, * G11_L2Even_HIP_hits_E;
  



  // Separation [0 - 250ns] and [250ns - beyond]
  // TH1F * ME0_el_hits_000_E,  * ME0_mu_hits_000_E,  * ME0_pi_hits_000_E,  * ME0_ka_hits_000_E,  * ME0_p_hits_000_E,  * ME0_n_hits_000_E,  * ME0_g_hits_000_E,  * ME0_N_hits_000_E, * ME0_OH_hits_000_E;
  // TH1F * ME0_el_hits_250_E,  * ME0_mu_hits_250_E,  * ME0_pi_hits_250_E,  * ME0_ka_hits_250_E,  * ME0_p_hits_250_E,  * ME0_n_hits_250_E,  * ME0_g_hits_250_E,  * ME0_N_hits_250_E, * ME0_OH_hits_250_E;
  // TH1F * G11_el_hits_000_E,  * G11_mu_hits_000_E,  * G11_pi_hits_000_E,  * G11_ka_hits_000_E,  * G11_p_hits_000_E,  * G11_n_hits_000_E,  * G11_g_hits_000_E,  * G11_N_hits_000_E, * G11_OH_hits_000_E;
  // TH1F * G11_el_hits_250_E,  * G11_mu_hits_250_E,  * G11_pi_hits_250_E,  * G11_ka_hits_250_E,  * G11_p_hits_250_E,  * G11_n_hits_250_E,  * G11_g_hits_250_E,  * G11_N_hits_250_E, * G11_OH_hits_250_E;
  // TH1F * G21_el_hits_000_E,  * G21_mu_hits_000_E,  * G21_pi_hits_000_E,  * G21_ka_hits_000_E,  * G21_p_hits_000_E,  * G21_n_hits_000_E,  * G21_g_hits_000_E,  * G21_N_hits_000_E, * G21_OH_hits_000_E;
  // TH1F * G21_el_hits_250_E,  * G21_mu_hits_250_E,  * G21_pi_hits_250_E,  * G21_ka_hits_250_E,  * G21_p_hits_250_E,  * G21_n_hits_250_E,  * G21_g_hits_250_E,  * G21_N_hits_250_E, * G21_OH_hits_250_E;

  // Separation [0 - 25ns] and [25ns - beyond]
  // TH1F * ME0_el_hits_00_E,  * ME0_mu_hits_00_E,  * ME0_pi_hits_00_E,  * ME0_ka_hits_00_E,  * ME0_p_hits_00_E,  * ME0_n_hits_00_E,  * ME0_g_hits_00_E,  * ME0_N_hits_00_E, * ME0_OH_hits_00_E;
  // TH1F * ME0_el_hits_25_E,  * ME0_mu_hits_25_E,  * ME0_pi_hits_25_E,  * ME0_ka_hits_25_E,  * ME0_p_hits_25_E,  * ME0_n_hits_25_E,  * ME0_g_hits_25_E,  * ME0_N_hits_25_E, * ME0_OH_hits_25_E;
  // TH1F * G11_el_hits_00_E,  * G11_mu_hits_00_E,  * G11_pi_hits_00_E,  * G11_ka_hits_00_E,  * G11_p_hits_00_E,  * G11_n_hits_00_E,  * G11_g_hits_00_E,  * G11_N_hits_00_E, * G11_OH_hits_00_E;
  // TH1F * G11_el_hits_25_E,  * G11_mu_hits_25_E,  * G11_pi_hits_25_E,  * G11_ka_hits_25_E,  * G11_p_hits_25_E,  * G11_n_hits_25_E,  * G11_g_hits_25_E,  * G11_N_hits_25_E, * G11_OH_hits_25_E;
  // TH1F * G21_el_hits_00_E,  * G21_mu_hits_00_E,  * G21_pi_hits_00_E,  * G21_ka_hits_00_E,  * G21_p_hits_00_E,  * G21_n_hits_00_E,  * G21_g_hits_00_E,  * G21_N_hits_00_E, * G21_OH_hits_00_E;
  // TH1F * G21_el_hits_25_E,  * G21_mu_hits_25_E,  * G21_pi_hits_25_E,  * G21_ka_hits_25_E,  * G21_p_hits_25_E,  * G21_n_hits_25_E,  * G21_g_hits_25_E,  * G21_N_hits_25_E, * G21_OH_hits_25_E;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 


  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_hits_process;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH1F * GEM_el_process,  * GEM_mu_process,  * GEM_pi_process,  * GEM_ka_process,  * GEM_p_process,  * GEM_n_process,  * GEM_g_process,  * GEM_N_process, * GEM_OH_process;
  TH1F * ME0_el_process,  * ME0_mu_process,  * ME0_pi_process,  * ME0_ka_process,  * ME0_p_process,  * ME0_n_process,  * ME0_g_process,  * ME0_N_process, * ME0_OH_process;
  TH1F * GEM_el_proc000,  * GEM_mu_proc000,  * GEM_pi_proc000,  * GEM_ka_proc000,  * GEM_p_proc000,  * GEM_n_proc000,  * GEM_g_proc000,  * GEM_N_proc000, * GEM_OH_proc000;
  TH1F * GEM_el_proc250,  * GEM_mu_proc250,  * GEM_pi_proc250,  * GEM_ka_proc250,  * GEM_p_proc250,  * GEM_n_proc250,  * GEM_g_proc250,  * GEM_N_proc250, * GEM_OH_proc250;
  TH1F * ME0_el_proc000,  * ME0_mu_proc000,  * ME0_pi_proc000,  * ME0_ka_proc000,  * ME0_p_proc000,  * ME0_n_proc000,  * ME0_g_proc000,  * ME0_N_proc000, * ME0_OH_proc000;
  TH1F * ME0_el_proc250,  * ME0_mu_proc250,  * ME0_pi_proc250,  * ME0_ka_proc250,  * ME0_p_proc250,  * ME0_n_proc250,  * ME0_g_proc250,  * ME0_N_proc250, * ME0_OH_proc250;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  std::vector< TH1F* > VTH1F_Muon_hits_process;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 




  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
  TDirectoryFile * TDir_Muon_hits_2Dplots;
  TDirectoryFile * TDir_Muon_hits_deposits;
  TDirectoryFile * TDir_Muon_hits_kinetics;
  TDirectoryFile * TDir_Muon_hits_timing;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  // 2D plot of Kinetic Energy of Particles vs TOF
  TH2F * GEM_el_hits,  * GEM_mu_hits,  * GEM_pi_hits,  * GEM_ka_hits,  * GEM_p_hits,  * GEM_n_hits,  * GEM_g_hits,  * GEM_N_hits, * GEM_OH_hits, * GEM_ha_hits;
  TH2F * ME0_el_hits,  * ME0_mu_hits,  * ME0_pi_hits,  * ME0_ka_hits,  * ME0_p_hits,  * ME0_n_hits,  * ME0_g_hits,  * ME0_N_hits, * ME0_OH_hits, * ME0_ha_hits;
  // 2D plot of Energy Deposit of Particles vs TOF
  TH2F * GEM_el_deposits,  * GEM_mu_deposits,  * GEM_pi_deposits,  * GEM_ka_deposits,  * GEM_p_deposits,  * GEM_n_deposits,  * GEM_g_deposits,  * GEM_N_deposits,  *GEM_OH_deposits;
  TH2F * ME0_el_deposits,  * ME0_mu_deposits,  * ME0_pi_deposits,  * ME0_ka_deposits,  * ME0_p_deposits,  * ME0_n_deposits,  * ME0_g_deposits,  * ME0_N_deposits,  *ME0_OH_deposits;
  // 2D plot of Kinetic Energy vs Energy Deposit
  TH2F * GEM_el_ekindep,  * GEM_mu_ekindep,  * GEM_pi_ekindep,  * GEM_ka_ekindep,  * GEM_p_ekindep,  * GEM_n_ekindep,  * GEM_g_ekindep,  * GEM_N_ekindep, *GEM_OH_ekindep;
  TH2F * ME0_el_ekindep,  * ME0_mu_ekindep,  * ME0_pi_ekindep,  * ME0_ka_ekindep,  * ME0_p_ekindep,  * ME0_n_ekindep,  * ME0_g_ekindep,  * ME0_N_ekindep, *ME0_OH_ekindep;

  // 1D plot of Kinetic Energy
  TH1F * GEM_el_kins,  * GEM_mu_kins,  * GEM_pi_kins,  * GEM_ka_kins,  * GEM_p_kins,  * GEM_n_kins,  * GEM_g_kins,  * GEM_N_kins, * GEM_OH_kins, * GEM_ha_kins, * GEM_All_kins, * GEM_HIP_kins;   
  TH1F * ME0_el_kins,  * ME0_mu_kins,  * ME0_pi_kins,  * ME0_ka_kins,  * ME0_p_kins,  * ME0_n_kins,  * ME0_g_kins,  * ME0_N_kins, * ME0_OH_kins, * ME0_ha_kins, * ME0_All_kins, * ME0_HIP_kins;
  TH1F * GEM_el_linkin,  * GEM_mu_linkin,  * GEM_pi_linkin,  * GEM_ka_linkin,  * GEM_p_linkin,  * GEM_n_linkin,  * GEM_g_linkin,  * GEM_N_linkin, * GEM_OH_linkin, * GEM_ha_linkin, * GEM_All_linkin, * GEM_HIP_linkin;   
  TH1F * ME0_el_linkin,  * ME0_mu_linkin,  * ME0_pi_linkin,  * ME0_ka_linkin,  * ME0_p_linkin,  * ME0_n_linkin,  * ME0_g_linkin,  * ME0_N_linkin, * ME0_OH_linkin, * ME0_ha_linkin, * ME0_All_linkin, * ME0_HIP_linkin;
  
  TH1F * G11_el_kins,    * G11_mu_kins,    * G11_pi_kins,    * G11_ka_kins,    * G11_p_kins,    * G11_n_kins,    * G11_g_kins,    * G11_N_kins,   * G11_ha_kins,   *G11_OH_kins,   * G11_All_kins,   *G11_HIP_kins;  
  TH1F * G21_el_kins,    * G21_mu_kins,    * G21_pi_kins,    * G21_ka_kins,    * G21_p_kins,    * G21_n_kins,    * G21_g_kins,    * G21_N_kins,   * G21_ha_kins,   *G21_OH_kins,   * G21_All_kins,   *G21_HIP_kins;  
  TH1F * G11_el_linkin,  * G11_mu_linkin,  * G11_pi_linkin,  * G11_ka_linkin,  * G11_p_linkin,  * G11_n_linkin,  * G11_g_linkin,  * G11_N_linkin, * G11_ha_linkin, *G11_OH_linkin, * G11_All_linkin, *G11_HIP_linkin;  
  TH1F * G21_el_linkin,  * G21_mu_linkin,  * G21_pi_linkin,  * G21_ka_linkin,  * G21_p_linkin,  * G21_n_linkin,  * G21_g_linkin,  * G21_N_linkin, * G21_ha_linkin, *G21_OH_linkin, * G21_All_linkin, *G21_HIP_linkin;  

  // 1D plot of Energy Deposit
  TH1F * ME0_el_deps,  * ME0_mu_deps,  * ME0_pi_deps,  * ME0_ka_deps,  * ME0_p_deps,  * ME0_n_deps,  * ME0_g_deps,  * ME0_N_deps, * ME0_ha_deps, *ME0_OH_deps, *ME0_All_deps, *ME0_HIP_deps;  
  TH1F * ME0_el_deps_000,  * ME0_mu_deps_000,  * ME0_pi_deps_000,  * ME0_ka_deps_000,  * ME0_p_deps_000,  * ME0_n_deps_000,  * ME0_g_deps_000,  * ME0_N_deps_000, * ME0_ha_deps_000, *ME0_OH_deps_000;  
  TH1F * ME0_el_deps_250,  * ME0_mu_deps_250,  * ME0_pi_deps_250,  * ME0_ka_deps_250,  * ME0_p_deps_250,  * ME0_n_deps_250,  * ME0_g_deps_250,  * ME0_N_deps_250, * ME0_ha_deps_250, *ME0_OH_deps_250;  
  TH1F * ME0_el_lindep,    * ME0_mu_lindep,    * ME0_pi_lindep,    * ME0_ka_lindep,    * ME0_p_lindep,    * ME0_n_lindep,    * ME0_g_lindep,    * ME0_N_lindep,   * ME0_ha_lindep,   *ME0_OH_lindep;  
  TH1F * ME0_All_lindep,   * ME0_HIP_lindep;

  TH1F * G11_el_deps,    * G11_mu_deps,    * G11_pi_deps,    * G11_ka_deps,    * G11_p_deps,    * G11_n_deps,    * G11_g_deps,    * G11_N_deps,   * G11_ha_deps,   *G11_OH_deps,   * G11_All_deps,   *G11_HIP_deps;  
  TH1F * G21_el_deps,    * G21_mu_deps,    * G21_pi_deps,    * G21_ka_deps,    * G21_p_deps,    * G21_n_deps,    * G21_g_deps,    * G21_N_deps,   * G21_ha_deps,   *G21_OH_deps,   * G21_All_deps,   *G21_HIP_deps;  
  TH1F * G11_el_lindep,  * G11_mu_lindep,  * G11_pi_lindep,  * G11_ka_lindep,  * G11_p_lindep,  * G11_n_lindep,  * G11_g_lindep,  * G11_N_lindep, * G11_ha_lindep, *G11_OH_lindep, * G11_All_lindep, *G11_HIP_lindep;  
  TH1F * G21_el_lindep,  * G21_mu_lindep,  * G21_pi_lindep,  * G21_ka_lindep,  * G21_p_lindep,  * G21_n_lindep,  * G21_g_lindep,  * G21_N_lindep, * G21_ha_lindep, *G21_OH_lindep, * G21_All_lindep, *G21_HIP_lindep;  

  TH2F * ME0_All_lindep_roll,  * G11_All_lindep_roll,  * G21_All_lindep_roll;
  TH1F * ME0_All_lindep_eta01, * ME0_All_lindep_eta02, * ME0_All_lindep_eta03, * ME0_All_lindep_eta04, * ME0_All_lindep_eta05, * ME0_All_lindep_eta06, * ME0_All_lindep_eta07, * ME0_All_lindep_eta08;
  TH1F * G11_All_lindep_eta01, * G11_All_lindep_eta02, * G11_All_lindep_eta03, * G11_All_lindep_eta04, * G11_All_lindep_eta05, * G11_All_lindep_eta06, * G11_All_lindep_eta07, * G11_All_lindep_eta08;
  TH1F * G21_All_lindep_eta01, * G21_All_lindep_eta02, * G21_All_lindep_eta03, * G21_All_lindep_eta04, * G21_All_lindep_eta05, * G21_All_lindep_eta06, * G21_All_lindep_eta07, * G21_All_lindep_eta08;

  // TH1F * G11_el_deps_000,  * G11_mu_deps_000,  * G11_pi_deps_000,  * G11_ka_deps_000,  * G11_p_deps_000,  * G11_n_deps_000,  * G11_g_deps_000,  * G11_N_deps_000, * G11_ha_deps_000, *G11_OH_deps_000;  
  // TH1F * G11_el_deps_250,  * G11_mu_deps_250,  * G11_pi_deps_250,  * G11_ka_deps_250,  * G11_p_deps_250,  * G11_n_deps_250,  * G11_g_deps_250,  * G11_N_deps_250, * G11_ha_deps_250, *G11_OH_deps_250;  
  // TH1F * G21_el_deps_000,  * G21_mu_deps_000,  * G21_pi_deps_000,  * G21_ka_deps_000,  * G21_p_deps_000,  * G21_n_deps_000,  * G21_g_deps_000,  * G21_N_deps_000, * G21_ha_deps_000, *G21_OH_deps_000;  
  // TH1F * G21_el_deps_250,  * G21_mu_deps_250,  * G21_pi_deps_250,  * G21_ka_deps_250,  * G21_p_deps_250,  * G21_n_deps_250,  * G21_g_deps_250,  * G21_N_deps_250, * G21_ha_deps_250, *G21_OH_deps_250;  
 
  TH1F * GEM_el_deps,      * GEM_mu_deps,      * GEM_pi_deps,  * GEM_ka_deps,  * GEM_p_deps,  * GEM_n_deps,  * GEM_g_deps,  * GEM_N_deps, * GEM_ha_deps, *GEM_OH_deps, *GEM_All_deps, *GEM_HIP_deps;  
  TH1F * GEM_el_deps_000,  * GEM_mu_deps_000,  * GEM_pi_deps_000,  * GEM_ka_deps_000,  * GEM_p_deps_000,  * GEM_n_deps_000,  * GEM_g_deps_000,  * GEM_N_deps_000, * GEM_ha_deps_000, *GEM_OH_deps_000;  
  TH1F * GEM_el_deps_250,  * GEM_mu_deps_250,  * GEM_pi_deps_250,  * GEM_ka_deps_250,  * GEM_p_deps_250,  * GEM_n_deps_250,  * GEM_g_deps_250,  * GEM_N_deps_250, * GEM_ha_deps_250, *GEM_OH_deps_250;  
  TH1F * GEM_el_lindep,    * GEM_mu_lindep,    * GEM_pi_lindep,    * GEM_ka_lindep,    * GEM_p_lindep,    * GEM_n_lindep,    * GEM_g_lindep,    * GEM_N_lindep,   * GEM_ha_lindep,   *GEM_OH_lindep;  
  TH1F * GEM_All_lindep,   * GEM_HIP_lindep;

  // 1D plot of TOF (log & lin)
  TH1F * GEM_el_tof,  * GEM_mu_tof,   * GEM_pi_tof,  * GEM_ka_tof,  * GEM_p_tof,  * GEM_n_tof,  * GEM_g_tof,  * GEM_N_tof,  * GEM_ha_tof,  * GEM_OH_tof;
  TH1F * ME0_el_tof,  * ME0_mu_tof,   * ME0_pi_tof,  * ME0_ka_tof,  * ME0_p_tof,  * ME0_n_tof,  * ME0_g_tof,  * ME0_N_tof,  * ME0_ha_tof,  * ME0_OH_tof;
  TH1F * GEM_el_time, * GEM_mu_time,  * GEM_pi_time, * GEM_ka_time, * GEM_p_time, * GEM_n_time, * GEM_g_time, * GEM_N_time, * GEM_ha_time, * GEM_OH_time;
  TH1F * ME0_el_time, * ME0_mu_time,  * ME0_pi_time, * ME0_ka_time, * ME0_p_time, * ME0_n_time, * ME0_g_time, * ME0_N_time, * ME0_ha_time, * ME0_OH_time;
  // Old 
  TH1F * GEM_hits_tof, * ME0_hits_tof, * GEM_hits_eta, * ME0_hits_eta, * GEM_hits_phi, * ME0_hits_phi, * GEM_hits_lin, * ME0_hits_lin;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TCanvas * Canvas_GEM_hits,       * Canvas_ME0_hits,       * Canvas_GEM_deposits,       * Canvas_ME0_deposits,       * Canvas_GEM_ekindep,       * Canvas_ME0_ekindep; 
  TCanvas * Canvas_GEM_1D_tof,     * Canvas_ME0_1D_tof,     * Canvas_GEM_1D_deps,        * Canvas_ME0_1D_deps;
  // TCanvas * Canvas_GEM_hits_fancy, * Canvas_ME0_hits_fancy, * Canvas_GEM_deposits_fancy, * Canvas_ME0_deposits_fancy;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 


  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TDirectoryFile * TDir_Muon_XY_RZ_views;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  TH2F * GEM_XY, *GEM_RZ, * ME0_XY, *ME0_RZ;
  TH2F * GEM_000ns_RZ, * ME0_000ns_RZ;
  TH2F * GEM_250ns_RZ, * ME0_250ns_RZ;
  TH2F * GEM_000ns_XY, * ME0_000ns_XY;
  TH2F * GEM_250ns_XY, * ME0_250ns_XY;
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------------- 

  // Nuclear codes are given as 10-digit numbers ±10LZZZAAAI. For a (hyper) nucleus consisting of n_p protons, n_n neutrons and n_Λ Λ’s, A=n_p+n_n+n_Λ gives the total 
  // baryon number,Z=n_p the total charge and L=n_Λ the total number of strange quarks. I gives the isomer level, with I= 0 corresponding to the ground state and I >0
  // to excitations, where states denoted m,n,p,q translate to I= 1−4. As examples, the deuteron is 1000010020 and 235U is 1000922350.
  // Examples:
  // - deuteron (Z = 1 A = 2) 1000010020 || 1414 su 100k || 8044 total
  // - tritium  (Z = 1 A = 3) 1000010030 ||  108 su 100k ||
  // - helium   (Z = 2 A = 3) 1000020030 ||   12 su 100k ||
  // -          (Z = 2 A = 4) 1000020040 ||   67 su 100k ||
  // - berillium(Z = 4 A = 8) 1000040080 ||    6 su 100k ||
  // - carbon   (Z = 6 A =12) 1000060120 || 1573 su 100k ||
  //            (Z = 6 A =13) 1000060130 ||   12 su 100k ||
  //            (Z = 6 A =14) 1000060140 ||    4 su 100k ||
  // - nitrogen (Z = 7 A =14) 1000070140 ||    0 su 100k ||
  //            (Z = 7 A =15) 1000070150 ||    2 su 100k ||
  // - oxygen   (Z = 8 A =16) 1000080160 || 1701 su 100k ||
  //            (Z = 8 A =17) 1000080170 ||    3 su 100k ||
  //            (Z = 8 A =18) 1000080180 ||    6 su 100k ||
  // - fluor    (Z = 9 A =19) 1000090190 || 2028 su 100k ||
  // - magnesium(Z =12 A =26) 1000120260 ||    2 su 100k ||
  // - aluminum (Z =13 A =27) 1000130270 ||    8 su 100k ||
  // - phospor  (Z =15 A =32) 1000150320 ||    2 su 100k ||
  // - sulfur   (Z =16 A =33) 1000160330 ||    1 su 100k ||
  //            (Z =16 A =35) 1000160350 ||    2 su 100k ||
  // - chlorine (Z =17 A =37) 1000170370 ||    2 su 100k ||
  // - argon    (Z =18 A =36) 1000180360 ||    2 su 100k ||
  // -          (Z =18 A =37) 1000180370 ||   85 su 100k ||
  // -          (Z =18 A =38) 1000180380 ||    2 su 100k ||
  // -          (Z =18 A =39) 1000180390 ||    5 su 100k ||
  // -          (Z =18 A =40) 1000180400 || 1010 su 100k ||

  TDirectoryFile * TDir_Muon_Nuclei;
  TH2F * ME0_Nuclei_A_Z;
  TH1F * ME0_Nuclei_A, * ME0_Nuclei_Z, * ME0_Nuclei_List;
  TH1F * GEM_HIP_id, * ME0_HIP_id;
};

//
// constants, enums and typedefs
//

int n_tof = 1100,  n1_tof = 1,  n2_tof = 12;
int m_tof = 110,   m1_tof = 1,  m2_tof = 12;
int n_time = 1000, n1_time = 0, n2_time = 1000;
int m_eta = 60;   double m1_eta =  0.0,  m2_eta = 3.0;
int m_phi = 36;   double m1_phi = -3.1415, m2_phi = 3.1415;
int m_lin = 1000, m1_lin = 0,  m2_lin = 100000000;
int n_hits = 20;  double n1_hits = 0.5, n2_hits = 20.5;
int n_lay  = 12;  double n1_lay  = 0.5, n2_lay  = 12.5;

int n_D   = 900,  n1_D   = -3, n2_D   = 6;
int n_E   = 900,  n1_E   = -4, n2_E   = 5;
// int n_F   = 90,   n1_F   = -3, n2_F   = 6;
int n_F   = 200,  n1_F   = -3, n2_F   = 6;
int n_R   = 90,   n1_R   = 60, n2_R   = 150;     // Radius ME0  :: 60 - 150 cm
int n_P   = 120,  n1_P   = 130,n2_P   = 250;     // Radius GE11 :: 130 - 250 cm
int n_Q   = 190,  n1_Q   = 130,n2_Q   = 320;     // Radius GE21 :: 130 - 320 cm
int n_N   = 8;    double n1_N  = 0.5, n2_N  = 8.5;  // Eta Partitions GE11 :: 8 
int n_N2  = 16;   double n1_N2 = 0.5, n2_N2 = 16.5; // Eta Partitions GE21 :: 16 

int n_Z   = 100,  n1_Z   =  1, n2_Z   = 100;
int n_A   = 200,  n1_A   = 1,  n2_A   = 200;
int n_DL  = 2500, n1_DL  = 0,  n2_DL  = 250;   // linear plot deposits: 0-250keV in steps of 100eV
int n_H   = 1000, n1_H   = 0,  n2_H   = 10000; // HIP deposits: 0-10MeV in steps of 5keV = steps of 167e-  
int n_EL  = 1000, n1_EL  = 0,  n2_EL  = 1000;  // linear kinetic energy 0 -1GeV in steps of 1 MeV
int n_CL  = 1000, n1_CL  = 0,  n2_CL  = 500;   // linear plot deposits: 0-500keV in steps of 500eV + 1 overflowbin

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
int n_pro  = 31; double n1_pro  = 0.5, n2_pro  = 31.5;

std::string evtcnt[7] = {"All", "ME0 Hits", "GE11 Hits", "GE21 Hits", "ME0 HIP", "GE11 HIP", "GE21 HIP"};

std::string allcat[27] = {"All", "Barrel", "W0_RB1", "W0_RB2", "W0_RB3", "W0_RB4", "W1_RB1", "W1_RB2", "W1_RB3", "W1_RB4", "W2_RB1", "W2_RB2", "W2_RB3", "W2_RB4", 
                                 "Endcap", "RE11", "RE12", "RE13", "RE21", "RE22", "RE23", "RE31", "RE32", "RE33", "RE41", "RE42", "RE43"};
std::string endcat[24] = {"RE12C", "RE12B", "RE12A", "RE13C", "RE13B", "RE13A", "RE22C", "RE22B", "RE22A", "RE23C", "RE23B", "RE23A", 
                          "RE32C", "RE32B", "RE32A", "RE33C", "RE33B", "RE33A", "RE42C", "RE42B", "RE42A", "RE43C", "RE43B", "RE43A"};
// std::string proc[15] = {"notDefined", "fCoulombScattering", "fIonisation", "fBremsstrahlung", "fPairProdByCharged", "fAnnihilation", "fAnnihilationToMuMu", "fAnnihilationToHadrons", 
//                         "fNuclearStopping", "notDefined", "fMultipleScattering", "fRayleigh", "fPhotoElectricEffect", "fComptonScattering", "fGammaConversion"};
std::string proc[31] = {"notDefined",
                        "CoulombScattering", "Ionisation",    "Bremsstrahlung",     "PairProdByCharged", "Annihilation",         "AnnihilationToMuMu", "AnnihilationToHadrons", 
                        "NuclearStopping",   "ElectronSuper", "MultipleScattering", "Rayleigh",          "PhotoElectricEffect",  "ComptonScattering",  "GammaConversion",
                        "GammaConvToMuMu",   "GammaSuper",    "Cherenkov",          "Scintillation",     "SynchrotronRadiation", "TransitionRadiation", "HadronElastic",
                        "HadronInelastic",   "HadronCapture", "HadronFission",      "HadronAtRest",      "HadronCEX",            "Decay",               "DecayWithSpin",
                        "DecayUnknown",      "Decay External"};
// :: Process Types :: 31 ::
// - - - - - - - - - - - - -
//   0 not Defined
//   1 Coulomb Scattering  2 Ionisation            3 Bremsstrahlung     4 PairProductionByCharged  5 Annihilation
//   6 AnnihilationToMuMu  7 AnnihilationToHad     8 NuclearStopping    9 ElectronSuper           10 MultipleScattering 
//  11 Rayleigh Scatter   12 PhotoElectricEffect  13 ComptonScattering 14 GammaConv               15 GammaConvToMuMu      16 GammaSuper
//  21 Cherenkov          22 Scintillation        23 Synchrotron Rad   24 Transition Radiation    
// 111 Hadron Elastic    121 Hadron Inelastic    131 Hadron Capture   141 Hadron Fission         151 Hadron At REst       161 HadronCEX
// 201 Decay             201 Decay with Spin     211 Decay Unknown    231 Decay External
// 

std::string nuclei[19] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "Other"};  
std::string particles[9] = {"el", "mu", "pi", "ka", "p", "n", "g", "N", "OH"};
//
// static data member definitions
//

//
// constructors and destructor
//
MyME0SimHitAnalyzer::MyME0SimHitAnalyzer(const edm::ParameterSet& iConfig)

{
  // now do what ever initialization is needed
  pdfFileNameBase = iConfig.getUntrackedParameter<std::string>("PdfFileNameBase");
  rootFileName    = iConfig.getUntrackedParameter<std::string>("RootFileName");
  bunchspacing    = iConfig.getUntrackedParameter<double>("BunchSpacing");
  comenergy       = iConfig.getUntrackedParameter<double>("COMEnergy");
  maxsimtime      = iConfig.getUntrackedParameter<double>("MaxSimTime");
  // gem_only_ge11   = iConfig.getUntrackedParameter<bool>("GEMOnlyGE11");
  edep_30eV       = iConfig.getUntrackedParameter<bool>("Edep30eV");
  edepMin         = iConfig.getUntrackedParameter<double>("MinHitEnergy");
  hipMin          = iConfig.getUntrackedParameter<double>("MinHIPEnergy");
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

  if(tech_debug) std::cout<<"[MyME0SimHitAnalyzer :: Constructor]"<<std::endl;
  
  EventCounter  = new TH1F("EventCounter", "EventCounter",   7, 0.5, 7.5);
  // HIPParticleId = new TH1F("HIPParticleId", "HIPParticleId", 9, 0.5, 9.5);


  TDir_Muon_hits_radius   = new TDirectoryFile("Muon_hits_radius",   "Muon_hits_radius");
  TDir_Muon_250ns_radius  = new TDirectoryFile("Muon_250ns_radius",  "Muon_250ns_radius");
  TDir_Muon_25ns_radius   = new TDirectoryFile("Muon_25ns_radius",   "Muon_25ns_radius");
  TDir_Muon_hits_etapart  = new TDirectoryFile("Muon_hits_etapart",  "Muon_hits_etapart");
  TDir_Muon_hits_process  = new TDirectoryFile("Muon_hits_process",  "Muon_hits_process");
  TDir_Muon_hits_2Dplots  = new TDirectoryFile("Muon_hits_2Dplots",  "Muon_hits_2Dplotss");
  TDir_Muon_hits_deposits = new TDirectoryFile("Muon_hits_deposits", "Muon_hits_deposits");
  TDir_Muon_hits_kinetics = new TDirectoryFile("Muon_hits_kinetics", "Muon_hits_kinematics");
  TDir_Muon_hits_timing   = new TDirectoryFile("Muon_hits_timing",   "Muon_hits_timing");
  TDir_Muon_XY_RZ_views   = new TDirectoryFile("Muon_XY_RZ_views",   "Muon_XY_RZ_views");
  TDir_Muon_Nuclei        = new TDirectoryFile("Muon_hits_Nuclei",   "Muon_hits_Nuclei");


  // Simhit vs R :: in future need to split between GE11 short & GE11 long
  G11_el_hits_R = new TH1F("G11_el_hits_R", "Simhit Radius :: Electrons", n_P, n1_P, n2_P);
  G11_mu_hits_R = new TH1F("G11_mu_hits_R", "Simhit Radius :: Muons",     n_P, n1_P, n2_P);
  G11_pi_hits_R = new TH1F("G11_pi_hits_R", "Simhit Radius :: Pions",     n_P, n1_P, n2_P);
  G11_ka_hits_R = new TH1F("G11_ka_hits_R", "Simhit Radius :: Kaons",     n_P, n1_P, n2_P);
  G11_p_hits_R  = new TH1F("G11_p_hits_R",  "Simhit Radius :: Protons",   n_P, n1_P, n2_P);
  G11_n_hits_R  = new TH1F("G11_n_hits_R",  "Simhit Radius :: Neutrons",  n_P, n1_P, n2_P);
  G11_g_hits_R  = new TH1F("G11_g_hits_R",  "Simhit Radius :: Photons",   n_P, n1_P, n2_P);
  G11_N_hits_R  = new TH1F("G11_N_hits_R",  "Simhit Radius :: Nuclei",    n_P, n1_P, n2_P);
  G11_OH_hits_R = new TH1F("G11_OH_hits_R", "Simhit Radius :: Other Hadrons", n_P, n1_P, n2_P);
  G11_All_hits_R= new TH1F("G11_All_hits_R","Simhit Radius :: All Particles", n_P, n1_P, n2_P);
  G11_HIP_hits_R= new TH1F("G11_HIP_hits_R","Simhit Radius :: Highly Ionizing", n_P, n1_P, n2_P);

  G21_el_hits_R = new TH1F("G21_el_hits_R", "Simhit Radius :: Electrons", n_Q, n1_Q, n2_Q);
  G21_mu_hits_R = new TH1F("G21_mu_hits_R", "Simhit Radius :: Muons",     n_Q, n1_Q, n2_Q);
  G21_pi_hits_R = new TH1F("G21_pi_hits_R", "Simhit Radius :: Pions",     n_Q, n1_Q, n2_Q);
  G21_ka_hits_R = new TH1F("G21_ka_hits_R", "Simhit Radius :: Kaons",     n_Q, n1_Q, n2_Q);
  G21_p_hits_R  = new TH1F("G21_p_hits_R",  "Simhit Radius :: Protons",   n_Q, n1_Q, n2_Q);
  G21_n_hits_R  = new TH1F("G21_n_hits_R",  "Simhit Radius :: Neutrons",  n_Q, n1_Q, n2_Q);
  G21_g_hits_R  = new TH1F("G21_g_hits_R",  "Simhit Radius :: Photons",   n_Q, n1_Q, n2_Q);
  G21_N_hits_R  = new TH1F("G21_N_hits_R",  "Simhit Radius :: Nuclei",    n_Q, n1_Q, n2_Q);
  G21_OH_hits_R = new TH1F("G21_OH_hits_R", "Simhit Radius :: Other Hadrons", n_Q, n1_Q, n2_Q);
  G21_All_hits_R= new TH1F("G21_All_hits_R","Simhit Radius :: All Particles", n_Q, n1_Q, n2_Q);
  G21_HIP_hits_R= new TH1F("G21_HIP_hits_R","Simhit Radius :: Highly Ionizing", n_Q, n1_Q, n2_Q);

  ME0_el_hits_R = new TH1F("ME0_el_hits_R", "Simhit Radius :: Electrons", n_R, n1_R, n2_R);
  ME0_mu_hits_R = new TH1F("ME0_mu_hits_R", "Simhit Radius :: Muons",     n_R, n1_R, n2_R);
  ME0_pi_hits_R = new TH1F("ME0_pi_hits_R", "Simhit Radius :: Pions",     n_R, n1_R, n2_R);
  ME0_ka_hits_R = new TH1F("ME0_ka_hits_R", "Simhit Radius :: Kaons",     n_R, n1_R, n2_R);
  ME0_p_hits_R  = new TH1F("ME0_p_hits_R",  "Simhit Radius :: Protons",   n_R, n1_R, n2_R);
  ME0_n_hits_R  = new TH1F("ME0_n_hits_R",  "Simhit Radius :: Neutrons",  n_R, n1_R, n2_R);
  ME0_g_hits_R  = new TH1F("ME0_g_hits_R",  "Simhit Radius :: Photons",   n_R, n1_R, n2_R);
  ME0_N_hits_R  = new TH1F("ME0_N_hits_R",  "Simhit Radius :: Nuclei",    n_R, n1_R, n2_R);
  ME0_OH_hits_R = new TH1F("ME0_OH_hits_R", "Simhit Radius :: Other Hadrons", n_R, n1_R, n2_R);
  ME0_All_hits_R= new TH1F("ME0_All_hits_R","Simhit Radius :: All Particles", n_R, n1_R, n2_R);
  ME0_HIP_hits_R= new TH1F("ME0_HIP_hits_R","Simhit Radius :: Highly Ionizing", n_R, n1_R, n2_R);

  ME0_el_hits_000_R = new TH1F("ME0_el_hits_000_R", "Simhit Radius :: Electrons :: t < 250 ns", n_R, n1_R, n2_R);
  ME0_mu_hits_000_R = new TH1F("ME0_mu_hits_000_R", "Simhit Radius :: Muons :: t < 250 ns",     n_R, n1_R, n2_R);
  ME0_pi_hits_000_R = new TH1F("ME0_pi_hits_000_R", "Simhit Radius :: Pions :: t < 250 ns",     n_R, n1_R, n2_R);
  ME0_ka_hits_000_R = new TH1F("ME0_ka_hits_000_R", "Simhit Radius :: Kaons :: t < 250 ns",     n_R, n1_R, n2_R);
  ME0_p_hits_000_R  = new TH1F("ME0_p_hits_000_R",  "Simhit Radius :: Protons :: t < 250 ns",   n_R, n1_R, n2_R);
  ME0_n_hits_000_R  = new TH1F("ME0_n_hits_000_R",  "Simhit Radius :: Neutrons :: t < 250 ns",  n_R, n1_R, n2_R);
  ME0_g_hits_000_R  = new TH1F("ME0_g_hits_000_R",  "Simhit Radius :: Photons :: t < 250 ns",   n_R, n1_R, n2_R);
  ME0_N_hits_000_R  = new TH1F("ME0_N_hits_000_R",  "Simhit Radius :: Nuclei :: t < 250 ns",    n_R, n1_R, n2_R);
  ME0_OH_hits_000_R = new TH1F("ME0_OH_hits_000_R", "Simhit Radius :: Other Hadrons :: t < 250 ns", n_R, n1_R, n2_R);
  ME0_All_hits_000_R= new TH1F("ME0_All_hits_000_R","Simhit Radius :: All Particles :: t < 250 ns", n_R, n1_R, n2_R);

  ME0_el_hits_250_R = new TH1F("ME0_el_hits_250_R", "Simhit Radius :: Electrons :: t > 250 ns", n_R, n1_R, n2_R);
  ME0_mu_hits_250_R = new TH1F("ME0_mu_hits_250_R", "Simhit Radius :: Muons :: t > 250 ns",     n_R, n1_R, n2_R);
  ME0_pi_hits_250_R = new TH1F("ME0_pi_hits_250_R", "Simhit Radius :: Pions :: t > 250 ns",     n_R, n1_R, n2_R);
  ME0_ka_hits_250_R = new TH1F("ME0_ka_hits_250_R", "Simhit Radius :: Kaons :: t > 250 ns",     n_R, n1_R, n2_R);
  ME0_p_hits_250_R  = new TH1F("ME0_p_hits_250_R",  "Simhit Radius :: Protons :: t > 250 ns",   n_R, n1_R, n2_R);
  ME0_n_hits_250_R  = new TH1F("ME0_n_hits_250_R",  "Simhit Radius :: Neutrons :: t > 250 ns",  n_R, n1_R, n2_R);
  ME0_g_hits_250_R  = new TH1F("ME0_g_hits_250_R",  "Simhit Radius :: Photons :: t > 250 ns",   n_R, n1_R, n2_R);
  ME0_N_hits_250_R  = new TH1F("ME0_N_hits_250_R",  "Simhit Radius :: Nuclei :: t > 250 ns",    n_R, n1_R, n2_R);
  ME0_OH_hits_250_R = new TH1F("ME0_OH_hits_250_R", "Simhit Radius :: Other Hadrons :: t > 250 ns", n_R, n1_R, n2_R);
  ME0_All_hits_250_R= new TH1F("ME0_All_hits_250_R","Simhit Radius :: All Particles :: t > 250 ns", n_R, n1_R, n2_R);

  G11_el_hits_000_R = new TH1F("G11_el_hits_000_R", "Simhit Radius :: Electrons :: t < 250 ns", n_P, n1_P, n2_P);
  G11_mu_hits_000_R = new TH1F("G11_mu_hits_000_R", "Simhit Radius :: Muons :: t < 250 ns",     n_P, n1_P, n2_P);
  G11_pi_hits_000_R = new TH1F("G11_pi_hits_000_R", "Simhit Radius :: Pions :: t < 250 ns",     n_P, n1_P, n2_P);
  G11_ka_hits_000_R = new TH1F("G11_ka_hits_000_R", "Simhit Radius :: Kaons :: t < 250 ns",     n_P, n1_P, n2_P);
  G11_p_hits_000_R  = new TH1F("G11_p_hits_000_R",  "Simhit Radius :: Protons :: t < 250 ns",   n_P, n1_P, n2_P);
  G11_n_hits_000_R  = new TH1F("G11_n_hits_000_R",  "Simhit Radius :: Neutrons :: t < 250 ns",  n_P, n1_P, n2_P);
  G11_g_hits_000_R  = new TH1F("G11_g_hits_000_R",  "Simhit Radius :: Photons :: t < 250 ns",   n_P, n1_P, n2_P);
  G11_N_hits_000_R  = new TH1F("G11_N_hits_000_R",  "Simhit Radius :: Nuclei :: t < 250 ns",    n_P, n1_P, n2_P);
  G11_OH_hits_000_R = new TH1F("G11_OH_hits_000_R", "Simhit Radius :: Other Hadrons :: t < 250 ns", n_P, n1_P, n2_P);
  G11_All_hits_000_R= new TH1F("G11_All_hits_000_R","Simhit Radius :: All Particles :: t < 250 ns", n_P, n1_P, n2_P);

  G11_el_hits_250_R = new TH1F("G11_el_hits_250_R", "Simhit Radius :: Electrons :: t > 250 ns", n_P, n1_P, n2_P);
  G11_mu_hits_250_R = new TH1F("G11_mu_hits_250_R", "Simhit Radius :: Muons :: t > 250 ns",     n_P, n1_P, n2_P);
  G11_pi_hits_250_R = new TH1F("G11_pi_hits_250_R", "Simhit Radius :: Pions :: t > 250 ns",     n_P, n1_P, n2_P);
  G11_ka_hits_250_R = new TH1F("G11_ka_hits_250_R", "Simhit Radius :: Kaons :: t > 250 ns",     n_P, n1_P, n2_P);
  G11_p_hits_250_R  = new TH1F("G11_p_hits_250_R",  "Simhit Radius :: Protons :: t > 250 ns",   n_P, n1_P, n2_P);
  G11_n_hits_250_R  = new TH1F("G11_n_hits_250_R",  "Simhit Radius :: Neutrons :: t > 250 ns",  n_P, n1_P, n2_P);
  G11_g_hits_250_R  = new TH1F("G11_g_hits_250_R",  "Simhit Radius :: Photons :: t > 250 ns",   n_P, n1_P, n2_P);
  G11_N_hits_250_R  = new TH1F("G11_N_hits_250_R",  "Simhit Radius :: Nuclei :: t > 250 ns",    n_P, n1_P, n2_P);
  G11_OH_hits_250_R = new TH1F("G11_OH_hits_250_R", "Simhit Radius :: Other Hadrons :: t > 250 ns", n_P, n1_P, n2_P);
  G11_All_hits_250_R= new TH1F("G11_All_hits_250_R","Simhit Radius :: All Particles :: t > 250 ns", n_P, n1_P, n2_P);

  G21_el_hits_000_R = new TH1F("G21_el_hits_000_R", "Simhit Radius :: Electrons :: t < 250 ns", n_Q, n1_Q, n2_Q);
  G21_mu_hits_000_R = new TH1F("G21_mu_hits_000_R", "Simhit Radius :: Muons :: t < 250 ns",     n_Q, n1_Q, n2_Q);
  G21_pi_hits_000_R = new TH1F("G21_pi_hits_000_R", "Simhit Radius :: Pions :: t < 250 ns",     n_Q, n1_Q, n2_Q);
  G21_ka_hits_000_R = new TH1F("G21_ka_hits_000_R", "Simhit Radius :: Kaons :: t < 250 ns",     n_Q, n1_Q, n2_Q);
  G21_p_hits_000_R  = new TH1F("G21_p_hits_000_R",  "Simhit Radius :: Protons :: t < 250 ns",   n_Q, n1_Q, n2_Q);
  G21_n_hits_000_R  = new TH1F("G21_n_hits_000_R",  "Simhit Radius :: Neutrons :: t < 250 ns",  n_Q, n1_Q, n2_Q);
  G21_g_hits_000_R  = new TH1F("G21_g_hits_000_R",  "Simhit Radius :: Photons :: t < 250 ns",   n_Q, n1_Q, n2_Q);
  G21_N_hits_000_R  = new TH1F("G21_N_hits_000_R",  "Simhit Radius :: Nuclei :: t < 250 ns",    n_Q, n1_Q, n2_Q);
  G21_OH_hits_000_R = new TH1F("G21_OH_hits_000_R", "Simhit Radius :: Other Hadrons :: t < 250 ns", n_Q, n1_Q, n2_Q);
  G21_All_hits_000_R= new TH1F("G21_All_hits_000_R","Simhit Radius :: All Particles :: t < 250 ns", n_Q, n1_Q, n2_Q);

  G21_el_hits_250_R = new TH1F("G21_el_hits_250_R", "Simhit Radius :: Electrons :: t > 250 ns", n_Q, n1_Q, n2_Q);
  G21_mu_hits_250_R = new TH1F("G21_mu_hits_250_R", "Simhit Radius :: Muons :: t > 250 ns",     n_Q, n1_Q, n2_Q);
  G21_pi_hits_250_R = new TH1F("G21_pi_hits_250_R", "Simhit Radius :: Pions :: t > 250 ns",     n_Q, n1_Q, n2_Q);
  G21_ka_hits_250_R = new TH1F("G21_ka_hits_250_R", "Simhit Radius :: Kaons :: t > 250 ns",     n_Q, n1_Q, n2_Q);
  G21_p_hits_250_R  = new TH1F("G21_p_hits_250_R",  "Simhit Radius :: Protons :: t > 250 ns",   n_Q, n1_Q, n2_Q);
  G21_n_hits_250_R  = new TH1F("G21_n_hits_250_R",  "Simhit Radius :: Neutrons :: t > 250 ns",  n_Q, n1_Q, n2_Q);
  G21_g_hits_250_R  = new TH1F("G21_g_hits_250_R",  "Simhit Radius :: Photons :: t > 250 ns",   n_Q, n1_Q, n2_Q);
  G21_N_hits_250_R  = new TH1F("G21_N_hits_250_R",  "Simhit Radius :: Nuclei :: t > 250 ns",    n_Q, n1_Q, n2_Q);
  G21_OH_hits_250_R = new TH1F("G21_OH_hits_250_R", "Simhit Radius :: Other Hadrons :: t > 250 ns", n_Q, n1_Q, n2_Q);
  G21_All_hits_250_R= new TH1F("G21_All_hits_250_R","Simhit Radius :: All Particles :: t > 250 ns", n_Q, n1_Q, n2_Q);

  // ME0_All_hits_00_R= new TH1F("ME0_All_hits_00_R","Simhit Radius :: All Particles :: t < 25 ns", n_R, n1_R, n2_R);
  // ME0_All_hits_25_R= new TH1F("ME0_All_hits_25_R","Simhit Radius :: All Particles :: t > 25 ns", n_R, n1_R, n2_R);
  // G11_All_hits_00_R= new TH1F("G11_All_hits_00_R","Simhit Radius :: All Particles :: t < 25 ns", n_P, n1_P, n2_P);
  // G11_All_hits_25_R= new TH1F("G11_All_hits_25_R","Simhit Radius :: All Particles :: t > 25 ns", n_P, n1_P, n2_P);
  // G21_All_hits_00_R= new TH1F("G21_All_hits_00_R","Simhit Radius :: All Particles :: t < 25 ns", n_Q, n1_Q, n2_Q);
  // G21_All_hits_25_R= new TH1F("G21_All_hits_25_R","Simhit Radius :: All Particles :: t > 25 ns", n_Q, n1_Q, n2_Q);

  ME0_el_hits_00_R = new TH1F("ME0_el_hits_00_R", "Simhit Radius :: Electrons :: t < 25 ns", n_R, n1_R, n2_R);
  ME0_mu_hits_00_R = new TH1F("ME0_mu_hits_00_R", "Simhit Radius :: Muons :: t < 25 ns",     n_R, n1_R, n2_R);
  ME0_pi_hits_00_R = new TH1F("ME0_pi_hits_00_R", "Simhit Radius :: Pions :: t < 25 ns",     n_R, n1_R, n2_R);
  ME0_ka_hits_00_R = new TH1F("ME0_ka_hits_00_R", "Simhit Radius :: Kaons :: t < 25 ns",     n_R, n1_R, n2_R);
  ME0_p_hits_00_R  = new TH1F("ME0_p_hits_00_R",  "Simhit Radius :: Protons :: t < 25 ns",   n_R, n1_R, n2_R);
  ME0_n_hits_00_R  = new TH1F("ME0_n_hits_00_R",  "Simhit Radius :: Neutrons :: t < 25 ns",  n_R, n1_R, n2_R);
  ME0_g_hits_00_R  = new TH1F("ME0_g_hits_00_R",  "Simhit Radius :: Photons :: t < 25 ns",   n_R, n1_R, n2_R);
  ME0_N_hits_00_R  = new TH1F("ME0_N_hits_00_R",  "Simhit Radius :: Nuclei :: t < 25 ns",    n_R, n1_R, n2_R);
  ME0_OH_hits_00_R = new TH1F("ME0_OH_hits_00_R", "Simhit Radius :: Other Hadrons :: t < 25 ns", n_R, n1_R, n2_R);
  ME0_All_hits_00_R= new TH1F("ME0_All_hits_00_R","Simhit Radius :: All Particles :: t < 25 ns", n_R, n1_R, n2_R);

  ME0_el_hits_25_R = new TH1F("ME0_el_hits_25_R", "Simhit Radius :: Electrons :: t > 25 ns", n_R, n1_R, n2_R);
  ME0_mu_hits_25_R = new TH1F("ME0_mu_hits_25_R", "Simhit Radius :: Muons :: t > 25 ns",     n_R, n1_R, n2_R);
  ME0_pi_hits_25_R = new TH1F("ME0_pi_hits_25_R", "Simhit Radius :: Pions :: t > 25 ns",     n_R, n1_R, n2_R);
  ME0_ka_hits_25_R = new TH1F("ME0_ka_hits_25_R", "Simhit Radius :: Kaons :: t > 25 ns",     n_R, n1_R, n2_R);
  ME0_p_hits_25_R  = new TH1F("ME0_p_hits_25_R",  "Simhit Radius :: Protons :: t > 25 ns",   n_R, n1_R, n2_R);
  ME0_n_hits_25_R  = new TH1F("ME0_n_hits_25_R",  "Simhit Radius :: Neutrons :: t > 25 ns",  n_R, n1_R, n2_R);
  ME0_g_hits_25_R  = new TH1F("ME0_g_hits_25_R",  "Simhit Radius :: Photons :: t > 25 ns",   n_R, n1_R, n2_R);
  ME0_N_hits_25_R  = new TH1F("ME0_N_hits_25_R",  "Simhit Radius :: Nuclei :: t > 25 ns",    n_R, n1_R, n2_R);
  ME0_OH_hits_25_R = new TH1F("ME0_OH_hits_25_R", "Simhit Radius :: Other Hadrons :: t > 25 ns", n_R, n1_R, n2_R);
  ME0_All_hits_25_R= new TH1F("ME0_All_hits_25_R","Simhit Radius :: All Particles :: t > 25 ns", n_R, n1_R, n2_R);

  G11_el_hits_00_R = new TH1F("G11_el_hits_00_R", "Simhit Radius :: Electrons :: t < 25 ns", n_P, n1_P, n2_P);
  G11_mu_hits_00_R = new TH1F("G11_mu_hits_00_R", "Simhit Radius :: Muons :: t < 25 ns",     n_P, n1_P, n2_P);
  G11_pi_hits_00_R = new TH1F("G11_pi_hits_00_R", "Simhit Radius :: Pions :: t < 25 ns",     n_P, n1_P, n2_P);
  G11_ka_hits_00_R = new TH1F("G11_ka_hits_00_R", "Simhit Radius :: Kaons :: t < 25 ns",     n_P, n1_P, n2_P);
  G11_p_hits_00_R  = new TH1F("G11_p_hits_00_R",  "Simhit Radius :: Protons :: t < 25 ns",   n_P, n1_P, n2_P);
  G11_n_hits_00_R  = new TH1F("G11_n_hits_00_R",  "Simhit Radius :: Neutrons :: t < 25 ns",  n_P, n1_P, n2_P);
  G11_g_hits_00_R  = new TH1F("G11_g_hits_00_R",  "Simhit Radius :: Photons :: t < 25 ns",   n_P, n1_P, n2_P);
  G11_N_hits_00_R  = new TH1F("G11_N_hits_00_R",  "Simhit Radius :: Nuclei :: t < 25 ns",    n_P, n1_P, n2_P);
  G11_OH_hits_00_R = new TH1F("G11_OH_hits_00_R", "Simhit Radius :: Other Hadrons :: t < 25 ns", n_P, n1_P, n2_P);
  G11_All_hits_00_R= new TH1F("G11_All_hits_00_R","Simhit Radius :: All Particles :: t < 25 ns", n_P, n1_P, n2_P);

  G11_el_hits_25_R = new TH1F("G11_el_hits_25_R", "Simhit Radius :: Electrons :: t > 25 ns", n_P, n1_P, n2_P);
  G11_mu_hits_25_R = new TH1F("G11_mu_hits_25_R", "Simhit Radius :: Muons :: t > 25 ns",     n_P, n1_P, n2_P);
  G11_pi_hits_25_R = new TH1F("G11_pi_hits_25_R", "Simhit Radius :: Pions :: t > 25 ns",     n_P, n1_P, n2_P);
  G11_ka_hits_25_R = new TH1F("G11_ka_hits_25_R", "Simhit Radius :: Kaons :: t > 25 ns",     n_P, n1_P, n2_P);
  G11_p_hits_25_R  = new TH1F("G11_p_hits_25_R",  "Simhit Radius :: Protons :: t > 25 ns",   n_P, n1_P, n2_P);
  G11_n_hits_25_R  = new TH1F("G11_n_hits_25_R",  "Simhit Radius :: Neutrons :: t > 25 ns",  n_P, n1_P, n2_P);
  G11_g_hits_25_R  = new TH1F("G11_g_hits_25_R",  "Simhit Radius :: Photons :: t > 25 ns",   n_P, n1_P, n2_P);
  G11_N_hits_25_R  = new TH1F("G11_N_hits_25_R",  "Simhit Radius :: Nuclei :: t > 25 ns",    n_P, n1_P, n2_P);
  G11_OH_hits_25_R = new TH1F("G11_OH_hits_25_R", "Simhit Radius :: Other Hadrons :: t > 25 ns", n_P, n1_P, n2_P);
  G11_All_hits_25_R= new TH1F("G11_All_hits_25_R","Simhit Radius :: All Particles :: t > 25 ns", n_P, n1_P, n2_P);

  G21_el_hits_00_R = new TH1F("G21_el_hits_00_R", "Simhit Radius :: Electrons :: t < 25 ns", n_Q, n1_Q, n2_Q);
  G21_mu_hits_00_R = new TH1F("G21_mu_hits_00_R", "Simhit Radius :: Muons :: t < 25 ns",     n_Q, n1_Q, n2_Q);
  G21_pi_hits_00_R = new TH1F("G21_pi_hits_00_R", "Simhit Radius :: Pions :: t < 25 ns",     n_Q, n1_Q, n2_Q);
  G21_ka_hits_00_R = new TH1F("G21_ka_hits_00_R", "Simhit Radius :: Kaons :: t < 25 ns",     n_Q, n1_Q, n2_Q);
  G21_p_hits_00_R  = new TH1F("G21_p_hits_00_R",  "Simhit Radius :: Protons :: t < 25 ns",   n_Q, n1_Q, n2_Q);
  G21_n_hits_00_R  = new TH1F("G21_n_hits_00_R",  "Simhit Radius :: Neutrons :: t < 25 ns",  n_Q, n1_Q, n2_Q);
  G21_g_hits_00_R  = new TH1F("G21_g_hits_00_R",  "Simhit Radius :: Photons :: t < 25 ns",   n_Q, n1_Q, n2_Q);
  G21_N_hits_00_R  = new TH1F("G21_N_hits_00_R",  "Simhit Radius :: Nuclei :: t < 25 ns",    n_Q, n1_Q, n2_Q);
  G21_OH_hits_00_R = new TH1F("G21_OH_hits_00_R", "Simhit Radius :: Other Hadrons :: t < 25 ns", n_Q, n1_Q, n2_Q);
  G21_All_hits_00_R= new TH1F("G21_All_hits_00_R","Simhit Radius :: All Particles :: t < 25 ns", n_Q, n1_Q, n2_Q);

  G21_el_hits_25_R = new TH1F("G21_el_hits_25_R", "Simhit Radius :: Electrons :: t > 25 ns", n_Q, n1_Q, n2_Q);
  G21_mu_hits_25_R = new TH1F("G21_mu_hits_25_R", "Simhit Radius :: Muons :: t > 25 ns",     n_Q, n1_Q, n2_Q);
  G21_pi_hits_25_R = new TH1F("G21_pi_hits_25_R", "Simhit Radius :: Pions :: t > 25 ns",     n_Q, n1_Q, n2_Q);
  G21_ka_hits_25_R = new TH1F("G21_ka_hits_25_R", "Simhit Radius :: Kaons :: t > 25 ns",     n_Q, n1_Q, n2_Q);
  G21_p_hits_25_R  = new TH1F("G21_p_hits_25_R",  "Simhit Radius :: Protons :: t > 25 ns",   n_Q, n1_Q, n2_Q);
  G21_n_hits_25_R  = new TH1F("G21_n_hits_25_R",  "Simhit Radius :: Neutrons :: t > 25 ns",  n_Q, n1_Q, n2_Q);
  G21_g_hits_25_R  = new TH1F("G21_g_hits_25_R",  "Simhit Radius :: Photons :: t > 25 ns",   n_Q, n1_Q, n2_Q);
  G21_N_hits_25_R  = new TH1F("G21_N_hits_25_R",  "Simhit Radius :: Nuclei :: t > 25 ns",    n_Q, n1_Q, n2_Q);
  G21_OH_hits_25_R = new TH1F("G21_OH_hits_25_R", "Simhit Radius :: Other Hadrons :: t > 25 ns", n_Q, n1_Q, n2_Q);
  G21_All_hits_25_R= new TH1F("G21_All_hits_25_R","Simhit Radius :: All Particles :: t > 25 ns", n_Q, n1_Q, n2_Q);

  // 2B updated :: change binnings to n_R (ME0) n_P (GE11) and n_Q (GE21)
  G11_L1_el_hits_R  = new TH1F("G11_L1_el_hits_R",  "Simhit Radius :: G1E1 L1 :: Electrons",       n_P, n1_P, n2_P);
  G11_L1_mu_hits_R  = new TH1F("G11_L1_mu_hits_R",  "Simhit Radius :: GE11 L1 :: Muons",           n_P, n1_P, n2_P);
  G11_L1_pi_hits_R  = new TH1F("G11_L1_pi_hits_R",  "Simhit Radius :: GE11 L1 :: Pions",           n_P, n1_P, n2_P);
  G11_L1_ka_hits_R  = new TH1F("G11_L1_ka_hits_R",  "Simhit Radius :: GE11 L1 :: Kaons",           n_P, n1_P, n2_P);
  G11_L1_p_hits_R   = new TH1F("G11_L1_p_hits_R",   "Simhit Radius :: GE11 L1 :: Protons",         n_P, n1_P, n2_P);
  G11_L1_n_hits_R   = new TH1F("G11_L1_n_hits_R",   "Simhit Radius :: GE11 L1 :: Neutrons",        n_P, n1_P, n2_P);
  G11_L1_g_hits_R   = new TH1F("G11_L1_g_hits_R",   "Simhit Radius :: GE11 L1 :: Photons",         n_P, n1_P, n2_P);
  G11_L1_N_hits_R   = new TH1F("G11_L1_N_hits_R",   "Simhit Radius :: GE11 L1 :: Nuclei",          n_P, n1_P, n2_P);
  G11_L1_OH_hits_R  = new TH1F("G11_L1_OH_hits_R",  "Simhit Radius :: GE11 L1 :: Other Hadrons",   n_P, n1_P, n2_P);
  G11_L1_All_hits_R = new TH1F("G11_L1_All_hits_R", "Simhit Radius :: GE11 L1 :: All Particles",   n_P, n1_P, n2_P);
  G11_L1_HIP_hits_R = new TH1F("G11_L1_HIP_hits_R", "Simhit Radius :: GE11 L1 :: Highly Ionizing", n_P, n1_P, n2_P);

  G11_L2_el_hits_R  = new TH1F("G11_L2_el_hits_R",  "Simhit Radius :: GE11 L2 :: Electrons",       n_P, n1_P, n2_P);
  G11_L2_mu_hits_R  = new TH1F("G11_L2_mu_hits_R",  "Simhit Radius :: GE11 L2 :: Muons",           n_P, n1_P, n2_P);
  G11_L2_pi_hits_R  = new TH1F("G11_L2_pi_hits_R",  "Simhit Radius :: GE11 L2 :: Pions",           n_P, n1_P, n2_P);
  G11_L2_ka_hits_R  = new TH1F("G11_L2_ka_hits_R",  "Simhit Radius :: GE11 L2 :: Kaons",           n_P, n1_P, n2_P);
  G11_L2_p_hits_R   = new TH1F("G11_L2_p_hits_R",   "Simhit Radius :: GE11 L2 :: Protons",         n_P, n1_P, n2_P);
  G11_L2_n_hits_R   = new TH1F("G11_L2_n_hits_R",   "Simhit Radius :: GE11 L2 :: Neutrons",        n_P, n1_P, n2_P);
  G11_L2_g_hits_R   = new TH1F("G11_L2_g_hits_R",   "Simhit Radius :: GE11 L2 :: Photons",         n_P, n1_P, n2_P);
  G11_L2_N_hits_R   = new TH1F("G11_L2_N_hits_R",   "Simhit Radius :: GE11 L2 :: Nuclei",          n_P, n1_P, n2_P);
  G11_L2_OH_hits_R  = new TH1F("G11_L2_OH_hits_R",  "Simhit Radius :: GE11 L2 :: Other Hadrons",   n_P, n1_P, n2_P);
  G11_L2_All_hits_R = new TH1F("G11_L2_All_hits_R", "Simhit Radius :: GE11 L2 :: All Particles",   n_P, n1_P, n2_P);
  G11_L2_HIP_hits_R = new TH1F("G11_L2_HIP_hits_R", "Simhit Radius :: GE11 L2 :: Highly Ionizing", n_P, n1_P, n2_P);

  G11_Od_el_hits_R  = new TH1F("G11_Od_el_hits_R",  "Simhit Radius :: GE11 Odd :: Electrons",       n_P, n1_P, n2_P);
  G11_Od_mu_hits_R  = new TH1F("G11_Od_mu_hits_R",  "Simhit Radius :: GE11 Odd :: Muons",           n_P, n1_P, n2_P);
  G11_Od_pi_hits_R  = new TH1F("G11_Od_pi_hits_R",  "Simhit Radius :: GE11 Odd :: Pions",           n_P, n1_P, n2_P);
  G11_Od_ka_hits_R  = new TH1F("G11_Od_ka_hits_R",  "Simhit Radius :: GE11 Odd :: Kaons",           n_P, n1_P, n2_P);
  G11_Od_p_hits_R   = new TH1F("G11_Od_p_hits_R",   "Simhit Radius :: GE11 Odd :: Protons",         n_P, n1_P, n2_P);
  G11_Od_n_hits_R   = new TH1F("G11_Od_n_hits_R",   "Simhit Radius :: GE11 Odd :: Neutrons",        n_P, n1_P, n2_P);
  G11_Od_g_hits_R   = new TH1F("G11_Od_g_hits_R",   "Simhit Radius :: GE11 Odd :: Photons",         n_P, n1_P, n2_P);
  G11_Od_N_hits_R   = new TH1F("G11_Od_N_hits_R",   "Simhit Radius :: GE11 Odd :: Nuclei",          n_P, n1_P, n2_P);
  G11_Od_OH_hits_R  = new TH1F("G11_Od_OH_hits_R",  "Simhit Radius :: GE11 Odd :: Other Hadrons",   n_P, n1_P, n2_P);
  G11_Od_All_hits_R = new TH1F("G11_Od_All_hits_R", "Simhit Radius :: GE11 Odd :: All Particles",   n_P, n1_P, n2_P);
  G11_Od_HIP_hits_R = new TH1F("G11_Od_HIP_hits_R", "Simhit Radius :: GE11 Odd :: Highly Ionizing", n_P, n1_P, n2_P);

  G11_Ev_el_hits_R  = new TH1F("G11_Ev_el_hits_R",  "Simhit Radius :: GE11 Even :: Electrons",       n_P, n1_P, n2_P);
  G11_Ev_mu_hits_R  = new TH1F("G11_Ev_mu_hits_R",  "Simhit Radius :: GE11 Even :: Muons",           n_P, n1_P, n2_P);
  G11_Ev_pi_hits_R  = new TH1F("G11_Ev_pi_hits_R",  "Simhit Radius :: GE11 Even :: Pions",           n_P, n1_P, n2_P);
  G11_Ev_ka_hits_R  = new TH1F("G11_Ev_ka_hits_R",  "Simhit Radius :: GE11 Even :: Kaons",           n_P, n1_P, n2_P);
  G11_Ev_p_hits_R   = new TH1F("G11_Ev_p_hits_R",   "Simhit Radius :: GE11 Even :: Protons",         n_P, n1_P, n2_P);
  G11_Ev_n_hits_R   = new TH1F("G11_Ev_n_hits_R",   "Simhit Radius :: GE11 Even :: Neutrons",        n_P, n1_P, n2_P);
  G11_Ev_g_hits_R   = new TH1F("G11_Ev_g_hits_R",   "Simhit Radius :: GE11 Even :: Photons",         n_P, n1_P, n2_P);
  G11_Ev_N_hits_R   = new TH1F("G11_Ev_N_hits_R",   "Simhit Radius :: GE11 Even :: Nuclei",          n_P, n1_P, n2_P);
  G11_Ev_OH_hits_R  = new TH1F("G11_Ev_OH_hits_R",  "Simhit Radius :: GE11 Even :: Other Hadrons",   n_P, n1_P, n2_P);
  G11_Ev_All_hits_R = new TH1F("G11_Ev_All_hits_R", "Simhit Radius :: GE11 Even :: All Particles",   n_P, n1_P, n2_P);
  G11_Ev_HIP_hits_R = new TH1F("G11_Ev_HIP_hits_R", "Simhit Radius :: GE11 Even :: Highly Ionizing", n_P, n1_P, n2_P);

  /*
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
  */

  // Simhit vs EtaPartition
  G11_el_hits_E = new TH1F("G11_el_hits_E", "Simhit Radius :: Electrons", n_N, n1_N, n2_N);
  G11_mu_hits_E = new TH1F("G11_mu_hits_E", "Simhit Radius :: Muons",     n_N, n1_N, n2_N);
  G11_pi_hits_E = new TH1F("G11_pi_hits_E", "Simhit Radius :: Pions",     n_N, n1_N, n2_N);
  G11_ka_hits_E = new TH1F("G11_ka_hits_E", "Simhit Radius :: Kaons",     n_N, n1_N, n2_N);
  G11_p_hits_E  = new TH1F("G11_p_hits_E",  "Simhit Radius :: Protons",   n_N, n1_N, n2_N);
  G11_n_hits_E  = new TH1F("G11_n_hits_E",  "Simhit Radius :: Neutrons",  n_N, n1_N, n2_N);
  G11_g_hits_E  = new TH1F("G11_g_hits_E",  "Simhit Radius :: Photons",   n_N, n1_N, n2_N);
  G11_N_hits_E  = new TH1F("G11_N_hits_E",  "Simhit Radius :: Nuclei",    n_N, n1_N, n2_N);
  G11_OH_hits_E = new TH1F("G11_OH_hits_E", "Simhit Radius :: Other Hadrons", n_N, n1_N, n2_N);
  G11_All_hits_E= new TH1F("G11_All_hits_E","Simhit Radius :: All Particles", n_N, n1_N, n2_N);
  G11_HIP_hits_E= new TH1F("G11_HIP_hits_E","Simhit Radius :: HIP Particles", n_N, n1_N, n2_N);

  G21_el_hits_E = new TH1F("G21_el_hits_E", "Simhit Radius :: Electrons", n_N2, n1_N2, n2_N2);
  G21_mu_hits_E = new TH1F("G21_mu_hits_E", "Simhit Radius :: Muons",     n_N2, n1_N2, n2_N2);
  G21_pi_hits_E = new TH1F("G21_pi_hits_E", "Simhit Radius :: Pions",     n_N2, n1_N2, n2_N2);
  G21_ka_hits_E = new TH1F("G21_ka_hits_E", "Simhit Radius :: Kaons",     n_N2, n1_N2, n2_N2);
  G21_p_hits_E  = new TH1F("G21_p_hits_E",  "Simhit Radius :: Protons",   n_N2, n1_N2, n2_N2);
  G21_n_hits_E  = new TH1F("G21_n_hits_E",  "Simhit Radius :: Neutrons",  n_N2, n1_N2, n2_N2);
  G21_g_hits_E  = new TH1F("G21_g_hits_E",  "Simhit Radius :: Photons",   n_N2, n1_N2, n2_N2);
  G21_N_hits_E  = new TH1F("G21_N_hits_E",  "Simhit Radius :: Nuclei",    n_N2, n1_N2, n2_N2);
  G21_OH_hits_E = new TH1F("G21_OH_hits_E", "Simhit Radius :: Other Hadrons", n_N2, n1_N2, n2_N2);
  G21_All_hits_E= new TH1F("G21_All_hits_E","Simhit Radius :: All Particles", n_N2, n1_N2, n2_N2);
  G21_HIP_hits_E= new TH1F("G21_HIP_hits_E","Simhit Radius :: Highly Ionising", n_N2, n1_N2, n2_N2);

  ME0_el_hits_E = new TH1F("ME0_el_hits_E", "Simhit Radius :: Electrons", n_N, n1_N, n2_N);
  ME0_mu_hits_E = new TH1F("ME0_mu_hits_E", "Simhit Radius :: Muons",     n_N, n1_N, n2_N);
  ME0_pi_hits_E = new TH1F("ME0_pi_hits_E", "Simhit Radius :: Pions",     n_N, n1_N, n2_N);
  ME0_ka_hits_E = new TH1F("ME0_ka_hits_E", "Simhit Radius :: Kaons",     n_N, n1_N, n2_N);
  ME0_p_hits_E  = new TH1F("ME0_p_hits_E",  "Simhit Radius :: Protons",   n_N, n1_N, n2_N);
  ME0_n_hits_E  = new TH1F("ME0_n_hits_E",  "Simhit Radius :: Neutrons",  n_N, n1_N, n2_N);
  ME0_g_hits_E  = new TH1F("ME0_g_hits_E",  "Simhit Radius :: Photons",   n_N, n1_N, n2_N);
  ME0_N_hits_E  = new TH1F("ME0_N_hits_E",  "Simhit Radius :: Nuclei",    n_N, n1_N, n2_N);
  ME0_OH_hits_E = new TH1F("ME0_OH_hits_E", "Simhit Radius :: Other Hadrons", n_N, n1_N, n2_N);
  ME0_All_hits_E= new TH1F("ME0_All_hits_E","Simhit Radius :: All Particles", n_N, n1_N, n2_N);
  ME0_HIP_hits_E= new TH1F("ME0_HIP_hits_E","Simhit Radius :: Highly Ionising", n_N, n1_N, n2_N);

  ME0_All_hits_00_E= new TH1F("ME0_All_hits_00_E","Simhit Radius :: All Particles :: t < 25 ns", n_N, n1_N, n2_N);
  ME0_All_hits_25_E= new TH1F("ME0_All_hits_25_E","Simhit Radius :: All Particles :: t > 25 ns", n_N, n1_N, n2_N);
  G11_All_hits_00_E= new TH1F("G11_All_hits_00_E","Simhit Radius :: All Particles :: t < 25 ns", n_N, n1_N, n2_N);
  G11_All_hits_25_E= new TH1F("G11_All_hits_25_E","Simhit Radius :: All Particles :: t > 25 ns", n_N, n1_N, n2_N);
  G21_All_hits_00_E= new TH1F("G21_All_hits_00_E","Simhit Radius :: All Particles :: t < 25 ns", n_N, n1_N, n2_N);
  G21_All_hits_25_E= new TH1F("G21_All_hits_25_E","Simhit Radius :: All Particles :: t > 25 ns", n_N, n1_N, n2_N);

  ME0_All_hits_000_E= new TH1F("ME0_All_hits_000_E","Simhit Radius :: All Particles :: t < 250 ns", n_N, n1_N, n2_N);
  ME0_All_hits_250_E= new TH1F("ME0_All_hits_250_E","Simhit Radius :: All Particles :: t > 250 ns", n_N, n1_N, n2_N);
  G11_All_hits_000_E= new TH1F("G11_All_hits_000_E","Simhit Radius :: All Particles :: t < 250 ns", n_N, n1_N, n2_N);
  G11_All_hits_250_E= new TH1F("G11_All_hits_250_E","Simhit Radius :: All Particles :: t > 250 ns", n_N, n1_N, n2_N);
  G21_All_hits_000_E= new TH1F("G21_All_hits_000_E","Simhit Radius :: All Particles :: t < 250 ns", n_N, n1_N, n2_N);
  G21_All_hits_250_E= new TH1F("G21_All_hits_250_E","Simhit Radius :: All Particles :: t > 250 ns", n_N, n1_N, n2_N);

  G11_L1Odd_el_hits_E = new TH1F("G11_el_hits_E", "Simhit Radius :: GE11 L1 Odd :: Electrons", n_N, n1_N, n2_N);
  G11_L1Odd_mu_hits_E = new TH1F("G11_mu_hits_E", "Simhit Radius :: GE11 L1 Odd :: Muons",     n_N, n1_N, n2_N);
  G11_L1Odd_pi_hits_E = new TH1F("G11_pi_hits_E", "Simhit Radius :: GE11 L1 Odd :: Pions",     n_N, n1_N, n2_N);
  G11_L1Odd_ka_hits_E = new TH1F("G11_ka_hits_E", "Simhit Radius :: GE11 L1 Odd :: Kaons",     n_N, n1_N, n2_N);
  G11_L1Odd_p_hits_E  = new TH1F("G11_p_hits_E",  "Simhit Radius :: GE11 L1 Odd :: Protons",   n_N, n1_N, n2_N);
  G11_L1Odd_n_hits_E  = new TH1F("G11_n_hits_E",  "Simhit Radius :: GE11 L1 Odd :: Neutrons",  n_N, n1_N, n2_N);
  G11_L1Odd_g_hits_E  = new TH1F("G11_g_hits_E",  "Simhit Radius :: GE11 L1 Odd :: Photons",   n_N, n1_N, n2_N);
  G11_L1Odd_N_hits_E  = new TH1F("G11_N_hits_E",  "Simhit Radius :: GE11 L1 Odd :: Nuclei",    n_N, n1_N, n2_N);
  G11_L1Odd_OH_hits_E = new TH1F("G11_OH_hits_E", "Simhit Radius :: GE11 L1 Odd :: Other Hadrons", n_N, n1_N, n2_N);
  G11_L1Odd_All_hits_E= new TH1F("G11_All_hits_E","Simhit Radius :: GE11 L1 Odd :: All Particles", n_N, n1_N, n2_N);
  G11_L1Odd_HIP_hits_E= new TH1F("G11_HIP_hits_E","Simhit Radius :: GE11 L1 Odd :: HIP Particles", n_N, n1_N, n2_N);

  G11_L1Even_el_hits_E = new TH1F("G11_el_hits_E", "Simhit Radius :: GE11 L1 Even :: Electrons", n_N, n1_N, n2_N);
  G11_L1Even_mu_hits_E = new TH1F("G11_mu_hits_E", "Simhit Radius :: GE11 L1 Even :: Muons",     n_N, n1_N, n2_N);
  G11_L1Even_pi_hits_E = new TH1F("G11_pi_hits_E", "Simhit Radius :: GE11 L1 Even :: Pions",     n_N, n1_N, n2_N);
  G11_L1Even_ka_hits_E = new TH1F("G11_ka_hits_E", "Simhit Radius :: GE11 L1 Even :: Kaons",     n_N, n1_N, n2_N);
  G11_L1Even_p_hits_E  = new TH1F("G11_p_hits_E",  "Simhit Radius :: GE11 L1 Even :: Protons",   n_N, n1_N, n2_N);
  G11_L1Even_n_hits_E  = new TH1F("G11_n_hits_E",  "Simhit Radius :: GE11 L1 Even :: Neutrons",  n_N, n1_N, n2_N);
  G11_L1Even_g_hits_E  = new TH1F("G11_g_hits_E",  "Simhit Radius :: GE11 L1 Even :: Photons",   n_N, n1_N, n2_N);
  G11_L1Even_N_hits_E  = new TH1F("G11_N_hits_E",  "Simhit Radius :: GE11 L1 Even :: Nuclei",    n_N, n1_N, n2_N);
  G11_L1Even_OH_hits_E = new TH1F("G11_OH_hits_E", "Simhit Radius :: GE11 L1 Even :: Other Hadrons", n_N, n1_N, n2_N);
  G11_L1Even_All_hits_E= new TH1F("G11_All_hits_E","Simhit Radius :: GE11 L1 Even :: All Particles", n_N, n1_N, n2_N);
  G11_L1Even_HIP_hits_E= new TH1F("G11_HIP_hits_E","Simhit Radius :: GE11 L1 Even :: HIP Particles", n_N, n1_N, n2_N);

  G11_L2Odd_el_hits_E = new TH1F("G11_el_hits_E", "Simhit Radius :: GE11 L2 Odd :: Electrons", n_N, n1_N, n2_N);
  G11_L2Odd_mu_hits_E = new TH1F("G11_mu_hits_E", "Simhit Radius :: GE11 L2 Odd :: Muons",     n_N, n1_N, n2_N);
  G11_L2Odd_pi_hits_E = new TH1F("G11_pi_hits_E", "Simhit Radius :: GE11 L2 Odd :: Pions",     n_N, n1_N, n2_N);
  G11_L2Odd_ka_hits_E = new TH1F("G11_ka_hits_E", "Simhit Radius :: GE11 L2 Odd :: Kaons",     n_N, n1_N, n2_N);
  G11_L2Odd_p_hits_E  = new TH1F("G11_p_hits_E",  "Simhit Radius :: GE11 L2 Odd :: Protons",   n_N, n1_N, n2_N);
  G11_L2Odd_n_hits_E  = new TH1F("G11_n_hits_E",  "Simhit Radius :: GE11 L2 Odd :: Neutrons",  n_N, n1_N, n2_N);
  G11_L2Odd_g_hits_E  = new TH1F("G11_g_hits_E",  "Simhit Radius :: GE11 L2 Odd :: Photons",   n_N, n1_N, n2_N);
  G11_L2Odd_N_hits_E  = new TH1F("G11_N_hits_E",  "Simhit Radius :: GE11 L2 Odd :: Nuclei",    n_N, n1_N, n2_N);
  G11_L2Odd_OH_hits_E = new TH1F("G11_OH_hits_E", "Simhit Radius :: GE11 L2 Odd :: Other Hadrons", n_N, n1_N, n2_N);
  G11_L2Odd_All_hits_E= new TH1F("G11_All_hits_E","Simhit Radius :: GE11 L2 Odd :: All Particles", n_N, n1_N, n2_N);
  G11_L2Odd_HIP_hits_E= new TH1F("G11_HIP_hits_E","Simhit Radius :: GE11 L2 Odd :: HIP Particles", n_N, n1_N, n2_N);

  G11_L2Even_el_hits_E = new TH1F("G11_el_hits_E", "Simhit Radius :: GE11 L2 Even :: Electrons", n_N, n1_N, n2_N);
  G11_L2Even_mu_hits_E = new TH1F("G11_mu_hits_E", "Simhit Radius :: GE11 L2 Even :: Muons",     n_N, n1_N, n2_N);
  G11_L2Even_pi_hits_E = new TH1F("G11_pi_hits_E", "Simhit Radius :: GE11 L2 Even :: Pions",     n_N, n1_N, n2_N);
  G11_L2Even_ka_hits_E = new TH1F("G11_ka_hits_E", "Simhit Radius :: GE11 L2 Even :: Kaons",     n_N, n1_N, n2_N);
  G11_L2Even_p_hits_E  = new TH1F("G11_p_hits_E",  "Simhit Radius :: GE11 L2 Even :: Protons",   n_N, n1_N, n2_N);
  G11_L2Even_n_hits_E  = new TH1F("G11_n_hits_E",  "Simhit Radius :: GE11 L2 Even :: Neutrons",  n_N, n1_N, n2_N);
  G11_L2Even_g_hits_E  = new TH1F("G11_g_hits_E",  "Simhit Radius :: GE11 L2 Even :: Photons",   n_N, n1_N, n2_N);
  G11_L2Even_N_hits_E  = new TH1F("G11_N_hits_E",  "Simhit Radius :: GE11 L2 Even :: Nuclei",    n_N, n1_N, n2_N);
  G11_L2Even_OH_hits_E = new TH1F("G11_OH_hits_E", "Simhit Radius :: GE11 L2 Even :: Other Hadrons", n_N, n1_N, n2_N);
  G11_L2Even_All_hits_E= new TH1F("G11_All_hits_E","Simhit Radius :: GE11 L2 Even :: All Particles", n_N, n1_N, n2_N);
  G11_L2Even_HIP_hits_E= new TH1F("G11_HIP_hits_E","Simhit Radius :: GE11 L2 Even :: HIP Particles", n_N, n1_N, n2_N);


  // Simhit Process Type
  GEM_el_process = new TH1F("GEM_el_process", "SimHit Process Type :: Electrons",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_el_process); 
  GEM_mu_process = new TH1F("GEM_mu_process", "SimHit Process Type :: Muons",         n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_mu_process); 
  GEM_pi_process = new TH1F("GEM_pi_process", "Simhit Process Type :: Pions",         n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_pi_process); 
  GEM_ka_process = new TH1F("GEM_ka_process", "Simhit Process Type :: Kaons",         n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_ka_process); 
  GEM_p_process  = new TH1F("GEM_p_process",  "Simhit Process Type :: Protons",       n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_p_process); 
  GEM_n_process  = new TH1F("GEM_n_process",  "Simhit Process Type :: Neutrons",      n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_n_process); 
  GEM_g_process  = new TH1F("GEM_g_process",  "Simhit Process Type :: Photons",       n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_g_process); 
  GEM_N_process  = new TH1F("GEM_N_process",  "Simhit Process Type :: Nuclei",        n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_N_process); 
  GEM_OH_process = new TH1F("GEM_OH_process", "Simhit Process Type :: Other Hadrons", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_OH_process); 

  ME0_el_process = new TH1F("ME0_el_process", "Simhit Process Type :: Electrons", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_el_process); 
  ME0_mu_process = new TH1F("ME0_mu_process", "Simhit Process Type :: Muons",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_mu_process); 
  ME0_pi_process = new TH1F("ME0_pi_process", "Simhit Process Type :: Pions",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_pi_process); 
  ME0_ka_process = new TH1F("ME0_ka_process", "Simhit Process Type :: Kaons",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_ka_process); 
  ME0_p_process  = new TH1F("ME0_p_process",  "Simhit Process Type :: Protons",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_p_process); 
  ME0_n_process  = new TH1F("ME0_n_process",  "Simhit Process Type :: Neutrons",  n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_n_process); 
  ME0_g_process  = new TH1F("ME0_g_process",  "Simhit Process Type :: Photons",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_g_process); 
  ME0_N_process  = new TH1F("ME0_N_process",  "Simhit Process Type :: Nuclei",    n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_N_process); 
  ME0_OH_process = new TH1F("ME0_OH_process", "Simhit Process Type :: Other Hadrons", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_OH_process); 

  ME0_el_proc000 = new TH1F("ME0_el_proc000", "Simhit Process Type :: Electrons :: t < 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_el_proc000); 
  ME0_mu_proc000 = new TH1F("ME0_mu_proc000", "Simhit Process Type :: Muons :: t < 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_mu_proc000); 
  ME0_pi_proc000 = new TH1F("ME0_pi_proc000", "Simhit Process Type :: Pions :: t < 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_pi_proc000); 
  ME0_ka_proc000 = new TH1F("ME0_ka_proc000", "Simhit Process Type :: Kaons :: t < 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_ka_proc000); 
  ME0_p_proc000  = new TH1F("ME0_p_proc000",  "Simhit Process Type :: Protons :: t < 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_p_proc000); 
  ME0_n_proc000  = new TH1F("ME0_n_proc000",  "Simhit Process Type :: Neutrons :: t < 250 ns",  n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_n_proc000); 
  ME0_g_proc000  = new TH1F("ME0_g_proc000",  "Simhit Process Type :: Photons :: t < 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_g_proc000); 
  ME0_N_proc000  = new TH1F("ME0_N_proc000",  "Simhit Process Type :: Nuclei :: t < 250 ns",    n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_N_proc000); 
  ME0_OH_proc000 = new TH1F("ME0_OH_proc000", "Simhit Process Type :: Other Hadrons :: t < 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_OH_proc000); 

  ME0_el_proc250 = new TH1F("ME0_el_proc250", "Simhit Process Type :: Electrons :: t > 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_el_proc250);
  ME0_mu_proc250 = new TH1F("ME0_mu_proc250", "Simhit Process Type :: Muons :: t > 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_mu_proc250); 
  ME0_pi_proc250 = new TH1F("ME0_pi_proc250", "Simhit Process Type :: Pions :: t > 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_pi_proc250); 
  ME0_ka_proc250 = new TH1F("ME0_ka_proc250", "Simhit Process Type :: Kaons :: t > 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_ka_proc250); 
  ME0_p_proc250  = new TH1F("ME0_p_proc250",  "Simhit Process Type :: Protons :: t > 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_p_proc250); 
  ME0_n_proc250  = new TH1F("ME0_n_proc250",  "Simhit Process Type :: Neutrons :: t > 250 ns",  n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_n_proc250); 
  ME0_g_proc250  = new TH1F("ME0_g_proc250",  "Simhit Process Type :: Photons :: t > 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_g_proc250); 
  ME0_N_proc250  = new TH1F("ME0_N_proc250",  "Simhit Process Type :: Nuclei :: t > 250 ns",    n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_N_proc250); 
  ME0_OH_proc250 = new TH1F("ME0_OH_proc250", "Simhit Process Type :: Other Hadrons :: t > 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(ME0_OH_proc250); 

  GEM_el_proc000 = new TH1F("GEM_el_proc000", "Simhit Process Type :: Electrons :: t < 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_el_proc000); 
  GEM_mu_proc000 = new TH1F("GEM_mu_proc000", "Simhit Process Type :: Muons :: t < 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_mu_proc000); 
  GEM_pi_proc000 = new TH1F("GEM_pi_proc000", "Simhit Process Type :: Pions :: t < 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_pi_proc000); 
  GEM_ka_proc000 = new TH1F("GEM_ka_proc000", "Simhit Process Type :: Kaons :: t < 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_ka_proc000); 
  GEM_p_proc000  = new TH1F("GEM_p_proc000",  "Simhit Process Type :: Protons :: t < 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_p_proc000); 
  GEM_n_proc000  = new TH1F("GEM_n_proc000",  "Simhit Process Type :: Neutrons :: t < 250 ns",  n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_n_proc000); 
  GEM_g_proc000  = new TH1F("GEM_g_proc000",  "Simhit Process Type :: Photons :: t < 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_g_proc000); 
  GEM_N_proc000  = new TH1F("GEM_N_proc000",  "Simhit Process Type :: Nuclei :: t < 250 ns",    n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_N_proc000); 
  GEM_OH_proc000 = new TH1F("GEM_OH_proc000", "Simhit Process Type :: Other Hadrons :: t < 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_OH_proc000); 

  GEM_el_proc250 = new TH1F("GEM_el_proc250", "Simhit Process Type :: Electrons :: t > 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_el_proc250); 
  GEM_mu_proc250 = new TH1F("GEM_mu_proc250", "Simhit Process Type :: Muons :: t > 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_mu_proc250); 
  GEM_pi_proc250 = new TH1F("GEM_pi_proc250", "Simhit Process Type :: Pions :: t > 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_pi_proc250); 
  GEM_ka_proc250 = new TH1F("GEM_ka_proc250", "Simhit Process Type :: Kaons :: t > 250 ns",     n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_ka_proc250); 
  GEM_p_proc250  = new TH1F("GEM_p_proc250",  "Simhit Process Type :: Protons :: t > 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_p_proc250); 
  GEM_n_proc250  = new TH1F("GEM_n_proc250",  "Simhit Process Type :: Neutrons :: t > 250 ns",  n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_n_proc250); 
  GEM_g_proc250  = new TH1F("GEM_g_proc250",  "Simhit Process Type :: Photons :: t > 250 ns",   n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_g_proc250); 
  GEM_N_proc250  = new TH1F("GEM_N_proc250",  "Simhit Process Type :: Nuclei :: t > 250 ns",    n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_N_proc250); 
  GEM_OH_proc250 = new TH1F("GEM_OH_proc250", "Simhit Process Type :: Other Hadrons :: t > 250 ns", n_pro, n1_pro, n2_pro);  VTH1F_Muon_hits_process.push_back(GEM_OH_proc250); 


  // Simhit Time vs Ekin
  GEM_el_hits = new TH2F("GEM_el_hits", "Simhit time vs E_{kin} :: GEM :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_mu_hits = new TH2F("GEM_mu_hits", "Simhit time vs E_{kin} :: GEM :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_pi_hits = new TH2F("GEM_pi_hits", "Simhit time vs E_{kin} :: GEM :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_ka_hits = new TH2F("GEM_ka_hits", "Simhit time vs E_{kin} :: GEM :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_p_hits  = new TH2F("GEM_p_hits",  "Simhit time vs E_{kin} :: GEM :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_n_hits  = new TH2F("GEM_n_hits",  "Simhit time vs E_{kin} :: GEM :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_g_hits  = new TH2F("GEM_g_hits",  "Simhit time vs E_{kin} :: GEM :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_N_hits  = new TH2F("GEM_N_hits",  "Simhit time vs E_{kin} :: GEM :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  GEM_OH_hits = new TH2F("GEM_OH_hits", "Simhit time vs E_{kin} :: GEM :: Other Hadrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  ME0_el_hits = new TH2F("ME0_el_hits", "Simhit time vs E_{kin} :: ME0 :: Electrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_mu_hits = new TH2F("ME0_mu_hits", "Simhit time vs E_{kin} :: ME0 :: Muons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_pi_hits = new TH2F("ME0_pi_hits", "Simhit time vs E_{kin} :: ME0 :: Pions",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_ka_hits = new TH2F("ME0_ka_hits", "Simhit time vs E_{kin} :: ME0 :: Kaons",     n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_p_hits  = new TH2F("ME0_p_hits",  "Simhit time vs E_{kin} :: ME0 :: Protons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_n_hits  = new TH2F("ME0_n_hits",  "Simhit time vs E_{kin} :: ME0 :: Neutrons",  n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_g_hits  = new TH2F("ME0_g_hits",  "Simhit time vs E_{kin} :: ME0 :: Photons",   n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_N_hits  = new TH2F("ME0_N_hits",  "Simhit time vs E_{kin} :: ME0 :: Nuclei",    n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);
  ME0_OH_hits = new TH2F("ME0_OH_hits", "Simhit time vs E_{kin} :: ME0 :: Other Hadrons", n_E, n1_E, n2_E, n_tof, n1_tof, n2_tof);

  // Ekin (1D)
  GEM_el_kins  = new TH1F("GEM_el_kins",  "E_{kin} :: GEM :: Electrons",  n_E, n1_E, n2_E);
  GEM_mu_kins  = new TH1F("GEM_mu_kins",  "E_{kin} :: GEM :: Muons",      n_E, n1_E, n2_E);
  GEM_ha_kins  = new TH1F("GEM_ha_kins",  "E_{kin} :: GEM :: Hadrons",    n_E, n1_E, n2_E);
  GEM_pi_kins  = new TH1F("GEM_pi_kins",  "E_{kin} :: GEM :: Pions",      n_E, n1_E, n2_E);
  GEM_ka_kins  = new TH1F("GEM_ka_kins",  "E_{kin} :: GEM :: Kaons",      n_E, n1_E, n2_E);
  GEM_p_kins   = new TH1F("GEM_p_kins",   "E_{kin} :: GEM :: Protons",    n_E, n1_E, n2_E);
  GEM_n_kins   = new TH1F("GEM_n_kins",   "E_{kin} :: GEM :: Neutrons",   n_E, n1_E, n2_E);
  GEM_g_kins   = new TH1F("GEM_g_kins",   "E_{kin} :: GEM :: Photons",    n_E, n1_E, n2_E);
  GEM_N_kins   = new TH1F("GEM_N_kins",   "E_{kin} :: GEM :: Nuclei",     n_E, n1_E, n2_E);
  GEM_OH_kins  = new TH1F("GEM_OH_kins",  "E_{kin} :: GEM :: Other Hadrons", n_E, n1_E, n2_E);
  GEM_All_kins = new TH1F("GEM_All_kins", "E_{kin} :: GEM :: All Particles", n_E, n1_E, n2_E);
  GEM_HIP_kins = new TH1F("GEM_HIP_kins", "E_{kin} :: GEM :: Highly Ionising Particles", n_E, n1_E, n2_E);

  ME0_el_kins  = new TH1F("ME0_el_kins",  "E_{kin} :: ME0 :: Electrons",  n_E, n1_E, n2_E);
  ME0_mu_kins  = new TH1F("ME0_mu_kins",  "E_{kin} :: ME0 :: Muons",      n_E, n1_E, n2_E);
  ME0_ha_kins  = new TH1F("ME0_ha_kins",  "E_{kin} :: ME0 :: Hadrons",    n_E, n1_E, n2_E);
  ME0_pi_kins  = new TH1F("ME0_pi_kins",  "E_{kin} :: ME0 :: Pions",      n_E, n1_E, n2_E);
  ME0_ka_kins  = new TH1F("ME0_ka_kins",  "E_{kin} :: ME0 :: Kaons",      n_E, n1_E, n2_E);
  ME0_p_kins   = new TH1F("ME0_p_kins",   "E_{kin} :: ME0 :: Protons",    n_E, n1_E, n2_E);
  ME0_n_kins   = new TH1F("ME0_n_kins",   "E_{kin} :: ME0 :: Neutrons",   n_E, n1_E, n2_E);
  ME0_g_kins   = new TH1F("ME0_g_kins",   "E_{kin} :: ME0 :: Photons",    n_E, n1_E, n2_E);
  ME0_N_kins   = new TH1F("ME0_N_kins",   "E_{kin} :: ME0 :: Nuclei",     n_E, n1_E, n2_E);
  ME0_OH_kins  = new TH1F("ME0_OH_kins",  "E_{kin} :: ME0 :: Other Hadrons", n_E, n1_E, n2_E);
  ME0_All_kins = new TH1F("ME0_All_kins",  "E_{kin} :: ME0 :: All Particles", n_E, n1_E, n2_E);
  ME0_HIP_kins = new TH1F("ME0_HIP_kins", "E_{kin} :: ME0 :: Highly Ionising Particles", n_E, n1_E, n2_E);

  G11_el_kins  = new TH1F("G11_el_kins",  "E_{kin} :: GE11 :: Electrons",  n_E, n1_E, n2_E);
  G11_mu_kins  = new TH1F("G11_mu_kins",  "E_{kin} :: GE11 :: Muons",      n_E, n1_E, n2_E);
  G11_ha_kins  = new TH1F("G11_ha_kins",  "E_{kin} :: GE11 :: Hadrons",    n_E, n1_E, n2_E);
  G11_pi_kins  = new TH1F("G11_pi_kins",  "E_{kin} :: GE11 :: Pions",      n_E, n1_E, n2_E);
  G11_ka_kins  = new TH1F("G11_ka_kins",  "E_{kin} :: GE11 :: Kaons",      n_E, n1_E, n2_E);
  G11_p_kins   = new TH1F("G11_p_kins",   "E_{kin} :: GE11 :: Protons",    n_E, n1_E, n2_E);
  G11_n_kins   = new TH1F("G11_n_kins",   "E_{kin} :: GE11 :: Neutrons",   n_E, n1_E, n2_E);
  G11_g_kins   = new TH1F("G11_g_kins",   "E_{kin} :: GE11 :: Photons",    n_E, n1_E, n2_E);
  G11_N_kins   = new TH1F("G11_N_kins",   "E_{kin} :: GE11 :: Nuclei",     n_E, n1_E, n2_E);
  G11_OH_kins  = new TH1F("G11_OH_kins",  "E_{kin} :: GE11 :: Other Hadrons", n_E, n1_E, n2_E);
  G11_All_kins = new TH1F("G11_All_kins", "E_{kin} :: GE11 :: All Particles", n_E, n1_E, n2_E);
  G11_HIP_kins = new TH1F("G11_HIP_kins", "E_{kin} :: GE11 :: Highly Ionising Particles", n_E, n1_E, n2_E);

  G21_el_kins  = new TH1F("G21_el_kins",  "E_{kin} :: GE21 :: Electrons",  n_E, n1_E, n2_E);
  G21_mu_kins  = new TH1F("G21_mu_kins",  "E_{kin} :: GE21 :: Muons",      n_E, n1_E, n2_E);
  G21_ha_kins  = new TH1F("G21_ha_kins",  "E_{kin} :: GE21 :: Hadrons",    n_E, n1_E, n2_E);
  G21_pi_kins  = new TH1F("G21_pi_kins",  "E_{kin} :: GE21 :: Pions",      n_E, n1_E, n2_E);
  G21_ka_kins  = new TH1F("G21_ka_kins",  "E_{kin} :: GE21 :: Kaons",      n_E, n1_E, n2_E);
  G21_p_kins   = new TH1F("G21_p_kins",   "E_{kin} :: GE21 :: Protons",    n_E, n1_E, n2_E);
  G21_n_kins   = new TH1F("G21_n_kins",   "E_{kin} :: GE21 :: Neutrons",   n_E, n1_E, n2_E);
  G21_g_kins   = new TH1F("G21_g_kins",   "E_{kin} :: GE21 :: Photons",    n_E, n1_E, n2_E);
  G21_N_kins   = new TH1F("G21_N_kins",   "E_{kin} :: GE21 :: Nuclei",     n_E, n1_E, n2_E);
  G21_OH_kins  = new TH1F("G21_OH_kins",  "E_{kin} :: GE21 :: Other Hadrons", n_E, n1_E, n2_E);
  G21_All_kins = new TH1F("G21_All_kins", "E_{kin} :: GE21 :: All Particles", n_E, n1_E, n2_E);
  G21_HIP_kins = new TH1F("G21_HIP_kins", "E_{kin} :: GE21 :: Highly Ionising Particles", n_E, n1_E, n2_E);

  GEM_el_linkin  = new TH1F("GEM_el_linkin",  "E_{kin} :: GEM :: Electrons",  n_EL, n1_EL, n2_EL);
  GEM_mu_linkin  = new TH1F("GEM_mu_linkin",  "E_{kin} :: GEM :: Muons",      n_EL, n1_EL, n2_EL);
  GEM_ha_linkin  = new TH1F("GEM_ha_linkin",  "E_{kin} :: GEM :: Hadrons",    n_EL, n1_EL, n2_EL);
  GEM_pi_linkin  = new TH1F("GEM_pi_linkin",  "E_{kin} :: GEM :: Pions",      n_EL, n1_EL, n2_EL);
  GEM_ka_linkin  = new TH1F("GEM_ka_linkin",  "E_{kin} :: GEM :: Kaons",      n_EL, n1_EL, n2_EL);
  GEM_p_linkin   = new TH1F("GEM_p_linkin",   "E_{kin} :: GEM :: Protons",    n_EL, n1_EL, n2_EL);
  GEM_n_linkin   = new TH1F("GEM_n_linkin",   "E_{kin} :: GEM :: Neutrons",   n_EL, n1_EL, n2_EL);
  GEM_g_linkin   = new TH1F("GEM_g_linkin",   "E_{kin} :: GEM :: Photons",    n_EL, n1_EL, n2_EL);
  GEM_N_linkin   = new TH1F("GEM_N_linkin",   "E_{kin} :: GEM :: Nuclei",     n_EL, n1_EL, n2_EL);
  GEM_OH_linkin  = new TH1F("GEM_OH_linkin",  "E_{kin} :: GEM :: Other Hadrons", n_EL, n1_EL, n2_EL);
  GEM_All_linkin = new TH1F("GEM_All_linkin", "E_{kin} :: GEM :: All Particles", n_EL, n1_EL, n2_EL);
  GEM_HIP_linkin = new TH1F("GEM_HIP_linkin", "E_{kin} :: GEM :: Highly Ionising Particles", n_EL, n1_EL, n2_EL);

  ME0_el_linkin  = new TH1F("ME0_el_linkin",  "E_{kin} :: ME0 :: Electrons",  n_EL, n1_EL, n2_EL);
  ME0_mu_linkin  = new TH1F("ME0_mu_linkin",  "E_{kin} :: ME0 :: Muons",      n_EL, n1_EL, n2_EL);
  ME0_ha_linkin  = new TH1F("ME0_ha_linkin",  "E_{kin} :: ME0 :: Hadrons",    n_EL, n1_EL, n2_EL);
  ME0_pi_linkin  = new TH1F("ME0_pi_linkin",  "E_{kin} :: ME0 :: Pions",      n_EL, n1_EL, n2_EL);
  ME0_ka_linkin  = new TH1F("ME0_ka_linkin",  "E_{kin} :: ME0 :: Kaons",      n_EL, n1_EL, n2_EL);
  ME0_p_linkin   = new TH1F("ME0_p_linkin",   "E_{kin} :: ME0 :: Protons",    n_EL, n1_EL, n2_EL);
  ME0_n_linkin   = new TH1F("ME0_n_linkin",   "E_{kin} :: ME0 :: Neutrons",   n_EL, n1_EL, n2_EL);
  ME0_g_linkin   = new TH1F("ME0_g_linkin",   "E_{kin} :: ME0 :: Photons",    n_EL, n1_EL, n2_EL);
  ME0_N_linkin   = new TH1F("ME0_N_linkin",   "E_{kin} :: ME0 :: Nuclei",     n_EL, n1_EL, n2_EL);
  ME0_OH_linkin  = new TH1F("ME0_OH_linkin",  "E_{kin} :: ME0 :: Other Hadrons", n_EL, n1_EL, n2_EL);
  ME0_All_linkin = new TH1F("ME0_All_linkin", "E_{kin} :: ME0 :: All Particles", n_EL, n1_EL, n2_EL);
  ME0_HIP_linkin = new TH1F("ME0_HIP_linkin", "E_{kin} :: ME0 :: Highly Ionising Particles", n_EL, n1_EL, n2_EL);

  G11_el_linkin  = new TH1F("G11_el_linkin",  "E_{kin} :: GE11 :: Electrons",  n_EL, n1_EL, n2_EL);
  G11_mu_linkin  = new TH1F("G11_mu_linkin",  "E_{kin} :: GE11 :: Muons",      n_EL, n1_EL, n2_EL);
  G11_ha_linkin  = new TH1F("G11_ha_linkin",  "E_{kin} :: GE11 :: Hadrons",    n_EL, n1_EL, n2_EL);
  G11_pi_linkin  = new TH1F("G11_pi_linkin",  "E_{kin} :: GE11 :: Pions",      n_EL, n1_EL, n2_EL);
  G11_ka_linkin  = new TH1F("G11_ka_linkin",  "E_{kin} :: GE11 :: Kaons",      n_EL, n1_EL, n2_EL);
  G11_p_linkin   = new TH1F("G11_p_linkin",   "E_{kin} :: GE11 :: Protons",    n_EL, n1_EL, n2_EL);
  G11_n_linkin   = new TH1F("G11_n_linkin",   "E_{kin} :: GE11 :: Neutrons",   n_EL, n1_EL, n2_EL);
  G11_g_linkin   = new TH1F("G11_g_linkin",   "E_{kin} :: GE11 :: Photons",    n_EL, n1_EL, n2_EL);
  G11_N_linkin   = new TH1F("G11_N_linkin",   "E_{kin} :: GE11 :: Nuclei",     n_EL, n1_EL, n2_EL);
  G11_OH_linkin  = new TH1F("G11_OH_linkin",  "E_{kin} :: GE11 :: Other Hadrons", n_EL, n1_EL, n2_EL);
  G11_All_linkin = new TH1F("G11_All_linkin", "E_{kin} :: GE11 :: All Particles", n_EL, n1_EL, n2_EL);
  G11_HIP_linkin = new TH1F("G11_HIP_linkin", "E_{kin} :: GE11 :: Highly Ionising Particles", n_EL, n1_EL, n2_EL);

  G21_el_linkin  = new TH1F("G21_el_linkin",  "E_{kin} :: GE21 :: Electrons",  n_EL, n1_EL, n2_EL);
  G21_mu_linkin  = new TH1F("G21_mu_linkin",  "E_{kin} :: GE21 :: Muons",      n_EL, n1_EL, n2_EL);
  G21_ha_linkin  = new TH1F("G21_ha_linkin",  "E_{kin} :: GE21 :: Hadrons",    n_EL, n1_EL, n2_EL);
  G21_pi_linkin  = new TH1F("G21_pi_linkin",  "E_{kin} :: GE21 :: Pions",      n_EL, n1_EL, n2_EL);
  G21_ka_linkin  = new TH1F("G21_ka_linkin",  "E_{kin} :: GE21 :: Kaons",      n_EL, n1_EL, n2_EL);
  G21_p_linkin   = new TH1F("G21_p_linkin",   "E_{kin} :: GE21 :: Protons",    n_EL, n1_EL, n2_EL);
  G21_n_linkin   = new TH1F("G21_n_linkin",   "E_{kin} :: GE21 :: Neutrons",   n_EL, n1_EL, n2_EL);
  G21_g_linkin   = new TH1F("G21_g_linkin",   "E_{kin} :: GE21 :: Photons",    n_EL, n1_EL, n2_EL);
  G21_N_linkin   = new TH1F("G21_N_linkin",   "E_{kin} :: GE21 :: Nuclei",     n_EL, n1_EL, n2_EL);
  G21_OH_linkin  = new TH1F("G21_OH_linkin",  "E_{kin} :: GE21 :: Other Hadrons", n_EL, n1_EL, n2_EL);
  G21_All_linkin = new TH1F("G21_All_linkin", "E_{kin} :: GE21 :: All Particles", n_EL, n1_EL, n2_EL);
  G21_HIP_linkin = new TH1F("G21_HIP_linkin", "E_{kin} :: GE21 :: Highly Ionising Particles", n_EL, n1_EL, n2_EL);

  // Simhit Time vs E deposit
  GEM_el_deposits = new TH2F("GEM_el_deposits", "Simhit time vs E_{deposit} :: GEM :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_mu_deposits = new TH2F("GEM_mu_deposits", "Simhit time vs E_{deposit} :: GEM :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_pi_deposits = new TH2F("GEM_pi_deposits", "Simhit time vs E_{deposit} :: GEM :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_ka_deposits = new TH2F("GEM_ka_deposits", "Simhit time vs E_{deposit} :: GEM :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_p_deposits  = new TH2F("GEM_p_deposits",  "Simhit time vs E_{deposit} :: GEM :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_n_deposits  = new TH2F("GEM_n_deposits",  "Simhit time vs E_{deposit} :: GEM :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_g_deposits  = new TH2F("GEM_g_deposits",  "Simhit time vs E_{deposit} :: GEM :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_N_deposits  = new TH2F("GEM_N_deposits",  "Simhit time vs E_{deposit} :: GEM :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  GEM_OH_deposits = new TH2F("GEM_OH_deposits", "Simhit time vs E_{deposit} :: GEM :: Other Hadrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  ME0_el_deposits = new TH2F("ME0_el_deposits", "Simhit time vs E_{deposit} :: ME0 :: Electrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_mu_deposits = new TH2F("ME0_mu_deposits", "Simhit time vs E_{deposit} :: ME0 :: Muons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_pi_deposits = new TH2F("ME0_pi_deposits", "Simhit time vs E_{deposit} :: ME0 :: Pions",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_ka_deposits = new TH2F("ME0_ka_deposits", "Simhit time vs E_{deposit} :: ME0 :: Kaons",     n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_p_deposits  = new TH2F("ME0_p_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Protons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_n_deposits  = new TH2F("ME0_n_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Neutrons",  n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_g_deposits  = new TH2F("ME0_g_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Photons",   n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_N_deposits  = new TH2F("ME0_N_deposits",  "Simhit time vs E_{deposit} :: ME0 :: Nuclei",    n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);
  ME0_OH_deposits = new TH2F("ME0_OH_deposits", "Simhit time vs E_{deposit} :: ME0 :: Other Hadrons", n_D, n1_D, n2_D, n_tof, n1_tof, n2_tof);

  // Deposit (1D)
  /*
  G11_el_deps  = new TH1F("G11_el_deps",  "E_{deposit} :: GE11 :: Electrons",  n_F, n1_F, n2_F);
  G11_mu_deps  = new TH1F("G11_mu_deps",  "E_{deposit} :: GE11 :: Muons",      n_F, n1_F, n2_F);
  G11_ha_deps  = new TH1F("G11_ha_deps",  "E_{deposit} :: GE11 :: Hadrons",    n_F, n1_F, n2_F);
  G11_pi_deps  = new TH1F("G11_pi_deps",  "E_{deposit} :: GE11 :: Pions",      n_F, n1_F, n2_F);
  G11_ka_deps  = new TH1F("G11_ka_deps",  "E_{deposit} :: GE11 :: Kaons",      n_F, n1_F, n2_F);
  G11_p_deps   = new TH1F("G11_p_deps",   "E_{deposit} :: GE11 :: Protons",    n_F, n1_F, n2_F);
  G11_n_deps   = new TH1F("G11_n_deps",   "E_{deposit} :: GE11 :: Neutrons",   n_F, n1_F, n2_F);
  G11_g_deps   = new TH1F("G11_g_deps",   "E_{deposit} :: GE11 :: Photons",    n_F, n1_F, n2_F);
  G11_N_deps   = new TH1F("G11_N_deps",   "E_{deposit} :: GE11 :: Nuclei",     n_F, n1_F, n2_F);
  G11_OH_deps  = new TH1F("G11_OH_deps",  "E_{deposit} :: GE11 :: Other Hadrons", n_F, n1_F, n2_F);
  G11_All_deps = new TH1F("G11_All_deps", "E_{deposit} :: GE11 :: All Particles", n_F, n1_F, n2_F);

  G21_el_deps  = new TH1F("G21_el_deps",  "E_{deposit} :: GE21 :: Electrons",  n_F, n1_F, n2_F);
  G21_mu_deps  = new TH1F("G21_mu_deps",  "E_{deposit} :: GE21 :: Muons",      n_F, n1_F, n2_F);
  G21_ha_deps  = new TH1F("G21_ha_deps",  "E_{deposit} :: GE21 :: Hadrons",    n_F, n1_F, n2_F);
  G21_pi_deps  = new TH1F("G21_pi_deps",  "E_{deposit} :: GE21 :: Pions",      n_F, n1_F, n2_F);
  G21_ka_deps  = new TH1F("G21_ka_deps",  "E_{deposit} :: GE21 :: Kaons",      n_F, n1_F, n2_F);
  G21_p_deps   = new TH1F("G21_p_deps",   "E_{deposit} :: GE21 :: Protons",    n_F, n1_F, n2_F);
  G21_n_deps   = new TH1F("G21_n_deps",   "E_{deposit} :: GE21 :: Neutrons",   n_F, n1_F, n2_F);
  G21_g_deps   = new TH1F("G21_g_deps",   "E_{deposit} :: GE21 :: Photons",    n_F, n1_F, n2_F);
  G21_N_deps   = new TH1F("G21_N_deps",   "E_{deposit} :: GE21 :: Nuclei",     n_F, n1_F, n2_F);
  G21_OH_deps  = new TH1F("G21_OH_deps",  "E_{deposit} :: GE21 :: Other Hadrons", n_F, n1_F, n2_F);
  G21_All_deps = new TH1F("G21_All_deps", "E_{deposit} :: GE21 :: All Particles", n_F, n1_F, n2_F);
  */

  GEM_el_deps  = new TH1F("GEM_el_deps",  "E_{deposit} :: GEM :: Electrons",  n_F, n1_F, n2_F);
  GEM_mu_deps  = new TH1F("GEM_mu_deps",  "E_{deposit} :: GEM :: Muons",      n_F, n1_F, n2_F);
  GEM_ha_deps  = new TH1F("GEM_ha_deps",  "E_{deposit} :: GEM :: Hadrons",    n_F, n1_F, n2_F);
  GEM_pi_deps  = new TH1F("GEM_pi_deps",  "E_{deposit} :: GEM :: Pions",      n_F, n1_F, n2_F);
  GEM_ka_deps  = new TH1F("GEM_ka_deps",  "E_{deposit} :: GEM :: Kaons",      n_F, n1_F, n2_F);
  GEM_p_deps   = new TH1F("GEM_p_deps",   "E_{deposit} :: GEM :: Protons",    n_F, n1_F, n2_F);
  GEM_n_deps   = new TH1F("GEM_n_deps",   "E_{deposit} :: GEM :: Neutrons",   n_F, n1_F, n2_F);
  GEM_g_deps   = new TH1F("GEM_g_deps",   "E_{deposit} :: GEM :: Photons",    n_F, n1_F, n2_F);
  GEM_N_deps   = new TH1F("GEM_N_deps",   "E_{deposit} :: GEM :: Nuclei",     n_F, n1_F, n2_F);
  GEM_OH_deps  = new TH1F("GEM_OH_deps",  "E_{deposit} :: GEM :: Other Hadrons", n_F, n1_F, n2_F);
  GEM_All_deps = new TH1F("GEM_All_deps", "E_{deposit} :: GEM :: All Particles", n_F, n1_F, n2_F);
  GEM_HIP_deps = new TH1F("GEM_HIP_deps", "E_{deposit} :: GEM :: Highly Ionizing", n_F, n1_F, n2_F);

  ME0_el_deps  = new TH1F("ME0_el_deps",  "E_{deposit} :: ME0 :: Electrons",  n_F, n1_F, n2_F);
  ME0_mu_deps  = new TH1F("ME0_mu_deps",  "E_{deposit} :: ME0 :: Muons",      n_F, n1_F, n2_F);
  ME0_ha_deps  = new TH1F("ME0_ha_deps",  "E_{deposit} :: ME0 :: Hadrons",    n_F, n1_F, n2_F);
  ME0_pi_deps  = new TH1F("ME0_pi_deps",  "E_{deposit} :: ME0 :: Pions",      n_F, n1_F, n2_F);
  ME0_ka_deps  = new TH1F("ME0_ka_deps",  "E_{deposit} :: ME0 :: Kaons",      n_F, n1_F, n2_F);
  ME0_p_deps   = new TH1F("ME0_p_deps",   "E_{deposit} :: ME0 :: Protons",    n_F, n1_F, n2_F);
  ME0_n_deps   = new TH1F("ME0_n_deps",   "E_{deposit} :: ME0 :: Neutrons",   n_F, n1_F, n2_F);
  ME0_g_deps   = new TH1F("ME0_g_deps",   "E_{deposit} :: ME0 :: Photons",    n_F, n1_F, n2_F);
  ME0_N_deps   = new TH1F("ME0_N_deps",   "E_{deposit} :: ME0 :: Nuclei",     n_F, n1_F, n2_F);
  ME0_OH_deps  = new TH1F("ME0_OH_deps",  "E_{deposit} :: ME0 :: Other Hadrons", n_F, n1_F, n2_F);
  ME0_All_deps = new TH1F("ME0_All_deps", "E_{deposit} :: ME0 :: All Particles", n_F, n1_F, n2_F);
  ME0_HIP_deps = new TH1F("ME0_HIP_deps", "E_{deposit} :: GEM :: Highly Ionizing", n_F, n1_F, n2_F);

  G11_el_deps  = new TH1F("G11_el_deps",  "E_{deposit} :: GE11 :: Electrons",  n_F, n1_F, n2_F);
  G11_mu_deps  = new TH1F("G11_mu_deps",  "E_{deposit} :: GE11 :: Muons",      n_F, n1_F, n2_F);
  G11_ha_deps  = new TH1F("G11_ha_deps",  "E_{deposit} :: GE11 :: Hadrons",    n_F, n1_F, n2_F);
  G11_pi_deps  = new TH1F("G11_pi_deps",  "E_{deposit} :: GE11 :: Pions",      n_F, n1_F, n2_F);
  G11_ka_deps  = new TH1F("G11_ka_deps",  "E_{deposit} :: GE11 :: Kaons",      n_F, n1_F, n2_F);
  G11_p_deps   = new TH1F("G11_p_deps",   "E_{deposit} :: GE11 :: Protons",    n_F, n1_F, n2_F);
  G11_n_deps   = new TH1F("G11_n_deps",   "E_{deposit} :: GE11 :: Neutrons",   n_F, n1_F, n2_F);
  G11_g_deps   = new TH1F("G11_g_deps",   "E_{deposit} :: GE11 :: Photons",    n_F, n1_F, n2_F);
  G11_N_deps   = new TH1F("G11_N_deps",   "E_{deposit} :: GE11 :: Nuclei",     n_F, n1_F, n2_F);
  G11_OH_deps  = new TH1F("G11_OH_deps",  "E_{deposit} :: GE11 :: Other Hadrons", n_F, n1_F, n2_F);
  G11_All_deps = new TH1F("G11_All_deps", "E_{deposit} :: GE11 :: All Particles", n_F, n1_F, n2_F);
  G11_HIP_deps = new TH1F("G11_HIP_deps", "E_{deposit} :: GE11 :: Highly Ionizing", n_F, n1_F, n2_F);

  G21_el_deps  = new TH1F("G21_el_deps",  "E_{deposit} :: GE21 :: Electrons",  n_F, n1_F, n2_F);
  G21_mu_deps  = new TH1F("G21_mu_deps",  "E_{deposit} :: GE21 :: Muons",      n_F, n1_F, n2_F);
  G21_ha_deps  = new TH1F("G21_ha_deps",  "E_{deposit} :: GE21 :: Hadrons",    n_F, n1_F, n2_F);
  G21_pi_deps  = new TH1F("G21_pi_deps",  "E_{deposit} :: GE21 :: Pions",      n_F, n1_F, n2_F);
  G21_ka_deps  = new TH1F("G21_ka_deps",  "E_{deposit} :: GE21 :: Kaons",      n_F, n1_F, n2_F);
  G21_p_deps   = new TH1F("G21_p_deps",   "E_{deposit} :: GE21 :: Protons",    n_F, n1_F, n2_F);
  G21_n_deps   = new TH1F("G21_n_deps",   "E_{deposit} :: GE21 :: Neutrons",   n_F, n1_F, n2_F);
  G21_g_deps   = new TH1F("G21_g_deps",   "E_{deposit} :: GE21 :: Photons",    n_F, n1_F, n2_F);
  G21_N_deps   = new TH1F("G21_N_deps",   "E_{deposit} :: GE21 :: Nuclei",     n_F, n1_F, n2_F);
  G21_OH_deps  = new TH1F("G21_OH_deps",  "E_{deposit} :: GE21 :: Other Hadrons", n_F, n1_F, n2_F);
  G21_All_deps = new TH1F("G21_All_deps", "E_{deposit} :: GE21 :: All Particles", n_F, n1_F, n2_F);
  G21_HIP_deps = new TH1F("G21_HIP_deps", "E_{deposit} :: GE21 :: Highly Ionizing", n_F, n1_F, n2_F);

  GEM_el_deps_000  = new TH1F("GEM_el_deps_000",  "E_{deposit} :: GEM :: Electrons :: t < 250 ns",  n_F, n1_F, n2_F);
  GEM_mu_deps_000  = new TH1F("GEM_mu_deps_000",  "E_{deposit} :: GEM :: Muons :: t < 250 ns",      n_F, n1_F, n2_F);
  GEM_ha_deps_000  = new TH1F("GEM_ha_deps_000",  "E_{deposit} :: GEM :: Hadrons :: t < 250 ns",    n_F, n1_F, n2_F);
  GEM_pi_deps_000  = new TH1F("GEM_pi_deps_000",  "E_{deposit} :: GEM :: Pions :: t < 250 ns",      n_F, n1_F, n2_F);
  GEM_ka_deps_000  = new TH1F("GEM_ka_deps_000",  "E_{deposit} :: GEM :: Kaons :: t < 250 ns",      n_F, n1_F, n2_F);
  GEM_p_deps_000   = new TH1F("GEM_p_deps_000",   "E_{deposit} :: GEM :: Protons :: t < 250 ns",    n_F, n1_F, n2_F);
  GEM_n_deps_000   = new TH1F("GEM_n_deps_000",   "E_{deposit} :: GEM :: Neutrons :: t < 250 ns",   n_F, n1_F, n2_F);
  GEM_g_deps_000   = new TH1F("GEM_g_deps_000",   "E_{deposit} :: GEM :: Photons :: t < 250 ns",    n_F, n1_F, n2_F);
  GEM_N_deps_000   = new TH1F("GEM_N_deps_000",   "E_{deposit} :: GEM :: Nuclei :: t < 250 ns",     n_F, n1_F, n2_F);
  GEM_OH_deps_000  = new TH1F("GEM_OH_deps_000",  "E_{deposit} :: GEM :: Other Hadrons :: t < 250 ns", n_F, n1_F, n2_F);

  ME0_el_deps_000  = new TH1F("ME0_el_deps_000",  "E_{deposit} :: ME0 :: Electrons :: t < 250 ns",  n_F, n1_F, n2_F);
  ME0_mu_deps_000  = new TH1F("ME0_mu_deps_000",  "E_{deposit} :: ME0 :: Muons :: t < 250 ns",      n_F, n1_F, n2_F);
  ME0_ha_deps_000  = new TH1F("ME0_ha_deps_000",  "E_{deposit} :: ME0 :: Hadrons :: t < 250 ns",    n_F, n1_F, n2_F);
  ME0_pi_deps_000  = new TH1F("ME0_pi_deps_000",  "E_{deposit} :: ME0 :: Pions :: t < 250 ns",      n_F, n1_F, n2_F);
  ME0_ka_deps_000  = new TH1F("ME0_ka_deps_000",  "E_{deposit} :: ME0 :: Kaons :: t < 250 ns",      n_F, n1_F, n2_F);
  ME0_p_deps_000   = new TH1F("ME0_p_deps_000",   "E_{deposit} :: ME0 :: Protons :: t < 250 ns",    n_F, n1_F, n2_F);
  ME0_n_deps_000   = new TH1F("ME0_n_deps_000",   "E_{deposit} :: ME0 :: Neutrons :: t < 250 ns",   n_F, n1_F, n2_F);
  ME0_g_deps_000   = new TH1F("ME0_g_deps_000",   "E_{deposit} :: ME0 :: Photons :: t < 250 ns",    n_F, n1_F, n2_F);
  ME0_N_deps_000   = new TH1F("ME0_N_deps_000",   "E_{deposit} :: ME0 :: Nuclei :: t < 250 ns",     n_F, n1_F, n2_F);
  ME0_OH_deps_000  = new TH1F("ME0_OH_deps_000",  "E_{deposit} :: ME0 :: Other Hadrons :: t < 250 ns", n_F, n1_F, n2_F);

  GEM_el_deps_250  = new TH1F("GEM_el_deps_250",  "E_{deposit} :: GEM :: Electrons :: t > 250 ns",  n_F, n1_F, n2_F);
  GEM_mu_deps_250  = new TH1F("GEM_mu_deps_250",  "E_{deposit} :: GEM :: Muons :: t > 250 ns",      n_F, n1_F, n2_F);
  GEM_ha_deps_250  = new TH1F("GEM_ha_deps_250",  "E_{deposit} :: GEM :: Hadrons :: t > 250 ns",    n_F, n1_F, n2_F);
  GEM_pi_deps_250  = new TH1F("GEM_pi_deps_250",  "E_{deposit} :: GEM :: Pions :: t > 250 ns",      n_F, n1_F, n2_F);
  GEM_ka_deps_250  = new TH1F("GEM_ka_deps_250",  "E_{deposit} :: GEM :: Kaons :: t > 250 ns",      n_F, n1_F, n2_F);
  GEM_p_deps_250   = new TH1F("GEM_p_deps_250",   "E_{deposit} :: GEM :: Protons :: t > 250 ns",    n_F, n1_F, n2_F);
  GEM_n_deps_250   = new TH1F("GEM_n_deps_250",   "E_{deposit} :: GEM :: Neutrons :: t > 250 ns",   n_F, n1_F, n2_F);
  GEM_g_deps_250   = new TH1F("GEM_g_deps_250",   "E_{deposit} :: GEM :: Photons :: t > 250 ns",    n_F, n1_F, n2_F);
  GEM_N_deps_250   = new TH1F("GEM_N_deps_250",   "E_{deposit} :: GEM :: Nuclei :: t > 250 ns",     n_F, n1_F, n2_F);
  GEM_OH_deps_250  = new TH1F("GEM_OH_deps_250",  "E_{deposit} :: GEM :: Other Hadrons :: t > 250 ns", n_F, n1_F, n2_F);

  ME0_el_deps_250  = new TH1F("ME0_el_deps_250",  "E_{deposit} :: ME0 :: Electrons :: t > 250 ns",  n_F, n1_F, n2_F);
  ME0_mu_deps_250  = new TH1F("ME0_mu_deps_250",  "E_{deposit} :: ME0 :: Muons :: t > 250 ns",      n_F, n1_F, n2_F);
  ME0_ha_deps_250  = new TH1F("ME0_ha_deps_250",  "E_{deposit} :: ME0 :: Hadrons :: t > 250 ns",    n_F, n1_F, n2_F);
  ME0_pi_deps_250  = new TH1F("ME0_pi_deps_250",  "E_{deposit} :: ME0 :: Pions :: t > 250 ns",      n_F, n1_F, n2_F);
  ME0_ka_deps_250  = new TH1F("ME0_ka_deps_250",  "E_{deposit} :: ME0 :: Kaons :: t > 250 ns",      n_F, n1_F, n2_F);
  ME0_p_deps_250   = new TH1F("ME0_p_deps_250",   "E_{deposit} :: ME0 :: Protons :: t > 250 ns",    n_F, n1_F, n2_F);
  ME0_n_deps_250   = new TH1F("ME0_n_deps_250",   "E_{deposit} :: ME0 :: Neutrons :: t > 250 ns",   n_F, n1_F, n2_F);
  ME0_g_deps_250   = new TH1F("ME0_g_deps_250",   "E_{deposit} :: ME0 :: Photons :: t > 250 ns",    n_F, n1_F, n2_F);
  ME0_N_deps_250   = new TH1F("ME0_N_deps_250",   "E_{deposit} :: ME0 :: Nuclei :: t > 250 ns",     n_F, n1_F, n2_F);
  ME0_OH_deps_250  = new TH1F("ME0_OH_deps_250",  "E_{deposit} :: ME0 :: Other Hadrons :: t > 250 ns", n_F, n1_F, n2_F);

  GEM_el_lindep  = new TH1F("GEM_el_lindep",  "E_{deposit} :: GEM :: Electrons",  n_DL, n1_DL, n2_DL);
  GEM_mu_lindep  = new TH1F("GEM_mu_lindep",  "E_{deposit} :: GEM :: Muons",      n_DL, n1_DL, n2_DL);
  GEM_ha_lindep  = new TH1F("GEM_ha_lindep",  "E_{deposit} :: GEM :: Hadrons",    n_DL, n1_DL, n2_DL);
  GEM_pi_lindep  = new TH1F("GEM_pi_lindep",  "E_{deposit} :: GEM :: Pions",      n_DL, n1_DL, n2_DL);
  GEM_ka_lindep  = new TH1F("GEM_ka_lindep",  "E_{deposit} :: GEM :: Kaons",      n_DL, n1_DL, n2_DL);
  GEM_p_lindep   = new TH1F("GEM_p_lindep",   "E_{deposit} :: GEM :: Protons",    n_DL, n1_DL, n2_DL);
  GEM_n_lindep   = new TH1F("GEM_n_lindep",   "E_{deposit} :: GEM :: Neutrons",   n_DL, n1_DL, n2_DL);
  GEM_g_lindep   = new TH1F("GEM_g_lindep",   "E_{deposit} :: GEM :: Photons",    n_DL, n1_DL, n2_DL);
  GEM_N_lindep   = new TH1F("GEM_N_lindep",   "E_{deposit} :: GEM :: Nuclei",     n_DL, n1_DL, n2_DL);
  GEM_OH_lindep  = new TH1F("GEM_OH_lindep",  "E_{deposit} :: GEM :: Other Hadrons", n_DL, n1_DL, n2_DL);
  GEM_All_lindep = new TH1F("GEM_All_lindep",  "E_{deposit} :: GEM :: All Particles", n_DL, n1_DL, n2_DL);
  GEM_HIP_lindep = new TH1F("GEM_HIP_lindep",  "E_{deposit} :: GEM :: Highly Ionising", n_H, n1_H, n2_H);

  ME0_el_lindep  = new TH1F("ME0_el_lindep",  "E_{deposit} :: ME0 :: Electrons",  n_DL, n1_DL, n2_DL);
  ME0_mu_lindep  = new TH1F("ME0_mu_lindep",  "E_{deposit} :: ME0 :: Muons",      n_DL, n1_DL, n2_DL);
  ME0_ha_lindep  = new TH1F("ME0_ha_lindep",  "E_{deposit} :: ME0 :: Hadrons",    n_DL, n1_DL, n2_DL);
  ME0_pi_lindep  = new TH1F("ME0_pi_lindep",  "E_{deposit} :: ME0 :: Pions",      n_DL, n1_DL, n2_DL);
  ME0_ka_lindep  = new TH1F("ME0_ka_lindep",  "E_{deposit} :: ME0 :: Kaons",      n_DL, n1_DL, n2_DL);
  ME0_p_lindep   = new TH1F("ME0_p_lindep",   "E_{deposit} :: ME0 :: Protons",    n_DL, n1_DL, n2_DL);
  ME0_n_lindep   = new TH1F("ME0_n_lindep",   "E_{deposit} :: ME0 :: Neutrons",   n_DL, n1_DL, n2_DL);
  ME0_g_lindep   = new TH1F("ME0_g_lindep",   "E_{deposit} :: ME0 :: Photons",    n_DL, n1_DL, n2_DL);
  ME0_N_lindep   = new TH1F("ME0_N_lindep",   "E_{deposit} :: ME0 :: Nuclei",     n_DL, n1_DL, n2_DL);
  ME0_OH_lindep  = new TH1F("ME0_OH_lindep",  "E_{deposit} :: ME0 :: Other Hadrons", n_DL, n1_DL, n2_DL);
  ME0_All_lindep = new TH1F("ME0_All_lindep",  "E_{deposit} :: ME0 :: All Particles", n_DL, n1_DL, n2_DL);
  ME0_HIP_lindep = new TH1F("ME0_HIP_lindep",  "E_{deposit} :: ME0 :: Highly Ionising", n_H, n1_H, n2_H);

  G11_el_lindep  = new TH1F("G11_el_lindep",  "E_{deposit} :: GE11 :: Electrons",  n_DL, n1_DL, n2_DL);
  G11_mu_lindep  = new TH1F("G11_mu_lindep",  "E_{deposit} :: GE11 :: Muons",      n_DL, n1_DL, n2_DL);
  G11_ha_lindep  = new TH1F("G11_ha_lindep",  "E_{deposit} :: GE11 :: Hadrons",    n_DL, n1_DL, n2_DL);
  G11_pi_lindep  = new TH1F("G11_pi_lindep",  "E_{deposit} :: GE11 :: Pions",      n_DL, n1_DL, n2_DL);
  G11_ka_lindep  = new TH1F("G11_ka_lindep",  "E_{deposit} :: GE11 :: Kaons",      n_DL, n1_DL, n2_DL);
  G11_p_lindep   = new TH1F("G11_p_lindep",   "E_{deposit} :: GE11 :: Protons",    n_DL, n1_DL, n2_DL);
  G11_n_lindep   = new TH1F("G11_n_lindep",   "E_{deposit} :: GE11 :: Neutrons",   n_DL, n1_DL, n2_DL);
  G11_g_lindep   = new TH1F("G11_g_lindep",   "E_{deposit} :: GE11 :: Photons",    n_DL, n1_DL, n2_DL);
  G11_N_lindep   = new TH1F("G11_N_lindep",   "E_{deposit} :: GE11 :: Nuclei",     n_DL, n1_DL, n2_DL);
  G11_OH_lindep  = new TH1F("G11_OH_lindep",  "E_{deposit} :: GE11 :: Other Hadrons", n_DL, n1_DL, n2_DL);
  G11_All_lindep = new TH1F("G11_All_lindep",  "E_{deposit} :: GE11 :: All Particles", n_DL, n1_DL, n2_DL);
  G11_HIP_lindep = new TH1F("G11_HIP_lindep",  "E_{deposit} :: GE11 :: Highly Ionising", n_H, n1_H, n2_H);

  G21_el_lindep  = new TH1F("G21_el_lindep",  "E_{deposit} :: GE21 :: Electrons",  n_DL, n1_DL, n2_DL);
  G21_mu_lindep  = new TH1F("G21_mu_lindep",  "E_{deposit} :: GE21 :: Muons",      n_DL, n1_DL, n2_DL);
  G21_ha_lindep  = new TH1F("G21_ha_lindep",  "E_{deposit} :: GE21 :: Hadrons",    n_DL, n1_DL, n2_DL);
  G21_pi_lindep  = new TH1F("G21_pi_lindep",  "E_{deposit} :: GE21 :: Pions",      n_DL, n1_DL, n2_DL);
  G21_ka_lindep  = new TH1F("G21_ka_lindep",  "E_{deposit} :: GE21 :: Kaons",      n_DL, n1_DL, n2_DL);
  G21_p_lindep   = new TH1F("G21_p_lindep",   "E_{deposit} :: GE21 :: Protons",    n_DL, n1_DL, n2_DL);
  G21_n_lindep   = new TH1F("G21_n_lindep",   "E_{deposit} :: GE21 :: Neutrons",   n_DL, n1_DL, n2_DL);
  G21_g_lindep   = new TH1F("G21_g_lindep",   "E_{deposit} :: GE21 :: Photons",    n_DL, n1_DL, n2_DL);
  G21_N_lindep   = new TH1F("G21_N_lindep",   "E_{deposit} :: GE21 :: Nuclei",     n_DL, n1_DL, n2_DL);
  G21_OH_lindep  = new TH1F("G21_OH_lindep",  "E_{deposit} :: GE21 :: Other Hadrons", n_DL, n1_DL, n2_DL);
  G21_All_lindep = new TH1F("G21_All_lindep",  "E_{deposit} :: GE21 :: All Particles", n_DL, n1_DL, n2_DL);
  G21_HIP_lindep = new TH1F("G21_HIP_lindep",  "E_{deposit} :: GE21 :: Highly Ionising", n_H, n1_H, n2_H);

  ME0_All_lindep_roll = new TH2F("ME0_All_lindep_roll",  "E_{deposit} vs EtaPart :: ME0  :: All Particles", n_CL, n1_CL, n2_CL, 8, 0.5, 8.5);
  G11_All_lindep_roll = new TH2F("G11_All_lindep_roll",  "E_{deposit} vs EtaPart :: GE11 :: All Particles", n_CL, n1_CL, n2_CL, 8, 0.5, 8.5);
  G21_All_lindep_roll = new TH2F("G21_All_lindep_roll",  "E_{deposit} vs EtaPart :: GE21 :: All Particles", n_CL, n1_CL, n2_CL, 8, 0.5, 8.5);

  ME0_All_lindep_eta01  = new TH1F("ME0_All_lindep_eta01",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 01", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta02  = new TH1F("ME0_All_lindep_eta02",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 02", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta03  = new TH1F("ME0_All_lindep_eta03",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 03", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta04  = new TH1F("ME0_All_lindep_eta04",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 04", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta05  = new TH1F("ME0_All_lindep_eta05",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 05", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta06  = new TH1F("ME0_All_lindep_eta06",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 06", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta07  = new TH1F("ME0_All_lindep_eta07",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 07", n_CL, n1_CL, n2_CL);
  ME0_All_lindep_eta08  = new TH1F("ME0_All_lindep_eta08",  "E_{deposit} :: ME0 :: All Particles :: EtaPart 08", n_CL, n1_CL, n2_CL);

  G11_All_lindep_eta01  = new TH1F("G11_All_lindep_eta01",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 01", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta02  = new TH1F("G11_All_lindep_eta02",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 02", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta03  = new TH1F("G11_All_lindep_eta03",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 03", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta04  = new TH1F("G11_All_lindep_eta04",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 04", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta05  = new TH1F("G11_All_lindep_eta05",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 05", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta06  = new TH1F("G11_All_lindep_eta06",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 06", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta07  = new TH1F("G11_All_lindep_eta07",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 07", n_CL, n1_CL, n2_CL);
  G11_All_lindep_eta08  = new TH1F("G11_All_lindep_eta08",  "E_{deposit} :: GE11 :: All Particles :: EtaPart 08", n_CL, n1_CL, n2_CL);

  G21_All_lindep_eta01  = new TH1F("G21_All_lindep_eta01",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 01", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta02  = new TH1F("G21_All_lindep_eta02",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 02", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta03  = new TH1F("G21_All_lindep_eta03",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 03", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta04  = new TH1F("G21_All_lindep_eta04",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 04", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta05  = new TH1F("G21_All_lindep_eta05",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 05", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta06  = new TH1F("G21_All_lindep_eta06",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 06", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta07  = new TH1F("G21_All_lindep_eta07",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 07", n_CL, n1_CL, n2_CL);
  G21_All_lindep_eta08  = new TH1F("G21_All_lindep_eta08",  "E_{deposit} :: GE21 :: All Particles :: EtaPart 08", n_CL, n1_CL, n2_CL);


  // Ekin vs E deposit
  GEM_el_ekindep = new TH2F("GEM_el_ekindep", "Simhit E_{deposit} vs E_{kin} :: GEM :: Electrons", n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_mu_ekindep = new TH2F("GEM_mu_ekindep", "Simhit E_{deposit} vs E_{kin} :: GEM :: Muons",     n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_pi_ekindep = new TH2F("GEM_pi_ekindep", "Simhit E_{deposit} vs E_{kin} :: GEM :: Pions",     n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_ka_ekindep = new TH2F("GEM_ka_ekindep", "Simhit E_{deposit} vs E_{kin} :: GEM :: Kaons",     n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_p_ekindep  = new TH2F("GEM_p_ekindep",  "Simhit E_{deposit} vs E_{kin} :: GEM :: Protons",   n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_n_ekindep  = new TH2F("GEM_n_ekindep",  "Simhit E_{deposit} vs E_{kin} :: GEM :: Neutrons",  n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_g_ekindep  = new TH2F("GEM_g_ekindep",  "Simhit E_{deposit} vs E_{kin} :: GEM :: Photons",   n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_N_ekindep  = new TH2F("GEM_N_ekindep",  "Simhit E_{deposit} vs E_{kin} :: GEM :: Nuclei",    n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  GEM_OH_ekindep = new TH2F("GEM_OH_ekindep", "Simhit E_{deposit} vs E_{kin} :: GEM :: Other Hadrons", n_E, n1_E, n2_E, n_D, n1_D, n2_D);

  ME0_el_ekindep = new TH2F("ME0_el_ekindep", "Simhit E_{deposit} vs E_{kin} :: ME0 :: Electrons", n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_mu_ekindep = new TH2F("ME0_mu_ekindep", "Simhit E_{deposit} vs E_{kin} :: ME0 :: Muons",     n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_pi_ekindep = new TH2F("ME0_pi_ekindep", "Simhit E_{deposit} vs E_{kin} :: ME0 :: Pions",     n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_ka_ekindep = new TH2F("ME0_ka_ekindep", "Simhit E_{deposit} vs E_{kin} :: ME0 :: Kaons",     n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_p_ekindep  = new TH2F("ME0_p_ekindep",  "Simhit E_{deposit} vs E_{kin} :: ME0 :: Protons",   n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_n_ekindep  = new TH2F("ME0_n_ekindep",  "Simhit E_{deposit} vs E_{kin} :: ME0 :: Neutrons",  n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_g_ekindep  = new TH2F("ME0_g_ekindep",  "Simhit E_{deposit} vs E_{kin} :: ME0 :: Photons",   n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_N_ekindep  = new TH2F("ME0_N_ekindep",  "Simhit E_{deposit} vs E_{kin} :: ME0 :: Nuclei",    n_E, n1_E, n2_E, n_D, n1_D, n2_D);
  ME0_OH_ekindep = new TH2F("ME0_OH_ekindep", "Simhit E_{deposit} vs E_{kin} :: ME0 :: Other Hadrons", n_E, n1_E, n2_E, n_D, n1_D, n2_D);

  // TOF (1D)
  GEM_el_tof  = new TH1F("GEM_el_tof",  "Time Of Flight :: GEM :: Electrons",  n_tof, n1_tof, n2_tof);
  GEM_mu_tof  = new TH1F("GEM_mu_tof",  "Time Of Flight :: GEM :: Muons",      n_tof, n1_tof, n2_tof);
  GEM_ha_tof  = new TH1F("GEM_ha_tof",  "Time Of Flight :: GEM :: Hadrons",    n_tof, n1_tof, n2_tof);
  GEM_pi_tof  = new TH1F("GEM_pi_tof",  "Time Of Flight :: GEM :: Pions",      n_tof, n1_tof, n2_tof);
  GEM_ka_tof  = new TH1F("GEM_ka_tof",  "Time Of Flight :: GEM :: Kaons",      n_tof, n1_tof, n2_tof);
  GEM_p_tof   = new TH1F("GEM_p_tof",   "Time Of Flight :: GEM :: Protons",    n_tof, n1_tof, n2_tof);
  GEM_n_tof   = new TH1F("GEM_n_tof",   "Time Of Flight :: GEM :: Neutrons",   n_tof, n1_tof, n2_tof);
  GEM_g_tof   = new TH1F("GEM_g_tof",   "Time Of Flight :: GEM :: Photons",    n_tof, n1_tof, n2_tof);
  GEM_N_tof   = new TH1F("GEM_N_tof",   "Time Of Flight :: GEM :: Nuclei",     n_tof, n1_tof, n2_tof);
  GEM_OH_tof  = new TH1F("GEM_OH_tof",  "Time Of Flight :: GEM :: Other Hadrons", n_tof, n1_tof, n2_tof);

  ME0_el_tof  = new TH1F("ME0_el_tof",  "Time Of Flight :: ME0 :: Electrons",  n_tof, n1_tof, n2_tof);
  ME0_mu_tof  = new TH1F("ME0_mu_tof",  "Time Of Flight :: ME0 :: Muons",      n_tof, n1_tof, n2_tof);
  ME0_ha_tof  = new TH1F("ME0_ha_tof",  "Time Of Flight :: ME0 :: Hadrons",    n_tof, n1_tof, n2_tof);
  ME0_pi_tof  = new TH1F("ME0_pi_tof",  "Time Of Flight :: ME0 :: Pions",      n_tof, n1_tof, n2_tof);
  ME0_ka_tof  = new TH1F("ME0_ka_tof",  "Time Of Flight :: ME0 :: Kaons",      n_tof, n1_tof, n2_tof);
  ME0_p_tof   = new TH1F("ME0_p_tof",   "Time Of Flight :: ME0 :: Protons",    n_tof, n1_tof, n2_tof);
  ME0_n_tof   = new TH1F("ME0_n_tof",   "Time Of Flight :: ME0 :: Neutrons",   n_tof, n1_tof, n2_tof);
  ME0_g_tof   = new TH1F("ME0_g_tof",   "Time Of Flight :: ME0 :: Photons",    n_tof, n1_tof, n2_tof);
  ME0_N_tof   = new TH1F("ME0_N_tof",   "Time Of Flight :: ME0 :: Nuclei",     n_tof, n1_tof, n2_tof);
  ME0_OH_tof  = new TH1F("ME0_OH_tof",  "Time Of Flight :: ME0 :: Other Hadrons", n_tof, n1_tof, n2_tof);

  // Time (1D)
  GEM_el_time  = new TH1F("GEM_el_time",  "Time Of Flight :: GEM :: Electrons",  n_time, n1_time, n2_time);
  GEM_mu_time  = new TH1F("GEM_mu_time",  "Time Of Flight :: GEM :: Muons",      n_time, n1_time, n2_time);
  GEM_ha_time  = new TH1F("GEM_ha_time",  "Time Of Flight :: GEM :: Hadrons",    n_time, n1_time, n2_time);
  GEM_pi_time  = new TH1F("GEM_pi_time",  "Time Of Flight :: GEM :: Pions",      n_time, n1_time, n2_time);
  GEM_ka_time  = new TH1F("GEM_ka_time",  "Time Of Flight :: GEM :: Kaons",      n_time, n1_time, n2_time);
  GEM_p_time   = new TH1F("GEM_p_time",   "Time Of Flight :: GEM :: Protons",    n_time, n1_time, n2_time);
  GEM_n_time   = new TH1F("GEM_n_time",   "Time Of Flight :: GEM :: Neutrons",   n_time, n1_time, n2_time);
  GEM_g_time   = new TH1F("GEM_g_time",   "Time Of Flight :: GEM :: Photons",    n_time, n1_time, n2_time);
  GEM_N_time   = new TH1F("GEM_N_time",   "Time Of Flight :: GEM :: Nuclei",     n_time, n1_time, n2_time);
  GEM_OH_time  = new TH1F("GEM_OH_time",  "Time Of Flight :: GEM :: Other Hadrons", n_time, n1_time, n2_time);

  ME0_el_time  = new TH1F("ME0_el_time",  "Time Of Flight :: ME0 :: Electrons",  n_time, n1_time, n2_time);
  ME0_mu_time  = new TH1F("ME0_mu_time",  "Time Of Flight :: ME0 :: Muons",      n_time, n1_time, n2_time);
  ME0_ha_time  = new TH1F("ME0_ha_time",  "Time Of Flight :: ME0 :: Hadrons",    n_time, n1_time, n2_time);
  ME0_pi_time  = new TH1F("ME0_pi_time",  "Time Of Flight :: ME0 :: Pions",      n_time, n1_time, n2_time);
  ME0_ka_time  = new TH1F("ME0_ka_time",  "Time Of Flight :: ME0 :: Kaons",      n_time, n1_time, n2_time);
  ME0_p_time   = new TH1F("ME0_p_time",   "Time Of Flight :: ME0 :: Protons",    n_time, n1_time, n2_time);
  ME0_n_time   = new TH1F("ME0_n_time",   "Time Of Flight :: ME0 :: Neutrons",   n_time, n1_time, n2_time);
  ME0_g_time   = new TH1F("ME0_g_time",   "Time Of Flight :: ME0 :: Photons",    n_time, n1_time, n2_time);
  ME0_N_time   = new TH1F("ME0_N_time",   "Time Of Flight :: ME0 :: Nuclei",     n_time, n1_time, n2_time);
  ME0_OH_time  = new TH1F("ME0_OH_time",  "Time Of Flight :: ME0 :: Other Hadrons", n_time, n1_time, n2_time);

  GEM_XY  = new TH2F("GEM_XY",  "Simhits in XY :: GEM", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  GEM_RZ  = new TH2F("GEM_RZ",  "Simhits in RZ :: GEM", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  ME0_XY  = new TH2F("ME0_XY",  "Simhits in XY :: ME0", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  ME0_RZ  = new TH2F("ME0_RZ",  "Simhits in RZ :: ME0", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);

  GEM_000ns_XY  = new TH2F("GEM_000ns_XY",  "Simhits with tof < 250ns in XY :: GEM", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  GEM_000ns_RZ  = new TH2F("GEM_000ns_RZ",  "Simhits with tof < 250ns in RZ :: GEM", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  ME0_000ns_XY  = new TH2F("ME0_000ns_XY",  "Simhits with tof < 250ns in XY :: ME0", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  ME0_000ns_RZ  = new TH2F("ME0_000ns_RZ",  "Simhits with tof < 250ns in RZ :: ME0", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  GEM_250ns_XY  = new TH2F("GEM_250ns_XY",  "Simhits with tof > 250ns in XY :: GEM", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  GEM_250ns_RZ  = new TH2F("GEM_250ns_RZ",  "Simhits with tof > 250ns in RZ :: GEM", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);
  ME0_250ns_XY  = new TH2F("ME0_250ns_XY",  "Simhits with tof > 250ns in XY :: ME0", n_xy_x, n_xy_x1, n_xy_x2, n_xy_y, n_xy_y1, n_xy_y2);
  ME0_250ns_RZ  = new TH2F("ME0_250ns_RZ",  "Simhits with tof > 250ns in RZ :: ME0", n_zr_z, n_zr_z1, n_zr_z2, n_zr_r, n_zr_r1, n_zr_r2);

  GEM_hits_tof  = new TH1F("GEM_hits_tof",  "Simhit time :: GEM",  m_tof, m1_tof, m2_tof);
  GEM_hits_eta  = new TH1F("GEM_hits_eta",  "Simhit time :: GEM",  m_eta, m1_eta, m2_eta);
  GEM_hits_phi  = new TH1F("GEM_hits_phi",  "Simhit time :: GEM",  m_phi, m1_phi, m2_phi);
  GEM_hits_lin  = new TH1F("GEM_hits_lin",  "Simhit time :: GEM",  m_lin, m1_lin, m2_lin);

  ME0_hits_tof  = new TH1F("ME0_hits_tof",  "Simhit time :: ME0",  m_tof, m1_tof, m2_tof);
  ME0_hits_eta  = new TH1F("ME0_hits_eta",  "Simhit time :: ME0",  m_eta, m1_eta, m2_eta);
  ME0_hits_phi  = new TH1F("ME0_hits_phi",  "Simhit time :: ME0",  m_phi, m1_phi, m2_phi);
  ME0_hits_lin  = new TH1F("ME0_hits_lin",  "Simhit time :: ME0",  m_lin, m1_lin, m2_lin);

  ME0_Nuclei_A_Z = new TH2F("ME0_Nuclei_A_Z", "Nuclei A vs Z :: ME0", n_Z, n1_Z, n2_Z, n_A, n1_A, n2_A);
  ME0_Nuclei_A   = new TH1F("ME0_Nuclei_A",   "Nuclei A dist :: ME0", n_A, n1_A, n2_A); 
  ME0_Nuclei_Z   = new TH1F("ME0_Nuclei_Z",   "Nuclei Z dist :: ME0", n_Z, n1_Z, n2_Z); 
  ME0_Nuclei_List= new TH1F("ME0_Nuclei_List","Nuclei List :: ME0", 19, 0.5, 19.5);
  ME0_HIP_id     = new TH1F("ME0_HIP_id",     "Highly Ionising particle ID :: ME0", 9, 0.5, 9.5);
  GEM_HIP_id     = new TH1F("GEM_HIP_id",     "Highly Ionising particle ID :: GEM", 9, 0.5, 9.5);
}


MyME0SimHitAnalyzer::~MyME0SimHitAnalyzer(){

  if(tech_debug) std::cout<<"[MyME0SimHitAnalyzer :: Destructor]"<<std::endl; 
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
  // First Plot for PDF File :: print empty Dummy 
  pdfFileName = pdfFileNameBase + ".pdf[";
  Dummy->Print(pdfFileName.c_str());
  // Name for next plot for PDF File
  pdfFileName = pdfFileNameBase + ".pdf";

  // Event Counter
  for(int i=0; i<7; ++i) {
    EventCounter->GetXaxis()->SetBinLabel(i+1, evtcnt[i].c_str());
  }
  EventCounter->Write();


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
  // M11_Od_el_hits_R->Write(); M11_Od_mu_hits_R->Write(); M11_Od_pi_hits_R->Write(); M11_Od_ka_hits_R->Write(); M11_Od_p_hits_R->Write();
  // M11_Od_n_hits_R->Write(); M11_Od_g_hits_R->Write(); M11_Od_N_hits_R->Write(); M11_Od_OH_hits_R->Write(); M11_Od_All_hits_R->Write(); M11_Od_HIP_hits_R->Write();
  // M11_Ev_el_hits_R->Write(); M11_Ev_mu_hits_R->Write(); M11_Ev_pi_hits_R->Write(); M11_Ev_ka_hits_R->Write(); M11_Ev_p_hits_R->Write();
  // M11_Ev_n_hits_R->Write(); M11_Ev_g_hits_R->Write(); M11_Ev_N_hits_R->Write(); M11_Ev_OH_hits_R->Write(); M11_Ev_All_hits_R->Write(); M11_Ev_HIP_hits_R->Write();
  outputfile->cd();
  // ---------------------------

  TDir_Muon_250ns_radius->cd();
  // ---------------------------
  ME0_el_hits_000_R->Write();  ME0_mu_hits_000_R->Write();  ME0_pi_hits_000_R->Write();  ME0_ka_hits_000_R->Write();  ME0_p_hits_000_R->Write();
  ME0_n_hits_000_R->Write();   ME0_g_hits_000_R->Write();   ME0_N_hits_000_R->Write();  ME0_OH_hits_000_R->Write();   ME0_All_hits_000_R->Write();
  ME0_el_hits_250_R->Write();  ME0_mu_hits_250_R->Write();  ME0_pi_hits_250_R->Write();  ME0_ka_hits_250_R->Write();  ME0_p_hits_250_R->Write();  
  ME0_n_hits_250_R->Write();   ME0_g_hits_250_R->Write();   ME0_N_hits_250_R->Write();   ME0_OH_hits_250_R->Write();  ME0_All_hits_250_R->Write();
  G11_el_hits_000_R->Write();  G11_mu_hits_000_R->Write();  G11_pi_hits_000_R->Write();  G11_ka_hits_000_R->Write();  G11_p_hits_000_R->Write();
  G11_n_hits_000_R->Write();   G11_g_hits_000_R->Write();   G11_N_hits_000_R->Write();   G11_OH_hits_000_R->Write();  G11_All_hits_000_R->Write();
  G11_el_hits_250_R->Write();  G11_mu_hits_250_R->Write();  G11_pi_hits_250_R->Write();  G11_ka_hits_250_R->Write();  G11_p_hits_250_R->Write();
  G11_n_hits_250_R->Write();   G11_g_hits_250_R->Write();   G11_N_hits_250_R->Write();   G11_OH_hits_250_R->Write();  G11_All_hits_250_R->Write();
  G21_el_hits_000_R->Write();  G21_mu_hits_000_R->Write();  G21_pi_hits_000_R->Write();  G21_ka_hits_000_R->Write();  G21_p_hits_000_R->Write();
  G21_n_hits_000_R->Write();   G21_g_hits_000_R->Write();   G21_N_hits_000_R->Write();   G21_OH_hits_000_R->Write();  G21_All_hits_000_R->Write();
  G21_el_hits_250_R->Write();  G21_mu_hits_250_R->Write();  G21_pi_hits_250_R->Write();  G21_ka_hits_250_R->Write();  G21_p_hits_250_R->Write();
  G21_n_hits_250_R->Write();   G21_g_hits_250_R->Write();   G21_N_hits_250_R->Write();   G21_OH_hits_250_R->Write();  G21_All_hits_250_R->Write();
  outputfile->cd();
  // ---------------------------

  TDir_Muon_25ns_radius->cd();
  // ---------------------------
  ME0_el_hits_00_R->Write();  ME0_mu_hits_00_R->Write();  ME0_pi_hits_00_R->Write();  ME0_ka_hits_00_R->Write();  ME0_p_hits_00_R->Write();
  ME0_n_hits_00_R->Write();   ME0_g_hits_00_R->Write();   ME0_N_hits_00_R->Write();   ME0_OH_hits_00_R->Write();  ME0_All_hits_00_R->Write();
  ME0_el_hits_25_R->Write();  ME0_mu_hits_25_R->Write();  ME0_pi_hits_25_R->Write();  ME0_ka_hits_25_R->Write();  ME0_p_hits_25_R->Write();  
  ME0_n_hits_25_R->Write();   ME0_g_hits_25_R->Write();   ME0_N_hits_25_R->Write();   ME0_OH_hits_25_R->Write();  ME0_All_hits_25_R->Write();
  G11_el_hits_00_R->Write();  G11_mu_hits_00_R->Write();  G11_pi_hits_00_R->Write();  G11_ka_hits_00_R->Write();  G11_p_hits_00_R->Write();
  G11_n_hits_00_R->Write();   G11_g_hits_00_R->Write();   G11_N_hits_00_R->Write();   G11_OH_hits_00_R->Write();  G11_All_hits_00_R->Write();
  G11_el_hits_25_R->Write();  G11_mu_hits_25_R->Write();  G11_pi_hits_25_R->Write();  G11_ka_hits_25_R->Write();  G11_p_hits_25_R->Write();
  G11_n_hits_25_R->Write();   G11_g_hits_25_R->Write();   G11_N_hits_25_R->Write();   G11_OH_hits_25_R->Write();  G11_All_hits_25_R->Write();
  G21_el_hits_00_R->Write();  G21_mu_hits_00_R->Write();  G21_pi_hits_00_R->Write();  G21_ka_hits_00_R->Write();  G21_p_hits_00_R->Write();
  G21_n_hits_00_R->Write();   G21_g_hits_00_R->Write();   G21_N_hits_00_R->Write();   G21_OH_hits_00_R->Write();  G21_All_hits_00_R->Write();
  G21_el_hits_25_R->Write();  G21_mu_hits_25_R->Write();  G21_pi_hits_25_R->Write();  G21_ka_hits_25_R->Write();  G21_p_hits_25_R->Write();
  G21_n_hits_25_R->Write();   G21_g_hits_25_R->Write();   G21_N_hits_25_R->Write();   G21_OH_hits_25_R->Write();  G21_All_hits_25_R->Write();
  outputfile->cd();
  // ---------------------------

  TDir_Muon_hits_etapart->cd();
  // ---------------------------
  ME0_el_hits_E->Write(); ME0_mu_hits_E->Write(); ME0_pi_hits_E->Write(); ME0_ka_hits_E->Write(); ME0_p_hits_E->Write();
  ME0_n_hits_E->Write();  ME0_g_hits_E->Write();  ME0_N_hits_E->Write();  ME0_OH_hits_E->Write(); ME0_All_hits_E->Write(); ME0_HIP_hits_E->Write();
  G11_el_hits_E->Write(); G11_mu_hits_E->Write(); G11_pi_hits_E->Write(); G11_ka_hits_E->Write(); G11_p_hits_E->Write();
  G11_n_hits_E->Write();  G11_g_hits_E->Write();  G11_N_hits_E->Write();  G11_OH_hits_E->Write(); G11_All_hits_E->Write(); G11_HIP_hits_E->Write();
  G21_el_hits_E->Write(); G21_mu_hits_E->Write(); G21_pi_hits_E->Write(); G21_ka_hits_E->Write(); G21_p_hits_E->Write();
  G21_n_hits_E->Write();  G21_g_hits_E->Write();  G21_N_hits_E->Write();  G21_OH_hits_E->Write(); G21_All_hits_E->Write(); G21_HIP_hits_E->Write();

  ME0_All_hits_000_E->Write(); ME0_All_hits_250_E->Write(); ME0_All_hits_00_E->Write(); ME0_All_hits_25_E->Write(); 
  G11_All_hits_000_E->Write(); G11_All_hits_250_E->Write(); G11_All_hits_00_E->Write(); G11_All_hits_25_E->Write(); 
  G21_All_hits_000_E->Write(); G21_All_hits_250_E->Write(); G21_All_hits_00_E->Write(); G21_All_hits_25_E->Write(); 

  G11_L1Odd_el_hits_E->Write(); G11_L1Odd_mu_hits_E->Write(); G11_L1Odd_pi_hits_E->Write(); G11_L1Odd_ka_hits_E->Write(); G11_L1Odd_p_hits_E->Write();
  G11_L1Odd_n_hits_E->Write(); G11_L1Odd_g_hits_E->Write(); G11_L1Odd_N_hits_E->Write(); G11_L1Odd_OH_hits_E->Write(); G11_L1Odd_All_hits_E->Write(); G11_L1Odd_HIP_hits_E->Write();
  G11_L1Even_el_hits_E->Write(); G11_L1Even_mu_hits_E->Write(); G11_L1Even_pi_hits_E->Write(); G11_L1Even_ka_hits_E->Write(); G11_L1Even_p_hits_E->Write();
  G11_L1Even_n_hits_E->Write(); G11_L1Even_g_hits_E->Write(); G11_L1Even_N_hits_E->Write(); G11_L1Even_OH_hits_E->Write(); G11_L1Even_All_hits_E->Write(); G11_L1Even_HIP_hits_E->Write();
  G11_L2Odd_el_hits_E->Write(); G11_L2Odd_mu_hits_E->Write(); G11_L2Odd_pi_hits_E->Write(); G11_L2Odd_ka_hits_E->Write(); G11_L2Odd_p_hits_E->Write();
  G11_L2Odd_n_hits_E->Write(); G11_L2Odd_g_hits_E->Write(); G11_L2Odd_N_hits_E->Write(); G11_L2Odd_OH_hits_E->Write(); G11_L2Odd_All_hits_E->Write(); G11_L2Odd_HIP_hits_E->Write();
  G11_L2Even_el_hits_E->Write(); G11_L2Even_mu_hits_E->Write(); G11_L2Even_pi_hits_E->Write(); G11_L2Even_ka_hits_E->Write(); G11_L2Even_p_hits_E->Write();
  G11_L2Even_n_hits_E->Write(); G11_L2Even_g_hits_E->Write(); G11_L2Even_N_hits_E->Write(); G11_L2Even_OH_hits_E->Write(); G11_L2Even_All_hits_E->Write(); G11_L2Even_HIP_hits_E->Write();


  outputfile->cd();
  // ---------------------------

  TDir_Muon_hits_process->cd();
  // ---------------------------
  for(int i=0; i<9*2*3; ++i) {
    // std::cout<<"Muon Hits Process Histogram number ="<<i<<std::endl;
    for(int m=0; m<31; ++m) {
      VTH1F_Muon_hits_process[i]->GetXaxis()->SetBinLabel(m+1, proc[m].c_str());
    }
    VTH1F_Muon_hits_process[i]->Write();
  }
  // ---------------------------
  outputfile->cd();


  // Center of Mass Energy Label
  std::stringstream comlabelss; comlabelss<<"CMS Simulation #sqrt{s} = "<<comenergy<<" TeV"; std::string comlabel = comlabelss.str();

  // Legends
  double l1_x1, l1_y1, l1_x2, l1_y2;
  l1_x1 = 0.65; l1_x2 = 0.85; l1_y1 = 0.65; l1_y2 = 0.85;
  TLegend *l1 = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  l1->SetLineColor(0);    l1->SetLineStyle(0);  l1->SetLineWidth(0);
  l1->SetFillColor(4000); l1->SetBorderSize(0); l1->SetNColumns(2);
  l1->AddEntry(GEM_el_hits, "e","p");
  l1->AddEntry(GEM_mu_hits, "#mu","p");
  l1->AddEntry(GEM_g_hits,  "#gamma","p");
  l1->AddEntry(GEM_pi_hits, "#pi^#pm","p");
  // l1->AddEntry(GEM_ka_hits, "K^#pm","p");
  l1->AddEntry(GEM_p_hits,  "p","p");
  l1->AddEntry(GEM_N_hits,  "Nuclei","p");

  TLegend *l2 = new TLegend(l1_x1, l1_y1,l1_x2,l1_y2,NULL,"brNDC");
  l2->SetLineColor(0);    l2->SetLineStyle(0);  l2->SetLineWidth(0);
  l2->SetFillColor(4000); l2->SetBorderSize(0); l2->SetNColumns(1);
  l2->AddEntry(GEM_el_deps, "e","l");
  l2->AddEntry(GEM_mu_deps, "#mu","l");
  l2->AddEntry(GEM_g_deps,  "#gamma","l");
  l2->AddEntry(GEM_ha_deps, "had","l");


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

  Canvas_GEM_hits  = new TCanvas("Canvas_GEM_hits",  "Simhit time vs E_{kin} :: GEM",   600, 600);
  Canvas_ME0_hits  = new TCanvas("Canvas_ME0_hits",  "Simhit time vs E_{kin} :: ME0",   600, 600);

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
  Canvas_GEM_hits->Print(pdfFileName.c_str());

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
  Canvas_ME0_hits->Print(pdfFileName.c_str());

  // To Be Deleted
  // Canvas_GEM_hits_fancy  = new TCanvas("Canvas_GEM_hits_fancy",  "Simhit time vs E_{kin} :: GEM",   600, 600);
  // Canvas_ME0_hits_fancy  = new TCanvas("Canvas_ME0_hits_fancy",  "Simhit time vs E_{kin} :: ME0",   600, 600);

  TDir_Muon_hits_2Dplots->cd();
  // --------------------------- 
  GEM_el_hits->Write(); 
  GEM_mu_hits->Write(); 
  GEM_pi_hits->Write(); 
  GEM_ka_hits->Write(); 
  GEM_p_hits->Write(); 
  GEM_n_hits->Write(); 
  GEM_g_hits->Write(); 
  GEM_N_hits->Write(); 
  GEM_OH_hits->Write(); 

  ME0_el_hits->Write();   ME0_mu_hits->Write();   ME0_pi_hits->Write();   ME0_ka_hits->Write();   ME0_p_hits->Write(); 
  ME0_n_hits->Write();   ME0_g_hits->Write();   ME0_N_hits->Write();   ME0_OH_hits->Write();   // --------------------------- 
  outputfile->cd();

  Canvas_GEM_deposits  = new TCanvas("Canvas_GEM_deposits",  "Simhit time vs E_{deposit} :: GEM",  600, 600);
  Canvas_ME0_deposits  = new TCanvas("Canvas_ME0_deposits",  "Simhit time vs E_{deposit} :: ME0",  600, 600);

  // SimHit time vs Energy Deposit
  // Combine histograms in single plot
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
  Canvas_GEM_deposits->Print(pdfFileName.c_str());

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
  Canvas_ME0_deposits->Print(pdfFileName.c_str());


  Canvas_GEM_ekindep  = new TCanvas("Canvas_GEM_ekindep",  "E_{deposit} vs E_{kin}:: GEM",  600, 600);
  Canvas_ME0_ekindep  = new TCanvas("Canvas_ME0_ekindep",  "E_{deposit} vs E_{kin}:: ME0",  600, 600);

  // SimHit time vs Energy Deposit
  // Combine histograms in single plot
  Canvas_GEM_ekindep->cd();
  GEM_el_ekindep->GetYaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  GEM_el_ekindep->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  GEM_el_ekindep->GetXaxis()->SetTitleOffset(1.2);
  GEM_el_ekindep->SetTitle("SimHit time vs E_{deposit} :: GEM");
  GEM_el_ekindep->SetMarkerStyle(7);  GEM_el_ekindep->SetMarkerColor(kBlack);   GEM_el_ekindep->SetMarkerSize(1);  GEM_el_ekindep->Draw("P");
  GEM_mu_ekindep->SetMarkerStyle(24); GEM_mu_ekindep->SetMarkerColor(kBlue);    GEM_mu_ekindep->SetMarkerSize(1);  GEM_mu_ekindep->Draw("PSame");
  GEM_pi_ekindep->SetMarkerStyle(33); GEM_pi_ekindep->SetMarkerColor(kGreen);   GEM_pi_ekindep->SetMarkerSize(1);  GEM_pi_ekindep->Draw("PSame");
  GEM_ka_ekindep->SetMarkerStyle(5);  GEM_ka_ekindep->SetMarkerColor(kOrange);  GEM_ka_ekindep->SetMarkerSize(1);  GEM_ka_ekindep->Draw("PSame");
  GEM_p_ekindep->SetMarkerStyle(26);  GEM_p_ekindep->SetMarkerColor(kMagenta);  GEM_p_ekindep->SetMarkerSize(1);   GEM_p_ekindep->Draw("PSame");
  GEM_n_ekindep->SetMarkerStyle(32);  GEM_n_ekindep->SetMarkerColor(kViolet);   GEM_n_ekindep->SetMarkerSize(1);   GEM_n_ekindep->Draw("PSame");
  GEM_g_ekindep->SetMarkerStyle(30);  GEM_g_ekindep->SetMarkerColor(kCyan);     GEM_g_ekindep->SetMarkerSize(1);   GEM_g_ekindep->Draw("PSame");
  GEM_N_ekindep->SetMarkerStyle(2);   GEM_N_ekindep->SetMarkerColor(kRed);      GEM_N_ekindep->SetMarkerSize(1);   GEM_N_ekindep->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"GEM");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  Canvas_GEM_ekindep->SetTicks(1,1);
  Canvas_GEM_ekindep->Write();
  Canvas_GEM_ekindep->Print(pdfFileName.c_str());

  Canvas_ME0_ekindep->cd();
  ME0_el_ekindep->GetYaxis()->SetTitle("^{10}log E_{deposit} (keV)");
  ME0_el_ekindep->GetXaxis()->SetTitle("^{10}log E_{kin} (MeV)");
  ME0_el_ekindep->GetXaxis()->SetTitleOffset(1.2);
  ME0_el_ekindep->SetTitle("SimHit time vs E_{deposit} :: ME0");
  ME0_el_ekindep->SetMarkerStyle(7);  ME0_el_ekindep->SetMarkerColor(kBlack);   ME0_el_ekindep->SetMarkerSize(1);  ME0_el_ekindep->Draw("P");
  ME0_mu_ekindep->SetMarkerStyle(24); ME0_mu_ekindep->SetMarkerColor(kBlue);    ME0_mu_ekindep->SetMarkerSize(1);  ME0_mu_ekindep->Draw("PSame");
  ME0_pi_ekindep->SetMarkerStyle(33); ME0_pi_ekindep->SetMarkerColor(kGreen);   ME0_pi_ekindep->SetMarkerSize(1);  ME0_pi_ekindep->Draw("PSame");
  ME0_ka_ekindep->SetMarkerStyle(5);  ME0_ka_ekindep->SetMarkerColor(kOrange);  ME0_ka_ekindep->SetMarkerSize(1);  ME0_ka_ekindep->Draw("PSame");
  ME0_p_ekindep->SetMarkerStyle(26);  ME0_p_ekindep->SetMarkerColor(kMagenta);  ME0_p_ekindep->SetMarkerSize(1);   ME0_p_ekindep->Draw("PSame");
  ME0_n_ekindep->SetMarkerStyle(32);  ME0_n_ekindep->SetMarkerColor(kViolet);   ME0_n_ekindep->SetMarkerSize(1);   ME0_n_ekindep->Draw("PSame");
  ME0_g_ekindep->SetMarkerStyle(30);  ME0_g_ekindep->SetMarkerColor(kCyan);     ME0_g_ekindep->SetMarkerSize(1);   ME0_g_ekindep->Draw("PSame");
  ME0_N_ekindep->SetMarkerStyle(2);   ME0_N_ekindep->SetMarkerColor(kRed);      ME0_N_ekindep->SetMarkerSize(1);   ME0_N_ekindep->Draw("PSame");
  line_250ns_deps->Draw("][same"); line_max_deps->Draw("][same");
  l1->Draw();
  latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  latex_cmslab.DrawLatex(0.125, 0.850,"ME0");
  latex_legend.DrawLatex(0.850, 0.850,"particle type:");
  Canvas_ME0_ekindep->SetTicks(1,1);
  Canvas_ME0_ekindep->Write();
  Canvas_ME0_ekindep->Print(pdfFileName.c_str());

  TDir_Muon_hits_2Dplots->cd();
  // --------------------------- 
  GEM_el_deposits->Write();   GEM_mu_deposits->Write();   GEM_pi_deposits->Write();   GEM_ka_deposits->Write();   GEM_p_deposits->Write();   
  GEM_n_deposits->Write();    GEM_g_deposits->Write();    GEM_N_deposits->Write();    GEM_OH_deposits->Write(); 

  ME0_el_deposits->Write();   ME0_mu_deposits->Write();  ME0_pi_deposits->Write();  ME0_ka_deposits->Write();   ME0_p_deposits->Write(); 
  ME0_n_deposits->Write();    ME0_g_deposits->Write();   ME0_N_deposits->Write();   ME0_OH_deposits->Write(); 
  GEM_el_ekindep->Write();    GEM_mu_ekindep->Write();   GEM_pi_ekindep->Write();   GEM_ka_ekindep->Write();   GEM_p_ekindep->Write(); 
  GEM_n_ekindep->Write();     GEM_g_ekindep->Write();    GEM_N_ekindep->Write();    GEM_OH_ekindep->Write(); 

  ME0_el_ekindep->Write();   ME0_mu_ekindep->Write();   ME0_pi_ekindep->Write();   ME0_ka_ekindep->Write();   ME0_p_ekindep->Write(); 
  ME0_n_ekindep->Write();    ME0_g_ekindep->Write();    ME0_N_ekindep->Write();    ME0_OH_ekindep->Write(); 
  // --------------------------- 
  outputfile->cd();

  Canvas_GEM_1D_deps  = new TCanvas("Canvas_GEM_1D_deps",  "E_{deposit} :: GEM",  600, 600);
  Canvas_ME0_1D_deps  = new TCanvas("Canvas_ME0_1D_deps",  "E_{deposit} :: ME0",  600, 600);

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
  Canvas_GEM_1D_deps->Print(pdfFileName.c_str());

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
  Canvas_ME0_1D_deps->Print(pdfFileName.c_str());

  Canvas_GEM_1D_tof  = new TCanvas("Canvas_GEM_1D_tof",  "Time Of Flight :: GEM",  600, 600);
  Canvas_ME0_1D_tof  = new TCanvas("Canvas_ME0_1D_tof",  "Time Of Flight :: ME0",  600, 600);

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
  Canvas_GEM_1D_tof->Print(pdfFileName.c_str());

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
  Canvas_ME0_1D_tof->Print(pdfFileName.c_str());

  outputfile->cd();

  TDir_Muon_XY_RZ_views->cd();
  // --------------------------- 
  GEM_XY->Write();
  GEM_RZ->Write();
  ME0_XY->Write();
  ME0_RZ->Write();
  GEM_000ns_XY->Write();
  GEM_000ns_RZ->Write();
  ME0_000ns_XY->Write();
  ME0_000ns_RZ->Write();
  GEM_250ns_XY->Write();
  GEM_250ns_RZ->Write();
  ME0_250ns_XY->Write();
  ME0_250ns_RZ->Write();
  // --------------------------- 
  outputfile->cd();

  TDir_Muon_hits_deposits->cd();
  // --------------------------- 
  GEM_hits_tof->Write();
  GEM_hits_eta->Write();
  GEM_hits_phi->Write();
  GEM_hits_lin->Write();

  ME0_hits_tof->Write();
  ME0_hits_eta->Write();
  ME0_hits_phi->Write();
  ME0_hits_lin->Write();

  GEM_el_deps->Write();
  GEM_mu_deps->Write();
  GEM_ha_deps->Write();
  GEM_pi_deps->Write();
  GEM_ka_deps->Write();
  GEM_p_deps->Write();
  GEM_n_deps->Write();
  GEM_g_deps->Write();
  GEM_N_deps->Write();
  GEM_OH_deps->Write();
  GEM_All_deps->Write();
  GEM_HIP_deps->Write();

  ME0_el_deps->Write();
  ME0_mu_deps->Write();
  ME0_ha_deps->Write();
  ME0_pi_deps->Write();
  ME0_ka_deps->Write();
  ME0_p_deps->Write();
  ME0_n_deps->Write();
  ME0_g_deps->Write();
  ME0_N_deps->Write();
  ME0_OH_deps->Write();
  ME0_All_deps->Write();
  ME0_HIP_deps->Write();

  G11_el_deps->Write();
  G11_mu_deps->Write();
  G11_ha_deps->Write();
  G11_pi_deps->Write();
  G11_ka_deps->Write();
  G11_p_deps->Write();
  G11_n_deps->Write();
  G11_g_deps->Write();
  G11_N_deps->Write();
  G11_OH_deps->Write();
  G11_All_deps->Write();
  G11_HIP_deps->Write();

  G21_el_deps->Write();
  G21_mu_deps->Write();
  G21_ha_deps->Write();
  G21_pi_deps->Write();
  G21_ka_deps->Write();
  G21_p_deps->Write();
  G21_n_deps->Write();
  G21_g_deps->Write();
  G21_N_deps->Write();
  G21_OH_deps->Write();
  G21_All_deps->Write();
  G21_HIP_deps->Write();

  GEM_el_deps_000->Write();
  GEM_mu_deps_000->Write();
  GEM_ha_deps_000->Write();
  GEM_pi_deps_000->Write();
  GEM_ka_deps_000->Write();
  GEM_p_deps_000->Write();
  GEM_n_deps_000->Write();
  GEM_g_deps_000->Write();
  GEM_N_deps_000->Write();
  GEM_OH_deps_000->Write();

  ME0_el_deps_000->Write();
  ME0_mu_deps_000->Write();
  ME0_ha_deps_000->Write();
  ME0_pi_deps_000->Write();
  ME0_ka_deps_000->Write();
  ME0_p_deps_000->Write();
  ME0_n_deps_000->Write();
  ME0_g_deps_000->Write();
  ME0_N_deps_000->Write();
  ME0_OH_deps_000->Write();

  GEM_el_deps_250->Write();
  GEM_mu_deps_250->Write();
  GEM_ha_deps_250->Write();
  GEM_pi_deps_250->Write();
  GEM_ka_deps_250->Write();
  GEM_p_deps_250->Write();
  GEM_n_deps_250->Write();
  GEM_g_deps_250->Write();
  GEM_N_deps_250->Write();
  GEM_OH_deps_250->Write();

  ME0_el_deps_250->Write();
  ME0_mu_deps_250->Write();
  ME0_ha_deps_250->Write();
  ME0_pi_deps_250->Write();
  ME0_ka_deps_250->Write();
  ME0_p_deps_250->Write();
  ME0_n_deps_250->Write();
  ME0_g_deps_250->Write();
  ME0_N_deps_250->Write();
  ME0_OH_deps_250->Write();

  GEM_el_lindep->Write();
  GEM_mu_lindep->Write();
  GEM_ha_lindep->Write();
  GEM_pi_lindep->Write();
  GEM_ka_lindep->Write();
  GEM_p_lindep->Write();
  GEM_n_lindep->Write();
  GEM_g_lindep->Write();
  GEM_N_lindep->Write();
  GEM_OH_lindep->Write();
  GEM_All_lindep->Write();
  GEM_HIP_lindep->Write();

  ME0_el_lindep->Write();
  ME0_mu_lindep->Write();
  ME0_ha_lindep->Write();
  ME0_pi_lindep->Write();
  ME0_ka_lindep->Write();
  ME0_p_lindep->Write();
  ME0_n_lindep->Write();
  ME0_g_lindep->Write();
  ME0_N_lindep->Write();
  ME0_OH_lindep->Write();
  ME0_All_lindep->Write();
  ME0_HIP_lindep->Write();

  G11_el_lindep->Write();
  G11_mu_lindep->Write();
  G11_ha_lindep->Write();
  G11_pi_lindep->Write();
  G11_ka_lindep->Write();
  G11_p_lindep->Write();
  G11_n_lindep->Write();
  G11_g_lindep->Write();
  G11_N_lindep->Write();
  G11_OH_lindep->Write();
  G11_All_lindep->Write();
  G11_HIP_lindep->Write();

  G21_el_lindep->Write();
  G21_mu_lindep->Write();
  G21_ha_lindep->Write();
  G21_pi_lindep->Write();
  G21_ka_lindep->Write();
  G21_p_lindep->Write();
  G21_n_lindep->Write();
  G21_g_lindep->Write();
  G21_N_lindep->Write();
  G21_OH_lindep->Write();
  G21_All_lindep->Write();
  G21_HIP_lindep->Write();

  ME0_All_lindep_roll->Write();
  G11_All_lindep_roll->Write();
  G21_All_lindep_roll->Write();

  ME0_All_lindep_eta01->Write();
  ME0_All_lindep_eta02->Write();
  ME0_All_lindep_eta03->Write();
  ME0_All_lindep_eta04->Write();
  ME0_All_lindep_eta05->Write();
  ME0_All_lindep_eta06->Write();
  ME0_All_lindep_eta07->Write();
  ME0_All_lindep_eta08->Write();

  G11_All_lindep_eta01->Write();
  G11_All_lindep_eta02->Write();
  G11_All_lindep_eta03->Write();
  G11_All_lindep_eta04->Write();
  G11_All_lindep_eta05->Write();
  G11_All_lindep_eta06->Write();
  G11_All_lindep_eta07->Write();
  ME0_All_lindep_eta08->Write();

  G21_All_lindep_eta01->Write();
  G21_All_lindep_eta02->Write();
  G21_All_lindep_eta03->Write();
  G21_All_lindep_eta04->Write();
  G21_All_lindep_eta05->Write();
  G21_All_lindep_eta06->Write();
  G21_All_lindep_eta07->Write();
  G21_All_lindep_eta08->Write();
  // ---------------------------
  outputfile->cd();

  TDir_Muon_hits_kinetics->cd();
  // ---------------------------
  GEM_el_kins->Write();
  GEM_mu_kins->Write();
  GEM_ha_kins->Write();
  GEM_pi_kins->Write();
  GEM_ka_kins->Write();
  GEM_p_kins->Write();
  GEM_n_kins->Write();
  GEM_g_kins->Write();
  GEM_N_kins->Write();
  GEM_OH_kins->Write();
  GEM_All_kins->Write();
  GEM_HIP_kins->Write();

  ME0_el_kins->Write();
  ME0_mu_kins->Write();
  ME0_ha_kins->Write();
  ME0_pi_kins->Write();
  ME0_ka_kins->Write();
  ME0_p_kins->Write();
  ME0_n_kins->Write();
  ME0_g_kins->Write();
  ME0_N_kins->Write();
  ME0_OH_kins->Write();
  ME0_All_kins->Write();
  ME0_HIP_kins->Write();

  GEM_el_kins->Write();
  GEM_mu_kins->Write();
  GEM_ha_kins->Write();
  GEM_pi_kins->Write();
  GEM_ka_kins->Write();
  GEM_p_kins->Write();
  GEM_n_kins->Write();
  GEM_g_kins->Write();
  GEM_N_kins->Write();
  GEM_OH_kins->Write();
  GEM_All_kins->Write();
  GEM_HIP_kins->Write();

  G11_el_kins->Write();
  G11_mu_kins->Write();
  G11_ha_kins->Write();
  G11_pi_kins->Write();
  G11_ka_kins->Write();
  G11_p_kins->Write();
  G11_n_kins->Write();
  G11_g_kins->Write();
  G11_N_kins->Write();
  G11_OH_kins->Write();
  G11_All_kins->Write();
  G11_HIP_kins->Write();

  G21_el_linkin->Write();
  G21_mu_linkin->Write();
  G21_ha_linkin->Write();
  G21_pi_linkin->Write();
  G21_ka_linkin->Write();
  G21_p_linkin->Write();
  G21_n_linkin->Write();
  G21_g_linkin->Write();
  G21_N_linkin->Write();
  G21_OH_linkin->Write();
  G21_All_linkin->Write();
  G21_HIP_linkin->Write();

  ME0_el_linkin->Write();
  ME0_mu_linkin->Write();
  ME0_ha_linkin->Write();
  ME0_pi_linkin->Write();
  ME0_ka_linkin->Write();
  ME0_p_linkin->Write();
  ME0_n_linkin->Write();
  ME0_g_linkin->Write();
  ME0_N_linkin->Write();
  ME0_OH_linkin->Write();
  ME0_All_linkin->Write();
  ME0_HIP_linkin->Write();

  G11_el_linkin->Write();
  G11_mu_linkin->Write();
  G11_ha_linkin->Write();
  G11_pi_linkin->Write();
  G11_ka_linkin->Write();
  G11_p_linkin->Write();
  G11_n_linkin->Write();
  G11_g_linkin->Write();
  G11_N_linkin->Write();
  G11_OH_linkin->Write();
  G11_All_linkin->Write();
  G11_HIP_linkin->Write();

  G21_el_linkin->Write();
  G21_mu_linkin->Write();
  G21_ha_linkin->Write();
  G21_pi_linkin->Write();
  G21_ka_linkin->Write();
  G21_p_linkin->Write();
  G21_n_linkin->Write();
  G21_g_linkin->Write();
  G21_N_linkin->Write();
  G21_OH_linkin->Write();
  G21_All_linkin->Write();
  G21_HIP_linkin->Write();
  // ---------------------------
  outputfile->cd();


  TDir_Muon_hits_timing->cd();
  // ---------------------------
  GEM_el_tof->Write();
  GEM_mu_tof->Write();
  GEM_ha_tof->Write();
  GEM_pi_tof->Write();
  GEM_ka_tof->Write();
  GEM_p_tof->Write();
  GEM_n_tof->Write();
  GEM_g_tof->Write();
  GEM_N_tof->Write();
  GEM_OH_tof->Write();

  ME0_el_tof->Write();
  ME0_mu_tof->Write();
  ME0_ha_tof->Write();
  ME0_pi_tof->Write();
  ME0_ka_tof->Write();
  ME0_p_tof->Write();
  ME0_n_tof->Write();
  ME0_g_tof->Write();
  ME0_N_tof->Write();
  ME0_OH_tof->Write();

  GEM_el_time->Write();
  GEM_mu_time->Write();
  GEM_ha_time->Write();
  GEM_pi_time->Write();
  GEM_ka_time->Write();
  GEM_p_time->Write();
  GEM_n_time->Write();
  GEM_g_time->Write();
  GEM_N_time->Write();
  GEM_OH_time->Write();

  ME0_el_time->Write();
  ME0_mu_time->Write();
  ME0_ha_time->Write();
  ME0_pi_time->Write();
  ME0_ka_time->Write();
  ME0_p_time->Write();
  ME0_n_time->Write();
  ME0_g_time->Write();
  ME0_N_time->Write();
  ME0_OH_time->Write();
  // --------------------------- 
  outputfile->cd();

  TDir_Muon_Nuclei->cd();
  ME0_Nuclei_A_Z->Write();
  ME0_Nuclei_A->Write();
  ME0_Nuclei_Z->Write();
  for(int i=0; i<19; ++i) {ME0_Nuclei_List->GetXaxis()->SetBinLabel(i+1, (nuclei[i]).c_str());}
  ME0_Nuclei_List->Write();
  for(int i=0; i<9; ++i) {ME0_HIP_id->GetXaxis()->SetBinLabel(i+1, (particles[i]).c_str()); GEM_HIP_id->GetXaxis()->SetBinLabel(i+1, (particles[i]).c_str());}
  ME0_HIP_id->Write();
  GEM_HIP_id->Write();
  outputfile->cd();

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
  // int pu = 24; // normalized to 5 * 10^33 cm^-2 s^-1

  // New Calculation
  // ---------------
  // average amount of hits per event = (hits / entries)
  // ( average amount of hits per event ) / area
  // amount of collisions in 1s = 11245 (revolution frequency) * 1380 (colliding bunches) * pu (no of collisions per crossing)
  // => Rate [Hz/cm^2] = ( average amount of hits per event ) / area * amount of collisions in 1s
  // => Rate [Hz/cm^2] = N /entries / area * 11245 * 1380 * 21

  // getLumi(pu, space, com)
  // for(int i=0; i<50; ++i) {
  // std::cout<<"Inst Lumi for "<<i<<" PU interactions at "<<bunchspacing<<"ns bunch spacing and at "<<comenergy<<"TeV is "<<getLumi(i,bunchspacing,comenergy)<<" * 10^34 cm^-2 s^-1 "<<std::endl;
  // InstLumi[i] = getLumi(i,bunchspacing,comenergy);
  // }
  // getPU(lumi, space, com)
  // for(int i=0; i<60; ++i) {
  // std::cout<<"PU for Inst Lumi "<<i*1.0/100<<" * 10^34 cm^-2 s^-1 at "<<bunchspacing<<"ns bunch spacing and at "<<comenergy<<"TeV is "<<getPU(i*1.0/100,bunchspacing,comenergy)<<" PU"<<std::endl;
  // }

  // double bunches;
  // if(bunchspacing==50) bunches = 1380;
  // if(bunchspacing==25) bunches = 2808;
  // double corr_fact = 1.0/entries * 11245 * bunches * pu;

  // Rates as function of pu / instantaneous luminosity
  // const int max_pu = 33;
  // for(int i=0; i<max_pu; ++i) {
  //  RPC_rates_Summary[0][i] =  (rpc_barrel_hits+rpc_endcap_hits) * 1.0/(rpc_barrel_area+rpc_endcap_area) * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
  //  RPC_rates_Summary[1][i] =  (rpc_barrel_hits)                 * 1.0/(rpc_barrel_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
  //  RPC_rates_Summary[2][i] =  (rpc_endcap_hits)                 * 1.0/(rpc_endcap_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
  //  RPC_uncer_Rate[0][i] =  sqrt(rpc_barrel_hits+rpc_endcap_hits) * 1.0/(rpc_barrel_area+rpc_endcap_area) * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
  //  RPC_uncer_Rate[1][i] =  sqrt(rpc_barrel_hits)                 * 1.0/(rpc_barrel_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
  //  RPC_uncer_Rate[2][i] =  sqrt(rpc_endcap_hits)                 * 1.0/(rpc_endcap_area)                 * 1.0/entries * 11245 * bunches * i * sg_corr_fact;
  //  RPC_uncer_Lumi[0][i] =  0;
  //  RPC_uncer_Lumi[1][i] =  0;
  //  RPC_uncer_Lumi[2][i] =  0;
  //}

  // gr_RPC_Rates_All    = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[0], RPC_uncer_Lumi[0], RPC_uncer_Rate[0]); 
  // gr_RPC_Rates_Barrel = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[1], RPC_uncer_Lumi[1], RPC_uncer_Rate[1]); 
  // gr_RPC_Rates_Endcap = new TGraphErrors(max_pu, InstLumi, RPC_rates_Summary[2], RPC_uncer_Lumi[2], RPC_uncer_Rate[2]); 

  // gr_RPC_Rates_All->SetMarkerStyle(21);    gr_RPC_Rates_All->SetMarkerColor(kCyan);     gr_RPC_Rates_All->SetLineColor(kCyan);
  // gr_RPC_Rates_Barrel->SetMarkerStyle(21); gr_RPC_Rates_Barrel->SetMarkerColor(kRed);   gr_RPC_Rates_Barrel->SetLineColor(kRed);
  // gr_RPC_Rates_Endcap->SetMarkerStyle(21); gr_RPC_Rates_Endcap->SetMarkerColor(kBlack); gr_RPC_Rates_Endcap->SetLineColor(kBlack);

  // l3_x1 = 0.20; l3_x2 = 0.45; l3_y1 = 0.65; l3_y2 = 0.85;
  // l4_x1 = 0.20; l4_x2 = 0.35; l4_y1 = 0.75; l4_y2 = 0.85;
  // TLegend *l3 = new TLegend(l3_x1,l3_y1,l3_x2,l3_y2,NULL,"brNDC");
  // l3->SetLineColor(0);    l3->SetLineStyle(0);  l3->SetLineWidth(0);
  // l3->SetFillColor(4000); l3->SetBorderSize(0); l3->SetNColumns(1);
  // l3->AddEntry(gr_RPC_Rates_Endcap, "Endcap","pl");
  // l3->AddEntry(gr_RPC_Rates_All,    "Barrel + Endcap","pl");
  // l3->AddEntry(gr_RPC_Rates_Barrel, "Barrel","pl");

  // Canvas_RPC_Rates = new TCanvas("Canvas_RPC_Rates", "Rates in RPC System", 600, 600);
  // gr_RPC_Rates_Endcap->GetXaxis()->SetTitle("Instantaneous Luminosity #times 10^{34} (cm^{-2}s^{-1})");
  // gr_RPC_Rates_Endcap->GetYaxis()->SetTitle("Rate (Hz/cm^{2})");
  // gr_RPC_Rates_Endcap->GetYaxis()->SetTitleOffset(1.30);
  // gr_RPC_Rates_Endcap->GetXaxis()->SetRangeUser(0.00,0.68);
  // gr_RPC_Rates_Endcap->GetYaxis()->SetRangeUser(0.00,6.00);
  // gr_RPC_Rates_Endcap->SetTitle("Rates in RPC System");
  // gr_RPC_Rates_Endcap->Draw("PA");
  // gr_RPC_Rates_Barrel->Draw("Psame");
  // gr_RPC_Rates_All->Draw("Psame");
  // l3->Draw("same");
  // latex_cmslab.DrawLatex(0.10, 0.925,comlabel.c_str());
  // Canvas_RPC_Rates->SetTicks(1,1);
  // Canvas_RPC_Rates->Write(); 
  // Canvas_RPC_Rates->Print(pdfFileName.c_str());

  // TDir_Muon_rates->cd();
  // ------------------- 
  // RPC_hits->Write();
  // RPC_area->Write();
  // RPC_rates->Write();

  // last plot for PDF File :: print empty Dummy
  pdfFileName = pdfFileNameBase + ".pdf]";
  Dummy->Print(pdfFileName.c_str());

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyME0SimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(tech_debug) std::cout<<"[MyME0SimHitAnalyzer :: Analyze]"<<std::endl;
  
  EventCounter->Fill(1); // count Events
  bool ME0hit = 0, GE11hit = 0, GE21hit = 0, ME0hip = 0, GE11hip = 0, GE21hip = 0;

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
  
  // ===============================
  //      Loop over the SimHits
  // ===============================
  
  for (std::vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); ++iHit) {
    DetId theDetUnitId((*iHit).detUnitId());
    DetId simdetid= DetId((*iHit).detUnitId());
    
    
    int pid            = (*iHit).particleType();
    int process        = (*iHit).processType();
    double time        = (*iHit).timeOfFlight();
    double energy      = (*iHit).momentumAtEntry().perp()*1000;     // MeV
    double deposit     = (*iHit).energyLoss()*1000000;              // keV
    double log_t       = log10((*iHit).timeOfFlight());
    double log_e       = log10((*iHit).momentumAtEntry().perp()*1000); // MeV
    double log_d       = log10((*iHit).energyLoss()*1000000);          // keV
    
    if(edep_30eV && deposit < edepMin) continue; // skip hits with Energy Deposit < 30.0 eV

    // Reprocess the process type
    switch(process) {
    case  21: // 21 -> 17
    case  22: // 22 -> 18
    case  23: // 23 -> 19
    case  24: process -=4;  break; //  24 -> 20
    case 111: process = 21; break; // 111 -> 21
    case 121: process = 22; break;
    case 131: process = 23; break;
    case 141: process = 24; break;
    case 151: process = 25; break;    //   1 Coulomb Scattering  2 Ionisation            3 Bremsstrahlung     4 PairProductionByCharged  5 Annihilation
    case 161: process = 26; break;    //   6 AnnihilationToMuMu  7 AnnihilationToHad     8 NuclearStopping    9 ElectronSuper           10 MultipleScattering
    case 201: process = 27; break;    //  11 Rayleigh Scatter   12 PhotoElectricEffect  13 ComptonScattering 14 GammaConv               15 GammaConvToMuMu      16 GammaSuper
    case 202: process = 28; break;    //  21 Cherenkov          22 Scintillation        23 Synchrotron Rad   24 Transition Radiation                                         
    case 211: process = 29; break;    // 111 Hadron Elastic    121 Hadron Inelastic    131 Hadron Capture   141 Hadron Fission         151 Hadron At REst       161 HadronCEX
    case 231: process = 30; break;    // 201 Decay             201 Decay with Spin     211 Decay Unknown    231 Decay External    
    }

    // Reprocess the particle id (for histogramming purposes)
    // species e - mu - pi - ka - p - n - g - N - OH
    // binno   1   2    3    4    5   6   7   8   9

    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::GEM){ // Only GEMs
      // GEM Geometry
      // ============
      GEMDetId rollId(theDetUnitId);
      const GEMEtaPartition* rollasociated = gemGeom->etaPartition(rollId);
      const BoundPlane & GEMSurface = rollasociated->surface(); 
      GlobalPoint GEMGlobalPoint = GEMSurface.toGlobal((*iHit).localPosition());
      // GlobalPoint GEMGlobalEntry = GEMSurface.toGlobal((*iHit).entryPoint());
      // GlobalPoint GEMGlobalExit  = GEMSurface.toGlobal((*iHit).exitPoint());
      
      // if(gem_only_ge11 && rollId.station()==2) continue; // exit loop if SimHit is in GE2/1 and the GE11 only mode was selected
      int station = rollId.station(); // 1 = GE11 | 2 = GE21
      int roll    = rollId.roll();    // 1-8      

      bool isGE11L1 = 0, isGE11L2 =0, isGE11Odd = 0, isGE11Even = 0;

      if(station != 0) {  
	if(phys_debug) {
	  std::cout<<"GEM SimHit in "<<std::setw(12)<<(int)rollId<<std::setw(24)<<rollId;
	  std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<GEMGlobalPoint.z();
	  std::cout<<" | r = "<<std::setw(12)<<GEMGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<GEMGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<GEMGlobalPoint.eta();
	  std::cout<<" | global position = "<<GEMGlobalPoint<<std::endl;
	}
	
	if(station==1) { GE11hit = 1; if(deposit > hipMin) { GE11hip = 1; } }
	if(station==2) { GE21hit = 1; if(deposit > hipMin) { GE21hip = 1; } }
	
	if(station==1) {
	  if(rollId.chamber()%2==1) isGE11Odd = 1;
	  if(rollId.chamber()%2==0) isGE11Even = 1;
	  if(rollId.layer()==1) isGE11L1 = 1;
	  if(rollId.layer()==2) isGE11L2 = 1;
	}

	double GEM_GlobalPoint_R = sqrt(pow(GEMGlobalPoint.x(),2)+pow(GEMGlobalPoint.y(),2));
	double gemr = GEM_GlobalPoint_R;
	GEM_XY->Fill(GEMGlobalPoint.x(), GEMGlobalPoint.y());             
	GEM_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));  
	
	bool InTime  = (*iHit).timeOfFlight()<250;
	bool OutTime = (*iHit).timeOfFlight()>250; 
	bool BXTime  = (*iHit).timeOfFlight()<25;
	bool OXTime  = (*iHit).timeOfFlight()>25; 
	
	if(InTime) {
	  GEM_000ns_XY->Fill(GEMGlobalPoint.x(), GEMGlobalPoint.y());
	  GEM_000ns_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));
	}
	if(OutTime) {
	  GEM_250ns_XY->Fill(GEMGlobalPoint.x(), GEMGlobalPoint.y());
	  GEM_250ns_RZ->Fill(fabs(GEMGlobalPoint.z()), fabs(GEM_GlobalPoint_R));
	}
	
	if(phys_debug) std::cout<<"SimHit from particle = "<<pid<<" created in process = "<<process<<" aka "<<proc[process]<<std::endl;
	
	// std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	// std::cout<<" 10 log (tof) = "<<log_t<<" [ns] and 10 log (E) = "<<log_d<<" [MeV]"<<std::endl;
	
	// All Particles :: Before any particle ID
	station==1?G11_All_hits_R->Fill(gemr):G21_All_hits_R->Fill(gemr);
	station==1?G11_All_hits_E->Fill(roll):G21_All_hits_E->Fill(roll);      
	GEM_All_deps->Fill(log_d); GEM_All_lindep->Fill(deposit);
	GEM_All_kins->Fill(log_e); GEM_All_linkin->Fill(energy);
	if(station==1) {G11_All_deps->Fill(log_d); G11_All_lindep->Fill(deposit); G11_All_kins->Fill(log_e); G11_All_linkin->Fill(energy);}
	else           {G21_All_deps->Fill(log_d); G21_All_lindep->Fill(deposit); G21_All_kins->Fill(log_e); G21_All_linkin->Fill(energy);}
	if(isGE11L1) G11_L1_All_hits_R->Fill(gemr); 
	if(isGE11L2) G11_L2_All_hits_R->Fill(gemr); 
	if(isGE11Odd)  G11_Od_All_hits_R->Fill(gemr); 
	if(isGE11Even) G11_Ev_All_hits_R->Fill(gemr); 
	if(deposit > hipMin) {
	  station==1?G11_HIP_hits_R->Fill(gemr):G21_HIP_hits_R->Fill(gemr);
	  station==1?G11_HIP_hits_E->Fill(roll):G21_HIP_hits_E->Fill(roll);
	  GEM_HIP_deps->Fill(log_d);  GEM_HIP_lindep->Fill(deposit);
	  GEM_HIP_kins->Fill(log_e); 	GEM_HIP_linkin->Fill(energy);
	  if(station==1) {G11_HIP_deps->Fill(log_d); G11_HIP_lindep->Fill(deposit); G11_HIP_kins->Fill(log_e); G11_HIP_linkin->Fill(energy);}
	  else           {G21_HIP_deps->Fill(log_d); G21_HIP_lindep->Fill(deposit); G21_HIP_kins->Fill(log_e); G21_HIP_linkin->Fill(energy);}
	  if(isGE11L1) G11_L1_HIP_hits_R->Fill(gemr); 
	  if(isGE11L2) G11_L2_HIP_hits_R->Fill(gemr); 
	  if(isGE11Odd)  G11_Od_HIP_hits_R->Fill(gemr); 
	  if(isGE11Even) G11_Ev_HIP_hits_R->Fill(gemr); 
	}
	if(InTime)  {station==1?G11_All_hits_000_R->Fill(gemr):G21_All_hits_000_R->Fill(gemr); station==1?G11_All_hits_000_E->Fill(roll):G21_All_hits_000_E->Fill(roll);}
	if(OutTime) {station==1?G11_All_hits_250_R->Fill(gemr):G21_All_hits_250_R->Fill(gemr); station==1?G11_All_hits_000_E->Fill(roll):G21_All_hits_250_E->Fill(roll);}
	if(BXTime)  {station==1?G11_All_hits_00_R->Fill(gemr):G21_All_hits_00_R->Fill(gemr);   station==1?G11_All_hits_00_E->Fill(roll):G21_All_hits_00_E->Fill(roll);}
	if(OXTime)  {station==1?G11_All_hits_25_R->Fill(gemr):G21_All_hits_25_R->Fill(gemr);   station==1?G11_All_hits_25_E->Fill(roll):G21_All_hits_25_E->Fill(roll);}
	
	if(station==1) {
	  G11_All_lindep_roll->Fill(deposit,roll);
	  switch(roll) {
	  case 1: G11_All_lindep_eta01->Fill(deposit); break;
	  case 2: G11_All_lindep_eta02->Fill(deposit); break;
	  case 3: G11_All_lindep_eta03->Fill(deposit); break;
	  case 4: G11_All_lindep_eta04->Fill(deposit); break;
	  case 5: G11_All_lindep_eta05->Fill(deposit); break;
	  case 6: G11_All_lindep_eta06->Fill(deposit); break;
	  case 7: G11_All_lindep_eta07->Fill(deposit); break;
	  case 8: G11_All_lindep_eta08->Fill(deposit); break;
	  }
	}
	if(station==2) {
	  G21_All_lindep_roll->Fill(deposit,roll);
	  switch(roll) {
	  case 1: G21_All_lindep_eta01->Fill(deposit); break;
	  case 2: G21_All_lindep_eta02->Fill(deposit); break;
	  case 3: G21_All_lindep_eta03->Fill(deposit); break;
	  case 4: G21_All_lindep_eta04->Fill(deposit); break;
	  case 5: G21_All_lindep_eta05->Fill(deposit); break;
	  case 6: G21_All_lindep_eta06->Fill(deposit); break;
	  case 7: G21_All_lindep_eta07->Fill(deposit); break;
	  case 8: G21_All_lindep_eta08->Fill(deposit); break;
	  }
	}
	
	
	if(tech_debug) std::cout<<"Distinction by Particle Type :: part 1"<<std::endl;
	
	// Distinction by Particle Type
	if(abs(pid)==22) { // PHOTON
	  GEM_g_hits->Fill(log_e,log_t); GEM_g_deposits->Fill(log_d,log_t); GEM_g_ekindep->Fill(log_e, log_d); GEM_g_kins->Fill(log_e);   GEM_g_process->Fill(process);
	  GEM_g_deps->Fill(log_d);       GEM_g_tof->Fill(log_t);            GEM_g_time->Fill(time);            GEM_g_lindep->Fill(deposit); GEM_g_linkin->Fill(energy);
	  station==1?G11_g_hits_R->Fill(gemr):G21_g_hits_R->Fill(gemr);
	  station==1?G11_g_hits_E->Fill(roll):G21_g_hits_E->Fill(roll);
	  if(InTime)  { GEM_g_deps_000->Fill(log_d); GEM_g_proc000->Fill(process); station==1?G11_g_hits_000_R->Fill(gemr):G21_g_hits_000_R->Fill(gemr);}
	  if(OutTime) { GEM_g_deps_250->Fill(log_d); GEM_g_proc250->Fill(process); station==1?G11_g_hits_250_R->Fill(gemr):G21_g_hits_250_R->Fill(gemr);}
	  if(deposit > hipMin) { GEM_HIP_id->Fill(7); }
	  if(station==1) {G11_g_deps->Fill(log_d); G11_g_lindep->Fill(deposit); G11_g_kins->Fill(log_e); G11_g_linkin->Fill(energy);}
	  else           {G21_g_deps->Fill(log_d); G21_g_lindep->Fill(deposit); G21_g_kins->Fill(log_e); G21_g_linkin->Fill(energy);}
	  if(isGE11L1)   { G11_L1_g_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_g_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_g_hits_E->Fill(roll); } }
	  if(isGE11L2)   { G11_L2_g_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_g_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_g_hits_E->Fill(roll); } }
	  if(isGE11Odd)  G11_Od_g_hits_R->Fill(gemr); 
	  if(isGE11Even) G11_Ev_g_hits_R->Fill(gemr); 
	}
	else if(abs(pid)==2212) { // PROTON
	  GEM_p_hits->Fill(log_e,log_t); GEM_p_deposits->Fill(log_d,log_t); GEM_p_ekindep->Fill(log_e, log_d); GEM_p_kins->Fill(log_e);   GEM_p_process->Fill(process);
	  GEM_p_deps->Fill(log_d);       GEM_p_tof->Fill(log_t);            GEM_p_time->Fill(time);            GEM_p_lindep->Fill(deposit); GEM_p_linkin->Fill(energy);
	  station==1?G11_p_hits_R->Fill(gemr):G21_p_hits_R->Fill(gemr);
	  station==1?G11_p_hits_E->Fill(roll):G21_p_hits_E->Fill(roll);
	  if(InTime)  { GEM_p_deps_000->Fill(log_d); GEM_p_proc000->Fill(process); station==1?G11_p_hits_000_R->Fill(gemr):G21_p_hits_000_R->Fill(gemr);}
	  if(OutTime) { GEM_p_deps_250->Fill(log_d); GEM_p_proc250->Fill(process); station==1?G11_p_hits_250_R->Fill(gemr):G21_p_hits_250_R->Fill(gemr);}
	  if(deposit > hipMin) { GEM_HIP_id->Fill(5); }
	  if(station==1) {G11_p_deps->Fill(log_d); G11_p_lindep->Fill(deposit); G11_p_kins->Fill(log_e); G11_p_linkin->Fill(energy);}
	  else           {G21_p_deps->Fill(log_d); G21_p_lindep->Fill(deposit); G21_p_kins->Fill(log_e); G21_p_linkin->Fill(energy);}
	  if(isGE11L1)   { G11_L1_p_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_p_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_p_hits_E->Fill(roll); } }
	  if(isGE11L2)   { G11_L2_p_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_p_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_p_hits_E->Fill(roll); } }
	  if(isGE11Odd)  G11_Od_p_hits_R->Fill(gemr); 
	  if(isGE11Even) G11_Ev_p_hits_R->Fill(gemr); 
	}
	else if(abs(pid)==2112) { // NEUTRON
	  GEM_n_hits->Fill(log_e,log_t); GEM_n_deposits->Fill(log_d,log_t); GEM_n_ekindep->Fill(log_e, log_d); GEM_n_kins->Fill(log_e);   GEM_n_process->Fill(process);
	  GEM_n_deps->Fill(log_d);       GEM_n_tof->Fill(log_t);            GEM_n_time->Fill(time);            GEM_n_lindep->Fill(deposit); GEM_n_linkin->Fill(energy);
	  station==1?G11_n_hits_R->Fill(gemr):G21_n_hits_R->Fill(gemr);
	  station==1?G11_n_hits_E->Fill(roll):G21_n_hits_E->Fill(roll);
	  if(InTime)  { GEM_n_deps_000->Fill(log_d); GEM_n_proc000->Fill(process); station==1?G11_n_hits_000_R->Fill(gemr):G21_n_hits_000_R->Fill(gemr);}
	  if(OutTime) { GEM_n_deps_250->Fill(log_d); GEM_n_proc250->Fill(process); station==1?G11_n_hits_000_R->Fill(gemr):G21_n_hits_000_R->Fill(gemr);}
	  if(deposit > hipMin) { GEM_HIP_id->Fill(6); }
	  if(station==1) {G11_n_deps->Fill(log_d); G11_n_lindep->Fill(deposit); G11_n_kins->Fill(log_e); G11_n_linkin->Fill(energy);}
	  else           {G21_n_deps->Fill(log_d); G21_n_lindep->Fill(deposit); G21_n_kins->Fill(log_e); G21_n_linkin->Fill(energy);}
	  if(isGE11L1)   { G11_L1_n_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_n_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_n_hits_E->Fill(roll); } }
	  if(isGE11L2)   { G11_L2_n_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_n_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_n_hits_E->Fill(roll); } }
	  if(isGE11Odd)  G11_Od_n_hits_R->Fill(gemr); 
	  if(isGE11Even) G11_Ev_n_hits_R->Fill(gemr); 
	}
	else if(abs(pid)>1E9)   { // NUCLEI
	  GEM_N_hits->Fill(log_e,log_t); GEM_N_deposits->Fill(log_d,log_t); GEM_N_ekindep->Fill(log_e, log_d); GEM_N_kins->Fill(log_e);   GEM_N_process->Fill(process);
	  GEM_N_deps->Fill(log_d);       GEM_N_tof->Fill(log_t);            GEM_N_time->Fill(time);            GEM_N_lindep->Fill(deposit); GEM_N_linkin->Fill(energy);
	  station==1?G11_N_hits_R->Fill(gemr):G21_N_hits_R->Fill(gemr);
	  station==1?G11_N_hits_E->Fill(roll):G21_N_hits_E->Fill(roll);
	  if(InTime)  { GEM_N_deps_000->Fill(log_d); GEM_N_proc000->Fill(process); station==1?G11_N_hits_000_R->Fill(gemr):G21_N_hits_000_R->Fill(gemr);}
	  if(OutTime) { GEM_N_deps_250->Fill(log_d); GEM_N_proc250->Fill(process); station==1?G11_N_hits_000_R->Fill(gemr):G21_N_hits_000_R->Fill(gemr);}
	  std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	  std::cout<<" 10 log (tof) = "<<log_t<<" [ns] and 10 log (E) = "<<log_d<<" [keV] catalogued as NUCLEI"<<std::endl;
	  if(deposit > hipMin) { GEM_HIP_id->Fill(8); }
	  if(station==1) {G11_N_deps->Fill(log_d); G11_N_lindep->Fill(deposit); G11_N_kins->Fill(log_e); G11_N_linkin->Fill(energy);}
	  else           {G21_N_deps->Fill(log_d); G21_N_lindep->Fill(deposit); G21_N_kins->Fill(log_e); G21_N_linkin->Fill(energy);}
	  if(isGE11L1)   { G11_L1_N_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_N_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_N_hits_E->Fill(roll); } }
	  if(isGE11L2)   { G11_L2_N_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_N_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_N_hits_E->Fill(roll); } }
	  if(isGE11Odd)  G11_Od_N_hits_R->Fill(gemr); 
	  if(isGE11Even) G11_Ev_N_hits_R->Fill(gemr); 
	}
	else {
	  
	  if(tech_debug) std::cout<<"Distinction by Particle Type :: part 2"<<std::endl;
	  
	  if(abs(pid)!=11 && abs(pid)!=13) {
	    if(tech_debug) std::cout<<"PID = "<<pid<<" --> Filling hadron histograms"<<std::endl;
	    GEM_ha_kins->Fill(log_e); GEM_ha_deps->Fill(log_d); GEM_ha_tof->Fill(log_t); GEM_ha_time->Fill(time); GEM_ha_lindep->Fill(deposit);
	    if(InTime)  {GEM_ha_deps_000->Fill(log_d);}
	    if(OutTime) {GEM_ha_deps_250->Fill(log_d);}
	  }
	  
	  switch (abs(pid)%1000) {
	    // Leptons
	  case 11:  
	    GEM_el_hits->Fill(log_e,log_t); GEM_el_deposits->Fill(log_d,log_t); GEM_el_ekindep->Fill(log_e, log_d); GEM_el_tof->Fill(log_t);       GEM_el_time->Fill(time);
	    GEM_el_kins->Fill(log_e);       GEM_el_deps->Fill(log_d);           GEM_el_lindep->Fill(deposit);  GEM_el_process->Fill(process);      GEM_el_linkin->Fill(energy);
	    station==1?G11_el_hits_R->Fill(gemr):G21_el_hits_R->Fill(gemr);
	    station==1?G11_el_hits_E->Fill(roll):G21_el_hits_E->Fill(roll);
	    if(InTime)  { GEM_el_deps_000->Fill(log_d); GEM_el_proc000->Fill(process); station==1?G11_el_hits_000_R->Fill(gemr):G21_el_hits_000_R->Fill(gemr);}
	    if(OutTime) { GEM_el_deps_250->Fill(log_d); GEM_el_proc250->Fill(process); station==1?G11_el_hits_250_R->Fill(gemr):G21_el_hits_250_R->Fill(gemr);}
	    if(deposit > hipMin) { GEM_HIP_id->Fill(1); }
	    if(station==1) {G11_el_deps->Fill(log_d); G11_el_lindep->Fill(deposit); G11_el_kins->Fill(log_e); G11_el_linkin->Fill(energy);}
	    else           {G21_el_deps->Fill(log_d); G21_el_lindep->Fill(deposit); G21_el_kins->Fill(log_e); G21_el_linkin->Fill(energy);}
	    if(isGE11L1)   { G11_L1_el_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_el_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_el_hits_E->Fill(roll); } }
	    if(isGE11L2)   { G11_L2_el_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_el_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_el_hits_E->Fill(roll); } } 
	    if(isGE11Odd)  G11_Od_el_hits_R->Fill(gemr); 
	    if(isGE11Even) G11_Ev_el_hits_R->Fill(gemr); 
	    break;
	  case 13:  
	    GEM_mu_hits->Fill(log_e,log_t); GEM_mu_deposits->Fill(log_d,log_t); GEM_mu_ekindep->Fill(log_e, log_d); GEM_mu_tof->Fill(log_t);       GEM_mu_time->Fill(time);
	    GEM_mu_kins->Fill(log_e);       GEM_mu_deps->Fill(log_d);           GEM_mu_lindep->Fill(deposit);  GEM_mu_process->Fill(process);      GEM_mu_linkin->Fill(energy);
	    station==1?G11_mu_hits_R->Fill(gemr):G21_mu_hits_R->Fill(gemr);
	    station==1?G11_mu_hits_E->Fill(roll):G21_mu_hits_E->Fill(roll);
	    if(InTime)  { GEM_mu_deps_000->Fill(log_d); GEM_mu_proc000->Fill(process); station==1?G11_mu_hits_000_R->Fill(gemr):G21_mu_hits_000_R->Fill(gemr);}
	    if(OutTime) { GEM_mu_deps_250->Fill(log_d); GEM_mu_proc250->Fill(process); station==1?G11_mu_hits_250_R->Fill(gemr):G21_mu_hits_250_R->Fill(gemr);}
	    if(deposit > hipMin) { GEM_HIP_id->Fill(2); }
	    if(station==1) {G11_mu_deps->Fill(log_d); G11_mu_lindep->Fill(deposit); G11_mu_kins->Fill(log_e); G11_mu_linkin->Fill(energy);}
	    else           {G21_mu_deps->Fill(log_d); G21_mu_lindep->Fill(deposit); G21_mu_kins->Fill(log_e); G21_mu_linkin->Fill(energy);}
	    if(isGE11L1)   { G11_L1_mu_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_mu_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_mu_hits_E->Fill(roll); } }
	    if(isGE11L2)   { G11_L2_mu_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_mu_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_mu_hits_E->Fill(roll); } }
	    if(isGE11Odd)  G11_Od_mu_hits_R->Fill(gemr); 
	    if(isGE11Even) G11_Ev_mu_hits_R->Fill(gemr); 
	    break;
	    // Pions
	  case 111: 
	  case 211: 
	    GEM_pi_hits->Fill(log_e,log_t); GEM_pi_deposits->Fill(log_d,log_t); GEM_pi_ekindep->Fill(log_e, log_d); GEM_pi_tof->Fill(log_t);       GEM_pi_time->Fill(time);
	    GEM_pi_kins->Fill(log_e);       GEM_pi_deps->Fill(log_d);           GEM_pi_lindep->Fill(deposit);  GEM_pi_process->Fill(process);      GEM_pi_linkin->Fill(energy);
	    station==1?G11_pi_hits_R->Fill(gemr):G21_pi_hits_R->Fill(gemr);
	    station==1?G11_pi_hits_E->Fill(roll):G21_pi_hits_E->Fill(roll);
	    if(InTime)  { GEM_pi_deps_000->Fill(log_d); GEM_pi_proc000->Fill(process); station==1?G11_pi_hits_000_R->Fill(gemr):G21_pi_hits_000_R->Fill(gemr);}
	    if(OutTime) { GEM_pi_deps_250->Fill(log_d); GEM_pi_proc250->Fill(process); station==1?G11_pi_hits_250_R->Fill(gemr):G21_pi_hits_250_R->Fill(gemr);}
	    if(deposit > hipMin) { GEM_HIP_id->Fill(3); }
	    if(station==1) {G11_pi_deps->Fill(log_d); G11_pi_lindep->Fill(deposit); G11_pi_kins->Fill(log_e); G11_pi_linkin->Fill(energy);}
	    else           {G21_pi_deps->Fill(log_d); G21_pi_lindep->Fill(deposit); G21_pi_kins->Fill(log_e); G21_pi_linkin->Fill(energy);}
	    if(isGE11L1)   { G11_L1_pi_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_pi_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_pi_hits_E->Fill(roll); } }
	    if(isGE11L2)   { G11_L2_pi_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_pi_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_pi_hits_E->Fill(roll); } }
	    if(isGE11Odd)  G11_Od_pi_hits_R->Fill(gemr); 
	    if(isGE11Even) G11_Ev_pi_hits_R->Fill(gemr); 
	    break;
	    // Kaons
	  case 130:
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
	    GEM_ka_hits->Fill(log_e,log_t); GEM_ka_deposits->Fill(log_d,log_t); GEM_ka_ekindep->Fill(log_e, log_d); GEM_ka_tof->Fill(log_t);       GEM_ka_time->Fill(time);
	    GEM_ka_kins->Fill(log_e);       GEM_ka_deps->Fill(log_d);           GEM_ka_lindep->Fill(deposit);  GEM_ka_process->Fill(process);      GEM_ka_linkin->Fill(energy);
	    station==1?G11_ka_hits_R->Fill(gemr):G21_ka_hits_R->Fill(gemr);
	    station==1?G11_ka_hits_E->Fill(roll):G21_ka_hits_E->Fill(roll);
	    if(InTime)  { GEM_ka_deps_000->Fill(log_d); GEM_ka_proc000->Fill(process); station==1?G11_ka_hits_000_R->Fill(gemr):G21_ka_hits_000_R->Fill(gemr);}
	    if(OutTime) { GEM_ka_deps_250->Fill(log_d); GEM_ka_proc250->Fill(process); station==1?G11_ka_hits_250_R->Fill(gemr):G21_ka_hits_250_R->Fill(gemr);}
	    if(deposit > hipMin) { GEM_HIP_id->Fill(4); }
	    if(station==1) {G11_ka_deps->Fill(log_d); G11_ka_lindep->Fill(deposit); G11_ka_kins->Fill(log_e); G11_ka_linkin->Fill(energy);}
	    else           {G21_ka_deps->Fill(log_d); G21_ka_lindep->Fill(deposit); G21_ka_kins->Fill(log_e); G21_ka_linkin->Fill(energy);}
	    if(isGE11L1)   { G11_L1_ka_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_ka_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_ka_hits_E->Fill(roll); } }
	    if(isGE11L2)   { G11_L2_ka_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_ka_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_ka_hits_E->Fill(roll); } }
	    if(isGE11Odd)  G11_Od_ka_hits_R->Fill(gemr); 
	    if(isGE11Even) G11_Ev_ka_hits_R->Fill(gemr); 
	    break;
	    // Other Hadrons
	  default:
	    GEM_OH_hits->Fill(log_e,log_t); GEM_OH_deposits->Fill(log_d,log_t); GEM_OH_ekindep->Fill(log_e, log_d); GEM_OH_tof->Fill(log_t);       GEM_OH_time->Fill(time); 
	    GEM_OH_kins->Fill(log_e);       GEM_OH_deps->Fill(log_d);           GEM_OH_lindep->Fill(deposit);  GEM_OH_process->Fill(process);      GEM_OH_linkin->Fill(energy);
	    station==1?G11_OH_hits_R->Fill(gemr):G21_OH_hits_R->Fill(gemr);
	    station==1?G11_OH_hits_E->Fill(roll):G21_OH_hits_E->Fill(roll);
	    if(InTime)  { GEM_OH_deps_000->Fill(log_d); GEM_OH_proc000->Fill(process); station==1?G11_OH_hits_000_R->Fill(gemr):G21_OH_hits_000_R->Fill(gemr);} 
	    if(OutTime) { GEM_OH_deps_250->Fill(log_d); GEM_OH_proc250->Fill(process); station==1?G11_OH_hits_250_R->Fill(gemr):G21_OH_hits_250_R->Fill(gemr);}
	    if(deposit > hipMin) { GEM_HIP_id->Fill(9); }
	    if(station==1) {G11_OH_deps->Fill(log_d); G11_OH_lindep->Fill(deposit); G11_OH_kins->Fill(log_e); G11_OH_linkin->Fill(energy);}
	    else           {G21_OH_deps->Fill(log_d); G21_OH_lindep->Fill(deposit); G21_OH_kins->Fill(log_e); G21_OH_linkin->Fill(energy);}
	    if(isGE11L1)   { G11_L1_OH_hits_R->Fill(gemr); if(isGE11Odd) { G11_L1Odd_OH_hits_E->Fill(roll); } if(isGE11Even) { G11_L1Even_OH_hits_E->Fill(roll); } }
	    if(isGE11L2)   { G11_L2_OH_hits_R->Fill(gemr); if(isGE11Odd) { G11_L2Odd_OH_hits_E->Fill(roll); } if(isGE11Even) { G11_L2Even_OH_hits_E->Fill(roll); } }
	    if(isGE11Odd)  G11_Od_OH_hits_R->Fill(gemr); 
	    if(isGE11Even) G11_Ev_OH_hits_R->Fill(gemr); 
	    std::cout<<"GEM :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	    std::cout<<" 10 log (tof) = "<<log_t<<" [ns] and 10 log (E) = "<<log_d<<" [keV] catalogued as OTHER HADRON"<<std::endl;
	    break;
	  } // close switch
	} // close else
	GEM_hits_tof->Fill(log_t);
	GEM_hits_eta->Fill(fabs(GEMGlobalPoint.eta()));
	GEM_hits_phi->Fill(GEMGlobalPoint.phi());
	GEM_hits_lin->Fill(time);     
      } //close GEM Station 1 or 2
      else { // station == 0 => ME0
      
        // if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::ME0){ // Only ME0

	// ME0 Geometry
	// ============
	// ME0DetId rollId(theDetUnitId);
	// const ME0EtaPartition* rollasociated = me0Geom->etaPartition(rollId);
	GEMDetId rollId(theDetUnitId);
	const GEMEtaPartition* rollasociated = gemGeom->etaPartition(rollId);
	const BoundPlane & ME0Surface = rollasociated->surface(); 
	GlobalPoint ME0GlobalPoint = ME0Surface.toGlobal((*iHit).localPosition());
	// GlobalPoint ME0GlobalEntry = ME0Surface.toGlobal((*iHit).entryPoint());
	// GlobalPoint ME0GlobalExit  = ME0Surface.toGlobal((*iHit).exitPoint());
	int roll    = rollId.roll();
	
	if(phys_debug) {
	  std::cout<<"ME0 SimHit in "<<std::setw(12)<<(int)rollId<<std::setw(24)<<rollId;
	  std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<ME0GlobalPoint.z();
	  std::cout<<" | r = "<<std::setw(12)<<ME0GlobalPoint.mag()<<" | phi = "<<std::setw(12)<<ME0GlobalPoint.phi()<<" | eta = "<<std::setw(12)<<ME0GlobalPoint.eta();
	  std::cout<<" | global position = "<<ME0GlobalPoint<<std::endl;
	}
	
	ME0hit = 1; if(deposit > hipMin) { ME0hip = 1; }
	
	double ME0_GlobalPoint_R = sqrt(pow(ME0GlobalPoint.x(),2)+pow(ME0GlobalPoint.y(),2));
	double me0r = ME0_GlobalPoint_R;
	ME0_XY->Fill(ME0GlobalPoint.x(), ME0GlobalPoint.y());             
	ME0_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));  
	bool InTime  = (*iHit).timeOfFlight()<250;
	bool OutTime = (*iHit).timeOfFlight()>250;
	bool BXTime  = (*iHit).timeOfFlight()<25;
	bool OXTime  = (*iHit).timeOfFlight()>25;
	if(InTime) { // considered prompt
	  ME0_000ns_XY->Fill(ME0GlobalPoint.x(), ME0GlobalPoint.y());
	  ME0_000ns_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));
	}
	if(OutTime) { // considered neutron background
	  ME0_250ns_XY->Fill(ME0GlobalPoint.x(), ME0GlobalPoint.y());
	  ME0_250ns_RZ->Fill(fabs(ME0GlobalPoint.z()), fabs(ME0_GlobalPoint_R));
	}
	
	if(phys_debug) std::cout<<"SimHit from particle = "<<pid<<" created in process = "<<process<<" aka "<<proc[process]<<std::endl;
	// std::cout<<"SimHit from particle = "<<pid<<" created in process = "<<process<<" aka "<<proc[process]<<std::endl;
	
	// std::cout<<"ME0 :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	// std::cout<<" 10 log (tof) = "<<log_t<<" [ns] and 10 log (E) = "<<log_d<<" [MeV]"<<std::endl;
	
	// All Particles
	ME0_All_hits_R->Fill(me0r); ME0_All_hits_E->Fill(roll);
	ME0_All_deps->Fill(log_d);  ME0_All_lindep->Fill(deposit);
	ME0_All_kins->Fill(log_e);  ME0_All_linkin->Fill(energy);
	
	ME0_All_lindep_roll->Fill(deposit,roll);
	switch(roll) {
	case 1: ME0_All_lindep_eta01->Fill(deposit); break;
	case 2: ME0_All_lindep_eta02->Fill(deposit); break;
	case 3: ME0_All_lindep_eta03->Fill(deposit); break;
	case 4: ME0_All_lindep_eta04->Fill(deposit); break;
	case 5: ME0_All_lindep_eta05->Fill(deposit); break;
	case 6: ME0_All_lindep_eta06->Fill(deposit); break;
	case 7: ME0_All_lindep_eta07->Fill(deposit); break;
	case 8: ME0_All_lindep_eta08->Fill(deposit); break;
	}
	
	if(deposit > hipMin) {
	  ME0_HIP_hits_R->Fill(me0r); ME0_HIP_hits_E->Fill(roll);
	  ME0_HIP_deps->Fill(log_d);  ME0_HIP_lindep->Fill(deposit);
	  ME0_HIP_kins->Fill(log_e);  ME0_HIP_linkin->Fill(energy);
	}
	if(InTime)  {ME0_All_hits_000_R->Fill(me0r); ME0_All_hits_000_E->Fill(roll);}
	if(OutTime) {ME0_All_hits_250_R->Fill(me0r); ME0_All_hits_250_E->Fill(roll);}
	if(BXTime)  {ME0_All_hits_00_R->Fill(me0r);  ME0_All_hits_00_E->Fill(roll);}
	if(OXTime)  {ME0_All_hits_25_R->Fill(me0r);  ME0_All_hits_25_E->Fill(roll);}
	
	if(tech_debug) std::cout<<"Distinction by Particle Type :: part 1"<<std::endl;
	
	if(abs(pid)==22)        { // PHOTON
	  ME0_g_hits->Fill(log_e,log_t); ME0_g_deposits->Fill(log_d,log_t); ME0_g_ekindep->Fill(log_e, log_d); ME0_g_kins->Fill(log_e);   ME0_g_process->Fill(process); ME0_g_linkin->Fill(energy);
	  ME0_g_deps->Fill(log_d);       ME0_g_tof->Fill(log_t);            ME0_g_time->Fill(time);            ME0_g_hits_R->Fill(me0r);  ME0_g_lindep->Fill(deposit);  ME0_g_hits_E->Fill(roll);
	  if(InTime)  { ME0_g_hits_000_R->Fill(me0r); ME0_g_deps_000->Fill(log_d);  ME0_g_proc000->Fill(process);}
	  if(OutTime) { ME0_g_hits_250_R->Fill(me0r); ME0_g_deps_250->Fill(log_d);  ME0_g_proc250->Fill(process);}
	  if(BXTime)  { ME0_g_hits_00_R->Fill(me0r);}
	  if(OXTime)  { ME0_g_hits_25_R->Fill(me0r);}
	  if(deposit > hipMin) { ME0_HIP_id->Fill(7); }
	}
	else if(abs(pid)==2212) { // PROTON
	  ME0_p_hits->Fill(log_e,log_t); ME0_p_deposits->Fill(log_d,log_t); ME0_p_ekindep->Fill(log_e, log_d); ME0_p_kins->Fill(log_e);   ME0_p_process->Fill(process); ME0_p_linkin->Fill(energy);
	  ME0_p_deps->Fill(log_d);       ME0_p_tof->Fill(log_t);            ME0_p_time->Fill(time);            ME0_p_hits_R->Fill(me0r);  ME0_p_lindep->Fill(deposit);  ME0_p_hits_E->Fill(roll);
	  if(InTime)  { ME0_p_hits_000_R->Fill(me0r); ME0_p_deps_000->Fill(log_d); ME0_p_proc000->Fill(process);}
	  if(OutTime) { ME0_p_hits_250_R->Fill(me0r); ME0_p_deps_250->Fill(log_d); ME0_p_proc250->Fill(process);}
	  if(BXTime)  { ME0_p_hits_00_R->Fill(me0r);}
	  if(OXTime)  { ME0_p_hits_25_R->Fill(me0r);}
	  if(deposit > hipMin) { ME0_HIP_id->Fill(5); }
	}
	else if(abs(pid)==2112) { // NEUTRON
	  ME0_n_hits->Fill(log_e,log_t); ME0_n_deposits->Fill(log_d,log_t); ME0_n_ekindep->Fill(log_e, log_d); ME0_n_kins->Fill(log_e);   ME0_n_process->Fill(process); ME0_n_linkin->Fill(energy);
	  ME0_n_deps->Fill(log_d);       ME0_n_tof->Fill(log_t);            ME0_n_time->Fill(time);            ME0_n_hits_R->Fill(me0r);  ME0_n_lindep->Fill(deposit);  ME0_n_hits_E->Fill(roll);
	  if(InTime)  { ME0_n_hits_000_R->Fill(me0r); ME0_n_deps_000->Fill(log_d); ME0_n_proc000->Fill(process);}
	  if(OutTime) { ME0_n_hits_250_R->Fill(me0r); ME0_n_deps_250->Fill(log_d); ME0_n_proc250->Fill(process);}
	  if(BXTime)  { ME0_n_hits_00_R->Fill(me0r);}
	  if(OXTime)  { ME0_n_hits_25_R->Fill(me0r);}
	  if(deposit > hipMin) { ME0_HIP_id->Fill(6); }
	}
	else if(abs(pid)>1E9)   { // NUCLEI
	  ME0_N_hits->Fill(log_e,log_t); ME0_N_deposits->Fill(log_d,log_t); ME0_N_ekindep->Fill(log_e, log_d); ME0_N_kins->Fill(log_e);   ME0_N_process->Fill(process); ME0_N_linkin->Fill(energy);
	  ME0_N_deps->Fill(log_d);       ME0_N_tof->Fill(log_t);            ME0_N_time->Fill(time);            ME0_N_hits_R->Fill(me0r);  ME0_N_lindep->Fill(deposit);  ME0_N_hits_E->Fill(roll);
	  if(InTime)  { ME0_N_hits_000_R->Fill(me0r); ME0_N_deps_000->Fill(log_d); ME0_N_proc000->Fill(process);}
	  if(OutTime) { ME0_N_hits_250_R->Fill(me0r); ME0_N_deps_250->Fill(log_d); ME0_N_proc250->Fill(process);}
	  if(BXTime)  { ME0_N_hits_00_R->Fill(me0r);}
	  if(OXTime)  { ME0_N_hits_25_R->Fill(me0r);}
	  if(deposit > hipMin) { ME0_HIP_id->Fill(8); }
	  std::cout<<"ME0 :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	  std::cout<<" 10 log (tof) = "<<log_t<<" [ns] and 10 log (E) = "<<log_d<<" [keV] catalogued as NUCLEI"<<std::endl;
	  // Nuclei code: 100ZZZAAA0
	  // - alpha    (Z = 2 A= 4) 1000020040
	  // - carbon   (Z = 6 A=12) 1000060120
	  int nuclei_Code = abs(pid)/10;
	  double nuclei_A = (nuclei_Code)%1000;
	  double nuclei_Z = (nuclei_Code/1000)%1000;
	  ME0_Nuclei_Z->Fill(nuclei_Z);
	  ME0_Nuclei_A->Fill(nuclei_A);
	  ME0_Nuclei_A_Z->Fill(nuclei_Z, nuclei_A);
	  ME0_Nuclei_List->Fill(nuclei_Z);
	}
	else {
	  
	  if(tech_debug) std::cout<<"Distinction by Particle Type :: part 2"<<std::endl;
	  
	  if(abs(pid)!=11 && abs(pid)!=13) {
	    if(tech_debug) std::cout<<"PID = "<<pid<<" --> Filling hadron histograms"<<std::endl;
	    ME0_ha_kins->Fill(log_e); ME0_ha_deps->Fill(log_d); ME0_ha_tof->Fill(log_t); ME0_ha_time->Fill(time);  ME0_ha_lindep->Fill(deposit);
	    if(InTime)  { ME0_ha_deps_000->Fill(log_d);}
	    if(OutTime) { ME0_ha_deps_250->Fill(log_d);}
	  }
	  
	  switch (abs(pid)%1000) {
	    // Leptons
	  case 11:  
	    ME0_el_hits->Fill(log_e,log_t); ME0_el_deposits->Fill(log_d,log_t); ME0_el_ekindep->Fill(log_e, log_d); ME0_el_tof->Fill(log_t);       ME0_el_time->Fill(time);        ME0_el_linkin->Fill(energy);
	    ME0_el_kins->Fill(log_e);       ME0_el_deps->Fill(log_d);           ME0_el_hits_R->Fill(me0r);          ME0_el_lindep->Fill(deposit);  ME0_el_process->Fill(process);  ME0_el_hits_E->Fill(roll);
	    if(InTime)  { ME0_el_hits_000_R->Fill(me0r); ME0_el_deps_000->Fill(log_d); ME0_el_proc000->Fill(process);}
	    if(OutTime) { ME0_el_hits_250_R->Fill(me0r); ME0_el_deps_250->Fill(log_d); ME0_el_proc250->Fill(process);}
	    if(BXTime)  { ME0_el_hits_00_R->Fill(me0r);}
	    if(OXTime)  { ME0_el_hits_25_R->Fill(me0r);}
	    if(deposit > hipMin) { ME0_HIP_id->Fill(1); }
	    break;
	  case 13:  
	    ME0_mu_hits->Fill(log_e,log_t); ME0_mu_deposits->Fill(log_d,log_t); ME0_mu_ekindep->Fill(log_e, log_d); ME0_mu_tof->Fill(log_t);       ME0_mu_time->Fill(time);        ME0_mu_linkin->Fill(energy); 
	    ME0_mu_kins->Fill(log_e);       ME0_mu_deps->Fill(log_d);           ME0_mu_hits_R->Fill(me0r);          ME0_mu_lindep->Fill(deposit);  ME0_mu_process->Fill(process);  ME0_mu_hits_E->Fill(roll);
	    if(InTime)  { ME0_mu_hits_000_R->Fill(me0r); ME0_mu_deps_000->Fill(log_d); ME0_mu_proc000->Fill(process);}
	    if(OutTime) { ME0_mu_hits_250_R->Fill(me0r); ME0_mu_deps_250->Fill(log_d); ME0_mu_proc250->Fill(process);}
	    if(BXTime)  { ME0_mu_hits_00_R->Fill(me0r);}
	    if(OXTime)  { ME0_mu_hits_25_R->Fill(me0r);}
	    if(deposit > hipMin) { ME0_HIP_id->Fill(2); }
	    break;
	    // Pions
	  case 111: 
	  case 211: 
	    ME0_pi_hits->Fill(log_e,log_t); ME0_pi_deposits->Fill(log_d,log_t); ME0_pi_ekindep->Fill(log_e, log_d); ME0_pi_tof->Fill(log_t);       ME0_pi_time->Fill(time);        ME0_pi_linkin->Fill(energy);
	    ME0_pi_kins->Fill(log_e);       ME0_pi_deps->Fill(log_d);           ME0_pi_hits_R->Fill(me0r);          ME0_pi_lindep->Fill(deposit);  ME0_pi_process->Fill(process);  ME0_pi_hits_E->Fill(roll);
	    if(InTime)  { ME0_pi_hits_000_R->Fill(me0r); ME0_pi_deps_000->Fill(log_d); ME0_mu_proc000->Fill(process);}
	    if(OutTime) { ME0_pi_hits_250_R->Fill(me0r); ME0_pi_deps_250->Fill(log_d); ME0_mu_proc250->Fill(process);}
	    if(BXTime)  { ME0_pi_hits_00_R->Fill(me0r);}
	    if(OXTime)  { ME0_pi_hits_25_R->Fill(me0r);}
	    if(deposit > hipMin) { ME0_HIP_id->Fill(3); }
	    break;
	    // Kaons
	  case 130:
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
	    ME0_ka_hits->Fill(log_e,log_t); ME0_ka_deposits->Fill(log_d,log_t); ME0_ka_ekindep->Fill(log_e, log_d); ME0_ka_tof->Fill(log_t); ME0_ka_time->Fill(time);      ME0_ka_linkin->Fill(energy);
	    ME0_ka_kins->Fill(log_e);       ME0_ka_deps->Fill(log_d);           ME0_ka_hits_R->Fill(me0r);  ME0_ka_lindep->Fill(deposit);  ME0_ka_process->Fill(process);  ME0_ka_hits_E->Fill(roll);
	    if(InTime)  { ME0_ka_hits_000_R->Fill(me0r); ME0_ka_deps_000->Fill(log_d); ME0_ka_proc000->Fill(process);}
	    if(OutTime) { ME0_ka_hits_250_R->Fill(me0r); ME0_ka_deps_250->Fill(log_d); ME0_ka_proc250->Fill(process);}
	    if(BXTime)  { ME0_ka_hits_00_R->Fill(me0r);}
	    if(OXTime)  { ME0_ka_hits_25_R->Fill(me0r);}
	    if(deposit > hipMin) { ME0_HIP_id->Fill(4); }
	    break;
	    // Other Hadrons
	  default:   
	    ME0_OH_hits->Fill(log_e,log_t); ME0_OH_deposits->Fill(log_d,log_t); ME0_OH_ekindep->Fill(log_e, log_d); ME0_OH_tof->Fill(log_t); ME0_OH_time->Fill(time);     ME0_OH_linkin->Fill(energy);
	    ME0_OH_kins->Fill(log_e);       ME0_OH_deps->Fill(log_d);           ME0_OH_hits_R->Fill(me0r);  ME0_OH_lindep->Fill(deposit);  ME0_OH_process->Fill(process); ME0_OH_hits_E->Fill(roll);
	    if(InTime)  { ME0_OH_hits_000_R->Fill(me0r); ME0_OH_deps_000->Fill(log_d); ME0_OH_proc000->Fill(process);}
	    if(OutTime) { ME0_OH_hits_250_R->Fill(me0r); ME0_OH_deps_250->Fill(log_d); ME0_OH_proc250->Fill(process);}
	    if(BXTime)  { ME0_OH_hits_00_R->Fill(me0r);}
	    if(OXTime)  { ME0_OH_hits_25_R->Fill(me0r);}
	    if(deposit > hipMin) { ME0_HIP_id->Fill(9); }
	    std::cout<<"ME0 :: SimHit from Particle id "<<pid<<" with time of flight "<<(*iHit).timeOfFlight()<<" [ns] and energy deposit "<<(*iHit).energyLoss()<< " [GeV]";
	    std::cout<<" 10 log (tof) = "<<log_t<<" [ns] and 10 log (E) = "<<log_d<<" [keV] catalogued as OTHER HADRONS"<<std::endl;
	    break;
	  } // close switch
	} // close else
	ME0_hits_tof->Fill(log_t);
	ME0_hits_eta->Fill(fabs(ME0GlobalPoint.eta()));
	ME0_hits_phi->Fill(ME0GlobalPoint.phi());
	ME0_hits_lin->Fill(time);     
      } // close ME0
    } // close GEM
  } // close simhitloop
    
  // Fill EventCounter
  if(ME0hit)  EventCounter->Fill(2);
  if(GE11hit) EventCounter->Fill(3);
  if(GE21hit) EventCounter->Fill(4);
  if(ME0hip)  EventCounter->Fill(5);
  if(GE11hip) EventCounter->Fill(6);
  if(GE21hip) EventCounter->Fill(7);
  
} // end analyze method

double
MyME0SimHitAnalyzer::getLumi(int pu, int space, int com)
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

double MyME0SimHitAnalyzer::getPU(double lumi, int space, int com)
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
MyME0SimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyME0SimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyME0SimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  // is not called anymore ...

  if(tech_debug) std::cout<<"[MyME0SimHitAnalyzer :: BeginRun]"<<std::endl;

  // iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  // iSetup.get<MuonGeometryRecord>().get(cscGeom);
  // iSetup.get<MuonGeometryRecord>().get(dtGeom);
  // iSetup.get<MuonGeometryRecord>().get(gemGeom);
  // iSetup.get<MuonGeometryRecord>().get(me0Geom);

  // const edm::ESTransientHandle<CSCGeometry> cscGeom = iSetup.getTransientHandle(cscGeom_Token);
  // const edm::ESTransientHandle<DTGeometry>  dtGeom  = iSetup.getTransientHandle(dtGeom_Token);
  // const edm::ESTransientHandle<GEMGeometry> gemGeom = iSetup.getTransientHandle(gemGeom_Token);
  // const edm::ESTransientHandle<RPCGeometry> rpcGeom = iSetup.getTransientHandle(rpcGeom_Token);

}


// ------------ method called when ending the processing of a run  ------------
void 
MyME0SimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{

}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyME0SimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyME0SimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyME0SimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyME0SimHitAnalyzer);

//  LocalWords:  str
