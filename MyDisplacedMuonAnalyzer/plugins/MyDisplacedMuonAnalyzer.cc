// -*- C++ -*-
//
// Package:    MyDisplacedMuonAnalyzer
// Class:      MyDisplacedMuonAnalyzer
// 
/**\class MyDisplacedMuonAnalyzer MyDisplacedMuonAnalyzer.cc MyAnalyzers/MyDisplacedMuonAnalyzer/plugins/MyDisplacedMuonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 17 Dec 2014 15:05:11 GMT
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

// ROOT include files
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
#include "TLorentzVector.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include "DataFormats/MuonDetId/interface/DTWireId.h"

#include "RecoMuon/TrackingTools/interface/SegmentsTrackAssociator.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class declaration
//

class MyDisplacedMuonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyDisplacedMuonAnalyzer(const edm::ParameterSet&);
      ~MyDisplacedMuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // bool myptsort(double*, double*);
  // bool myetasort(double*, double*);

      // ----------member data ---------------------------
  edm::InputTag STAMuLabel1, STAMuLabel2, STAMuLabel3, STAMuLabel4;
 
  std::string rootFileName;
  bool physDebug, techDebug, genDebug, mcTruthMatch, recoDRMatch;

  const reco::GenParticle     *g, *b, *p, *q, *r, *m, *n;

  TFile * outputfile;
  TDirectoryFile * GenPhotons, * GenMuons, * SimMuons, * SAMuons, * RSAMuons, * SAMuons_Vtx, * RSAMuons_Vtx;

  // GEN Dark Photon Histograms
  TH1F * Gen_All_DarkPhot_pt, * Gen_All_DarkPhot_eta, * Gen_All_DarkPhot_phi, * Gen_All_DarkPhot_lxy, * Gen_All_DarkPhot_lz, * Gen_All_DarkPhot_m, * Gen_All_DarkPhot_mt;
  TH1F * Gen_DarkPhot_pt1_pt, * Gen_DarkPhot_pt1_eta, * Gen_DarkPhot_pt1_phi, * Gen_DarkPhot_pt1_lxy, * Gen_DarkPhot_pt1_lz, * Gen_DarkPhot_pt1_m, * Gen_DarkPhot_pt1_mt;
  TH1F * Gen_DarkPhot_pt2_pt, * Gen_DarkPhot_pt2_eta, * Gen_DarkPhot_pt2_phi, * Gen_DarkPhot_pt2_lxy, * Gen_DarkPhot_pt2_lz, * Gen_DarkPhot_pt2_m, * Gen_DarkPhot_pt2_mt;
  // GEN DiMuon Histograms
  TH1F * Gen_All_DiMuon_pt, * Gen_All_DiMuon_eta, * Gen_All_DiMuon_phi, * Gen_All_DiMuon_lxy, * Gen_All_DiMuon_lz, * Gen_All_DiMuon_m, * Gen_All_DiMuon_dR;
  // GEN Muon Histograms
  TH1F * Gen_All_Muon_pt, * Gen_All_Muon_eta, * Gen_All_Muon_phi, * Gen_All_Muon_lxy, * Gen_All_Muon_lz;
  TH1F * Gen_Muon_pt1_pt, * Gen_Muon_pt1_eta, * Gen_Muon_pt1_phi, * Gen_Muon_pt1_lxy, * Gen_Muon_pt1_lz;
  TH1F * Gen_Muon_pt2_pt, * Gen_Muon_pt2_eta, * Gen_Muon_pt2_phi, * Gen_Muon_pt2_lxy, * Gen_Muon_pt2_lz;
  TH1F * Gen_Muon_pt3_pt, * Gen_Muon_pt3_eta, * Gen_Muon_pt3_phi, * Gen_Muon_pt3_lxy, * Gen_Muon_pt3_lz;
  TH1F * Gen_Muon_pt4_pt, * Gen_Muon_pt4_eta, * Gen_Muon_pt4_phi, * Gen_Muon_pt4_lxy, * Gen_Muon_pt4_lz;
  TH1F * Gen_Muon_eta1_pt, * Gen_Muon_eta1_eta, * Gen_Muon_eta1_phi, * Gen_Muon_eta1_lxy, * Gen_Muon_eta1_lz;
  TH1F * Gen_Muon_eta2_pt, * Gen_Muon_eta2_eta, * Gen_Muon_eta2_phi, * Gen_Muon_eta2_lxy, * Gen_Muon_eta2_lz;
  TH1F * Gen_Muon_eta3_pt, * Gen_Muon_eta3_eta, * Gen_Muon_eta3_phi, * Gen_Muon_eta3_lxy, * Gen_Muon_eta3_lz;
  TH1F * Gen_Muon_eta4_pt, * Gen_Muon_eta4_eta, * Gen_Muon_eta4_phi, * Gen_Muon_eta4_lxy, * Gen_Muon_eta4_lz;
  // SIM DiMuon Histograms
  TH1F * Sim_All_DiMuon_pt, * Sim_All_DiMuon_eta, * Sim_All_DiMuon_phi, * Sim_All_DiMuon_lxy, * Sim_All_DiMuon_lz, * Sim_All_DiMuon_m, * Sim_All_DiMuon_dR;
  // SIM Muon Histograms
  TH1F * Sim_All_Muon_pt, * Sim_All_Muon_eta, * Sim_All_Muon_phi, * Sim_All_Muon_lxy, * Sim_All_Muon_lz;
  TH1F * Sim_Muon_pt1_pt, * Sim_Muon_pt1_eta, * Sim_Muon_pt1_phi, * Sim_Muon_pt1_lxy, * Sim_Muon_pt1_lz;
  TH1F * Sim_Muon_pt2_pt, * Sim_Muon_pt2_eta, * Sim_Muon_pt2_phi, * Sim_Muon_pt2_lxy, * Sim_Muon_pt2_lz;
  TH1F * Sim_Muon_pt3_pt, * Sim_Muon_pt3_eta, * Sim_Muon_pt3_phi, * Sim_Muon_pt3_lxy, * Sim_Muon_pt3_lz;
  TH1F * Sim_Muon_pt4_pt, * Sim_Muon_pt4_eta, * Sim_Muon_pt4_phi, * Sim_Muon_pt4_lxy, * Sim_Muon_pt4_lz;
  // RSA DiMuon Histograms
  TH1F * RSA_All_DiMuon_pt, * RSA_All_DiMuon_eta, * RSA_All_DiMuon_phi, * RSA_All_DiMuon_lxy, * RSA_All_DiMuon_lz, * RSA_All_DiMuon_m, * RSA_All_DiMuon_dR;
  // RSA Muon Histograms
  TH1F * RSA_All_Muon_pt, * RSA_All_Muon_eta, * RSA_All_Muon_phi, * RSA_All_Muon_lxy, * RSA_All_Muon_lz;
  TH1F * RSA_Muon_pt1_pt, * RSA_Muon_pt1_eta, * RSA_Muon_pt1_phi, * RSA_Muon_pt1_lxy, * RSA_Muon_pt1_lz;
  TH1F * RSA_Muon_pt2_pt, * RSA_Muon_pt2_eta, * RSA_Muon_pt2_phi, * RSA_Muon_pt2_lxy, * RSA_Muon_pt2_lz;
  TH1F * RSA_Muon_pt3_pt, * RSA_Muon_pt3_eta, * RSA_Muon_pt3_phi, * RSA_Muon_pt3_lxy, * RSA_Muon_pt3_lz;
  TH1F * RSA_Muon_pt4_pt, * RSA_Muon_pt4_eta, * RSA_Muon_pt4_phi, * RSA_Muon_pt4_lxy, * RSA_Muon_pt4_lz;

  // Counters
  int tot_gen, tot_sim, tot_rsa;
  int tot_gen_ev, tot_sim_ev, tot_rsa_ev;
  // int tot_sim_matched, tot_rsa_matched;
};

//
// constants, enums and typedefs
//
double etabinning[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4};
double phibinning[] = {-3.2288, -2.7052, -2.1816, -1.6580, -1.1344, -0.6108, -0.0872, +0.4364, 0.960, 1.483, 2.007, 2.531, 3.0544};
double delta_rel_pt = 0.1;
// double delta_R = 0.1;
double delta_R_gensim = 0.1;
double delta_R_simrec = 0.15;
std::string staMuonLabels[] = {"Stand Alone Muon (SA)", "Refitted Stand Alone Muon (RSA)", "Stand Alone Muon updated at Vtx (SAVtx)", "Refitted Stand Alone Muon updated at Vtx (RSAVtx)"}; 
//
// static data member definitions
//

//
// constructors and destructor
//
MyDisplacedMuonAnalyzer::MyDisplacedMuonAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName   = iConfig.getUntrackedParameter<std::string>("RootFileName");
  physDebug      = iConfig.getUntrackedParameter<bool>("PhysicsDebug");
  techDebug      = iConfig.getUntrackedParameter<bool>("TechnicDebug");
  genDebug       = iConfig.getUntrackedParameter<bool>("GenPartDebug");
  mcTruthMatch   = iConfig.getUntrackedParameter<bool>("MCTruthMatching");
  recoDRMatch    = iConfig.getUntrackedParameter<bool>("RecoDRMatching");
  STAMuLabel1 = iConfig.getParameter<edm::InputTag>("StandAloneTrackCollectionLabel1");
  STAMuLabel2 = iConfig.getParameter<edm::InputTag>("StandAloneTrackCollectionLabel2");
  STAMuLabel3 = iConfig.getParameter<edm::InputTag>("StandAloneTrackCollectionLabel3");
  STAMuLabel4 = iConfig.getParameter<edm::InputTag>("StandAloneTrackCollectionLabel4");

  if(techDebug) std::cout<<"[MyDisplacedMuonAnalyzer :: Constructor]"<<std::endl;
  outputfile      = new TFile(rootFileName.c_str(), "RECREATE" );

  GenPhotons      = new TDirectoryFile("GeneratorDarkPhotons",                 "GenPhotons");
  GenMuons        = new TDirectoryFile("GeneratorMuons",                       "GenMuons");
  SimMuons        = new TDirectoryFile("SimTrackMuons",                        "SimMuons");
  SAMuons         = new TDirectoryFile("StandAloneMuons",                      "SAMuons");
  RSAMuons        = new TDirectoryFile("RefittedStandAloneMuons",              "RSAMuons");
  SAMuons_Vtx     = new TDirectoryFile("StandAloneMuons_UpdatedAtVtx",         "SAMuonsVtx");
  RSAMuons_Vtx    = new TDirectoryFile("RefittedStandAloneMuons_UpdatedAtVtx", "RSAMuonsVtx");

  // Counters
  tot_gen  = 0; tot_sim = 0; tot_rsa = 0;
  // tot_sim_matched = 0; tot_rsa_matched = 0;


  // GEN Histograms
  Gen_All_DarkPhot_pt  = new TH1F("Gen_All_DarkPhot_pt",  "Gen_All_DarkPhot_pt",   100,  0.0, 500.00);
  Gen_All_DarkPhot_eta = new TH1F("Gen_All_DarkPhot_eta", "Gen_All_DarkPhot_eta",   48,  0.0,   2.40); 
  Gen_All_DarkPhot_phi = new TH1F("Gen_All_DarkPhot_phi", "Gen_All_DarkPhot_phi",   12, -3.14,  3.14); 
  Gen_All_DarkPhot_lxy = new TH1F("Gen_All_DarkPhot_lxy", "Gen_All_DarkPhot_lxy",  160,  0.0,   8.00);
  Gen_All_DarkPhot_lz  = new TH1F("Gen_All_DarkPhot_lz",  "Gen_All_DarkPhot_lz",   150,  0.0,  15.00);
  Gen_All_DarkPhot_m   = new TH1F("Gen_All_DarkPhot_m",   "Gen_All_DarkPhot_m",      3,  0.3,   0.50);
  Gen_All_DarkPhot_mt  = new TH1F("Gen_All_DarkPhot_mt",  "Gen_All_DarkPhot_mt",   100,  0.0, 100.00);
  Gen_DarkPhot_pt1_pt  = new TH1F("Gen_DarkPhot_pt1_pt",  "Gen_DarkPhot_pt1_pt",   100,  0.0, 500.00);
  Gen_DarkPhot_pt1_eta = new TH1F("Gen_DarkPhot_pt1_eta", "Gen_DarkPhot_pt1_eta",   48,  0.0,   2.40); 
  Gen_DarkPhot_pt1_phi = new TH1F("Gen_DarkPhot_pt1_phi", "Gen_DarkPhot_pt1_phi",   12, -3.14,  3.14); 
  Gen_DarkPhot_pt1_lxy = new TH1F("Gen_DarkPhot_pt1_lxy", "Gen_DarkPhot_pt1_lxy",  160,  0.0,   8.00);
  Gen_DarkPhot_pt1_lz  = new TH1F("Gen_DarkPhot_pt1_lz",  "Gen_DarkPhot_pt1_lz",   150,  0.0,  15.00);
  Gen_DarkPhot_pt1_m   = new TH1F("Gen_DarkPhot_pt1_m",   "Gen_DarkPhot_pt1_m",      3,  0.3,   0.50);
  Gen_DarkPhot_pt1_mt  = new TH1F("Gen_DarkPhot_pt1_mt",  "Gen_DarkPhot_pt1_mt",   100,  0.0, 100.00);
  Gen_DarkPhot_pt2_pt  = new TH1F("Gen_DarkPhot_pt2_pt",  "Gen_DarkPhot_pt2_pt",   100,  0.0, 500.00);
  Gen_DarkPhot_pt2_eta = new TH1F("Gen_DarkPhot_pt2_eta", "Gen_DarkPhot_pt2_eta",   48,  0.0,   2.40); 
  Gen_DarkPhot_pt2_phi = new TH1F("Gen_DarkPhot_pt2_phi", "Gen_DarkPhot_pt2_phi",   12, -3.14,  3.14); 
  Gen_DarkPhot_pt2_lxy = new TH1F("Gen_DarkPhot_pt2_lxy", "Gen_DarkPhot_pt2_lxy",  160,  0.0,   8.00);
  Gen_DarkPhot_pt2_lz  = new TH1F("Gen_DarkPhot_pt2_lz",  "Gen_DarkPhot_pt2_lz",   150,  0.0,  15.00);
  Gen_DarkPhot_pt2_m   = new TH1F("Gen_DarkPhot_pt2_m",   "Gen_DarkPhot_pt2_m",      3,  0.3,   0.50);
  Gen_DarkPhot_pt2_mt  = new TH1F("Gen_DarkPhot_pt2_mt",  "Gen_DarkPhot_pt2_mt",   100,  0.0, 100.00);

  Gen_All_DiMuon_pt  = new TH1F("Gen_All_DiMuon_pt",  "Gen_All_DiMuon_pt",   100,  0.0, 500.00);
  Gen_All_DiMuon_eta = new TH1F("Gen_All_DiMuon_eta", "Gen_All_DiMuon_eta",   48,  0.0,   2.40);
  Gen_All_DiMuon_phi = new TH1F("Gen_All_DiMuon_phi", "Gen_All_DiMuon_phi",   12, -3.14,  3.14);
  Gen_All_DiMuon_lxy = new TH1F("Gen_All_DiMuon_lxy", "Gen_All_DiMuon_lxy",  160,  0.0,   8.00);
  Gen_All_DiMuon_lz  = new TH1F("Gen_All_DiMuon_lz",  "Gen_All_DiMuon_lz",   150,  0.0,  15.00);
  Gen_All_DiMuon_m   = new TH1F("Gen_All_DiMuon_m",   "Gen_All_DiMuon_m",      3,  0.3,   0.50);
  Gen_All_DiMuon_dR  = new TH1F("Gen_All_DiMuon_dR",  "Gen_All_DiMuon_m",    100,  0.0,   5.00);

  Gen_All_Muon_pt  = new TH1F("Gen_All_Muon_pt",  "Gen_All_Muon_pt",   100,  0.0, 500.00);
  Gen_All_Muon_eta = new TH1F("Gen_All_Muon_eta", "Gen_All_Muon_eta",   48,  0.0,   2.40); 
  Gen_All_Muon_phi = new TH1F("Gen_All_Muon_phi", "Gen_All_Muon_phi",   12, -3.14,  3.14); 
  Gen_All_Muon_lxy = new TH1F("Gen_All_Muon_lxy", "Gen_All_Muon_lxy",  160,  0.0,   8.00);
  Gen_All_Muon_lz  = new TH1F("Gen_All_Muon_lz",  "Gen_All_Muon_lz",   150,  0.0,  15.00);
  Gen_Muon_pt1_pt  = new TH1F("Gen_Muon_pt1_pt",  "Gen_Muon_pt1_pt",   100,  0.0, 500.00);
  Gen_Muon_pt1_eta = new TH1F("Gen_Muon_pt1_eta", "Gen_Muon_pt1_eta",   48,  0.0,   2.40); 
  Gen_Muon_pt1_phi = new TH1F("Gen_Muon_pt1_phi", "Gen_Muon_pt1_phi",   12, -3.14,  3.14); 
  Gen_Muon_pt1_lxy = new TH1F("Gen_Muon_pt1_lxy", "Gen_Muon_pt1_lxy",  160,  0.0,   8.00);
  Gen_Muon_pt1_lz  = new TH1F("Gen_Muon_pt1_lz",  "Gen_Muon_pt1_lz",   150,  0.0,  15.00);
  Gen_Muon_pt2_pt  = new TH1F("Gen_Muon_pt2_pt",  "Gen_Muon_pt2_pt",   100,  0.0, 500.00);
  Gen_Muon_pt2_eta = new TH1F("Gen_Muon_pt2_eta", "Gen_Muon_pt2_eta",   48,  0.0,   2.40); 
  Gen_Muon_pt2_phi = new TH1F("Gen_Muon_pt2_phi", "Gen_Muon_pt2_phi",   12, -3.14,  3.14); 
  Gen_Muon_pt2_lxy = new TH1F("Gen_Muon_pt2_lxy", "Gen_Muon_pt2_lxy",  160,  0.0,   8.00);
  Gen_Muon_pt2_lz  = new TH1F("Gen_Muon_pt2_lz",  "Gen_Muon_pt2_lz",   150,  0.0,  15.00);
  Gen_Muon_pt3_pt  = new TH1F("Gen_Muon_pt3_pt",  "Gen_Muon_pt3_pt",   100,  0.0, 500.00);
  Gen_Muon_pt3_eta = new TH1F("Gen_Muon_pt3_eta", "Gen_Muon_pt3_eta",   48,  0.0,   2.40); 
  Gen_Muon_pt3_phi = new TH1F("Gen_Muon_pt3_phi", "Gen_Muon_pt3_phi",   12, -3.14,  3.14); 
  Gen_Muon_pt3_lxy = new TH1F("Gen_Muon_pt3_lxy", "Gen_Muon_pt3_lxy",  160,  0.0,   8.00);
  Gen_Muon_pt3_lz  = new TH1F("Gen_Muon_pt3_lz",  "Gen_Muon_pt3_lz",   150,  0.0,  15.00);
  Gen_Muon_pt4_pt  = new TH1F("Gen_Muon_pt4_pt",  "Gen_Muon_pt4_pt",   100,  0.0, 500.00);
  Gen_Muon_pt4_eta = new TH1F("Gen_Muon_pt4_eta", "Gen_Muon_pt4_eta",   48,  0.0,   2.40); 
  Gen_Muon_pt4_phi = new TH1F("Gen_Muon_pt4_phi", "Gen_Muon_pt4_phi",   12, -3.14,  3.14); 
  Gen_Muon_pt4_lxy = new TH1F("Gen_Muon_pt4_lxy", "Gen_Muon_pt4_lxy",  160,  0.0,   8.00);
  Gen_Muon_pt4_lz  = new TH1F("Gen_Muon_pt4_lz",  "Gen_Muon_pt4_lz",   150,  0.0,  15.00);
  Gen_Muon_eta1_pt  = new TH1F("Gen_Muon_eta1_pt",  "Gen_Muon_eta1_pt",   100,  0.0, 500.00);
  Gen_Muon_eta1_eta = new TH1F("Gen_Muon_eta1_eta", "Gen_Muon_eta1_eta",   48,  0.0,   2.40); 
  Gen_Muon_eta1_phi = new TH1F("Gen_Muon_eta1_phi", "Gen_Muon_eta1_phi",   12, -3.14,  3.14); 
  Gen_Muon_eta1_lxy = new TH1F("Gen_Muon_eta1_lxy", "Gen_Muon_eta1_lxy",  160,  0.0,   8.00);
  Gen_Muon_eta1_lz  = new TH1F("Gen_Muon_eta1_lz",  "Gen_Muon_eta1_lz",   150,  0.0,  15.00);
  Gen_Muon_eta2_pt  = new TH1F("Gen_Muon_eta2_pt",  "Gen_Muon_eta2_pt",   100,  0.0, 500.00);
  Gen_Muon_eta2_eta = new TH1F("Gen_Muon_eta2_eta", "Gen_Muon_eta2_eta",   48,  0.0,   2.40); 
  Gen_Muon_eta2_phi = new TH1F("Gen_Muon_eta2_phi", "Gen_Muon_eta2_phi",   12, -3.14,  3.14); 
  Gen_Muon_eta2_lxy = new TH1F("Gen_Muon_eta2_lxy", "Gen_Muon_eta2_lxy",  160,  0.0,   8.00);
  Gen_Muon_eta2_lz  = new TH1F("Gen_Muon_eta2_lz",  "Gen_Muon_eta2_lz",   150,  0.0,  15.00);
  Gen_Muon_eta3_pt  = new TH1F("Gen_Muon_eta3_pt",  "Gen_Muon_eta3_pt",   100,  0.0, 500.00);
  Gen_Muon_eta3_eta = new TH1F("Gen_Muon_eta3_eta", "Gen_Muon_eta3_eta",   48,  0.0,   2.40); 
  Gen_Muon_eta3_phi = new TH1F("Gen_Muon_eta3_phi", "Gen_Muon_eta3_phi",   12, -3.14,  3.14); 
  Gen_Muon_eta3_lxy = new TH1F("Gen_Muon_eta3_lxy", "Gen_Muon_eta3_lxy",  160,  0.0,   8.00);
  Gen_Muon_eta3_lz  = new TH1F("Gen_Muon_eta3_lz",  "Gen_Muon_eta3_lz",   150,  0.0,  15.00);
  Gen_Muon_eta4_pt  = new TH1F("Gen_Muon_eta4_pt",  "Gen_Muon_eta4_pt",   100,  0.0, 500.00);
  Gen_Muon_eta4_eta = new TH1F("Gen_Muon_eta4_eta", "Gen_Muon_eta4_eta",   48,  0.0,   2.40); 
  Gen_Muon_eta4_phi = new TH1F("Gen_Muon_eta4_phi", "Gen_Muon_eta4_phi",   12, -3.14,  3.14); 
  Gen_Muon_eta4_lxy = new TH1F("Gen_Muon_eta4_lxy", "Gen_Muon_eta4_lxy",  160,  0.0,   8.00);
  Gen_Muon_eta4_lz  = new TH1F("Gen_Muon_eta4_lz",  "Gen_Muon_eta4_lz",   150,  0.0,  15.00);

  Sim_All_DiMuon_pt  = new TH1F("Sim_All_DiMuon_pt",  "Sim_All_DiMuon_pt",   100,  0.0, 500.00);
  Sim_All_DiMuon_eta = new TH1F("Sim_All_DiMuon_eta", "Sim_All_DiMuon_eta",   48,  0.0,   2.40);
  Sim_All_DiMuon_phi = new TH1F("Sim_All_DiMuon_phi", "Sim_All_DiMuon_phi",   12, -3.14,  3.14);
  Sim_All_DiMuon_lxy = new TH1F("Sim_All_DiMuon_lxy", "Sim_All_DiMuon_lxy",  160,  0.0,   8.00);
  Sim_All_DiMuon_lz  = new TH1F("Sim_All_DiMuon_lz",  "Sim_All_DiMuon_lz",   150,  0.0,  15.00);
  Sim_All_DiMuon_m   = new TH1F("Sim_All_DiMuon_m",   "Sim_All_DiMuon_m",      3,  0.3,   0.50);
  Sim_All_DiMuon_dR  = new TH1F("Sim_All_DiMuon_dR",  "Sim_All_DiMuon_m",    100,  0.0,   5.00);

  Sim_All_Muon_pt  = new TH1F("Sim_All_Muon_pt",  "Sim_All_Muon_pt",   100,  0.0, 500.00);
  Sim_All_Muon_eta = new TH1F("Sim_All_Muon_eta", "Sim_All_Muon_eta",   48,  0.0,   2.40); 
  Sim_All_Muon_phi = new TH1F("Sim_All_Muon_phi", "Sim_All_Muon_phi",   12, -3.14,  3.14); 
  Sim_All_Muon_lxy = new TH1F("Sim_All_Muon_lxy", "Sim_All_Muon_lxy",  160,  0.0,   8.00);
  Sim_All_Muon_lz  = new TH1F("Sim_All_Muon_lz",  "Sim_All_Muon_lz",   150,  0.0,  15.00);
  Sim_Muon_pt1_pt  = new TH1F("Sim_Muon_pt1_pt",  "Sim_Muon_pt1_pt",   100,  0.0, 500.00);
  Sim_Muon_pt1_eta = new TH1F("Sim_Muon_pt1_eta", "Sim_Muon_pt1_eta",   48,  0.0,   2.40); 
  Sim_Muon_pt1_phi = new TH1F("Sim_Muon_pt1_phi", "Sim_Muon_pt1_phi",   12, -3.14,  3.14); 
  Sim_Muon_pt1_lxy = new TH1F("Sim_Muon_pt1_lxy", "Sim_Muon_pt1_lxy",  160,  0.0,   8.00);
  Sim_Muon_pt1_lz  = new TH1F("Sim_Muon_pt1_lz",  "Sim_Muon_pt1_lz",   150,  0.0,  15.00);
  Sim_Muon_pt2_pt  = new TH1F("Sim_Muon_pt2_pt",  "Sim_Muon_pt2_pt",   100,  0.0, 500.00);
  Sim_Muon_pt2_eta = new TH1F("Sim_Muon_pt2_eta", "Sim_Muon_pt2_eta",   48,  0.0,   2.40); 
  Sim_Muon_pt2_phi = new TH1F("Sim_Muon_pt2_phi", "Sim_Muon_pt2_phi",   12, -3.14,  3.14); 
  Sim_Muon_pt2_lxy = new TH1F("Sim_Muon_pt2_lxy", "Sim_Muon_pt2_lxy",  160,  0.0,   8.00);
  Sim_Muon_pt2_lz  = new TH1F("Sim_Muon_pt2_lz",  "Sim_Muon_pt2_lz",   150,  0.0,  15.00);
  Sim_Muon_pt3_pt  = new TH1F("Sim_Muon_pt3_pt",  "Sim_Muon_pt3_pt",   100,  0.0, 500.00);
  Sim_Muon_pt3_eta = new TH1F("Sim_Muon_pt3_eta", "Sim_Muon_pt3_eta",   48,  0.0,   2.40); 
  Sim_Muon_pt3_phi = new TH1F("Sim_Muon_pt3_phi", "Sim_Muon_pt3_phi",   12, -3.14,  3.14); 
  Sim_Muon_pt3_lxy = new TH1F("Sim_Muon_pt3_lxy", "Sim_Muon_pt3_lxy",  160,  0.0,   8.00);
  Sim_Muon_pt3_lz  = new TH1F("Sim_Muon_pt3_lz",  "Sim_Muon_pt3_lz",   150,  0.0,  15.00);
  Sim_Muon_pt4_pt  = new TH1F("Sim_Muon_pt4_pt",  "Sim_Muon_pt4_pt",   100,  0.0, 500.00);
  Sim_Muon_pt4_eta = new TH1F("Sim_Muon_pt4_eta", "Sim_Muon_pt4_eta",   48,  0.0,   2.40); 
  Sim_Muon_pt4_phi = new TH1F("Sim_Muon_pt4_phi", "Sim_Muon_pt4_phi",   12, -3.14,  3.14); 
  Sim_Muon_pt4_lxy = new TH1F("Sim_Muon_pt4_lxy", "Sim_Muon_pt4_lxy",  160,  0.0,   8.00);
  Sim_Muon_pt4_lz  = new TH1F("Sim_Muon_pt4_lz",  "Sim_Muon_pt4_lz",   150,  0.0,  15.00);

  RSA_All_DiMuon_pt  = new TH1F("RSA_All_DiMuon_pt",  "RSA_All_DiMuon_pt",   100,  0.0, 500.00);
  RSA_All_DiMuon_eta = new TH1F("RSA_All_DiMuon_eta", "RSA_All_DiMuon_eta",   48,  0.0,   2.40);
  RSA_All_DiMuon_phi = new TH1F("RSA_All_DiMuon_phi", "RSA_All_DiMuon_phi",   12, -3.14,  3.14);
  RSA_All_DiMuon_lxy = new TH1F("RSA_All_DiMuon_lxy", "RSA_All_DiMuon_lxy",  160,  0.0,   8.00);
  RSA_All_DiMuon_lz  = new TH1F("RSA_All_DiMuon_lz",  "RSA_All_DiMuon_lz",   150,  0.0,  15.00);
  RSA_All_DiMuon_m   = new TH1F("RSA_All_DiMuon_m",   "RSA_All_DiMuon_m",      3,  0.3,   0.50);
  RSA_All_DiMuon_dR  = new TH1F("RSA_All_DiMuon_dR",  "RSA_All_DiMuon_m",    100,  0.0,   5.00);

  RSA_All_Muon_pt  = new TH1F("RSA_All_Muon_pt",  "RSA_All_Muon_pt",   100,  0.0, 500.00);
  RSA_All_Muon_eta = new TH1F("RSA_All_Muon_eta", "RSA_All_Muon_eta",   48,  0.0,   2.40); 
  RSA_All_Muon_phi = new TH1F("RSA_All_Muon_phi", "RSA_All_Muon_phi",   12, -3.14,  3.14); 
  RSA_All_Muon_lxy = new TH1F("RSA_All_Muon_lxy", "RSA_All_Muon_lxy",  160,  0.0,   8.00);
  RSA_All_Muon_lz  = new TH1F("RSA_All_Muon_lz",  "RSA_All_Muon_lz",   150,  0.0,  15.00);
  RSA_Muon_pt1_pt  = new TH1F("RSA_Muon_pt1_pt",  "RSA_Muon_pt1_pt",   100,  0.0, 500.00);
  RSA_Muon_pt1_eta = new TH1F("RSA_Muon_pt1_eta", "RSA_Muon_pt1_eta",   48,  0.0,   2.40); 
  RSA_Muon_pt1_phi = new TH1F("RSA_Muon_pt1_phi", "RSA_Muon_pt1_phi",   12, -3.14,  3.14); 
  RSA_Muon_pt1_lxy = new TH1F("RSA_Muon_pt1_lxy", "RSA_Muon_pt1_lxy",  160,  0.0,   8.00);
  RSA_Muon_pt1_lz  = new TH1F("RSA_Muon_pt1_lz",  "RSA_Muon_pt1_lz",   150,  0.0,  15.00);
  RSA_Muon_pt2_pt  = new TH1F("RSA_Muon_pt2_pt",  "RSA_Muon_pt2_pt",   100,  0.0, 500.00);
  RSA_Muon_pt2_eta = new TH1F("RSA_Muon_pt2_eta", "RSA_Muon_pt2_eta",   48,  0.0,   2.40); 
  RSA_Muon_pt2_phi = new TH1F("RSA_Muon_pt2_phi", "RSA_Muon_pt2_phi",   12, -3.14,  3.14); 
  RSA_Muon_pt2_lxy = new TH1F("RSA_Muon_pt2_lxy", "RSA_Muon_pt2_lxy",  160,  0.0,   8.00);
  RSA_Muon_pt2_lz  = new TH1F("RSA_Muon_pt2_lz",  "RSA_Muon_pt2_lz",   150,  0.0,  15.00);
  RSA_Muon_pt3_pt  = new TH1F("RSA_Muon_pt3_pt",  "RSA_Muon_pt3_pt",   100,  0.0, 500.00);
  RSA_Muon_pt3_eta = new TH1F("RSA_Muon_pt3_eta", "RSA_Muon_pt3_eta",   48,  0.0,   2.40); 
  RSA_Muon_pt3_phi = new TH1F("RSA_Muon_pt3_phi", "RSA_Muon_pt3_phi",   12, -3.14,  3.14); 
  RSA_Muon_pt3_lxy = new TH1F("RSA_Muon_pt3_lxy", "RSA_Muon_pt3_lxy",  160,  0.0,   8.00);
  RSA_Muon_pt3_lz  = new TH1F("RSA_Muon_pt3_lz",  "RSA_Muon_pt3_lz",   150,  0.0,  15.00);
  RSA_Muon_pt4_pt  = new TH1F("RSA_Muon_pt4_pt",  "RSA_Muon_pt4_pt",   100,  0.0, 500.00);
  RSA_Muon_pt4_eta = new TH1F("RSA_Muon_pt4_eta", "RSA_Muon_pt4_eta",   48,  0.0,   2.40); 
  RSA_Muon_pt4_phi = new TH1F("RSA_Muon_pt4_phi", "RSA_Muon_pt4_phi",   12, -3.14,  3.14); 
  RSA_Muon_pt4_lxy = new TH1F("RSA_Muon_pt4_lxy", "RSA_Muon_pt4_lxy",  160,  0.0,   8.00);
  RSA_Muon_pt4_lz  = new TH1F("RSA_Muon_pt4_lz",  "RSA_Muon_pt4_lz",   150,  0.0,  15.00);
}


MyDisplacedMuonAnalyzer::~MyDisplacedMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  // Counters
  std::cout<<"\n"<<std::endl;
  std::cout<<"========================================================"<<std::endl;
  std::cout<<"= Analyzed Gen Muons from Dark Photon Decay: "<<std::setw(9)<<tot_gen<<" ="<<std::endl;
  std::cout<<"= Analyzed Sim Muons from Dark Photon Decay: "<<std::setw(9)<<tot_sim<<" ="<<std::endl;
  std::cout<<"= Analyzed RSA Muons from Dark Photon Decay: "<<std::setw(9)<<tot_rsa<<" ="<<std::endl;
  std::cout<<"--------------------------------------------------------"<<std::endl;
  std::cout<<"=== Sim Matching Efficiency:                 "<<std::setw(9)<<tot_sim*1.0/tot_gen<<" ="<<std::endl;
  std::cout<<"=== RSA Efficiency (x) Matching Efficiency:  "<<std::setw(9)<<tot_rsa*1.0/tot_gen<<" ="<<std::endl;
  std::cout<<"========================================================"<<std::endl;
  std::cout<<"\n"<<std::endl;

  // Output File
  outputfile->cd();
  GenPhotons->cd();
  Gen_All_DarkPhot_pt->Write(); Gen_All_DarkPhot_eta->Write(); Gen_All_DarkPhot_phi->Write(); Gen_All_DarkPhot_lxy->Write(); Gen_All_DarkPhot_lz->Write(); Gen_All_DarkPhot_m->Write(); Gen_All_DarkPhot_mt->Write();
  Gen_DarkPhot_pt1_pt->Write(); Gen_DarkPhot_pt1_eta->Write(); Gen_DarkPhot_pt1_phi->Write(); Gen_DarkPhot_pt1_lxy->Write(); Gen_DarkPhot_pt1_lz->Write(); Gen_DarkPhot_pt1_m->Write(); Gen_DarkPhot_pt1_mt->Write();
  Gen_DarkPhot_pt2_pt->Write(); Gen_DarkPhot_pt2_eta->Write(); Gen_DarkPhot_pt2_phi->Write(); Gen_DarkPhot_pt2_lxy->Write(); Gen_DarkPhot_pt2_lz->Write(); Gen_DarkPhot_pt2_m->Write(); Gen_DarkPhot_pt2_mt->Write();
  GenMuons->cd();
  Gen_All_DiMuon_pt->Write();  Gen_All_DiMuon_eta->Write();  Gen_All_DiMuon_phi->Write();  Gen_All_DiMuon_lxy->Write();  Gen_All_DiMuon_lz->Write();  Gen_All_DiMuon_m->Write();  Gen_All_DiMuon_dR->Write();
  Gen_All_Muon_pt->Write(); Gen_All_Muon_eta->Write(); Gen_All_Muon_phi->Write(); Gen_All_Muon_lxy->Write(); Gen_All_Muon_lz->Write();
  Gen_Muon_pt1_pt->Write(); Gen_Muon_pt1_eta->Write(); Gen_Muon_pt1_phi->Write(); Gen_Muon_pt1_lxy->Write(); Gen_Muon_pt1_lz->Write();
  Gen_Muon_pt2_pt->Write(); Gen_Muon_pt2_eta->Write(); Gen_Muon_pt2_phi->Write(); Gen_Muon_pt2_lxy->Write(); Gen_Muon_pt2_lz->Write();
  Gen_Muon_pt3_pt->Write(); Gen_Muon_pt3_eta->Write(); Gen_Muon_pt3_phi->Write(); Gen_Muon_pt3_lxy->Write(); Gen_Muon_pt3_lz->Write();
  Gen_Muon_pt4_pt->Write(); Gen_Muon_pt4_eta->Write(); Gen_Muon_pt4_phi->Write(); Gen_Muon_pt4_lxy->Write(); Gen_Muon_pt4_lz->Write();
  Gen_Muon_eta1_pt->Write(); Gen_Muon_eta1_eta->Write(); Gen_Muon_eta1_phi->Write(); Gen_Muon_eta1_lxy->Write(); Gen_Muon_eta1_lz->Write();
  Gen_Muon_eta2_pt->Write(); Gen_Muon_eta2_eta->Write(); Gen_Muon_eta2_phi->Write(); Gen_Muon_eta2_lxy->Write(); Gen_Muon_eta2_lz->Write();
  Gen_Muon_eta3_pt->Write(); Gen_Muon_eta3_eta->Write(); Gen_Muon_eta3_phi->Write(); Gen_Muon_eta3_lxy->Write(); Gen_Muon_eta3_lz->Write();
  Gen_Muon_eta4_pt->Write(); Gen_Muon_eta4_eta->Write(); Gen_Muon_eta4_phi->Write(); Gen_Muon_eta4_lxy->Write(); Gen_Muon_eta4_lz->Write();
  SimMuons->cd();
  Sim_All_DiMuon_pt->Write();  Sim_All_DiMuon_eta->Write();  Sim_All_DiMuon_phi->Write();  Sim_All_DiMuon_lxy->Write();  Sim_All_DiMuon_lz->Write();  Sim_All_DiMuon_m->Write();  Sim_All_DiMuon_dR->Write();
  Sim_All_Muon_pt->Write(); Sim_All_Muon_eta->Write(); Sim_All_Muon_phi->Write(); Sim_All_Muon_lxy->Write(); Sim_All_Muon_lz->Write();
  Sim_Muon_pt1_pt->Write(); Sim_Muon_pt1_eta->Write(); Sim_Muon_pt1_phi->Write(); Sim_Muon_pt1_lxy->Write(); Sim_Muon_pt1_lz->Write();
  Sim_Muon_pt2_pt->Write(); Sim_Muon_pt2_eta->Write(); Sim_Muon_pt2_phi->Write(); Sim_Muon_pt2_lxy->Write(); Sim_Muon_pt2_lz->Write();
  Sim_Muon_pt3_pt->Write(); Sim_Muon_pt3_eta->Write(); Sim_Muon_pt3_phi->Write(); Sim_Muon_pt3_lxy->Write(); Sim_Muon_pt3_lz->Write();
  Sim_Muon_pt4_pt->Write(); Sim_Muon_pt4_eta->Write(); Sim_Muon_pt4_phi->Write(); Sim_Muon_pt4_lxy->Write(); Sim_Muon_pt4_lz->Write();
  RSAMuons->cd();
  RSA_All_DiMuon_pt->Write();  RSA_All_DiMuon_eta->Write();  RSA_All_DiMuon_phi->Write();  RSA_All_DiMuon_lxy->Write();  RSA_All_DiMuon_lz->Write();  RSA_All_DiMuon_m->Write();  RSA_All_DiMuon_dR->Write();
  RSA_All_Muon_pt->Write(); RSA_All_Muon_eta->Write(); RSA_All_Muon_phi->Write(); RSA_All_Muon_lxy->Write(); RSA_All_Muon_lz->Write();
  RSA_Muon_pt1_pt->Write(); RSA_Muon_pt1_eta->Write(); RSA_Muon_pt1_phi->Write(); RSA_Muon_pt1_lxy->Write(); RSA_Muon_pt1_lz->Write();
  RSA_Muon_pt2_pt->Write(); RSA_Muon_pt2_eta->Write(); RSA_Muon_pt2_phi->Write(); RSA_Muon_pt2_lxy->Write(); RSA_Muon_pt2_lz->Write();
  RSA_Muon_pt3_pt->Write(); RSA_Muon_pt3_eta->Write(); RSA_Muon_pt3_phi->Write(); RSA_Muon_pt3_lxy->Write(); RSA_Muon_pt3_lz->Write();
  RSA_Muon_pt4_pt->Write(); RSA_Muon_pt4_eta->Write(); RSA_Muon_pt4_phi->Write(); RSA_Muon_pt4_lxy->Write(); RSA_Muon_pt4_lz->Write();
  //Rechits_All->Write();
  outputfile->cd();
}


//
// member functions
//
/*
bool MyDisplacedMuonAnalyzer::myptsort(double * a, double * b)
{
  return (*a) > (*b); 
}
bool MyDisplacedMuonAnalyzer::myetasort(double * a, double * b)
{
  return (*(a+1)) > (*(b+1));
}
*/


// ------------ method called for each event  ------------
void
MyDisplacedMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // std::vector<reco::TrackCollection> staTracksCollection;
  edm::Handle<reco::TrackCollection> staTracks_SAMuon;
  iEvent.getByLabel(STAMuLabel1, staTracks_SAMuon);
  // staTracksCollection.push_back(staTracks_SAMuon); 
  if(techDebug) std::cout<<"Reconstructed  SA Muon Tracks:       " <<staTracks_SAMuon->size()<< std::endl;
  edm::Handle<reco::TrackCollection> staTracks_RSAMuon;
  iEvent.getByLabel(STAMuLabel2, staTracks_RSAMuon);
  // staTracksCollection.push_back(staTracks_RSAMuon); 
  if(techDebug) std::cout<<"Reconstructed RSA Muon Tracks:       " <<staTracks_RSAMuon->size()<< std::endl;
  edm::Handle<reco::TrackCollection> staTracks_SAMuonVtx;
  iEvent.getByLabel(STAMuLabel3, staTracks_SAMuonVtx);
  // staTracksCollection.push_back(staTracks_SAMuonVtx); 
  if(techDebug) std::cout<<"Reconstructed  SA Muon (Vtx) Tracks: " <<staTracks_SAMuonVtx->size()<< std::endl;
  edm::Handle<reco::TrackCollection> staTracks_RSAMuonVtx;
  iEvent.getByLabel(STAMuLabel4, staTracks_RSAMuonVtx);
  // staTracksCollection.push_back(staTracks_RSAMuonVtx); 
  if(techDebug) std::cout<<"Reconstructed RSA Muon (Vtx) Tracks: " <<staTracks_RSAMuonVtx->size()<< std::endl;

  // std::vector<reco::TrackCollection>::const_iterator staTrackCollection;
  reco::TrackCollection::const_iterator staTrack;

  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  const int n_prop = 15;
  const int n_prop_reco = 16;
  double dark_me_properties[2][n_prop];        // 2 Dark Photons         // Save for each phot PT, ETA, PHI, E, M, Q, Boolean Match                    vx, vy, vz, Lxy, Lz, Dark Photon Cand [0,1]  
  double muon_me_properties[4][n_prop];        // 4 Matrix Element Muons // Save for each muon PT, ETA, PHI, E, M, Q, Boolean Matched to sim track,    vx, vy, vz, Lxy, Lz, Dark Photon Cand [0,1]  Muon Numbering [0,1,2,3], dR(1,2)
  double muon_tk_properties[4][n_prop];        // 4 SimTrack Muons       // Save for each muon PT, ETA, PHI, E, M, Q, Boolean Matched to gen particle,                                              Muon Numbering [0,1,2,3], dR(1,2)
  // double muon_sa_properties[4][n_prop_rco];      // 4 SimTrack Muons  // Save for each muon PT, ETA, PHI, E, M, Q, Boolean Matched to gen particle,                                              Muon Numbering [0,1,2,3], dR(1,2)
  double muon_rsa_properties[4][n_prop_reco];       // 4 SimTrack Muons  // Save for each muon PT, ETA, PHI, E, M, Q, Boolean Matched to gen particle,                                              Muon Numbering [0,1,2,3], dR(1,2), dR(S,R)
  // double muon_savtx_properties[4][n_prop_reco];  // 4 SimTrack Muons  // Save for each muon PT, ETA, PHI, E, M, Q, Boolean Matched to gen particle,                                              Muon Numbering [0,1,2,3], dR(1,2). dR(S,R)
  // double muon_rsavtx_properties[4][n_prop_reco]; // 4 SimTrack Muons  // Save for each muon PT, ETA, PHI, E, M, Q, Boolean Matched to gen particle,                                              Muon Numbering [0,1,2,3], dR(1,2), dR(S,R)

  double dimu_me_properties[2][n_prop];    // 2 DiMuons
  double dimu_tk_properties[2][n_prop];    // 2 DiMuons
  double dimu_rsa_properties[2][n_prop];   // 2 DiMuons

  // double vimu_me_properties[1][n_prop_ext];    // 1 FourMuon             // Save for each Fourmu pair in addition Delta R. Delta Eta, Delta Phi, Invariant Mass

  for(int i=0; i<2; ++i) { for(int j=0; j<n_prop; ++j) { dark_me_properties[i][j]=-1; dimu_me_properties[i][j]=-1;} dark_me_properties[i][6]=false; dimu_me_properties[i][6]=false; }
  for(int i=0; i<4; ++i) { for(int j=0; j<n_prop; ++j) { muon_me_properties[i][j]=-1; muon_tk_properties[i][j]=-1;  muon_rsa_properties[i][j]=-1; }  muon_me_properties[i][6]=false; muon_tk_properties[i][6]=false; muon_rsa_properties[i][6]=false; }
  for(int i=0; i<2; ++i) { for(int j=0; j<n_prop; ++j) { dimu_tk_properties[i][j]=-1; }  dimu_tk_properties[i][6]=false; }
  for(int i=0; i<4; ++i) { for(int j=0; j<n_prop; ++j) { muon_rsa_properties[i][j]=-1; } muon_rsa_properties[i][6]=false; }
  for(int i=0; i<2; ++i) { for(int j=0; j<n_prop; ++j) { dimu_rsa_properties[i][j]=-1; } dimu_rsa_properties[i][6]=false; }



  std::vector< std::vector< double > > dark_me_prop_vec;
  std::vector< std::vector< double > > muon_me_prop_vec;
  std::vector< std::vector< double > > dimu_me_prop_vec;
  std::vector< std::vector< double > > muon_tk_prop_vec;
  // std::vector< std::vector< double > > dimu_tk_prop_vec;
  std::vector< std::vector< double > > muon_rsa_prop_vec;
  // std::vector< std::vector< double > > dimu_rsa_prop_vec;

  tot_gen_ev  = 0; tot_sim_ev = 0; tot_rsa_ev = 0;

  // ===================================================================
  // ===      MC Truth Matching                                      ===
  // ===================================================================  
  // Get the Muon SimTracks / Muon GenParticles here && match to H(yy)-decay
  // Do a pt match between SimTrack and Muon GenParticles from H(yy)-decay  
  // Then do dR match between SimTrack and STA Muon               
  if(mcTruthMatch) {
    // 1) GEN-LEVEL :: Save Matrix Element Muon Properties (PT, ETA, PHI)
    edm::Handle<reco::GenParticleCollection>      genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    bool HbosonFound = 0;
    if(genDebug)  std::cout<<"\n "<<std::endl;
    if(genDebug)  std::cout<<"=== Analysis of GEN Particles :: ============   Generated Particles: "<<genParticles->size()<<std::endl;
    if(genDebug)  std::cout<<"=============================================================================================================================================================================================="<<std::endl;
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if(HbosonFound) continue;
      if (g->pdgId() == 25) {
	// All Status 3 Particles [ = All Particles of the Matrix Element Calculation]
	if(genDebug) {
	  std::cout<<"Higgs Boson:                  id = "<<std::showpos<<std::setw(8)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(10)<<g->eta()<<" | phi = "<<std::setw(10)<<g->phi();
	  std::cout<<" | pt = "<<std::setw(10)<<g->pt()<<" GeV/c";  // " | et = "<<std::setw(10)<<g->et()<<" GeV | E = "<<g->energy()<<" GeV ";
	  std::cout<<" | mt = "<<std::setw(10)<<g->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<g->mass()<<" GeV/c^2 | Daughters: "<<(g->daughterRefVector()).size();
	  std::cout<<" | vtx = ["<<std::setw(10)<<g->vx()<<", "<<std::setw(10)<<g->vy()<<", "<<std::setw(10)<<g->vz()<<"]"<<std::endl;
	}
	HbosonFound = 1;
	reco::GenParticle::daughters d = g->daughterRefVector();
	int dark_me_count = 0;
	int muon_me_count = 0;
	for (reco::GenParticle::daughters::const_iterator it_d = d.begin(); it_d != d.end(); ++it_d) {
	  if((*it_d)->status()==3 && fabs((*it_d)->pdgId())>3000000) {
	    n = &(*(*it_d));
	    if(genDebug) {
	      std::cout<<"|--> Neutralino N1:           id = "<<std::showpos<<std::setw(8)<<n->pdgId()<<" | st = "<<std::setw(2)<<n->status()<<" | eta = "<<std::setw(10)<<n->eta()<<" | phi = "<<std::setw(10)<<n->phi();
	      std::cout<<" | pt = "<<std::setw(10)<<n->pt()<<" GeV/c"; // " | et = "<<std::setw(10)<<n->et()<<" GeV | E = "<<n->energy()<<" GeV ";
	      std::cout<<" | mt = "<<std::setw(10)<<n->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<n->mass()<<" GeV/c^2 | Daughters: "<<(n->daughterRefVector()).size();
	      std::cout<<" | vtx = ["<<std::setw(10)<<n->vx()<<", "<<std::setw(10)<<n->vy()<<", "<<std::setw(10)<<n->vz()<<"]"<<std::endl;
	    }
	    reco::GenParticle::daughters e = n->daughterRefVector();
	    for (reco::GenParticle::daughters::const_iterator it_e = e.begin(); it_e != e.end(); ++it_e) {
	      if(fabs((*it_e)->pdgId())==3000022) {
		p = &(*(*it_e));
		if(genDebug) {
		  std::cout<<"     |--> Dark Photon         id = "<<std::showpos<<std::setw(8)<<p->pdgId()<<" | st = "<<std::setw(2)<<p->status()<<" | eta = "<<std::setw(10)<<p->eta()<<" | phi = "<<std::setw(10)<<p->phi();
		  std::cout<<" | pt = "<<std::setw(10)<<p->pt()<<" GeV/c"; // " | et = "<<std::setw(10)<<p->et()<<" GeV | E = "<<p->energy()<<" GeV ";
		  std::cout<<" | mt = "<<std::setw(10)<<p->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<p->mass()<<" GeV/c^2 | Daughters: "<<(p->daughterRefVector()).size();
		  std::cout<<" | vtx = ["<<std::setw(10)<<p->vx()<<", "<<std::setw(10)<<p->vy()<<", "<<std::setw(10)<<p->vz()<<"]"<<std::endl;
		}
		dark_me_properties[dark_me_count][0] = p->pt();
		dark_me_properties[dark_me_count][1] = p->eta();
		dark_me_properties[dark_me_count][2] = p->phi();
		dark_me_properties[dark_me_count][3] = p->energy();
		dark_me_properties[dark_me_count][4] = p->mass();
		dark_me_properties[dark_me_count][5] = p->charge();
		dark_me_properties[dark_me_count][6] = false;
		dark_me_properties[dark_me_count][7] = p->vx();
		dark_me_properties[dark_me_count][8] = p->vy();
		dark_me_properties[dark_me_count][9] = p->vz();
		dark_me_properties[dark_me_count][12] = dark_me_count; // This is Dark Photon number X
		++dark_me_count;
		// math::XYZTLorentzVector v4 = p->LorentzVector();
	      }
	      else {
		q = &(*(*it_e));
		if(genDebug) {
		  std::cout<<"     |--> Neutralino N0       id = "<<std::showpos<<std::setw(8)<<q->pdgId()<<" | st = "<<std::setw(2)<<q->status()<<" | eta = "<<std::setw(10)<<q->eta()<<" | phi = "<<std::setw(10)<<q->phi();
		  std::cout<<" | pt = "<<std::setw(10)<<q->pt()<<" GeV/c"; // " | et = "<<std::setw(10)<<q->et()<<" GeV | E = "<<q->energy()<<" GeV ";
		  std::cout<<" | mt = "<<std::setw(10)<<q->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<q->mass()<<" GeV/c^2 | Daughters: "<<(q->daughterRefVector()).size();
		  std::cout<<" | vtx = ["<<std::setw(10)<<q->vx()<<", "<<std::setw(10)<<q->vy()<<", "<<std::setw(10)<<q->vz()<<"]"<<std::endl;
		}
	      }
	      reco::GenParticle::daughters f = (*it_e)->daughterRefVector();
	      for (reco::GenParticle::daughters::const_iterator it_f = f.begin(); it_f != f.end(); ++it_f) {
		if((*it_f)->status()==3 && fabs((*it_f)->pdgId())==13) {
		  m = &(*(*it_f));
		  if(genDebug) std::cout<<"          |--> Muon:          id = "<<std::showpos<<std::setw(8)<<m->pdgId()<<" | st = "<<std::setw(2)<<m->status()<<" | eta = "<<std::setw(10)<<m->eta()<<" | phi = "<<std::setw(10)<<m->phi();
		  if(genDebug) std::cout<<" | pt = "<<std::setw(10)<<m->pt()<<" GeV/c"; // " | et = "<<std::setw(10)<<m->et()<<" GeV | E = "<<m->energy()<<" GeV ";
		  if(genDebug) std::cout<<" | mt = "<<std::setw(10)<<m->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<m->mass()<<" GeV/c^2 | Daughters: "<<(m->daughterRefVector()).size();
		  if(genDebug) std::cout<<" | vtx = ["<<std::setw(10)<<m->vx()<<", "<<std::setw(10)<<m->vy()<<", "<<std::setw(10)<<m->vz()<<"]"<<std::endl;
		  muon_me_properties[muon_me_count][0] = m->pt();
		  muon_me_properties[muon_me_count][1] = m->eta();
		  muon_me_properties[muon_me_count][2] = m->phi();
		  muon_me_properties[muon_me_count][3] = m->energy();
		  muon_me_properties[muon_me_count][4] = m->mass();
		  muon_me_properties[muon_me_count][5] = m->charge();
		  muon_me_properties[muon_me_count][6] = false;
		  muon_me_properties[muon_me_count][7] = m->vx();
		  muon_me_properties[muon_me_count][8] = m->vy();
		  muon_me_properties[muon_me_count][9] = m->vz();
		  muon_me_properties[muon_me_count][10] = sqrt(pow(p->vx()-m->vx(),2)+pow(p->vy()-m->vy(),2));
		  muon_me_properties[muon_me_count][11] = sqrt(pow(p->vz()-m->vz(),2));
		  dark_me_properties[dark_me_count-1][10] = sqrt(pow(p->vx()-m->vx(),2)+pow(p->vy()-m->vy(),2)); // Fill the Dark Photon Array with this information
		  dark_me_properties[dark_me_count-1][11] = sqrt(pow(p->vz()-m->vz(),2));                        // Fill the Dark Photon Array with this information
		  muon_me_properties[muon_me_count][12] = dark_me_count-1; // To which Dark Photon belongs this muon
		  muon_me_properties[muon_me_count][13] = muon_me_count;   // This is muon number X
		  ++muon_me_count; 
		  // count events within acceptance
		  if(fabs(m->eta()) < 2.4) ++tot_gen; ++tot_gen_ev;
		  // math::XYZTLorentzVector v4 = m->LorentzVector();
		}
		else if(fabs((*it_f)->pdgId())>3000000) {
		  r = &(*(*it_f));
		  if(fabs(r->pdgId())==3000022) {
		    if(genDebug) std::cout<<"          |--> Dark Photon:   id = "<<std::showpos<<std::setw(7)<<r->pdgId()<<" | st = "<<std::setw(2)<<r->status()<<" | eta = "<<std::setw(10)<<r->eta()<<" | phi = "<<std::setw(10)<<r->phi();
		    if(genDebug) std::cout<<" | pt = "<<std::setw(10)<<r->pt()<<" GeV/c"; // " | et = "<<std::setw(10)<<r->et()<<" GeV | E = "<<r->energy()<<" GeV ";
		    if(genDebug) std::cout<<" | mt = "<<std::setw(10)<<r->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<r->mass()<<" GeV/c^2 | Daughters: "<<(r->daughterRefVector()).size();
		    if(genDebug) std::cout<<" | vtx = ["<<std::setw(10)<<r->vx()<<", "<<std::setw(10)<<r->vy()<<", "<<std::setw(10)<<r->vz()<<"]"<<std::endl;
		  }
		  else {
		    if(genDebug) std::cout<<"          |--> Neutralino N0  id = "<<std::showpos<<std::setw(8)<<r->pdgId()<<" | st = "<<std::setw(2)<<r->status()<<" | eta = "<<std::setw(10)<<r->eta()<<" | phi = "<<std::setw(10)<<r->phi();
		    if(genDebug) std::cout<<" | pt = "<<std::setw(10)<<r->pt()<<" GeV/c"; // " | et = "<<std::setw(10)<<r->et()<<" GeV | E = "<<r->energy()<<" GeV ";
		    if(genDebug) std::cout<<" | mt = "<<std::setw(10)<<r->mt()<<" GeV/c^2 | m = "<<std::setw(10)<<r->mass()<<" GeV/c^2 | Daughters: "<<(r->daughterRefVector()).size();
		    if(genDebug) std::cout<<" | vtx = ["<<std::setw(10)<<r->vx()<<", "<<std::setw(10)<<r->vy()<<", "<<std::setw(10)<<r->vz()<<"]"<<std::endl;
		  }
		}
		else {}
	      }
	    }
	  }
	}
      }
    }
    // Muon Pair Properties
    int dimu_me_count = 0;
    for(int i=0; i<4; ++i) {
      for(int j=0; j<i; ++j) {
	// std::cout<<"Dimuon: comparing muons [i,j] = ["<<i<<","<<j<<"] :: ";
	if(dimu_me_count==2) continue;
	if(muon_me_properties[i][12]==muon_me_properties[j][12] && muon_me_properties[i][12]!=dimu_me_properties[dimu_me_count][12] && !dimu_me_properties[dimu_me_count][6]) {
	  // std::cout<<"inside if :: DP of Muon i = "<<muon_me_properties[i][12]<<" == DP of Muon j = "<<muon_me_properties[j][12]<<" DP of DiMu = "<<dimu_me_properties[dimu_me_count][12]<<" Matched = "<<dimu_me_properties[dimu_me_count][6]<<std::endl;
	  TLorentzVector Mu1, Mu2, DiMu;
          // Mu1.SetPtEtaPhiM
          // Mu1.SetPxPyPzE
	  Mu1.SetPtEtaPhiE(muon_me_properties[i][0],muon_me_properties[i][1],muon_me_properties[i][2],muon_me_properties[i][3]);
	  Mu2.SetPtEtaPhiE(muon_me_properties[j][0],muon_me_properties[j][1],muon_me_properties[j][2],muon_me_properties[j][3]);
	  DiMu = Mu1+Mu2;
	  dimu_me_properties[dimu_me_count][0] = DiMu.Pt();
	  dimu_me_properties[dimu_me_count][1] = DiMu.Eta();
	  dimu_me_properties[dimu_me_count][2] = DiMu.Phi();
	  dimu_me_properties[dimu_me_count][3] = DiMu.E();
	  dimu_me_properties[dimu_me_count][4] = DiMu.M();
	  dimu_me_properties[dimu_me_count][5] = muon_me_properties[i][5]+muon_me_properties[j][5];
	  dimu_me_properties[dimu_me_count][6] = true;
	  for(int k=7; k<n_prop-1; ++k) { dimu_me_properties[dimu_me_count][k] = muon_me_properties[i][k]; }
	  // number
	  dimu_me_properties[dimu_me_count][13] = dimu_me_count;
	  // Additional Info
	  dimu_me_properties[dimu_me_count][14] = deltaR(muon_me_properties[i][1],muon_me_properties[i][2], muon_me_properties[j][1], muon_me_properties[j][2]);
	  muon_me_properties[i][14] = dimu_me_properties[dimu_me_count][14];
	  muon_me_properties[j][14] = dimu_me_properties[dimu_me_count][14]; 
	  // std::cout<<"saving delta R :: diMuon ["<<dimu_me_count<<"] dR = "<<dimu_me_properties[dimu_me_count][14]<<" --> saving to Muon ["<<i<<"] = "<<muon_me_properties[i][14]<<" and Muon ["<<j<<"] = "<<muon_me_properties[j][14]<<std::endl;
	  // Counter
	  ++dimu_me_count;
	}
	else {} // std::cout<<""<<std::endl; }
      }
    }
    if(genDebug)  {
      // std::cout<<"\n"<<std::endl;
      std::cout<<"--- Dark Photons found: ----------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      for(int i=0; i<2; ++i) {
	std::cout<<"--- Gamm "<<i<<" pt = "<<std::setw(11)<<dark_me_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<dark_me_properties[i][1]<<" | phi = "<<std::setw(10)<<dark_me_properties[i][2];
	// std::cout<<" E = "<<std::setw(11)<<dark_me_properties[i][3]<<" GeV | M = "<<std::setw(11)<<dark_me_properties[i][4]<<" GeV/c2";
	std::cout<<" | q = "<<dark_me_properties[i][5];
	std::cout<<" | vtx = ["<<std::setw(10)<<dark_me_properties[i][7]<<", "<<std::setw(10)<<dark_me_properties[i][8]<<", "<<std::setw(10)<<dark_me_properties[i][9]<<"]";
	std::cout<<" | Lxy = "<<std::setw(10)<<dark_me_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<dark_me_properties[i][11]<<" cm";
	std::cout<<" | Dark Photon no "<<dark_me_properties[i][12]<<std::endl;
      }
      std::cout<<"--- Muons found: -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      for(int i=0; i<4; ++i) {
	std::cout<<"--- Muon "<<i<<" pt = "<<std::setw(11)<<muon_me_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<muon_me_properties[i][1];
	std::cout<<" | phi = "<<std::setw(10)<<muon_me_properties[i][2]<<" | q = "<<muon_me_properties[i][5];
	std::cout<<" | vtx = ["<<std::setw(10)<<muon_me_properties[i][7]<<", "<<std::setw(10)<<muon_me_properties[i][8]<<", "<<std::setw(10)<<muon_me_properties[i][9]<<"]";
	std::cout<<" | Lxy = "<<std::setw(10)<<muon_me_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<muon_me_properties[i][11]<<" cm";
	std::cout<<" | Dark Photon no "<<muon_me_properties[i][12]<<" | Muon no "<<muon_me_properties[i][13]<<" | Delta R = "<<std::setw(10)<<muon_me_properties[i][14]<<std::endl;
      }
      std::cout<<"--- Muon Pair Properties: --------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      for(int i=0; i<2; ++i) {
	std::cout<<"--- DiMu "<<i<<" pt = "<<std::setw(11)<<dimu_me_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<dimu_me_properties[i][1];
	std::cout<<" | phi = "<<std::setw(10)<<dimu_me_properties[i][2]<<" | q = "<<dimu_me_properties[i][5];
	std::cout<<" | vtx = ["<<std::setw(10)<<dimu_me_properties[i][7]<<", "<<std::setw(10)<<dimu_me_properties[i][8]<<", "<<std::setw(10)<<dimu_me_properties[i][9]<<"]";
	std::cout<<" | Lxy = "<<std::setw(10)<<dimu_me_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<dimu_me_properties[i][11]<<" cm";
	std::cout<<" | Dark Photon no "<<dimu_me_properties[i][12]<<" | Dimu no "<<dimu_me_properties[i][13]<<" | Delta R = "<<std::setw(10)<<dimu_me_properties[i][14]<<std::endl;
      }
      std::cout<<"=============================================================================================================================================================================================="<<std::endl;
      std::cout<<"\n"<<std::endl;
    }


    // Change from Arrays to Vector Format ...
    for(int i=0; i<2; ++i) {
      std::vector<double> dark_me_prop_temp;
      for(int j=0; j<n_prop; ++j) { dark_me_prop_temp.push_back(dark_me_properties[i][j]); }
      dark_me_prop_vec.push_back(dark_me_prop_temp);
    }
    for(int i=0; i<4; ++i) {
      std::vector<double> muon_me_prop_temp;
      for(int j=0; j<n_prop; ++j) { muon_me_prop_temp.push_back(muon_me_properties[i][j]); }
      muon_me_prop_vec.push_back(muon_me_prop_temp);
    }
    // Sort the Vectors here based on pt (first element in vector) ...
    std::sort(dark_me_prop_vec.begin(), dark_me_prop_vec.end(),
	      [](const std::vector<double>& a, const std::vector<double>& b) {
		return a[0] > b[0];
	      });
    std::sort(muon_me_prop_vec.begin(), muon_me_prop_vec.end(),
	      [](const std::vector<double>& a, const std::vector<double>& b) {
		return a[0] > b[0];
	      });
    // Little Test to see whether PT sorting worked
    if(techDebug)  {
      std::cout<<"--- Muons PT sorted --------------------------"<<std::endl;
      for (std::vector< std::vector<double> >::const_iterator iMuon = muon_me_prop_vec.begin(); iMuon != muon_me_prop_vec.end(); ++iMuon) {
	std::cout<<"--- Muon :: pt = "<<std::setw(10)<<(*iMuon)[0]<<" eta = "<<std::setw(10)<<(*iMuon)[1]<<" phi = "<<std::setw(10)<<(*iMuon)[2]<<" matched to Dark Photon no "<<(*iMuon)[12]<<" Muon no "<<(*iMuon)[13]<<std::endl;
      }
      std::cout<<"----------------------------------------------"<<std::endl;
      std::cout<<"\n"<<std::endl;
    }
    // Fill Histograms here ...
    // First Dark Photon Histograms
    for(int i=0; i<2; ++i) {
      Gen_All_DarkPhot_pt->Fill(dark_me_prop_vec[i][0]);   Gen_All_DarkPhot_eta->Fill(fabs(dark_me_prop_vec[i][1]));  Gen_All_DarkPhot_phi->Fill(dark_me_prop_vec[i][2]); 
      Gen_All_DarkPhot_lxy->Fill(dark_me_prop_vec[i][10]); Gen_All_DarkPhot_lz->Fill(dark_me_prop_vec[i][11]);        Gen_All_DarkPhot_m->Fill(dark_me_prop_vec[i][3]);  Gen_All_DarkPhot_mt->Fill(dark_me_prop_vec[i][4]);
    }
    Gen_DarkPhot_pt1_pt->Fill(dark_me_prop_vec[0][0]);   Gen_DarkPhot_pt1_eta->Fill(fabs(dark_me_prop_vec[0][1]));  Gen_DarkPhot_pt1_phi->Fill(dark_me_prop_vec[0][2]); 
    Gen_DarkPhot_pt1_lxy->Fill(dark_me_prop_vec[0][10]); Gen_DarkPhot_pt1_lz->Fill(dark_me_prop_vec[0][11]);        Gen_DarkPhot_pt1_m->Fill(dark_me_prop_vec[0][3]);    Gen_DarkPhot_pt1_mt->Fill(dark_me_prop_vec[0][4]);
    Gen_DarkPhot_pt2_pt->Fill(dark_me_prop_vec[1][0]);   Gen_DarkPhot_pt2_eta->Fill(fabs(dark_me_prop_vec[1][1]));  Gen_DarkPhot_pt2_phi->Fill(dark_me_prop_vec[1][2]); 
    Gen_DarkPhot_pt2_lxy->Fill(dark_me_prop_vec[1][10]); Gen_DarkPhot_pt2_lz->Fill(dark_me_prop_vec[1][11]);        Gen_DarkPhot_pt2_m->Fill(dark_me_prop_vec[1][3]);    Gen_DarkPhot_pt2_mt->Fill(dark_me_prop_vec[1][4]);
    // Then DiMuon Histograms
    for(int i=0; i<2; ++i) {
      Gen_All_DiMuon_pt->Fill(dimu_me_properties[i][0]);   Gen_All_DiMuon_eta->Fill(dimu_me_properties[i][1]);  Gen_All_DiMuon_phi->Fill(dimu_me_properties[i][2]);  
      Gen_All_DiMuon_lxy->Fill(dimu_me_properties[i][10]);  Gen_All_DiMuon_lz->Fill(dimu_me_properties[i][11]);  Gen_All_DiMuon_m->Fill(dimu_me_properties[i][3]);  Gen_All_DiMuon_dR->Fill(dimu_me_properties[i][14]);
    }
    // Then Muon Histograms 
    for(int i=0; i<4; ++i) {
      Gen_All_Muon_pt->Fill(muon_me_prop_vec[i][0]);   Gen_All_Muon_eta->Fill(fabs(muon_me_prop_vec[i][1]));       Gen_All_Muon_phi->Fill(muon_me_prop_vec[i][2]);
      Gen_All_Muon_lxy->Fill(muon_me_prop_vec[i][10]); Gen_All_Muon_lz->Fill(muon_me_prop_vec[i][11]); 
    }
    Gen_Muon_pt1_pt->Fill(muon_me_prop_vec[0][0]);   Gen_Muon_pt1_eta->Fill(fabs(muon_me_prop_vec[0][1])); Gen_Muon_pt1_phi->Fill(muon_me_prop_vec[0][2]); 
    Gen_Muon_pt1_lxy->Fill(muon_me_prop_vec[0][10]); Gen_Muon_pt1_lz->Fill(muon_me_prop_vec[0][11]); 
    Gen_Muon_pt2_pt->Fill(muon_me_prop_vec[1][0]);   Gen_Muon_pt2_eta->Fill(fabs(muon_me_prop_vec[1][1])); Gen_Muon_pt2_phi->Fill(muon_me_prop_vec[1][2]); 
    Gen_Muon_pt2_lxy->Fill(muon_me_prop_vec[1][10]); Gen_Muon_pt2_lz->Fill(muon_me_prop_vec[1][11]);
    Gen_Muon_pt3_pt->Fill(muon_me_prop_vec[2][0]);   Gen_Muon_pt3_eta->Fill(fabs(muon_me_prop_vec[2][1])); Gen_Muon_pt3_phi->Fill(muon_me_prop_vec[2][2]); 
    Gen_Muon_pt3_lxy->Fill(muon_me_prop_vec[2][10]); Gen_Muon_pt3_lz->Fill(muon_me_prop_vec[2][11]); 
    Gen_Muon_pt4_pt->Fill(muon_me_prop_vec[3][0]);   Gen_Muon_pt4_eta->Fill(fabs(muon_me_prop_vec[3][1])); Gen_Muon_pt4_phi->Fill(muon_me_prop_vec[3][2]); 
    Gen_Muon_pt4_lxy->Fill(muon_me_prop_vec[3][10]); Gen_Muon_pt4_lz->Fill(muon_me_prop_vec[3][11]);      

    // Sort the Vectors here based on eta (second element in vector) ...
    // std::cout<<"Sorting based on ETA"<<std::endl;
    std::sort(muon_me_prop_vec.begin(), muon_me_prop_vec.end(),
	      [](const std::vector<double>& a, const std::vector<double>& b) {
		return fabs(a[1]) < fabs(b[1]);
	      });
    // Little test to see whether ETA sorting worked ...
    if(techDebug)  {
      std::cout<<"--- Muons ABS ETA sorted ---------------------"<<std::endl;
      for (std::vector< std::vector<double> >::const_iterator iMuon = muon_me_prop_vec.begin(); iMuon != muon_me_prop_vec.end(); ++iMuon) {
	std::cout<<"--- Muon :: pt = "<<std::setw(10)<<(*iMuon)[0]<<" eta = "<<std::setw(10)<<(*iMuon)[1]<<" phi = "<<std::setw(10)<<(*iMuon)[2]<<" matched to Dark Photon no "<<(*iMuon)[12]<<" Muon no "<<(*iMuon)[13]<<std::endl;
      }
      std::cout<<"----------------------------------------------"<<std::endl;
      std::cout<<"\n"<<std::endl;
    } 
    // Now we are ready to fill the histograms   
    Gen_Muon_eta1_pt->Fill(muon_me_prop_vec[0][0]);   Gen_Muon_eta1_eta->Fill(fabs(muon_me_prop_vec[0][1])); Gen_Muon_eta1_phi->Fill(muon_me_prop_vec[0][2]); 
    Gen_Muon_eta1_lxy->Fill(muon_me_prop_vec[0][10]); Gen_Muon_eta1_lz->Fill(muon_me_prop_vec[0][11]); 
    Gen_Muon_eta2_pt->Fill(muon_me_prop_vec[1][0]);   Gen_Muon_eta2_eta->Fill(fabs(muon_me_prop_vec[1][1])); Gen_Muon_eta2_phi->Fill(muon_me_prop_vec[1][2]); 
    Gen_Muon_eta2_lxy->Fill(muon_me_prop_vec[1][10]); Gen_Muon_eta2_lz->Fill(muon_me_prop_vec[1][11]); 
    Gen_Muon_eta3_pt->Fill(muon_me_prop_vec[2][0]);   Gen_Muon_eta3_eta->Fill(fabs(muon_me_prop_vec[2][1])); Gen_Muon_eta3_phi->Fill(muon_me_prop_vec[2][2]); 
    Gen_Muon_eta3_lxy->Fill(muon_me_prop_vec[2][10]); Gen_Muon_eta3_lz->Fill(muon_me_prop_vec[2][11]); 
    Gen_Muon_eta4_pt->Fill(muon_me_prop_vec[3][0]);   Gen_Muon_eta4_eta->Fill(fabs(muon_me_prop_vec[3][1])); Gen_Muon_eta4_phi->Fill(muon_me_prop_vec[3][2]); 
    Gen_Muon_eta4_lxy->Fill(muon_me_prop_vec[3][10]); Gen_Muon_eta4_lz->Fill(muon_me_prop_vec[3][11]);      

    // 2) SIM-LEVEL :: Are there SimTracks related to the long living dark photons?
    std::vector<SimTrack> theSimTracks;
    edm::Handle<edm::SimTrackContainer> SimTk;
    iEvent.getByLabel("g4SimHits",SimTk);
    theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
    if(genDebug)  {
      // std::cout<<"\n"<<std::endl;
      std::cout<<"=== Analysis of Sim Tracks    :: ============   Simulated Tracks: "<<SimTk->size()<<std::endl;
      std::cout<<"=============================================================================================================================================================================================="<<std::endl;
    }
    // for (std::vector<SimTrack>::const_iterator iTrack = theSimTracks.begin(); iTrack != theSimTracks.end(); ++iTrack) {
    // SimTrack simtrack = (*iTrack);
    // Print First 25 SimTracks in Collection
    // if(simtrack.trackId() < 25) {
    // if(genDebug) {
    // std::cout<<"SimTrack Found: track = "<<std::setw(2)<<simtrack.trackId()<<" id = "<<std::setw(4)<<simtrack.type()<<" pt = "<<std::setw(10)<<simtrack.momentum().pt();
    // std::cout<<" eta = "<<std::setw(10)<<simtrack.momentum().eta()<<" phi = "<<std::setw(10)<<simtrack.momentum().phi()<<" vtx = "<<std::setw(3)<<simtrack.vertIndex()<<std::endl;
    // }
    // }
    // else continue;
    // }
    // if(genDebug)  std::cout<<"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    // 3) SIM-LEVEL :: Match Muon SimTracks to GenParticle Muons
    int muon_tk_count = 0;
    // for(int i=0; i<4; ++i) { muon_tk_properties[i][6] = false; }
    for (std::vector<SimTrack>::const_iterator iTrack = theSimTracks.begin(); iTrack != theSimTracks.end(); ++iTrack) {
      SimTrack simtrack = (*iTrack);
      if(fabs(simtrack.type()) != 13) continue;
      else {
        if(genDebug) {
	  std::cout<<"Muon SimTrack Found: track = "<<std::setw(6)<<simtrack.trackId()<<" id = "<<std::setw(4)<<simtrack.type()<<" pt = "<<std::setw(11)<<simtrack.momentum().pt();
	  std::cout<<" eta = "<<std::setw(11)<<simtrack.momentum().eta()<<" phi = "<<std::setw(11)<<simtrack.momentum().phi()<<" sim vtx = "<<std::setw(10)<<simtrack.vertIndex()<<" || ";
	}
        double simtrack_eta = simtrack.momentum().eta();
        double simtrack_phi = simtrack.momentum().phi();
	double delta_rel_pt_mu[] = {-10, -10, -10, -10};
	// Loop over 4 (or less) ME muon gen particles
	for(int i=0; i<4; ++i) {
	  if(muon_me_properties[i][6]) {
	    if(i==3) {if(genDebug) std::cout<<" ---> NOT MATCHED"<<std::endl;}
	    continue; // if GenParticle is already matched, stop the loop
	  }
	  // Calculate Relative Delta Pt between this simtrack and the ith saved Gen Particle
	  delta_rel_pt_mu[i] = fabs(simtrack.momentum().pt()-muon_me_properties[i][0])/muon_me_properties[i][0];
	  // std::cout<<"Delta Pt between this Sim Track and Gen Particle "<<i<<" = "<<delta_rel_pt_mu[i]<<std::endl;
	  // Try to match the Sim Track and the Gen Particle based on Delta Pt
	  if(!muon_me_properties[i][6] && delta_rel_pt_mu[i] > 0 && delta_rel_pt_mu[i]<delta_rel_pt) {
	    muon_me_properties[i][6] = true;
	    muon_tk_properties[muon_tk_count][0] = simtrack.momentum().pt();
	    muon_tk_properties[muon_tk_count][1] = simtrack.momentum().eta();
	    muon_tk_properties[muon_tk_count][2] = simtrack.momentum().phi();
	    muon_tk_properties[muon_tk_count][3] = simtrack.momentum().E();
	    muon_tk_properties[muon_tk_count][4] = simtrack.momentum().M();
	    muon_tk_properties[muon_tk_count][5] = simtrack.charge();
	    muon_tk_properties[muon_tk_count][6] = true;
	    if(genDebug) std::cout<<" ---> MATCHED (based on Delta PT)"<<std::endl;
	    // Matched to Muon i ==> write properties also in SimTrack array
            for(int j=7; j<n_prop; ++j) { muon_tk_properties[muon_tk_count][j] = muon_me_properties[i][j]; } // i or muon_tk_count ???
	  }
	  // I found some cases of severe bremsstrahlung leading to Pt mismatch ==> Try Delta R matching for these                                    
	  else if(!muon_me_properties[i][6] && deltaR(simtrack_eta, simtrack_phi, muon_me_properties[i][1], muon_me_properties[i][2]) < delta_R_gensim) {
	    muon_me_properties[i][6] = true;
	    muon_tk_properties[muon_tk_count][0] = simtrack.momentum().pt();
	    muon_tk_properties[muon_tk_count][1] = simtrack.momentum().eta();
	    muon_tk_properties[muon_tk_count][2] = simtrack.momentum().phi();
	    muon_tk_properties[muon_tk_count][3] = simtrack.momentum().E();
	    muon_tk_properties[muon_tk_count][4] = simtrack.momentum().M();
	    muon_tk_properties[muon_tk_count][5] = simtrack.charge();
	    muon_tk_properties[muon_tk_count][6] = true;
	    if(genDebug) std::cout<<" ---> MATCHED (based on Delta R)"<<std::endl;
	    // Matched to Muon i ==> write properties also in SimTrack array
            for(int j=7; j<n_prop; ++j) { muon_tk_properties[muon_tk_count][j] = muon_me_properties[i][j]; } // i or muon_tk_count ???
	  }
	  else {
	    if(genDebug) std::cout<<" ---> NOT MATCHED"<<std::endl;
	  }
	  // If Sim Track is matched, stop the loop
	  if(muon_tk_properties[muon_tk_count][6] == true) {
	    ++muon_tk_count;
	    if(fabs(simtrack.momentum().eta()) < 2.4) ++tot_sim; ++tot_sim_ev;
	    break;
	  }
	}
      }
    }
    // Muon Pair Properties
    int dimu_tk_count = 0;
    for(int i=0; i<4; ++i) {
      for(int j=0; j<i; ++j) {
	// std::cout<<"Dimuon: comparing muons [i,j] = ["<<i<<","<<j<<"] :: "<<std::endl;
	if(dimu_tk_count==2) continue;
	// check whether Dark Photon Candidate (saved in [12]) is same;
	// check whether the Dark Photon Candidate (saved in [12]) is not yet filled for the dimuon tk candidate
	// check whether the dimuon tk candidate is not yet matched
	if(muon_tk_properties[i][12]==muon_tk_properties[j][12] && muon_tk_properties[i][12]!=dimu_tk_properties[dimu_tk_count][12] && !dimu_tk_properties[dimu_tk_count][6]) {
          // std::cout<<"inside if :: DP of Muon i = "<<muon_me_properties[i][12]<<" == DP of Muon j = "<<muon_me_properties[j][12]<<" DP of DiMu = "<<dimu_me_properties[dimu_me_count][12]<<" Matched = "<<dimu_me_properties[dimu_me_count][6]<<std::endl;
	  TLorentzVector Mu1, Mu2, DiMu;
          // Mu1.SetPtEtaPhiM
          // Mu1.SetPxPyPzE
	  Mu1.SetPtEtaPhiE(muon_tk_properties[i][0],muon_tk_properties[i][1],muon_tk_properties[i][2],muon_tk_properties[i][3]);
	  Mu2.SetPtEtaPhiE(muon_tk_properties[j][0],muon_tk_properties[j][1],muon_tk_properties[j][2],muon_tk_properties[j][3]);
	  DiMu = Mu1+Mu2;
	  dimu_tk_properties[dimu_tk_count][0] = DiMu.Pt();
	  dimu_tk_properties[dimu_tk_count][1] = DiMu.Eta();
	  dimu_tk_properties[dimu_tk_count][2] = DiMu.Phi();
	  dimu_tk_properties[dimu_tk_count][3] = DiMu.E();
	  dimu_tk_properties[dimu_tk_count][4] = DiMu.M();
	  dimu_tk_properties[dimu_tk_count][5] = muon_tk_properties[i][5]+muon_tk_properties[j][5];
	  dimu_tk_properties[dimu_tk_count][6] = true;
	  for(int k=7; k<n_prop-1; ++k) { dimu_tk_properties[dimu_tk_count][k] = muon_tk_properties[i][k]; }
	  // number
	  dimu_tk_properties[dimu_tk_count][13] = dimu_tk_count;
	  // Additional Info
	  dimu_tk_properties[dimu_tk_count][14] = deltaR(muon_tk_properties[i][1],muon_tk_properties[i][2], muon_tk_properties[j][1], muon_tk_properties[j][2]);
	  muon_tk_properties[i][14] = dimu_tk_properties[dimu_tk_count][14];
	  muon_tk_properties[j][14] = dimu_tk_properties[dimu_tk_count][14]; 
	  // std::cout<<"saving delta R :: diMuon ["<<dimu_tk_count<<"] dR = "<<dimu_tk_properties[dimu_tk_count][14]<<" --> saving to Muon ["<<i<<"] = "<<muon_tk_properties[i][14]<<" and Muon ["<<j<<"] = "<<muon_tk_properties[j][14]<<std::endl;
	  // Counter
	  ++dimu_tk_count;
	}
	else {} // std::cout<<""<<std::endl; }
      }
    }
    if(genDebug) {
      // std::cout<<"\n"<<std::endl; 
      std::cout<<"--- Muons found: -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      for(int i=0; i<4; ++i) {
	std::cout<<"--- Muon "<<i<<" pt = "<<std::setw(11)<<muon_tk_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<muon_tk_properties[i][1];
	std::cout<<" | phi = "<<std::setw(10)<<muon_tk_properties[i][2]<<" | q = "<<muon_tk_properties[i][5];
        std::cout<<" | vtx = ["<<std::setw(10)<<muon_tk_properties[i][7]<<", "<<std::setw(10)<<muon_tk_properties[i][8]<<", "<<std::setw(10)<<muon_tk_properties[i][9]<<"]";
	std::cout<<" | Lxy = "<<std::setw(10)<<muon_tk_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<muon_tk_properties[i][11]<<" cm";
	std::cout<<" | Dark Photon no "<<muon_tk_properties[i][12]<<" | Muon no "<<muon_tk_properties[i][13]<<" | dR(s1,s2) = "<<std::setw(10)<<muon_tk_properties[i][14]<<std::endl;
      }
      std::cout<<"--- Muon Pair Properties: --------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      for(int i=0; i<2; ++i) {
	std::cout<<"--- DiMu "<<i<<" pt = "<<std::setw(11)<<dimu_tk_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<dimu_tk_properties[i][1];
	std::cout<<" | phi = "<<std::setw(10)<<dimu_tk_properties[i][2]<<" | q = "<<dimu_tk_properties[i][5];
	std::cout<<" | vtx = ["<<std::setw(10)<<dimu_tk_properties[i][7]<<", "<<std::setw(10)<<dimu_tk_properties[i][8]<<", "<<std::setw(10)<<dimu_tk_properties[i][9]<<"]";
	std::cout<<" | Lxy = "<<std::setw(10)<<dimu_tk_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<dimu_tk_properties[i][11]<<" cm";
	std::cout<<" | Dark Photon no "<<dimu_tk_properties[i][12]<<" | Dimu no "<<dimu_tk_properties[i][13]<<" | dR(s1,s2) = "<<std::setw(10)<<dimu_tk_properties[i][14]<<std::endl;
      }      
      std::cout<<"=============================================================================================================================================================================================="<<std::endl;
      std::cout<<"\n"<<std::endl;
    }

    // Change from Arrays to Vector Format ...
    for(int i=0; i<4; ++i) {
      std::vector<double> muon_tk_prop_temp;
      for(int j=0; j<n_prop; ++j) { muon_tk_prop_temp.push_back(muon_tk_properties[i][j]); }
      muon_tk_prop_vec.push_back(muon_tk_prop_temp);
    }
    // Sort the Vectors here based on pt (first element in vector) ...
    std::sort(muon_tk_prop_vec.begin(), muon_tk_prop_vec.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
                return a[0] > b[0];
              });
    // Now we are ready to fill the histograms
    for(int i=0; i<2; ++i) {
      Sim_All_DiMuon_pt->Fill(dimu_tk_properties[i][0]);   Sim_All_DiMuon_eta->Fill(dimu_tk_properties[i][1]);  Sim_All_DiMuon_phi->Fill(dimu_tk_properties[i][2]);
      Sim_All_DiMuon_lxy->Fill(dimu_tk_properties[i][10]); Sim_All_DiMuon_lz->Fill(muon_tk_prop_vec[i][11]);    Sim_All_DiMuon_m->Fill(dimu_tk_properties[i][3]);   Sim_All_DiMuon_dR->Fill(dimu_tk_properties[i][14]);
    }
    for(int i=0; i<4; ++i) {
      Sim_All_Muon_pt->Fill(muon_tk_prop_vec[i][0]);   Sim_All_Muon_eta->Fill(muon_tk_prop_vec[i][1]); Sim_All_Muon_phi->Fill(muon_tk_prop_vec[i][2]);
      Sim_All_Muon_lxy->Fill(muon_tk_prop_vec[i][10]); Sim_All_Muon_lz->Fill(muon_tk_prop_vec[i][11]);
    }
    Sim_Muon_pt1_pt->Fill(muon_tk_prop_vec[0][0]);   Sim_Muon_pt1_eta->Fill(muon_tk_prop_vec[0][1]); Sim_Muon_pt1_phi->Fill(muon_tk_prop_vec[0][2]);
    Sim_Muon_pt1_lxy->Fill(muon_tk_prop_vec[0][10]); Sim_Muon_pt1_lz->Fill(muon_tk_prop_vec[0][11]);
    Sim_Muon_pt2_pt->Fill(muon_tk_prop_vec[1][0]);   Sim_Muon_pt2_eta->Fill(muon_tk_prop_vec[1][1]); Sim_Muon_pt2_phi->Fill(muon_tk_prop_vec[1][2]);
    Sim_Muon_pt2_lxy->Fill(muon_tk_prop_vec[1][10]); Sim_Muon_pt2_lz->Fill(muon_tk_prop_vec[1][11]);
    Sim_Muon_pt3_pt->Fill(muon_tk_prop_vec[2][0]);   Sim_Muon_pt3_eta->Fill(muon_tk_prop_vec[2][1]); Sim_Muon_pt3_phi->Fill(muon_tk_prop_vec[2][2]);
    Sim_Muon_pt3_lxy->Fill(muon_tk_prop_vec[2][10]); Sim_Muon_pt3_lz->Fill(muon_tk_prop_vec[2][11]);
    Sim_Muon_pt4_pt->Fill(muon_tk_prop_vec[3][0]);   Sim_Muon_pt4_eta->Fill(muon_tk_prop_vec[3][1]); Sim_Muon_pt4_phi->Fill(muon_tk_prop_vec[3][2]);
    Sim_Muon_pt4_lxy->Fill(muon_tk_prop_vec[3][10]); Sim_Muon_pt4_lz->Fill(muon_tk_prop_vec[3][11]);

    /*
    if(genDebug)  std::cout<<"=============================================================================================================================================================================================="<<std::endl;
    if(genDebug)  std::cout<<"\n"<<std::endl;
    // 3) SIM-LEVEL :: Match Muon SimTracks to SimVertices
    std::vector<SimVertex> theSimVertices;
    edm::Handle<edm::SimVertexContainer> SimVtx;
    iEvent.getByLabel("g4SimHits",SimVtx);
    if(genDebug)  std::cout<<"=== Analysis of Sim Vertices  :: ============   Simulated Vertices: "<<SimVtx->size()<<std::endl;
    if(genDebug)  std::cout<<"=============================================================================================================================================================================================="<<std::endl;
    */
  }



  // ===================================================================
  // ===      Loop over Stand Alone Muon Collection                  ===
  // ===        ---> and use the MC Truth Matching                   ===
  // ===        ---> or Tight Muon ID and pt cut                     ===
  // ===================================================================

  // reset the Matched value for the SimTrackMuons --> from now on this is used to determine the efficiency of RSA muons ...
  for(int i=0; i<4; ++i) { muon_tk_properties[i][6] = false; }

  if(genDebug) {
    // std::cout<<"\n"<<std::endl;
    std::cout<<"=== Analysis of Refitted Stand Alone Muon Tracks    :: ============   RSA Muon Tracks: "<<staTracks_RSAMuon->size()<<std::endl;
    std::cout<<"=============================================================================================================================================================================================="<<std::endl;
  }
  // loop over the stand alone muon collection
  int muon_rsa_count = 0;
  double muon_rsa_prop_temp[8][4][n_prop_reco]; 
  std::vector< std::vector< double > > indexRSAMuon; // Keep here the delta R for all SimTrack Muon and Stand Alone Muon combination < SimTrack Muon index, StandAlone Muon index, Delta R > 

  for(int i=0; i<4; ++i) { for (int j=0; j<8; ++j) { for(int k=0; k<n_prop_reco; ++k) { muon_rsa_prop_temp[j][i][k]=-1;}}}

  for (staTrack = staTracks_RSAMuon->begin(); staTrack != staTracks_RSAMuon->end(); ++staTrack) {
    reco::TransientTrack track(*staTrack,&*theMGField,theTrackingGeometry);
    if(genDebug) {
      std::cout<<"RSA Muon Track Found: p = "<<std::setw(10)<<track.impactPointTSCP().momentum().mag()<<" pt = "<<std::setw(11)<<track.impactPointTSCP().momentum().perp();
      std::cout<<" eta = "<<std::setw(10)<<track.impactPointTSCP().momentum().eta()<<" phi = "<<std::setw(10)<<track.impactPointTSCP().momentum().phi();
      std::cout<<" q = "<<track.charge()<<" chi2 = "<<std::setw(4)<<track.chi2()<<" with "<<std::setw(2)<<staTrack->recHitsSize()<<" rechits || ";
    }
    double track_eta = track.impactPointTSCP().momentum().eta();
    double track_phi = track.impactPointTSCP().momentum().phi();

    // Match RecoMuon to closest by Sim Track Muon
    // Traditional Matching :: dR < 0.15
    // -------------------------------------------
    if(!recoDRMatch) {
      // Loop over 4 (or less) SimTrack Muons
      for(int i=0; i<tot_sim_ev; ++i) {
	if(muon_rsa_properties[i][6]) {
	  if(i==3) {if(genDebug) std::cout<<" ---> NOT MATCHED"<<std::endl;}
	  continue; // if SimTrack is already matched, stop the loop                                                                                                                                                                                     
	}
	if(genDebug) std::cout<<"Comparing to SimTrack "<<i<<" --> delta R = "<<deltaR(track_eta, track_phi, muon_tk_properties[i][1], muon_tk_properties[i][2]);
	if(!muon_tk_properties[i][6] && deltaR(track_eta, track_phi, muon_tk_properties[i][1], muon_tk_properties[i][2]) < delta_R_simrec) {
	  muon_tk_properties[i][6] = true;
	  muon_rsa_properties[muon_rsa_count][0] = track.impactPointTSCP().momentum().perp();
	  muon_rsa_properties[muon_rsa_count][1] = track.impactPointTSCP().momentum().eta();
	  muon_rsa_properties[muon_rsa_count][2] = track.impactPointTSCP().momentum().phi();
	  // muon_rsa_properties[muon_rsa_count][3] = track.LorentzVector().E();
	  // muon_rsa_properties[muon_rsa_count][4] = track.LorentzVector().M();
	  muon_rsa_properties[muon_rsa_count][5] = track.charge();
	  muon_rsa_properties[muon_rsa_count][6] = true;
	  if(genDebug) std::cout<<" ---> MATCHED (based on Delta R)"<<std::endl;
	  // Matched to Muon i ==> write properties also in SimTrack array                                                                                                                                                                                  
	  for(int j=7; j<n_prop; ++j) { muon_rsa_properties[muon_rsa_count][j] = muon_tk_properties[i][j]; } 
	}
	else {
	  if(genDebug) std::cout<<" ---> NOT MATCHED"<<std::endl;
	}
	// If Stand Alone Muon Track is matched, stop the loop                                                                                                                                            
	if(muon_rsa_properties[muon_rsa_count][6] == true) {
	  ++muon_rsa_count;
	  ++tot_rsa;
	  break;
	}
      }
    }
    // ---------------------------------------------------

    // Match RecoMuon to closest by Sim Track Muon :: pt 1
    // ---------------------------------------------------
    if(recoDRMatch) {
      if(genDebug) std::cout<<""<<std::endl;
      for(int i=0; i<tot_sim_ev; ++i) {
	muon_rsa_prop_temp[muon_rsa_count][i][0] = track.impactPointTSCP().momentum().perp();
        muon_rsa_prop_temp[muon_rsa_count][i][1] = track.impactPointTSCP().momentum().eta();
	muon_rsa_prop_temp[muon_rsa_count][i][2] = track.impactPointTSCP().momentum().phi();
	// muon_rsa_prop_temp[muon_rsa_count][i][3] = track.LorentzVector().E();
	// muon_rsa_prop_temp[muon_rsa_count][i][4] = track.LorentzVector().M();
	muon_rsa_prop_temp[muon_rsa_count][i][5] = track.charge();
	muon_rsa_prop_temp[muon_rsa_count][i][6] = true;
	for(int j=7; j<n_prop; ++j) { muon_rsa_prop_temp[muon_rsa_count][i][j] = muon_tk_properties[i][j]; } 
	double dR = deltaR(track_eta, track_phi, muon_tk_properties[i][1], muon_tk_properties[i][2]);
	muon_rsa_prop_temp[muon_rsa_count][i][15] = dR;

	std::vector<double> tempindexvect;
	if(techDebug) std::cout<<" Added to indexRSAMuon vector :: ["<<muon_tk_properties[i][13]<<","<<muon_rsa_count<<","<<dR<<"]"<<std::endl;
	tempindexvect.push_back(muon_tk_properties[i][13]); tempindexvect.push_back(muon_rsa_count); tempindexvect.push_back(dR);
	indexRSAMuon.push_back(tempindexvect);
      }
      ++muon_rsa_count;
    }
  }

  // Match RecoMuon to closest by Sim Track Muon :: pt 2
  // ---------------------------------------------------
  if(recoDRMatch) {
    // Sort indexRSAMuon vector based on dR
    // Sorting needs to be done only once
    if(techDebug) {
      std::cout<<"Vectors are not sorted :: ";
      for (std::vector< std::vector<double> >::const_iterator iMuon = indexRSAMuon.begin(); iMuon != indexRSAMuon.end(); ++iMuon) {
	std::cout<<"["<<(*iMuon)[0]<<","<<(*iMuon)[1]<<","<<(*iMuon)[2]<<"], ";
      }
      std::cout<<""<<std::endl;
    }
    // Sort the Vectors here based on dR (3rd element in vector) ...
    std::sort(indexRSAMuon.begin(), indexRSAMuon.end(),
	      [](const std::vector<double>& a, const std::vector<double>& b) {
		return a[2] < b[2];
	      });
    if(techDebug) {
      std::cout<<"Vectors are sorted now :: ";
      for (std::vector< std::vector<double> >::const_iterator iMuon = indexRSAMuon.begin(); iMuon != indexRSAMuon.end(); ++iMuon) {
	std::cout<<"["<<(*iMuon)[0]<<","<<(*iMuon)[1]<<","<<(*iMuon)[2]<<"], ";
      }
      std::cout<<""<<std::endl;
    }
    // Run this as long as the vector indexRSAMuon contains elements ...
    while(indexRSAMuon.size()>0) {
      // Take element with smalles dR and save this combination
      int sim_muon_index = (*indexRSAMuon.begin())[0];
      int rsa_muon_index = (*indexRSAMuon.begin())[1];
      if(techDebug) std::cout<<"Fill muon_rsa_properties vector at position :: sim_muon_index = "<<sim_muon_index<<" rsa_muon_index = "<<rsa_muon_index<<std::endl;
      for(int j=0; j<n_prop; ++j) {
	muon_rsa_properties[sim_muon_index][j] = muon_rsa_prop_temp[rsa_muon_index][sim_muon_index][j];
      }
      // muon_rsa_properties[sim_muon_index][14] = (*indexRSAMuon.begin())[2]; // instead of keeping the Delta R between 2 muons at sim level, i prefer to store the Delta R between the matched sim and reco muon
      // Remove all other matches with the same SimTrack Muon and RSA Muon in the vector with the indices
      std::vector<int> pos_to_be_removed; int pos_counter = 0;
      if(techDebug) std::cout<<" Size of Vector before deletion :: "<<indexRSAMuon.size()<<std::endl;
      for (std::vector< std::vector<double> >::const_iterator iMuon = indexRSAMuon.begin(); iMuon != indexRSAMuon.end(); ++iMuon) {
	if((*iMuon)[0]==sim_muon_index || (*iMuon)[1]==rsa_muon_index) { 
	  if(techDebug) std::cout<<"["<<(*iMuon)[0]<<","<<(*iMuon)[1]<<","<<(*iMuon)[2]<<"] at pos = "<<pos_counter<<" scheduled for removal"<<std::endl; 
	  pos_to_be_removed.push_back(pos_counter); }
	++pos_counter;
      }
      // tricky :: sort this vector and start removing elements with the highest position
      // imagine you want to remove the elements at position 2 and 3 in a vector
      // if you remove first the element at position 2, position 3 becomes position 2 and position 4 becomes position 3
      // then removing position 3 makes you actually remove the element that was originally at position 4
      std::sort(pos_to_be_removed.begin(), pos_to_be_removed.end(),
		[](const int a, const int b) {
		  return a > b;
		});
      // now remove the elements, starting with the element at highest position
      for (std::vector< int >::const_iterator iPos = pos_to_be_removed.begin(); iPos != pos_to_be_removed.end(); ++iPos) {
	int pos = *iPos;
	indexRSAMuon.erase(indexRSAMuon.begin()+pos);
      }
      if(techDebug) {
	std::cout<<" Size of Vector after deletion :: "<<indexRSAMuon.size()<<" Elements :: ";
	for (std::vector< std::vector<double> >::const_iterator iMuon = indexRSAMuon.begin(); iMuon != indexRSAMuon.end(); ++iMuon) {
	  std::cout<<"["<<(*iMuon)[0]<<","<<(*iMuon)[1]<<","<<(*iMuon)[2]<<"], ";
	}
	std::cout<<""<<std::endl;
      }
    }
  }
  // ---------------------------------------------------

  // Muon Pair Properties
  int dimu_rsa_count = 0;
  for(int i=0; i<4; ++i) {
    for(int j=0; j<i; ++j) {
      if(dimu_rsa_count==2) continue;
      if(muon_rsa_properties[i][12]==muon_rsa_properties[j][12] && muon_rsa_properties[i][12]!=dimu_rsa_properties[dimu_rsa_count][12] && !dimu_rsa_properties[dimu_rsa_count][6]) {
	TLorentzVector Mu1, Mu2, DiMu;
	// Mu1.SetPtEtaPhiM
	// Mu1.SetPxPyPzE
	Mu1.SetPtEtaPhiE(muon_rsa_properties[i][0],muon_rsa_properties[i][1],muon_rsa_properties[i][2],muon_rsa_properties[i][3]);
	Mu2.SetPtEtaPhiE(muon_rsa_properties[j][0],muon_rsa_properties[j][1],muon_rsa_properties[j][2],muon_rsa_properties[j][3]);
	DiMu = Mu1+Mu2;
	dimu_rsa_properties[dimu_rsa_count][0] = DiMu.Pt();
	dimu_rsa_properties[dimu_rsa_count][1] = DiMu.Eta();
	dimu_rsa_properties[dimu_rsa_count][2] = DiMu.Phi();
	dimu_rsa_properties[dimu_rsa_count][3] = DiMu.E();
	dimu_rsa_properties[dimu_rsa_count][4] = DiMu.M();
	dimu_rsa_properties[dimu_rsa_count][5] = muon_rsa_properties[i][5]+muon_rsa_properties[j][5];
	dimu_rsa_properties[dimu_rsa_count][6] = true;
	for(int k=7; k<n_prop-1; ++k) { dimu_rsa_properties[dimu_rsa_count][k] = muon_rsa_properties[i][k]; }
	// number
	dimu_rsa_properties[dimu_rsa_count][13] = dimu_rsa_count;
	// Additional Info
	dimu_rsa_properties[dimu_rsa_count][14] = deltaR(muon_rsa_properties[i][1],muon_rsa_properties[i][2], muon_rsa_properties[j][1], muon_rsa_properties[j][2]);
	muon_rsa_properties[i][14] = dimu_rsa_properties[dimu_rsa_count][14];
	muon_rsa_properties[j][14] = dimu_rsa_properties[dimu_rsa_count][14]; 
	// Counter
	++dimu_rsa_count;
      }
      else { std::cout<<""<<std::endl; }
    }
  }
  if(genDebug)  {
    // std::cout<<"\n"<<std::endl;
    std::cout<<"--- Muons found: -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    for(int i=0; i<4; ++i) {
      std::cout<<"--- Muon "<<i<<" pt = "<<std::setw(11)<<muon_rsa_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<muon_rsa_properties[i][1];
      std::cout<<" | phi = "<<std::setw(10)<<muon_rsa_properties[i][2]<<" | q = "<<muon_rsa_properties[i][5];
      std::cout<<" | vtx = ["<<std::setw(10)<<muon_rsa_properties[i][7]<<", "<<std::setw(10)<<muon_rsa_properties[i][8]<<", "<<std::setw(10)<<muon_rsa_properties[i][9]<<"]";
      std::cout<<" | Lxy = "<<std::setw(10)<<muon_rsa_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<muon_rsa_properties[i][11]<<" cm";
      std::cout<<" | Dark Photon no "<<muon_rsa_properties[i][12]<<" | Muon no "<<muon_rsa_properties[i][13]<< " | dR(r1,r2) = "<<std::setw(10)<<muon_rsa_properties[i][14]<<std::endl;
    }
    std::cout<<"--- Muon Pair Properties: --------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    for(int i=0; i<2; ++i) {
      std::cout<<"--- DiMu "<<i<<" pt = "<<std::setw(11)<<dimu_rsa_properties[i][0]<<" GeV/c | eta = "<<std::setw(10)<<dimu_rsa_properties[i][1];
      std::cout<<" | phi = "<<std::setw(10)<<dimu_rsa_properties[i][2]<<" | q = "<<dimu_rsa_properties[i][5];
      std::cout<<" | vtx = ["<<std::setw(10)<<dimu_rsa_properties[i][7]<<", "<<std::setw(10)<<dimu_rsa_properties[i][8]<<", "<<std::setw(10)<<dimu_rsa_properties[i][9]<<"]";
      std::cout<<" | Lxy = "<<std::setw(10)<<dimu_rsa_properties[i][10]<<" cm | Lz = "<<std::setw(10)<<dimu_rsa_properties[i][11]<<" cm";
      std::cout<<" | Dark Photon no "<<dimu_rsa_properties[i][12]<<" | Dimu no "<<dimu_rsa_properties[i][13]<<" | dR(r1,r2) = "<<std::setw(10)<<dimu_rsa_properties[i][14]<<std::endl;
    }      
    std::cout<<"=============================================================================================================================================================================================="<<std::endl;
    std::cout<<"\n"<<std::endl;
  }

  // Change from Arrays to Vector Format ...
  for(int i=0; i<4; ++i) {
    std::vector<double> muon_rsa_prop_temp;
    for(int j=0; j<n_prop; ++j) { muon_rsa_prop_temp.push_back(muon_rsa_properties[i][j]); }
    muon_rsa_prop_vec.push_back(muon_rsa_prop_temp);
  }

  // Sort the Vectors here based on pt (first element in vector) ...
  std::sort(muon_rsa_prop_vec.begin(), muon_rsa_prop_vec.end(),
	    [](const std::vector<double>& a, const std::vector<double>& b) {
	      return a[0] > b[0];
	    });

  // Now we are ready to fill the histograms
  for(int i=0; i<2; ++i) {
    RSA_All_DiMuon_pt->Fill(dimu_rsa_properties[i][0]);   RSA_All_DiMuon_eta->Fill(dimu_rsa_properties[i][1]);  RSA_All_DiMuon_phi->Fill(dimu_rsa_properties[i][2]);
    RSA_All_DiMuon_lxy->Fill(dimu_rsa_properties[i][10]); RSA_All_DiMuon_lz->Fill(muon_rsa_prop_vec[i][11]);    RSA_All_DiMuon_m->Fill(dimu_rsa_properties[i][3]);   RSA_All_DiMuon_dR->Fill(dimu_rsa_properties[i][14]);
  }
  for(int i=0; i<4; ++i) {
    RSA_All_Muon_pt->Fill(muon_rsa_prop_vec[i][0]);   RSA_All_Muon_eta->Fill(muon_rsa_prop_vec[i][1]); RSA_All_Muon_phi->Fill(muon_rsa_prop_vec[i][2]);
    RSA_All_Muon_lxy->Fill(muon_rsa_prop_vec[i][10]); RSA_All_Muon_lz->Fill(muon_rsa_prop_vec[i][11]);
  }
  RSA_Muon_pt1_pt->Fill(muon_rsa_prop_vec[0][0]);   RSA_Muon_pt1_eta->Fill(muon_rsa_prop_vec[0][1]); RSA_Muon_pt1_phi->Fill(muon_rsa_prop_vec[0][2]);
  RSA_Muon_pt1_lxy->Fill(muon_rsa_prop_vec[0][10]); RSA_Muon_pt1_lz->Fill(muon_rsa_prop_vec[0][11]);
  RSA_Muon_pt2_pt->Fill(muon_rsa_prop_vec[1][0]);   RSA_Muon_pt2_eta->Fill(muon_rsa_prop_vec[1][1]); RSA_Muon_pt2_phi->Fill(muon_rsa_prop_vec[1][2]);
  RSA_Muon_pt2_lxy->Fill(muon_rsa_prop_vec[1][10]); RSA_Muon_pt2_lz->Fill(muon_rsa_prop_vec[1][11]);
  RSA_Muon_pt3_pt->Fill(muon_rsa_prop_vec[2][0]);   RSA_Muon_pt3_eta->Fill(muon_rsa_prop_vec[2][1]); RSA_Muon_pt3_phi->Fill(muon_rsa_prop_vec[2][2]);
  RSA_Muon_pt3_lxy->Fill(muon_rsa_prop_vec[2][10]); RSA_Muon_pt3_lz->Fill(muon_rsa_prop_vec[2][11]);
  RSA_Muon_pt4_pt->Fill(muon_rsa_prop_vec[3][0]);   RSA_Muon_pt4_eta->Fill(muon_rsa_prop_vec[3][1]); RSA_Muon_pt4_phi->Fill(muon_rsa_prop_vec[3][2]);
  RSA_Muon_pt4_lxy->Fill(muon_rsa_prop_vec[3][10]); RSA_Muon_pt4_lz->Fill(muon_rsa_prop_vec[3][11]);
  
  RSA_All_DiMuon_pt->Fill(dimu_rsa_properties[0][0]);     RSA_All_DiMuon_pt->Fill(dimu_rsa_properties[1][0]);
  RSA_All_DiMuon_eta->Fill(dimu_rsa_properties[0][1]);    RSA_All_DiMuon_eta->Fill(dimu_rsa_properties[1][1]);
  RSA_All_DiMuon_phi->Fill(dimu_rsa_properties[0][2]);    RSA_All_DiMuon_phi->Fill(dimu_rsa_properties[1][2]);
  RSA_All_DiMuon_lxy->Fill(dimu_rsa_properties[0][10]);   RSA_All_DiMuon_lxy->Fill(dimu_rsa_properties[1][10]);
  RSA_All_DiMuon_lz->Fill(muon_rsa_prop_vec[0][11]);      RSA_All_DiMuon_lz->Fill(muon_rsa_prop_vec[1][11]);
  RSA_All_DiMuon_m->Fill(dimu_rsa_properties[0][3]);      RSA_All_DiMuon_m->Fill(dimu_rsa_properties[1][3]);
  RSA_All_DiMuon_dR->Fill(dimu_rsa_properties[0][14]);    RSA_All_DiMuon_dR->Fill(dimu_rsa_properties[1][14]);
  
  // Clear all vectors
  /*
  dark_me_prop_vec.clear();
  muon_me_prop_vec.clear();
  dimu_me_prop_vec.clear();
  muon_tk_prop_vec.clear();
  // dimu_tk_prop_vec.clear();
  muon_rsa_prop_vec.clear();
  // dimu_rsa_prop_vec.clear();
  indexRSAMuon.clear();
  */
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyDisplacedMuonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyDisplacedMuonAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MyDisplacedMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyDisplacedMuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyDisplacedMuonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyDisplacedMuonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyDisplacedMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyDisplacedMuonAnalyzer);
