// -*- C++ -*-
//
// Package:    MyStandAloneMuonAnalyzer
// Class:      MyStandAloneMuonAnalyzer
// 
/**\class MyStandAloneMuonAnalyzer MyStandAloneMuonAnalyzer.cc MyAnalyzers/MyStandAloneMuonAnalyzer/plugins/MyStandAloneMuonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen
//         Created:  Thu, 05 Jun 2014 17:16:14 GMT
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

#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//
// class declaration
//

class MuonServiceProxy;

class MyStandAloneMuonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyStandAloneMuonAnalyzer(const edm::ParameterSet&);
      ~MyStandAloneMuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::InputTag STAMuLabel, GLBMuLabel;
  // std::string muonLabel1, muonLabel2;
  std::string rootFileName;
  bool physDebug, techDebug, genDebug, mcTruthMatch, wantTightId, wantLooseId;

  const reco::GenParticle     *g, *b, *d1, *d2;

  TFile * outputfile;
  TDirectoryFile * AllHits, * SegmentsAndHits;  

  TH1F * Muon_All, * Muon_All_MyReBin, * Muon_All_Reduced;

  TH1F * Rechits_All,         * Rechits_DT,          * Rechits_CSC,           * Rechits_GEM,           * Rechits_RPC;
  TH1F * Rechits_All_Hits,    * Rechits_DT_Hits,     * Rechits_CSC_Hits,      * Rechits_GEM_Hits;
  TH1F * Rechits_All_PhiHits, * Rechits_All_EtaHits, * Rechits_DT_Hits_EtaSL, * Rechits_DT_Hits_PhiSL; 


  TH1F * Rechits_All_Eta_1D,         * Rechits_DT_Eta_1D,            * Rechits_CSC_Eta_1D,           * Rechits_GEM_Eta_1D,           * Rechits_RPC_Eta_1D;
  TH1F * Rechits_All_Hits_Eta_1D,    * Rechits_DT_Hits_Eta_1D,       * Rechits_CSC_Hits_Eta_1D,      * Rechits_GEM_Hits_Eta_1D;
  TH1F * Rechits_All_PhiHits_Eta_1D, * Rechits_All_EtaHits_Eta_1D,   * Rechits_DT_Hits_EtaSL_Eta_1D, * Rechits_DT_Hits_PhiSL_Eta_1D; 

  TH2F * Rechits_All_Eta_2D, * Rechits_DT_Eta_2D, * Rechits_CSC_Eta_2D, * Rechits_RPC_Eta_2D, * Rechits_GEM_Eta_2D;
  TH2F * Rechits_All_Eta_2D_MyReBin, * Rechits_DT_Eta_2D_MyReBin, * Rechits_CSC_Eta_2D_MyReBin, * Rechits_RPC_Eta_2D_MyReBin, * Rechits_GEM_Eta_2D_MyReBin;
  TH2F * Rechits_All_Eta_2D_Reduced, * Rechits_DT_Eta_2D_Reduced, * Rechits_CSC_Eta_2D_Reduced, * Rechits_RPC_Eta_2D_Reduced, * Rechits_GEM_Eta_2D_Reduced;
  // New request Daniele
  TH2F * Rechits_Seg_Eta_2D_Reduced, * Rechits_Hit_Eta_2D_Reduced;


  TH1F * StandAloneMuon_PT, * StandAloneMuon_Phi, * StandAloneMuon_Eta, * StandAloneMuon_Hit;
  TH1F * StandAloneMuon_1p8To2p5_2Hit_PT, * StandAloneMuon_1p8To2p5_2Hit_Phi, * StandAloneMuon_1p8To2p5_2Hit_Eta, * StandAloneMuon_1p8To2p5_2Hit_Hit, * StandAloneMuon_1p8To2p5_2Hit_Sta;
  TH1F * Event_NumStaMuons, * Event_NumStaMuons_1p8To2p5_2Hit;
  TH1F * GlobalMuon_PT;

  // switch this off in CMSSW_6_2_0_SLHCX
  // also not necessary in CMSSW_7_X_Y if the Muon POG stuff is not used 
  // Find the segments associated to the track
  // SegmentsTrackAssociator* theSegmentsAssociator;
  // MuonServiceProxy * theService;
};

//
// constants, enums and typedefs
//
//     reduced ETA binning (looking at holes in detector acceptance) :: 15 bins
double reduced[]    = {-2.5, -1.7, -1.6, -1.1, -0.9, -0.8, -0.3, -0.2, 0.2, 0.3, 0.8, 0.9, 1.1, 1.6, 1.7, 2.5};
//     reduced ETA binning (used by Muon POG) :: 15 bins
double pogbinning[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4};
//     realistic PHI binning (10 degree sectors, 0 is in the middle of chamber 1) but for statistics lets take 3 chambers together ==> 12 30* sectors
double realist[] = {-3.2288, -2.7052, -2.1816, -1.6580, -1.1344, -0.6108, -0.0872, +0.4364, 0.960, 1.483, 2.007, 2.531, 3.0544};
//     names of the ENDCAP Muon Stations 
std::string station[] = {"ME1", "ME2", "ME3", "ME4", "RE1", "RE2", "RE3", "RE4"};
//     DELTA REL PT and DELTA R matching cones
double delta_rel_pt = 0.1;
double delta_R      = 0.1;

//
// static data member definitions
//

//
// constructors and destructor
//
MyStandAloneMuonAnalyzer::MyStandAloneMuonAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName   = iConfig.getUntrackedParameter<std::string>("RootFileName");
  physDebug      = iConfig.getUntrackedParameter<bool>("PhysicsDebug");
  techDebug      = iConfig.getUntrackedParameter<bool>("TechnicDebug");
  genDebug       = iConfig.getUntrackedParameter<bool>("GenPartDebug");
  mcTruthMatch   = iConfig.getUntrackedParameter<bool>("MCTruthMatching");
  wantTightId    = iConfig.getUntrackedParameter<bool>("TightID");
  wantLooseId    = iConfig.getUntrackedParameter<bool>("LooseID");
  // muonLabel1     = iConfig.getUntrackedParameter<std::string>("MuonLabel1");
  // muonLabel2     = iConfig.getUntrackedParameter<std::string>("MuonLabel2");
  STAMuLabel = iConfig.getParameter<edm::InputTag>("StandAloneTrackCollectionLabel");
  GLBMuLabel = iConfig.getParameter<edm::InputTag>("GlobalTrackCollectionLabel");

  // switch this off in CMSSW_6_2_0_SLHCX
  // also not necessary in CMSSW_7_X_Y if the Muon POG stuff is not used
  // const edm::ParameterSet SegmentsTrackAssociatorParameters = iConfig.getParameter<edm::ParameterSet>("SegmentsTrackAssociatorParameters");
  // theService = new MuonServiceProxy(iConfig.getParameter<edm::ParameterSet>("ServiceParameters"));

  if(techDebug) std::cout<<"[MyStandAloneMuonAnalyzer :: Constructor]"<<std::endl;
  outputfile      = new TFile(rootFileName.c_str(), "RECREATE" );

  AllHits         = new TDirectoryFile("AllHits",         "AllHits");
  SegmentsAndHits = new TDirectoryFile("SegmentsAndHits", "SegmentsAndHits");

  Muon_All         = new TH1F("Muon_All", "Stand Alone Muons :: All", 50, -2.5, 2.5);
  Muon_All_MyReBin = new TH1F("Muon_All_MyReBin", "Stand Alone Muons :: All :: Rebinned #eta Range", 15, reduced);
  Muon_All_Reduced = new TH1F("Muon_All_Reduced", "Stand Alone Muons :: All :: Rebinned #eta Range", 15, pogbinning);

  StandAloneMuon_PT  = new TH1F("StandAloneMuon_PT",  "Stand Alone Muon momentum", 100, 0, 500);
  StandAloneMuon_Phi = new TH1F("StandAloneMuon_Phi", "Stand Alone Muon azimuthal angle", 12, realist);
  StandAloneMuon_Eta = new TH1F("StandAloneMuon_Eta", "Stand Alone Muon pseudorapidity",  100, -2.5, 2.5);
  StandAloneMuon_Hit = new TH1F("StandAloneMuon_Hit", "Number of Muon Hits (Segments/RPC Hits) of the Stand Alone Muon", 16, -0.5, 15.5);

  StandAloneMuon_1p8To2p5_2Hit_PT  = new TH1F("StandAloneMuon_1p8To2p5_2Hit_PT",  "Stand Alone Muon momentum [1.7 < eta < 2.5 && 2 hits]", 100, 0, 500);
  StandAloneMuon_1p8To2p5_2Hit_Phi = new TH1F("StandAloneMuon_1p8To2p5_2Hit_Phi", "Stand Alone Muon azimuthal angle [1.7 < eta < 2.5 && 2 hits]", 12, realist);
  StandAloneMuon_1p8To2p5_2Hit_Eta = new TH1F("StandAloneMuon_1p8To2p5_2Hit_Eta", "Stand Alone Muon pseudorapidity [1.7 < eta < 2.5 && 2 hits]",  100, -2.5, 2.5);
  StandAloneMuon_1p8To2p5_2Hit_Hit = new TH1F("StandAloneMuon_1p8To2p5_2Hit_Hit", "Number of Muon Hits (Segments/RPC Hits) of the Stand Alone Muon [1.7 < eta < 2.5 && 2 hits]", 16, -0.5, 15.5);
  StandAloneMuon_1p8To2p5_2Hit_Sta = new TH1F("StandAloneMuon_1p8To2p5_2Hit_Sta", "Station of the Muon Hits (Segments/RPC Hits)  of the Stand Alone Muon [1.7 < eta < 2.5 && 2 hits]", 8, 0.5, 8.5);

  Event_NumStaMuons                = new TH1F("Event_NumStaMuons",               "Number of Stand Alone Muons Per Event", 11, -0.5, 10.5);
  Event_NumStaMuons_1p8To2p5_2Hit  = new TH1F("Event_NumStaMuons_1p8To2p5_2Hit", "Number of Stand Alone Muons Per Event [1.7 < eta < 2.5 && 2 hits]", 11, -0.5, 10.5);
  GlobalMuon_PT                    = new TH1F("GlobalMuon_PT",                   "Global Muon momentum", 100, 0, 500);


  // Plots for segments + RPC hits
  // -----------------------------
  Rechits_All = new TH1F("Rechits_All", "RecHits :: All (Segments/RPC Hits)", 50, -2.5, 2.5);
  Rechits_DT  = new TH1F("Rechits_DT",  "RecHits :: DT  (Segments)", 50, -2.5, 2.5);
  Rechits_CSC = new TH1F("Rechits_CSC", "RecHits :: CSC (Segments)", 50, -2.5, 2.5);
  Rechits_RPC = new TH1F("Rechits_RPC", "RecHits :: RPC (RPC Hits)", 50, -2.5, 2.5);
  Rechits_GEM = new TH1F("Rechits_GEM", "RecHits :: GEM (GEM Hits)", 50, -2.5, 2.5);

  Rechits_All_Eta_1D = new TH1F("Rechits_All_Eta_1D", "Average amount of RecHits :: All (Segments/RPC Hits)", 50, -2.5, 2.5);
  Rechits_DT_Eta_1D  = new TH1F("Rechits_DT_Eta_1D",  "Average amount of RecHits :: DT  (Segments)",  50, -2.5, 2.5);
  Rechits_CSC_Eta_1D = new TH1F("Rechits_CSC_Eta_1D", "Average amount of RecHits :: CSC (Segments)", 50, -2.5, 2.5);
  Rechits_RPC_Eta_1D = new TH1F("Rechits_RPC_Eta_1D", "Average amount of RecHits :: RPC (RPC Hits)", 50, -2.5, 2.5);
  Rechits_GEM_Eta_1D = new TH1F("Rechits_RPC_Eta_1D", "Average amount of RecHits :: RPC (RPC Hits)", 50, -2.5, 2.5);


  // Plots for single hits
  // ---------------------
  Rechits_All_Hits = new TH1F("Rechits_All_Hits", "RecHits :: All (Hits)", 50, -2.5, 2.5);
  Rechits_DT_Hits  = new TH1F("Rechits_DT_Hits",  "RecHits :: DT  (Hits)", 50, -2.5, 2.5);
  Rechits_CSC_Hits = new TH1F("Rechits_CSC_Hits", "RecHits :: CSC (Hits)", 50, -2.5, 2.5);
  Rechits_GEM_Hits = new TH1F("Rechits_GEM_Hits", "RecHits :: GEM (Hits)", 50, -2.5, 2.5);

  Rechits_All_PhiHits    = new TH1F("Rechits_All_PhiHits", "RecHits :: All (Hits) :: #phi", 50, -2.5, 2.5);
  Rechits_All_EtaHits    = new TH1F("Rechits_All_EtaHits", "RecHits :: All (Hits) :: #eta", 50, -2.5, 2.5);
  Rechits_DT_Hits_PhiSL  = new TH1F("Rechits_DT_Hits_Phi", "RecHits :: DT :: #phi",  50, -2.5, 2.5);
  Rechits_DT_Hits_EtaSL  = new TH1F("Rechits_DT_Hits_Eta", "RecHits :: DT :: #eta",  50, -2.5, 2.5);

  Rechits_All_Hits_Eta_1D = new TH1F("Rechits_All_Hits_Eta_1D", "Average amount of RecHits :: All (Hits)", 50, -2.5, 2.5);
  Rechits_DT_Hits_Eta_1D  = new TH1F("Rechits_DT_Hits_Eta_1D",  "Average amount of RecHits :: DT  (Hits)", 50, -2.5, 2.5);
  Rechits_CSC_Hits_Eta_1D = new TH1F("Rechits_CSC_Hits_Eta_1D", "Average amount of RecHits :: CSC (Hits)", 50, -2.5, 2.5);
  Rechits_GEM_Hits_Eta_1D = new TH1F("Rechits_Eta_Hits_Eta_1D", "Average amount of RecHits :: Eta (Hits)", 50, -2.5, 2.5);

  Rechits_All_PhiHits_Eta_1D    = new TH1F("Rechits_All_PhiHits_Eta_1D",    "Average amount of RecHits :: All (Hits) :: #phi", 50, -2.5, 2.5);
  Rechits_All_EtaHits_Eta_1D    = new TH1F("Rechits_All_EtaHits_Eta_1D",    "Average amount of RecHits :: All (Hits) :: #eta", 50, -2.5, 2.5);
  Rechits_DT_Hits_PhiSL_Eta_1D  = new TH1F("Rechits_DT_Hits_PhiSL_Eta_1D",  "Average amount of RecHits :: DT ::Phi",  50, -2.5, 2.5);
  Rechits_DT_Hits_EtaSL_Eta_1D  = new TH1F("Rechits_DT_Hits_EtaSL_Eta_1D",  "Average amount of RecHits :: DT ::Eta",  50, -2.5, 2.5);

  // 2D Plots are only for Segments + RPC Hits
  // -----------------------------------------
  // Nothing changed here
  Rechits_All_Eta_2D = new TH2F("Rechits_All_Eta_2D", "Amount of RecHits :: All", 50, -2.5, 2.5, 101, -0.5, 100.5);
  Rechits_DT_Eta_2D  = new TH2F("Rechits_DT_Eta_2D",  "Amount of RecHits :: DT",  50, -2.5, 2.5, 101, -0.5, 100.5);
  Rechits_CSC_Eta_2D = new TH2F("Rechits_CSC_Eta_2D", "Amount of RecHits :: CSC", 50, -2.5, 2.5, 101, -0.5, 100.5);
  Rechits_RPC_Eta_2D = new TH2F("Rechits_RPC_Eta_2D", "Amount of RecHits :: RPC", 50, -2.5, 2.5, 101, -0.5, 100.5);
  Rechits_GEM_Eta_2D = new TH2F("Rechits_GEM_Eta_2D", "Amount of RecHits :: GEM", 50, -2.5, 2.5, 101, -0.5, 100.5);

  // My Reduced Binning
  Rechits_All_Eta_2D_MyReBin = new TH2F("Rechits_All_Eta_2D_MyReBin", "Amount of RecHits :: All :: Rebinned #eta Range", 15, reduced, 101, -0.5, 100.5);
  Rechits_DT_Eta_2D_MyReBin  = new TH2F("Rechits_DT_Eta_2D_MyReBin",  "Amount of RecHits :: DT :: Rebinned #eta Range",  15, reduced, 101, -0.5, 100.5);
  Rechits_CSC_Eta_2D_MyReBin = new TH2F("Rechits_CSC_Eta_2D_MyReBin", "Amount of RecHits :: CSC :: Rebinned #eta Range", 15, reduced, 101, -0.5, 100.5);
  Rechits_RPC_Eta_2D_MyReBin = new TH2F("Rechits_RPC_Eta_2D_MyReBin", "Amount of RecHits :: RPC :: Rebinned #eta Range", 15, reduced, 101, -0.5, 100.5);
  Rechits_GEM_Eta_2D_MyReBin = new TH2F("Rechits_GEM_Eta_2D_MyReBin", "Amount of RecHits :: GEM :: Rebinned #eta Range", 15, reduced, 101, -0.5, 100.5);

  // POG Reduced Binning
  Rechits_All_Eta_2D_Reduced = new TH2F("Rechits_All_Eta_2D_Reduced", "Amount of RecHits :: All :: Rebinned #eta Range", 15, pogbinning, 101, -0.5, 100.5);
  Rechits_DT_Eta_2D_Reduced  = new TH2F("Rechits_DT_Eta_2D_Reduced",  "Amount of RecHits :: DT :: Rebinned #eta Range",  15, pogbinning, 101, -0.5, 100.5);
  Rechits_CSC_Eta_2D_Reduced = new TH2F("Rechits_CSC_Eta_2D_Reduced", "Amount of RecHits :: CSC :: Rebinned #eta Range", 15, pogbinning, 101, -0.5, 100.5);
  Rechits_RPC_Eta_2D_Reduced = new TH2F("Rechits_RPC_Eta_2D_Reduced", "Amount of RecHits :: RPC :: Rebinned #eta Range", 15, pogbinning, 101, -0.5, 100.5);
  Rechits_GEM_Eta_2D_Reduced = new TH2F("Rechits_GEM_Eta_2D_Reduced", "Amount of RecHits :: GEM :: Rebinned #eta Range", 15, pogbinning, 101, -0.5, 100.5);

  Rechits_Seg_Eta_2D_Reduced = new TH2F("Rechits_Seg_Eta_2D_Reduced", "Amount of RecHits :: All DT/CSC Segments :: Rebinned #eta Range", 15, pogbinning, 101, -0.5, 100.5);
  Rechits_Hit_Eta_2D_Reduced = new TH2F("Rechits_Hit_Eta_2D_Reduced", "Amount of RecHits :: All RPC/GEM Hits    :: Rebinned #eta Range", 15, pogbinning, 101, -0.5, 100.5);

  // switch this off in CMSSW_6_2_0_SLHCX
  // also not necessary in CMSSW_7_X_Y if the Muon POG stuff is not used 
  // edm::ConsumesCollector iC = consumesCollector();
  // theSegmentsAssociator = new SegmentsTrackAssociator(SegmentsTrackAssociatorParameters,iC);
}


MyStandAloneMuonAnalyzer::~MyStandAloneMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  Muon_All->Write();
  Muon_All_MyReBin->Write();
  Muon_All_Reduced->Write();

  StandAloneMuon_PT->Write();
  StandAloneMuon_Phi->Write();
  StandAloneMuon_Eta->Write();
  StandAloneMuon_Hit->Write();

  StandAloneMuon_1p8To2p5_2Hit_PT->Write();
  StandAloneMuon_1p8To2p5_2Hit_Phi->Write();
  StandAloneMuon_1p8To2p5_2Hit_Eta->Write();
  StandAloneMuon_1p8To2p5_2Hit_Hit->Write();

  for(int i=0; i<8; ++i) { StandAloneMuon_1p8To2p5_2Hit_Sta->GetXaxis()->SetBinLabel(i+1,(station[i]).c_str()); }
  StandAloneMuon_1p8To2p5_2Hit_Sta->Write();

  Event_NumStaMuons->Write();
  Event_NumStaMuons_1p8To2p5_2Hit->Write();
  GlobalMuon_PT->Write();

  SegmentsAndHits->cd();
  Rechits_All->Write();
  Rechits_DT->Write();
  Rechits_CSC->Write();
  Rechits_RPC->Write();
  Rechits_GEM->Write();
  outputfile->cd();

  AllHits->cd();
  Rechits_All_Hits->Write();
  Rechits_DT_Hits->Write();
  Rechits_CSC_Hits->Write();

  Rechits_All_PhiHits->Write();
  Rechits_All_EtaHits->Write();
  Rechits_DT_Hits_EtaSL->Write();
  Rechits_DT_Hits_PhiSL->Write();
  outputfile->cd();


  for(int i=0; i<50; ++i) {

    int num1  = Rechits_All->GetBinContent(i+1); 
    int num2  = Rechits_DT->GetBinContent(i+1); 
    int num3  = Rechits_CSC->GetBinContent(i+1); 
    int num4  = Rechits_RPC->GetBinContent(i+1); 
    int num5  = Rechits_GEM->GetBinContent(i+1); 

    int denom = Muon_All->GetBinContent(i+1);

    double ave1=0.0, ave2=0.0, ave3=0.0, ave4=0.0, ave5=0.0;
    double err1=0.0, err2=0.0, err3=0.0, err4=0.0, err5=0.0;

    if(denom>0) {
      ave1 = 1.0*num1/denom; err1 = sqrt(num1)/denom;
      ave2 = 1.0*num2/denom; err2 = sqrt(num2)/denom;
      ave3 = 1.0*num3/denom; err3 = sqrt(num3)/denom;
      ave4 = 1.0*num4/denom; err4 = sqrt(num4)/denom;
      ave5 = 1.0*num5/denom; err5 = sqrt(num5)/denom;
    }  
    Rechits_All_Eta_1D->SetBinContent(i+1, ave1);
    Rechits_DT_Eta_1D->SetBinContent(i+1,  ave2);
    Rechits_CSC_Eta_1D->SetBinContent(i+1, ave3);
    Rechits_RPC_Eta_1D->SetBinContent(i+1, ave4);
    Rechits_GEM_Eta_1D->SetBinContent(i+1, ave5);

    Rechits_All_Eta_1D->SetBinError(i+1, err1);
    Rechits_DT_Eta_1D->SetBinError(i+1,  err2);
    Rechits_CSC_Eta_1D->SetBinError(i+1, err3);
    Rechits_RPC_Eta_1D->SetBinError(i+1, err4);
    Rechits_GEM_Eta_1D->SetBinError(i+1, err5);
  }

  SegmentsAndHits->cd();
  Rechits_All_Eta_1D->Write();
  Rechits_DT_Eta_1D->Write();
  Rechits_CSC_Eta_1D->Write();
  Rechits_RPC_Eta_1D->Write();
  Rechits_GEM_Eta_1D->Write();

  Rechits_All_Eta_2D->Write();
  Rechits_DT_Eta_2D->Write();
  Rechits_CSC_Eta_2D->Write();
  Rechits_RPC_Eta_2D->Write();
  Rechits_GEM_Eta_2D->Write();

  Rechits_All_Eta_2D_MyReBin->Write();
  Rechits_DT_Eta_2D_MyReBin->Write();
  Rechits_CSC_Eta_2D_MyReBin->Write();
  Rechits_RPC_Eta_2D_MyReBin->Write();
  Rechits_GEM_Eta_2D_MyReBin->Write();

  Rechits_All_Eta_2D_Reduced->Write();
  Rechits_DT_Eta_2D_Reduced->Write();
  Rechits_CSC_Eta_2D_Reduced->Write();
  Rechits_RPC_Eta_2D_Reduced->Write();
  Rechits_GEM_Eta_2D_Reduced->Write();

  Rechits_Seg_Eta_2D_Reduced->Write();
  Rechits_Hit_Eta_2D_Reduced->Write();

  outputfile->cd();

  for(int i=0; i<50; ++i) {

    int num1  = Rechits_All_Hits->GetBinContent(i+1); 
    int num2  = Rechits_DT_Hits->GetBinContent(i+1); 
    int num3  = Rechits_CSC_Hits->GetBinContent(i+1); 

    int num4  = Rechits_All_PhiHits->GetBinContent(i+1); 
    int num5  = Rechits_All_EtaHits->GetBinContent(i+1); 
    int num6  = Rechits_DT_Hits_EtaSL->GetBinContent(i+1); 
    int num7  = Rechits_DT_Hits_PhiSL->GetBinContent(i+1); 

    int denom = Muon_All->GetBinContent(i+1);

    double ave1=0.0, ave2=0.0, ave3=0.0, ave4=0.0, ave5 = 0.0, ave6 =0.0, ave7 = 0.0;
    double err1=0.0, err2=0.0, err3=0.0, err4=0.0, err5 = 0.0, err6 =0.0, err7 = 0.0;

    if(denom>0) {
      ave1 = 1.0*num1/denom; err1 = sqrt(num1)/denom;
      ave2 = 1.0*num2/denom; err2 = sqrt(num2)/denom;
      ave3 = 1.0*num3/denom; err3 = sqrt(num3)/denom;
      ave4 = 1.0*num4/denom; err4 = sqrt(num4)/denom;
      ave5 = 1.0*num5/denom; err5 = sqrt(num5)/denom;
      ave6 = 1.0*num6/denom; err6 = sqrt(num6)/denom;
      ave7 = 1.0*num7/denom; err7 = sqrt(num7)/denom;
    }  
    Rechits_All_Hits_Eta_1D->SetBinContent(i+1, ave1);
    Rechits_DT_Hits_Eta_1D->SetBinContent(i+1,  ave2);
    Rechits_CSC_Hits_Eta_1D->SetBinContent(i+1, ave3);

    Rechits_All_Hits_Eta_1D->SetBinError(i+1, err1);
    Rechits_DT_Hits_Eta_1D->SetBinError(i+1,  err2);
    Rechits_CSC_Hits_Eta_1D->SetBinError(i+1, err3);

    Rechits_All_PhiHits_Eta_1D->SetBinContent(i+1,  ave4);
    Rechits_All_EtaHits_Eta_1D->SetBinContent(i+1,  ave5);
    Rechits_DT_Hits_EtaSL_Eta_1D->SetBinContent(i+1,  ave6);
    Rechits_DT_Hits_PhiSL_Eta_1D->SetBinContent(i+1,  ave7);

    Rechits_All_PhiHits_Eta_1D->SetBinError(i+1,  err4);
    Rechits_All_EtaHits_Eta_1D->SetBinError(i+1,  err5);
    Rechits_DT_Hits_EtaSL_Eta_1D->SetBinError(i+1,  err6);
    Rechits_DT_Hits_PhiSL_Eta_1D->SetBinError(i+1,  err7);
  }

  AllHits->cd();
  Rechits_All_Hits_Eta_1D->Write();
  Rechits_DT_Hits_Eta_1D->Write();
  Rechits_CSC_Hits_Eta_1D->Write();
  Rechits_GEM_Hits_Eta_1D->Write();

  Rechits_All_PhiHits_Eta_1D->Write();
  Rechits_All_EtaHits_Eta_1D->Write();
  Rechits_DT_Hits_EtaSL_Eta_1D->Write();
  Rechits_DT_Hits_PhiSL_Eta_1D->Write();
  outputfile->cd();


}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyStandAloneMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // switch this off in CMSSW_6_2_0_SLHCX
  // also not necessary in CMSSW_7_X_Y if the Muon POG stuff is not used 
  // theService->update(iSetup);

  edm::Handle<reco::TrackCollection> staTracks;
  iEvent.getByLabel(STAMuLabel, staTracks);

  edm::Handle<reco::TrackCollection> glbTracks;
  iEvent.getByLabel(GLBMuLabel, glbTracks);

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByLabel("muons", recoMuons);

  edm::Handle<reco::VertexCollection> recoVertices;
  iEvent.getByLabel("offlinePrimaryVertices", recoVertices);

  reco::TrackCollection::const_iterator staTrack;
  if(techDebug) std::cout<<"Reconstructed STA Muon Tracks: " <<staTracks->size()<< std::endl;
  reco::TrackCollection::const_iterator glbTrack;
  if(techDebug) std::cout<<"Reconstructed GLO Muon Tracks: " <<glbTracks->size()<<std::endl;
  reco::MuonCollection::const_iterator recoMuon;
  if(techDebug) std::cout<<"Reconstructed Muons: "<<recoMuons->size()<<std::endl;
  reco::VertexCollection::const_iterator recoVertex;
  if(techDebug) std::cout<<"Reconstructed Vertices: "<<recoVertices->size()<<std::endl;

  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  int stamuons = 0;
  int stamuons_eta_hits = 0;

  double muon_me_properties[2][3]; // 2 Matrix Element Muons // Save for each muon PT, ETA, PHI
  double muon_tk_properties[2][3]; // 2 SimTrack Muons       // Save for each muon PT, ETA, PHI

  // ===================================================================
  // ===      MC Truth Matching                                      === 
  // ===================================================================
  // Get the Muon SimTracks / Muon GenParticles here && match to Z-decay
  // Do a pt match between SimTrack and Muon GenParticles from Z-decay
  // Then do dR match between SimTrack and STA / Glb Muon
  if(mcTruthMatch) {
    // 1) GEN-LEVEL :: Save Matrix Element Muon Properties (PT, ETA, PHI)
    edm::Handle<reco::GenParticleCollection>      genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    bool ZbosonFound = 0;
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      if(ZbosonFound) continue;
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if (g->pdgId() == 23) {
	if(genDebug) std::cout<<"ME Z boson: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = "<<std::setw(9)<<g->phi();
	if(genDebug) std::cout<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	if(genDebug) std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	reco::GenParticle::daughters d = g->daughterRefVector();
	int muon_me_count = 0;
	for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
	  if((*it_d)->status()==3 && fabs((*it_d)->pdgId())==13) {
	    b = &(*(*it_d));
	    if(genDebug) std::cout<<"ME Muon: id = "<<std::setw(2)<<b->pdgId()<<" | st = "<<std::setw(2)<<b->status()<<" | eta = "<<std::setw(9)<<b->eta()<<" | phi = "<<std::setw(9)<<b->phi();
	    if(genDebug) std::cout<<" | pt = "<<std::setw(9)<<b->pt()<<" GeV/c | et = "<<std::setw(9)<<b->et()<<" GeV | Energy = "<<b->energy()<<" GeV ";
	    if(genDebug) std::cout<<" | transverse mass = "<<std::setw(9)<<b->mt()<<" GeV/c^2"<<std::endl;
	    muon_me_properties[muon_me_count][0] = b->pt();
	    muon_me_properties[muon_me_count][1] = b->eta();
	    muon_me_properties[muon_me_count][2] = b->phi();
	    ++muon_me_count;
	  }
	}
      }
    }
    // 2) SIM-LEVEL :: Match Muon SimTracks to GenParticle Muons
    int muon_tk_count = 0;
    bool mu1_me_matched = false, mu2_me_matched = false;
    std::vector<SimTrack> theSimTracks;
    edm::Handle<edm::SimTrackContainer> SimTk;
    iEvent.getByLabel("g4SimHits",SimTk);
    theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
    for (std::vector<SimTrack>::const_iterator iTrack = theSimTracks.begin(); iTrack != theSimTracks.end(); ++iTrack) {
      SimTrack simtrack = (*iTrack);
      // if(fabs(simtrack.type()) != 13 || simtrack.momentum().vect().perp() < 95.0) continue;
      if(fabs(simtrack.type()) != 13) continue;
      else {
	if(genDebug) std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<simtrack.momentum().pt()<<" eta = "<<simtrack.momentum().eta()<<" phi = "<<simtrack.momentum().phi();
	// Compare to first Muon
	double delta_rel_pt_mu1 = fabs(simtrack.momentum().pt()-muon_me_properties[0][0])/muon_me_properties[0][0];
	double delta_rel_pt_mu2 = fabs(simtrack.momentum().pt()-muon_me_properties[1][0])/muon_me_properties[1][0];
	double simtrack_eta = simtrack.momentum().eta();
	double simtrack_phi = simtrack.momentum().phi();
	if(genDebug) std::cout<<" Rel PT diff mu1 = "<<delta_rel_pt_mu1<<" Rel PT diff mu2 = "<<delta_rel_pt_mu2;
	if(!mu1_me_matched && delta_rel_pt_mu1<delta_rel_pt) {
	  mu1_me_matched = true;
	  muon_tk_properties[muon_tk_count][0] = simtrack.momentum().pt();
	  muon_tk_properties[muon_tk_count][1] = simtrack.momentum().eta();
	  muon_tk_properties[muon_tk_count][2] = simtrack.momentum().phi();
	  ++muon_tk_count;
	  if(genDebug) std::cout<<" ---> MATCHED"<<std::endl;
	}
	else if(!mu2_me_matched && delta_rel_pt_mu2<delta_rel_pt) {
	  mu2_me_matched = true;
	  muon_tk_properties[muon_tk_count][0] = simtrack.momentum().pt();
	  muon_tk_properties[muon_tk_count][1] = simtrack.momentum().eta();
	  muon_tk_properties[muon_tk_count][2] = simtrack.momentum().phi();
	  ++muon_tk_count;
	  if(genDebug) std::cout<<" ---> MATCHED"<<std::endl;
	}
	// I found some cases of severe bremsstrahlung leading to Pt mismatch
	// Try Delta R matching for these
	else if(!mu1_me_matched && deltaR(simtrack_eta, simtrack_phi, muon_me_properties[0][1], muon_me_properties[0][2]) < delta_R) {
	  mu1_me_matched = true;
	  muon_tk_properties[muon_tk_count][0] = simtrack.momentum().pt();
	  muon_tk_properties[muon_tk_count][1] = simtrack.momentum().eta();
	  muon_tk_properties[muon_tk_count][2] = simtrack.momentum().phi();
	  ++muon_tk_count;
	  if(genDebug) std::cout<<" ---> MATCHED"<<std::endl;
	}
	else if(!mu2_me_matched && deltaR(simtrack_eta, simtrack_phi, muon_me_properties[1][1], muon_me_properties[1][2]) < delta_R) {
	  mu2_me_matched = true;
	  muon_tk_properties[muon_tk_count][0] = simtrack.momentum().pt();
	  muon_tk_properties[muon_tk_count][1] = simtrack.momentum().eta();
	  muon_tk_properties[muon_tk_count][2] = simtrack.momentum().phi();
	  ++muon_tk_count;
	  if(genDebug) std::cout<<" ---> MATCHED"<<std::endl;
	}
	else {
	  if(genDebug) std::cout<<" ---> NOT MATCHED"<<std::endl;
	}
      }
    }
  }


  bool mu1_tk_matched = false, mu2_tk_matched = false;
  // ===================================================================
  // ===      Loop over Stand Alone Muon TRACK Collection            === 
  // ===        ---> and use the MC Truth Matching                   ===
  // ===        ---> or Tight Muon ID and pt cut                     ===
  // ===================================================================
  // Quick Loop just for Printout
  if(physDebug) {
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack) {
      reco::TransientTrack track(*staTrack,&*theMGField,theTrackingGeometry); 
      std::cout<<"TRACKS :: Stand Alone Muon :: ";
      std::cout<<" Track pt: "<<std::setw(9)<<track.impactPointTSCP().momentum().perp();
      std::cout<<" eta: "     <<std::setw(9)<<track.impactPointTSCP().momentum().eta();
      std::cout<<" phi: "     <<std::setw(9)<<track.impactPointTSCP().momentum().phi();
      std::cout<<" p: "       <<std::setw(9)<<track.impactPointTSCP().momentum().mag();
      std::cout<<" chi2: "    <<std::setw(9)<<track.chi2();
      std::cout<<" with "     <<std::setw(2)<<staTrack->recHitsSize()<<" rechits"<<std::endl;
    }
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    for (recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); ++recoMuon) {
      if(recoMuon->isStandAloneMuon()) {
	std::cout<<" MUONS :: Stand Alone Muon :: Muon pt = "<<std::setw(9)<<recoMuon->pt()<<" eta = "<<std::setw(6)<<recoMuon->eta()<<" phi= "<<std::setw(6)<<recoMuon->phi();
	if(recoMuon->time().direction()<0)  std::cout<<" direction = "<<("OutsideIn")<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" +/- "<<recoMuon->time().timeAtIpInOutErr<<" ns";
	if(recoMuon->time().direction()>0)  std::cout<<" direction = "<<("InsideOut")<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" +/- "<<recoMuon->time().timeAtIpInOutErr<<" ns";
	if(recoMuon->time().direction()==0) std::cout<<" direction = "<<("Undefined")<<" time at IP = "<<("Undefined")<<" ns";
	std::cout<<std::endl;
      }
    }
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    for (recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); ++recoMuon) {
      if(recoMuon->isGlobalMuon()) {
	std::cout<<" MUONS ::      Global Muon :: Muon pt = "<<std::setw(9)<<recoMuon->pt()<<" eta = "<<std::setw(6)<<recoMuon->eta()<<" phi= "<<std::setw(6)<<recoMuon->phi();
	if(recoMuon->time().direction()<0)  std::cout<<" direction = "<<("OutsideIn")<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" +/- "<<recoMuon->time().timeAtIpInOutErr<<" ns";
	if(recoMuon->time().direction()>0)  std::cout<<" direction = "<<("InsideOut")<<" time at IP = "<<recoMuon->time().timeAtIpInOut<<" +/- "<<recoMuon->time().timeAtIpInOutErr<<" ns";
	if(recoMuon->time().direction()==0) std::cout<<" direction = "<<("Undefined")<<" time at IP = "<<("Undefined")<<" ns";
	std::cout<<std::endl;
      }
    }
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    std::cout<<"\n\n\n"<<std::endl;
  }
  // Starts the real stuff
  for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack) {
    reco::TransientTrack track(*staTrack,&*theMGField,theTrackingGeometry); 

    if(physDebug) {
      std::cout<<"--------------------------------------------------------------"<<std::endl;
      std::cout<<"--- Stand Alone Muon Track :: ";
      std::cout<<" pT: "  <<std::setw(9)<<track.impactPointTSCP().momentum().perp();
      std::cout<<" eta: " <<std::setw(9)<<track.impactPointTSCP().momentum().eta();
      std::cout<<" phi: " <<std::setw(9)<<track.impactPointTSCP().momentum().phi();
      std::cout<<" p: "   <<std::setw(9)<<track.impactPointTSCP().momentum().mag();
      std::cout<<" chi2: "<<std::setw(9)<<track.chi2();
      std::cout<<" with " <<std::setw(2)<<staTrack->recHitsSize()<<" rechits"<<std::endl;
      std::cout<<"--------------------------------------------------------------"<<std::endl;
    }
    double track_eta = track.impactPointTSCP().momentum().eta();
    double track_phi = track.impactPointTSCP().momentum().phi();

    bool matched = false;
    if(mcTruthMatch) { // do matching to SimTrack Muons
      if(physDebug) {
	std::cout<<" Track eta = "<<track_eta<<" phi = "<<track_phi<<" SimTrk eta = "<<muon_tk_properties[0][1]<<" phi = "<<muon_tk_properties[0][2];
	std::cout<<" => delta R (mu1) = "<<deltaR(track_eta, track_phi, muon_tk_properties[0][1], muon_tk_properties[0][2])<<std::endl;
	std::cout<<" Track eta = "<<track_eta<<" phi = "<<track_phi<<" SimTrk eta = "<<muon_tk_properties[1][1]<<" phi = "<<muon_tk_properties[1][2];
	std::cout<<" => delta R (mu2) = "<<deltaR(track_eta, track_phi, muon_tk_properties[1][1], muon_tk_properties[1][2])<<std::endl;
      }
      if(!mu1_tk_matched && deltaR(track_eta, track_phi, muon_tk_properties[0][1], muon_tk_properties[0][2]) < delta_R) {
	mu1_tk_matched = true; matched = true;
      }
      else if(!mu2_tk_matched && deltaR(track_eta, track_phi, muon_tk_properties[1][1], muon_tk_properties[1][2]) < delta_R ) { 
	mu2_tk_matched = true; matched = true;
      }
      else {
	matched = false;
      }
    }
    else { // no matching, use Tight ID and pt > 10 GeV/c cut
      // Loop over reco::Muon collection and do dR matching 
      // (I do not have time now to look up how the track reference works)
      const reco::Muon * myMuon = &(*recoMuons->begin());
      bool foundMyMuon = false;
      for (recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); ++recoMuon) {
	if(deltaR(track_eta, track_phi, recoMuon->eta(), recoMuon->phi())< delta_R && recoMuon->isStandAloneMuon()) {
	  myMuon = &(*recoMuon); foundMyMuon = true; 
	  if(physDebug) {
	    std::cout<<" Reco Muon :: pt = "<<std::setw(9)<<recoMuon->pt()<<" eta = "<<std::setw(5)<<recoMuon->eta()<<" phi= "<<std::setw(5)<<recoMuon->phi();
	    std::cout<<" is StandAloneMuon = "<<recoMuon->isStandAloneMuon()<<" is GlobalMuon = "<<recoMuon->isGlobalMuon();
	    std::cout<<" is TrackerMuon = "<<recoMuon->isTrackerMuon()<<" is ParticleFlowMuon = "<<recoMuon->isPFMuon()<<std::endl;
	  }
	  if(physDebug) std::cout<<" ---> this Reco (STA) Muon is dR matched to the StandAlone Muon"<<std::endl;
	}
	else {
	  if(techDebug) {
	    std::cout<<" Reco Muon :: pt = "<<std::setw(9)<<recoMuon->pt()<<" eta = "<<std::setw(5)<<recoMuon->eta()<<" phi= "<<std::setw(5)<<recoMuon->phi();
	    std::cout<<" is StandAloneMuon = "<<recoMuon->isStandAloneMuon()<<" is GlobalMuon = "<<recoMuon->isGlobalMuon();
	    std::cout<<" is TrackerMuon = "<<recoMuon->isTrackerMuon()<<" is ParticleFlowMuon = "<<recoMuon->isPFMuon()<<std::endl;
	  }
	}
      }
      if(foundMyMuon==false) { std::cout<<" ---> No Reco (STA) Muon found that is dR matched to the StandAlone Muon"<<std::endl; }
      if(foundMyMuon==true) {
	if(physDebug) {
	  std::cout<<" ---> this Reco Muon is a Tight Muon = "<<muon::isTightMuon(*myMuon, *(recoVertices->begin()));
	  std::cout<<" and StandAlone Pt = "<<track.impactPointTSCP().momentum().perp()<<std::endl;
	  std::cout<<"--------------------------------------------------------------"<<std::endl;
	}
      }
      if (wantTightId && foundMyMuon && muon::isTightMuon(*myMuon, *(recoVertices->begin())) && track.impactPointTSCP().momentum().perp() > 10) {
	matched = true;
      }
      else if(wantLooseId && foundMyMuon && muon::isLooseMuon(*myMuon) && track.impactPointTSCP().momentum().perp() > 10) {
	matched = true;
      }
      else {
	matched = false;
      }
    }



    // ===================================================================
    // ===      Now Matching or Identification is done                 ===     
    // ===      from here continue with the real program               ===
    // ===================================================================

    if(!matched) continue;
    if(physDebug) {
      std::cout<<"      ---   RECHITS of the MATCHED MUON"<<std::endl;
      std::cout<<"      --------------------------------------------------------------"<<std::endl;
    }

    StandAloneMuon_PT->Fill(track.impactPointTSCP().momentum().perp());
    double phi = 0.0;
    if(track.impactPointTSCP().momentum().phi()>3.0544) phi = track.impactPointTSCP().momentum().phi()-2*3.141592; 
    else phi = track.impactPointTSCP().momentum().phi();
    StandAloneMuon_Phi->Fill(phi);
    StandAloneMuon_Eta->Fill(track.impactPointTSCP().momentum().eta());
    ++stamuons; // count how many stand alone muons in this event

    trackingRecHit_iterator rhbegin = staTrack->recHitsBegin();
    trackingRecHit_iterator rhend = staTrack->recHitsEnd();

    double eta = track.impactPointTSCP().momentum().eta();
    Muon_All->Fill(eta);
    Muon_All_MyReBin->Fill(eta);
    Muon_All_Reduced->Fill(eta);

    int All_Rechits   = 0;
    int RPC_Rechits   = 0;
    int GEM_Rechits   = 0;
    int CSC_Rechits   = 0;
    int DT_Rechits    = 0;

    int muonhits = 0;
    int station_fired_sel[8] = {0,0,0,0,0,0,0,0};
    bool gemstation_1_fired = false, gemstation_2_fired =false, gemstation_3_fired = false;

    for(trackingRecHit_iterator recHit = rhbegin; recHit != rhend; ++recHit) {
      const GeomDet* geomDet = theTrackingGeometry->idToDet((*recHit)->geographicalId());
      double r = geomDet->surface().position().perp();
      double z = geomDet->toGlobal((*recHit)->localPosition()).z();

      DetId detid = DetId((*recHit)->geographicalId());

      // Muon Hits
      // ---------
      if(detid.det()==DetId::Muon && (detid.subdetId()== MuonSubdetId::RPC || detid.subdetId()== MuonSubdetId::CSC || detid.subdetId()== MuonSubdetId::DT || detid.subdetId()== MuonSubdetId::GEM)) {
	++muonhits;
      }

      // RPC Hits
      // --------
      if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::RPC) {
	++All_Rechits; ++RPC_Rechits; Rechits_All->Fill(eta); Rechits_RPC->Fill(eta);
	Rechits_All_Hits->Fill(eta);
	Rechits_All_EtaHits->Fill(eta); Rechits_All_PhiHits->Fill(eta);
	DetId idRivHit = (*recHit)->geographicalId();
	RPCDetId rollId(idRivHit.rawId());
	if(physDebug) std::cout<<"      RPC RecHit at "<<"r: "<<std::setw(9)<<r<<" cm"<<" z: "<<std::setw(9)<<z<<" cm in DetId "<<idRivHit.rawId()<<" = "<<rollId<<std::endl;
	// which RPC stations fired for a muon in 1.7 < eta < 2.5
	if(fabs(track.impactPointTSCP().momentum().eta()) > 1.7 && fabs(track.impactPointTSCP().momentum().eta()) < 2.5) {
	  station_fired_sel[3+rollId.station()] = 1;
	}
      }

      // GEM Hits
      // --------
      if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::GEM) {

	DetId idRivHit = (*recHit)->geographicalId();
	GEMDetId rollId(idRivHit.rawId());

	if(physDebug) std::cout<<"      GEM RecHit at "<<"r: "<<std::setw(9)<<r<<" cm"<<" z: "<<std::setw(9)<<z<<" cm in DetId "<<idRivHit.rawId()<<" = "<<rollId<<std::endl;
	// if(physDebug) std::cout<<"GEM DetId = "<<rollId<<std::endl;

	// count only one hit per station and skip station 2
	// GE21 short will disappear in the future 
	// GE11 = Station 1 | GE21 short = Station 2 | GE21 long = Station 3
	if(rollId.station() == 1 && !gemstation_1_fired) {
	  ++All_Rechits; ++GEM_Rechits; Rechits_All->Fill(eta); Rechits_GEM->Fill(eta); gemstation_1_fired = true;
	  //std::cout<<"Specific GEM RecHit --- Station 1 --- at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm | GEM Det Id = "<<rollId<<std::endl;
	}
	if(rollId.station() == 2 && !gemstation_2_fired) {
	  gemstation_2_fired = true; // skip station 2 = GE21 short
	  // std::cout<<"Specific GEM RecHit --- Station 2 --- at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm | GEM Det Id = "<<rollId<<std::endl;
	} 
	if(rollId.station() == 3 && !gemstation_3_fired) {
	  ++All_Rechits; ++GEM_Rechits; Rechits_All->Fill(eta); Rechits_GEM->Fill(eta); gemstation_3_fired = true;
	  // std::cout<<"Specific GEM RecHit --- Station 3 --- at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm | GEM Det Id = "<<rollId<<std::endl;
	}
	// now count all hits (except station 2)  
	if(rollId.station() != 2) {
	  Rechits_All_Hits->Fill(eta);
	  Rechits_GEM_Hits->Fill(eta);    Rechits_All_Hits->Fill(eta);
	  Rechits_All_EtaHits->Fill(eta); Rechits_All_PhiHits->Fill(eta);
	}
	// which GEM stations fired for a muon in 1.7 < eta < 2.5
	if(fabs(track.impactPointTSCP().momentum().eta()) > 1.7 && fabs(track.impactPointTSCP().momentum().eta()) < 2.5) {
	  station_fired_sel[3+rollId.station()] = 1;
	}
      }

      // DT Hits 
      // -------
      if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::DT) {
	++All_Rechits; ++DT_Rechits; Rechits_All->Fill(eta); Rechits_DT->Fill(eta);
	// DetId & Printout
	DetId idRivHit = (*recHit)->geographicalId();
	DTChamberId    chamberId(idRivHit.rawId()); DTSuperLayerId sulayerId(idRivHit.rawId()); DTLayerId layerId(idRivHit.rawId()); DTWireId wireId(idRivHit.rawId());
	if(physDebug) std::cout<<"      DT RecHit at "<<"r: "<<std::setw(9)<<r<<" cm"<<" z: "<<std::setw(9)<<z<<" cm in DetId "<<chamberId.rawId()<<" = "<<chamberId<<std::endl;
	// This is actually a Segment ... Try to access different rechits of Segment
	if(techDebug) std::cout<<"DT Rec Hit in Raw ID = "<<idRivHit.rawId()<<" Chamber ID = "<<chamberId<<" SL ID = "<<sulayerId<<" Layer ID = "<<layerId<<" Wire ID = "<<wireId<<std::endl;
	std::vector< const TrackingRecHit * > DTRecHitsL1 = (*recHit)->recHits();	
	if(techDebug) std::cout<<"Number of constituting rechits = "<<DTRecHitsL1.size()<<std::endl;
	std::vector< const TrackingRecHit * >::const_iterator recHitl1;
	for(recHitl1 = DTRecHitsL1.begin(); recHitl1 != DTRecHitsL1.end(); ++recHitl1) {
	  DetId detidl1 = DetId((*recHitl1)->geographicalId());
	  DTChamberId    chamberIdL1(detidl1.rawId()); DTSuperLayerId sulayerIdL1(detidl1.rawId()); DTLayerId layerIdL1(detidl1.rawId()); DTWireId wireIdL1(detidl1.rawId());
	  if(techDebug) std::cout<<"     |--> DT Rec Hit in Raw ID = "<<detidl1.rawId()<<" Chamber ID = "<<chamberIdL1<<" SL ID = "<<sulayerIdL1<<" Layer ID = "<<layerIdL1<<" Wire ID = "<<wireIdL1<<std::endl;
	  std::vector< const TrackingRecHit * > DTRecHitsL2 = (*recHitl1)->recHits();
	  if(techDebug) std::cout<<"     |--> Number of constituting rechits = "<<DTRecHitsL2.size()<<std::endl;
	  std::vector< const TrackingRecHit * >::const_iterator recHitl2;
	  for(recHitl2 = DTRecHitsL2.begin(); recHitl2 != DTRecHitsL2.end(); ++recHitl2) {
	    DetId detidl2 = DetId((*recHitl2)->geographicalId());
	    DTChamberId    chamberIdL2(detidl2.rawId()); DTSuperLayerId sulayerIdL2(detidl2.rawId()); DTLayerId layerIdL2(detidl2.rawId()); DTWireId wireIdL2(detidl2.rawId());
	    if(techDebug) std::cout<<"          |--> DT Rec Hit in Raw ID = "<<detidl2.rawId()<<" Chamber ID = "<<chamberIdL2<<" SL ID = "<<sulayerIdL2<<" Layer ID = "<<layerIdL2<<" Wire ID = "<<wireId<<std::endl;
	    // std::vector< const TrackingRecHit * > DTRecHitsL3 = (*recHitl2)->recHits();
	    // std::cout<<"          |--> Number of constituting rechits = "<<DTRecHitsL3.size()<<std::endl;
	    if(techDebug) std::cout<<"Super Layer Number = "<<sulayerIdL2.superLayer();
	    if(sulayerIdL2.superLayer()==0)      { if(techDebug) { std::cout<<" chamber"<<std::endl; } }
	    else if(sulayerIdL2.superLayer()==2) { Rechits_DT_Hits_EtaSL->Fill(eta); Rechits_All_EtaHits->Fill(eta); if(techDebug) { std::cout<<" theta-SL"<<std::endl; } }
	    else                                 { Rechits_DT_Hits_PhiSL->Fill(eta); Rechits_All_PhiHits->Fill(eta); if(techDebug) { std::cout<<" phi-SL"<<std::endl;   } }
	    Rechits_DT_Hits->Fill(eta);
            Rechits_All_Hits->Fill(eta);
	  }
	}
      }

      // CSC Hits
      // --------
      if(detid.det()==DetId::Muon && detid.subdetId()== MuonSubdetId::CSC) {
	++All_Rechits; ++CSC_Rechits; Rechits_All->Fill(eta); Rechits_CSC->Fill(eta);
	// DetId & Printout
	DetId idRivHit = (*recHit)->geographicalId();
        CSCDetId chamberId(idRivHit.rawId());
	if(physDebug) std::cout<<"      CSC RecHit at "<<"r: "<<std::setw(9)<<r<<" cm"<<" z: "<<std::setw(9)<<z<<" cm in DetId "<<chamberId.rawId()<<" = "<<chamberId<<std::endl;
        // This is actually a Segment ... Try to access different rechits of Segment
	std::vector< const TrackingRecHit * > CSCRecHitsL1 = (*recHit)->recHits();	
	if(techDebug) std::cout<<"Number of constituting rechits = "<<CSCRecHitsL1.size()<<std::endl;
	std::vector< const TrackingRecHit * >::const_iterator recHitl1;
	for(recHitl1 = CSCRecHitsL1.begin(); recHitl1 != CSCRecHitsL1.end(); ++recHitl1) {
	  Rechits_CSC_Hits->Fill(eta);
	  Rechits_All_Hits->Fill(eta);
	  Rechits_All_EtaHits->Fill(eta); Rechits_All_PhiHits->Fill(eta);
	}
	if(fabs(track.impactPointTSCP().momentum().eta()) > 1.7 && fabs(track.impactPointTSCP().momentum().eta()) < 2.5) {
          station_fired_sel[chamberId.station()-1] = 1;
        }
      }

      // Tracker Hits
      // ------------
      if(detid.det()==DetId::Tracker) {
	if(physDebug) std::cout<<"      Tracker RecHit at "<<"r: "<< r <<" cm"<<" z: "<<z<<" cm"<<std::endl;
      }
    } // end loop rechits
    if(physDebug) std::cout<<"\n"<<std::endl;


    Rechits_All_Eta_2D->Fill(eta, All_Rechits);
    Rechits_RPC_Eta_2D->Fill(eta, RPC_Rechits);
    Rechits_GEM_Eta_2D->Fill(eta, GEM_Rechits);
    Rechits_CSC_Eta_2D->Fill(eta, CSC_Rechits);
    Rechits_DT_Eta_2D->Fill(eta, DT_Rechits);

    Rechits_All_Eta_2D_MyReBin->Fill(eta, All_Rechits);
    Rechits_RPC_Eta_2D_MyReBin->Fill(eta, RPC_Rechits);
    Rechits_GEM_Eta_2D_MyReBin->Fill(eta, GEM_Rechits);
    Rechits_CSC_Eta_2D_MyReBin->Fill(eta, CSC_Rechits);
    Rechits_DT_Eta_2D_MyReBin->Fill(eta, DT_Rechits);

    Rechits_All_Eta_2D_Reduced->Fill(eta, All_Rechits);
    Rechits_Seg_Eta_2D_Reduced->Fill(eta, DT_Rechits+CSC_Rechits);
    Rechits_Hit_Eta_2D_Reduced->Fill(eta, RPC_Rechits+GEM_Rechits);
    Rechits_RPC_Eta_2D_Reduced->Fill(eta, RPC_Rechits);
    Rechits_GEM_Eta_2D_Reduced->Fill(eta, GEM_Rechits);
    Rechits_CSC_Eta_2D_Reduced->Fill(eta, CSC_Rechits);
    Rechits_DT_Eta_2D_Reduced->Fill(eta, DT_Rechits);



    StandAloneMuon_Hit->Fill(muonhits);
    if(muonhits < 3 && fabs(track.impactPointTSCP().momentum().eta()) > 1.7 && fabs(track.impactPointTSCP().momentum().eta()) < 2.5) {
      ++stamuons_eta_hits; // count how many stand alone muons in this phasespace in this event
      StandAloneMuon_1p8To2p5_2Hit_PT->Fill(track.impactPointTSCP().momentum().perp());
      phi = 0.0;
      if(track.impactPointTSCP().momentum().phi()>3.0544) phi = track.impactPointTSCP().momentum().phi()-2*3.141592;
      else phi = track.impactPointTSCP().momentum().phi();
      StandAloneMuon_1p8To2p5_2Hit_Phi->Fill(phi);
      StandAloneMuon_1p8To2p5_2Hit_Eta->Fill(track.impactPointTSCP().momentum().eta());
      StandAloneMuon_1p8To2p5_2Hit_Hit->Fill(muonhits);
      std::cout<<"Muon has only 2 hits && Forward :: Stations fired = [ME1-4,RE1-4] = [";
      for(int i=0; i<8; ++i) {
	if(station_fired_sel[i]==1) {
	  // this way the labels on the x-axis are in a random order :-(
	  // StandAloneMuon_1p8To2p5_2Hit_Sta->Fill((station[i]).c_str(),1); 
	  StandAloneMuon_1p8To2p5_2Hit_Sta->Fill(i+1,1);

	}
	std::cout<<station_fired_sel[i]<<", "; 
      }
      std::cout<<"]"<<std::endl;
    }
    if(muonhits < 2 && fabs(track.impactPointTSCP().momentum().eta()) > 1.7 && fabs(track.impactPointTSCP().momentum().eta()) < 2.5) {
      std::cout<<"??? Muon Reconstructed with only 1 hit ???"<<std::endl;
    }
  } // end loop stand alone muons

  Event_NumStaMuons->Fill(stamuons);
  Event_NumStaMuons_1p8To2p5_2Hit->Fill(stamuons_eta_hits);

  if(physDebug) {std::cout<<"=============================================================="<<std::endl; std::cout<<"\n\n\n"<<std::endl;}

  // Code from: DQMOffline/Muon/src/SegmentTrackAnalyzer.cc
  // GLOBAL MUONS
  for (glbTrack = glbTracks->begin(); glbTrack!=glbTracks->end(); ++glbTrack) {
   
    reco::TransientTrack track(*glbTrack,&*theMGField,theTrackingGeometry); 
    if(physDebug) {
      std::cout<<"--------------------------------------------------------------"<<std::endl;
      std::cout<<"--- Global Muon Track :: ";
      std::cout<<" p: "<<track.impactPointTSCP().momentum().mag();
      std::cout<<" pT: "<<track.impactPointTSCP().momentum().perp();
      std::cout<<" eta: "<<track.impactPointTSCP().momentum().eta();
      std::cout<<" chi2: "<<track.chi2();
      std::cout<<" with "<<glbTrack->recHitsSize()<<" rechits"<<std::endl;
      std::cout<<"--------------------------------------------------------------"<<std::endl;
    }
    GlobalMuon_PT->Fill(track.impactPointTSCP().momentum().perp());
    /*
    MuonTransientTrackingRecHit::MuonRecHitContainer segments = theSegmentsAssociator->associate(iEvent, iSetup, *glbTrack );
    if(techDebug) std::cout<<"Global Track :: segments associated = "<<segments.size()<<std::endl;

    for (MuonTransientTrackingRecHit::MuonRecHitContainer::const_iterator segment=segments.begin(); segment!=segments.end(); segment++) {
      DetId id = (*segment)->geographicalId();
      // hits from DT segments
      if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
        const DTRecSegment4D *seg4D = dynamic_cast<const DTRecSegment4D*>((*segment)->hit());
        if((*seg4D).hasPhi())
	  if(techDebug) std::cout<<"DT Segment :: phi hits = "<<(*seg4D).phiSegment()->specificRecHits().size()<<std::endl;
        if((*seg4D).hasZed())
	  if(techDebug) std::cout<<"DT Segment :: z hits = "<<(*seg4D).zSegment()->specificRecHits().size()<<std::endl;
      }
      // hits from CSC segments
      if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) {
	if(techDebug) std::cout<<"CSC Segment :: 2D hits = "<<(*segment)->recHits().size()<<std::endl;
      }
    }
    */
  }
  /*
  // Code from: DQMOffline/Muon/src/SegmentTrackAnalyzer.cc
  // STAND ALONE MUONS
  for (staTrack = staTracks->begin(); staTrack!=staTracks->end(); ++staTrack) {
    
    MuonTransientTrackingRecHit::MuonRecHitContainer segments = theSegmentsAssociator->associate(iEvent, iSetup, *staTrack );
    if(techDebug) std::cout<<"Stand Alone Track :: segments associated = "<<segments.size()<<std::endl;

    for (MuonTransientTrackingRecHit::MuonRecHitContainer::const_iterator segment=segments.begin(); segment!=segments.end(); segment++) {
      DetId id = (*segment)->geographicalId();
      // hits from DT segments
      if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
        const DTRecSegment4D *seg4D = dynamic_cast<const DTRecSegment4D*>((*segment)->hit());
        if((*seg4D).hasPhi())
	  if(techDebug) std::cout<<"DT Segment :: phi hits = "<<(*seg4D).phiSegment()->specificRecHits().size()<<std::endl;
        if((*seg4D).hasZed())
	  if(techDebug) std::cout<<"DT Segment :: z hits = "<<(*seg4D).zSegment()->specificRecHits().size()<<std::endl;
      }
      // hits from CSC segments
      if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) {
	std::cout<<"CSC Segment :: 2D hits = "<<(*segment)->recHits().size()<<std::endl;
      }
    }
  }
  if(physDebug) std::cout<<"=============================================================="<<std::endl;
  */

  // Special :: Analyze the DT Segments of W+0/S11
  /*
  edm::Handle<DTRecSegment4DCollection> dtSegmentCollection;
  iEvent.getByLabel("dt4DSegments", dtSegmentCollection);
  if(physDebug) std::cout<<"Size of DT 4D Segment Collection :: "<<dtSegmentCollection->size()<<std::endl;
  DTRecSegment4DCollection::const_iterator segmentDT;
  for (segmentDT = dtSegmentCollection->begin(); segmentDT != dtSegmentCollection->end(); ++segmentDT){
    // const GeomDet* geomDet = theTrackingGeometry->idToDet((*segmentDT).geographicalId());
    DetId idSeg = (*segmentDT).geographicalId();
    DTChamberId chamberId(idSeg.rawId());
    if(physDebug && (chamberId.wheel()==0 && chamberId.sector()==11)) {
      std::cout << " DT 4D Segment in DetId "<<idSeg.rawId()<<" = "<<chamberId;
      // std::cout << " global pos: " << geomDet->toGlobal((*segmentDT).localPosition()) << "global dir: "<<geomDet->toGlobal((*segmentDT).localDirection());
      std::cout << " local pos: " << (*segmentDT).localPosition() << "local dir: "<<(*segmentDT).localDirection()<< std::endl;

      // Searching for RecHits inside DT Segment:
      DTSuperLayerId sulayerId(idSeg.rawId());
      if(physDebug) std::cout<<"DT Rec Hit in Raw ID = "<<idSeg.rawId()<<" Chamber ID = "<<chamberId<<std::endl;
      std::vector< const TrackingRecHit * > DTRecHitsL1 = (*segmentDT).recHits();
      if(physDebug) std::cout<<"Number of constituting rechits = "<<DTRecHitsL1.size()<<std::endl;
      // First Hierarchical Layer :: Phi and Theta SuperLayers
      std::vector< const TrackingRecHit * >::const_iterator recHitl1;
      for(recHitl1 = DTRecHitsL1.begin(); recHitl1 != DTRecHitsL1.end(); ++recHitl1) {
	DetId detidl1 = DetId((*recHitl1)->geographicalId());
	DTChamberId    chamberIdL1(detidl1.rawId()); DTSuperLayerId sulayerIdL1(detidl1.rawId());
	if(physDebug) std::cout<<"     |--> DT Rec Hit in Raw ID = "<<detidl1.rawId()<<" SuperLayer ID = "<<sulayerIdL1<<std::endl;
	std::vector< const TrackingRecHit * > DTRecHitsL2 = (*recHitl1)->recHits();
	if(physDebug) std::cout<<"     |--> Number of constituting rechits = "<<DTRecHitsL2.size()<<std::endl;
	// Second Hierarchical Layer :: Layers (1-4) Theta SuperLayer and (1-8) for Combined Phi SuperLayer
	std::vector< const TrackingRecHit * >::const_iterator recHitl2;
	for(recHitl2 = DTRecHitsL2.begin(); recHitl2 != DTRecHitsL2.end(); ++recHitl2) {
	  DetId detidl2 = DetId((*recHitl2)->geographicalId());
	  DTChamberId    chamberIdL2(detidl2.rawId()); DTSuperLayerId sulayerIdL2(detidl2.rawId()); DTLayerId layerIdL2(detidl2.rawId());
	  if(physDebug) std::cout<<"          |--> DT Rec Hit in Raw ID = "<<detidl2.rawId()<<" Layer ID = "<<layerIdL2<<std::endl;
	  std::vector< const TrackingRecHit * > DTRecHitsL3 = (*recHitl2)->recHits();
	  // if(physDebug) std::cout<<"          |--> Number of constituting rechits = "<<DTRecHitsL3.size()<<std::endl;
	  // Third Hierarchical Layer :: Wires inside a Layer
	  // Information unfortunately not saved ...
	  / *
	  std::vector< const TrackingRecHit * >::const_iterator recHitl3;
	  for(recHitl3 = DTRecHitsL3.begin(); recHitl3 != DTRecHitsL3.end(); ++recHitl3) {
	    DetId detidl3 = DetId((*recHitl3)->geographicalId());
	    DTChamberId    chamberIdL3(detidl3.rawId()); DTSuperLayerId sulayerIdL3(detidl3.rawId()); DTLayerId layerIdL3(detidl3.rawId()); DTWireId wireIdL3(detidl3.rawId());
	    if(physDebug) std::cout<<"               |--> DT Rec Hit in Raw ID = "<<detidl3.rawId()<<" Wire ID = "<<wireIdL3<<std::endl;
	    // if(physDebug) std::cout<<"Super Layer Number = "<<sulayerIdL2.superLayer();
	    // if(sulayerIdL2.superLayer()==0)      { if(physDebug) { std::cout<<" chamber"<<std::endl; } }
	    // else if(sulayerIdL2.superLayer()==2) { if(physDebug) { std::cout<<" theta-SL"<<std::endl; } }
	    // else                                 { if(physDebug) { std::cout<<" phi-SL"<<std::endl;   } }
	  }
	  * /
	}
      }
    }
  }
  */

}


// ------------ method called once each job just before starting event loop  ------------
void 
MyStandAloneMuonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyStandAloneMuonAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MyStandAloneMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyStandAloneMuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyStandAloneMuonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyStandAloneMuonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyStandAloneMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyStandAloneMuonAnalyzer);
