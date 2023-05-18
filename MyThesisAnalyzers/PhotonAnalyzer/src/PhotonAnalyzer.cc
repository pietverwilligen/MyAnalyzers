// -*- C++ -*-
//
// Package:    PhotonAnalyzer
// Class:      PhotonAnalyzer
// 
/**\class PhotonAnalyzer PhotonAnalyzer.cc MyAnalyzers/PhotonAnalyzer/src/PhotonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen
//         Created:  Fri Nov  12 12:03:20 CEST 2010
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <math.h>

// root include files
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
#include "TStopwatch.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include <DataFormats/EgammaReco/interface/ElectronSeed.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/Photon.h>
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/METReco/interface/MET.h>
#include "MyAnalyzers/PhotonAnalyzer/interface/SHistContainer.h"
#include "MyAnalyzers/PhotonAnalyzer/interface/SGraphContainer.h"


#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

//
// class declaration
//

class PhotonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PhotonAnalyzer(const edm::ParameterSet&);
      ~PhotonAnalyzer();


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);

  // could be optimized by using pointers to vectors instead of passing the vector as a whole

  void Selection     (const edm::Event&, const edm::EventSetup&, bool&, const pat::Photon*&, std::vector<const pat::Jet*>&, double, int, int&);
  void RA2Selection  (const edm::Event&, const edm::EventSetup&, const pat::Photon*&, std::vector<const pat::Jet*>&, double, int);
  void FillHistograms(const edm::Event&, const edm::EventSetup&, int, const pat::Photon*&, std::vector<const pat::Jet*>&, double);
  void DedicatedPlots(const edm::Event&, const edm::EventSetup&, const reco::GenParticle*&);
  void Fill          (TH1*, int, double, double, double);
  void Fill          (TH1*, int, double, double, double, double);

  void RA2Book       (TH1F*, bool, bool, double, int);
  void RA2Book1      (TH1F*, bool, bool, double);
  void RA2Book2      (TH1F*, double, double);

  // ----------member data ---------------------------

  // WRITE
  TFile * outputfile;
  TDirectoryFile *RECO, *RA2SEL, *ExtraPlots; 
  TDirectoryFile *RECO_Photon, *RECO_Event, *RECO_Jets;
  int n_PhotonDir;
  TDirectoryFile *RECO_Egamma, *RECO_RA2Sel_Egamma, *RECO_Electron_Egamma, *RECO_Decay_Egamma, *RECO_Direct_Egamma, *RECO_Fragment_Egamma, *RECO_Barrel_Egamma, *RECO_Endcap_Egamma, *RECO_2Jets_Egamma, *RECO_3Jets_Egamma, *RECO_0_Egamma, *RECO_130_Egamma, *RECO_260_Egamma;
  std::vector<SHistContainer> histcontainers, photoncontainers;
  TH1F * EventCounters, * RA2SelectionHisto, *AdditionalPrediction, *Book2_2J_150, *Book2_3J_150, *ErrorFlags;
  TH1F *JetMult_All, *PT_All, *HT_All, *MHT_All, *MHT_J1_All, *MHT_J2_All, *MHT_J3_All, *Book_All, *JetMult_PTC, *PT_PTC, *HT_PTC, *MHT_PTC, *MHT_J1_PTC, *MHT_J2_PTC, *MHT_J3_PTC, *Book_PTC;
  TH1F *JetMult_AJC, *PT_AJC, *HT_AJC, *MHT_AJC, *MHT_J1_AJC, *MHT_J2_AJC, *MHT_J3_AJC, *Book_AJC, *JetMult_AHC, *PT_AHC, *HT_AHC, *MHT_AHC, *MHT_J1_AHC, *MHT_J2_AHC, *MHT_J3_AHC, *Book_AHC;
  TH1F *JetMult_AAC, *PT_AAC, *HT_AAC, *MHT_AAC, *MHT_J1_AAC, *MHT_J2_AAC, *MHT_J3_AAC, *Book_AAC, *JetMult_AMC, *PT_AMC, *HT_AMC, *MHT_AMC, *MHT_J1_AMC, *MHT_J2_AMC, *MHT_J3_AMC, *Book_AMC;
  TH1F *JetMult_HHS, *PT_HHS, *HT_HHS, *MHT_HHS, *MHT_J1_HHS, *MHT_J2_HHS, *MHT_J3_HHS, *Book_HHS, *JetMult_HMS, *PT_HMS, *HT_HMS, *MHT_HMS, *MHT_J1_HMS, *MHT_J2_HMS, *MHT_J3_HMS, *Book_HMS;
  TH1F *PH_JCl_All, *PH_JCl_PTC, *PH_JCl_AJC, *PH_JCl_AHC, *PH_JCl_AAC, *PH_JCl_AMC, *PH_JCl_HHS, *PH_JCl_HMS;
  TH1F *PH_RJC_All, *PH_RJC_PTC, *PH_RJC_AJC, *PH_RJC_AHC, *PH_RJC_AAC, *PH_RJC_AMC, *PH_RJC_HHS, *PH_RJC_HMS;
  TH1F *PH_MHT_All, *PH_MHT_PTC, *PH_MHT_AJC, *PH_MHT_AHC, *PH_MHT_AAC, *PH_MHT_AMC, *PH_MHT_HHS, *PH_MHT_HMS;
  TH1F *Pt_All_BaS, *Pt_All_HHS, *Pt_All_HMS,*Pt_R05_BaS, *Pt_R05_HHS, *Pt_R05_HMS,*Pt_R09_BaS, *Pt_R09_HHS, *Pt_R09_HMS;
  TH1F *Book1_All, *Book1_AJC, *Book1_AHC, *Book1_AAC, *Book1_AMC, *Book1_PTC, *Book1_HHS, *Book1_HMS;
  TH1F *Book2_All, *Book2_AJC, *Book2_AHC, *Book2_AAC, *Book2_AMC, *Book2_PTC, *Book2_HHS, *Book2_HMS;

  // READ
  std::vector<int> Debug;
  bool Data, SpikeCleaning, EGMIsolation;
  std::string RootFileName, GenParticles, GenJets, RecoPhotons, RecoJets;
  int AmountOfJets;
  // INTERNAL USE
  double egm_trk_prop, egm_ecal_prop, egm_hcal_prop;
  const reco::GenParticle     *g, *g1, *b;
  const pat::Photon           *p, *p1, *p2;
  const reco::GenJet          *gj, *gj1, *gj2, *gj3;
  const pat::Jet              *r, *r1, *r2, *r3;
  TStopwatch timer, globaltimer, isolatetimer;
  std::vector<int> Count, RA2Count, AddCount;
  // DEBUG STREAM
  bool print_tss, print_dss, print_rss, print_lss, print_fss;
  std::stringstream tss, dss, rss, lss, fss;  // stringstreams: tss = Timing Debug Stream | dss = Matching Debug Stream | rss = RA2Selection Debug Stream | lss = Luminosity Debug Stream | fss = Fill Debug Stream
  int IdIsoLepton;
  int Identify;
  bool PhotonFound;
  double min_r_pi, min_allr_pi;
  // Debug
  bool PasBaseline;
  bool inDR0104;
  bool FailEGM;
};

//
// constants, enums and typedefs
//
double pi = 3.1415926535;

// Photon ID cuts:
double photon_ptcut = 100.0;
double photon_etacut = 2.50;
double egm_trk  = 2.0; 
double egm_ecal = 4.2; 
double egm_hcal = 2.2;
double photon_he = 0.05;
double sigma_barrel = 0.01;
double sigma_endcap = 0.03;
double sigma_spike = 0.001; // sigma_IEtaIEta and sigma_IPhiIPhi spike cut: p_sIEIE > 0.001 p_sIPIP > 0.001 (manual)
double sc_spike = 0.95;     // Swiss Cross and E2 over E9 spike cut:        E1/E4 < 0.95 (38X rereco: ok) E2/E9 < 0.95 (manual)
double st_spike = -3.5;     // Seed Time spike cut (ns):                    seedTime > -3.5 ns (manual)
// ECAL Geom def
double EBEE_bord = 1.479;

double EB_accept = 1.379;
double EE_accept = 1.579;
// double EB_accept = 1.4442;  // use method: isEB()  --> implemented for supercluster Eta
// double EE_accept = 1.566;   // use method: isEE()

// Conversion cuts:                                                                                                                                                                                              
double dphi_cut = 0.2;
double dcottheta_cut = 0.3;
double chi2_cut = 0.0005;

// TH1F Settings
int n_eta  = 30;  double n1_eta  = -3.0,  n2_eta  = 3.0;
int n_phi  = 36;  double n1_phi   = -3.14, n2_phi  = 3.14;
int n_et   = 200; double n1_et   = 0.0,   n2_et   = 1000;
int n_ht   = 500; double n1_ht   = 0.0,   n2_ht   = 2500;
int n_iso  = 32;  double n1_iso  = -1.0,  n2_iso  = 15;
int n_mult = 11;  double n1_mult = -0.5,  n2_mult = 10.5;
int n_pe   = 15;  double n1_pe   = 0.0,   n2_pe   = 3.0;
int n_he   = 20;  double n1_he   = 0.00,  n2_he   = 0.10;
int n_sb   = 30;  double n1_sb   = 0.000, n2_sb   = 0.030;
int n_se   = 30;  double n1_se   = 0.015, n2_se   = 0.045;
int n_s    = 60;  double n1_s    = 0.000, n2_s    = 0.060;
int n_ps   = 2;   double n1_ps   = -0.5,  n2_ps   = 1.5;
int n_cv   = 3;   double n1_cv   = -0.5,  n2_cv   = 3.5;
int n_trk  = 31;  double n1_trk  = 0.00,  n2_trk  = 30.0;

int n_pd   = 31;  double n1_pd   = 0.00,  n2_pd   = 3.14;
int n_rdd  = 50;  double n1_rdd  = 0.00,  n2_rdd  = 0.50;
int n_rd   = 50;  double n1_rd   = 0.00,  n2_rd   = 5.00;

int n_t1   = 100; double n1_t1   = 0.00,  n2_t1   = 0.10;
int n_t2   = 10;  double n1_t2   = 0.00,  n2_t2   = 0.10;

int n_sT   = 100; double n1_sT   = -25,   n2_sT   = 25;
int n_RF   = 17;  double n1_RF   = -0.5,  n2_RF   = 15.5;
int n_sS   = 6;   double n1_sS   = -0.5,  n2_sS   = 5.5;
int n_SC   = 60;  double n1_SC   = 0.0,   n2_SC   = 1.2;

int n_bk   = 13;  double n1_bk   = 0.5,   n2_bk   = 13.5;
int n_b1   = 4;   double n1_b1   = 0.5,   n2_b1   = 4.5;
int n_b2   = 16;  double n1_b2   = 0.5,   n2_b2   = 16.5;
int n_or   = 9;   double n1_or   = 0.5,   n2_or   = 9.5;

const char * PhotonContainNameArray [] = { "RECO_Egamma", "RECO_Direct_Egamma", "RECO_Electron_Egamma", "RECO_Fragment_Egamma", "RECO_Decay_Egamma", 
                                           "RECO_Barrel_Egamma", "RECO_Endcap_Egamma","RECO_2Jets_Egamma", "RECO_3Jets_Egamma", "RECO_RA2Sel_Egamma", 
					   "RECO_0_Egamma", "RECO_130_Egamma", "RECO_260_Egamma"};
std::string RecoFlagString [] = {"kGood", "kPoorReco", "kOutOfTime", "kFaultyHardware", "kNoisy", "kPoorCalib", "kSaturated", "kLeadingEdgeRecovered", "kNeighboursRecovered", "kTowerRecovered", "kFake", "kFakeNeighbours", "kDead", "kKilled", "kTPSaturated", "kL1SpikeFlag", "kUnknown"};
std::string severityString [] = {"kGood", "kProblematic", "kRecovered", "kTime", "kWeird", "kBad"};

//
// static data member definitions
//

//
// constructors and destructor
//
PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig)
{
  timer.Start();
  IdIsoLepton = 0;
  Identify = 0;
  PhotonFound = 0;
  min_r_pi = 0; min_allr_pi = 0;
  // debug
  PasBaseline = 0;
  inDR0104 = 0;
  FailEGM = 0;
  // read parameters from config file
  Debug         = iConfig.getParameter< std::vector<int> >("Debug");
  Data          = iConfig.getParameter<bool>("Data");
  RootFileName  = iConfig.getParameter<std::string>("RootFileName");
  GenParticles  = iConfig.getParameter<std::string>("GenParticles");
  GenJets       = iConfig.getParameter<std::string>("GenJets");
  RecoPhotons   = iConfig.getParameter<std::string>("RecoPhotons");
  RecoJets      = iConfig.getParameter<std::string>("RecoJets");
  SpikeCleaning = iConfig.getParameter<bool>("SpikeCleaning");
  EGMIsolation  = iConfig.getParameter<bool>("EGMIsolation");
  AmountOfJets  = iConfig.getParameter<int>("AmountOfJets");
  EGMIsolation = true;

  if(EGMIsolation) {
    egm_trk_prop = 0.000;
    egm_ecal_prop = 0.000;
    egm_hcal_prop = 0.000;
  }
  else {
    egm_trk_prop = 0.001;
    egm_ecal_prop = 0.003;
    egm_hcal_prop = 0.001;
  }
   //now do what ever initialization is needed
  outputfile = new TFile(RootFileName.c_str(), "RECREATE" );
  RECO                  = (TDirectoryFile*) outputfile->mkdir("RECO","RECO");
  RA2SEL                = (TDirectoryFile*) outputfile->mkdir("RA2SEL","RA2SEL");
  ExtraPlots            = (TDirectoryFile*) outputfile->mkdir("ExtraPlots","ExtraPlots");

  SHistContainer RECO_hist; histcontainers.push_back(RECO_hist);
  RECO_Photon           = (TDirectoryFile*) RECO->mkdir("RECO_Photon", "RECO_Photon");
  n_PhotonDir           = 13;    // Keep up to date !!!
  RECO_Egamma           = (TDirectoryFile*) RECO_Photon->mkdir("01_Egamma",           "01_Egamma");
  RECO_Direct_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("02_Direct_Egamma",    "02_Direct_Egamma");
  RECO_Electron_Egamma  = (TDirectoryFile*) RECO_Photon->mkdir("03_Electron_Egamma",  "03_Electron_Egamma"); 
  RECO_Fragment_Egamma  = (TDirectoryFile*) RECO_Photon->mkdir("04_Fragment_Egamma",  "04_Fragment_Egamma");
  RECO_Decay_Egamma     = (TDirectoryFile*) RECO_Photon->mkdir("05_Decay_Egamma",     "05_Decay_Egamma");
  RECO_Barrel_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("06_Barrel_Egamma",    "06_Barrel_Egamma");
  RECO_Endcap_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("07_Endcap_Egamma",    "07_Endcap_Egamma");    
  RECO_2Jets_Egamma     = (TDirectoryFile*) RECO_Photon->mkdir("08_2Jets_Egamma",     "08_2Jets_Egamma");
  RECO_3Jets_Egamma     = (TDirectoryFile*) RECO_Photon->mkdir("09_3Jets_Egamma",     "09_3Jets_Egamma");
  RECO_RA2Sel_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("10_RA2Sel_Egamma",    "10_RA2Sel_Egamma");
  RECO_0_Egamma         = (TDirectoryFile*) RECO_Photon->mkdir("11_Bump_0_Egamma",    "11_Bump_0_Egamma");
  RECO_130_Egamma       = (TDirectoryFile*) RECO_Photon->mkdir("12_Bump_130_Egamma",  "12_Bump_130_Egamma");
  RECO_260_Egamma       = (TDirectoryFile*) RECO_Photon->mkdir("13_Bump_260_Egamma",  "13_Bump_260_Egamma");
  SHistContainer RECO_Egamma_c, RECO_Direct_Egamma_c, RECO_Electron_Egamma_c, RECO_Fragment_Egamma_c, RECO_Decay_Egamma_c, 
    RECO_Barrel_Egamma_c, RECO_Endcap_Egamma_c, RECO_2Jets_Egamma_c, RECO_3Jets_Egamma_c, RECO_RA2Sel_Egamma_c, RECO_Bump_0_Egamma_c, RECO_Bump_130_Egamma_c, RECO_Bump_260_Egamma_c;
  SHistContainer PhotonContainerArray [] = { RECO_Egamma_c, RECO_Direct_Egamma_c, RECO_Electron_Egamma_c, RECO_Fragment_Egamma_c, RECO_Decay_Egamma_c, 
					     RECO_Barrel_Egamma_c, RECO_Endcap_Egamma_c, RECO_2Jets_Egamma_c, RECO_3Jets_Egamma_c, RECO_RA2Sel_Egamma_c, 
					     RECO_Bump_0_Egamma_c, RECO_Bump_130_Egamma_c, RECO_Bump_260_Egamma_c};
  for(int i=0; i<n_PhotonDir; ++i) { photoncontainers.push_back(PhotonContainerArray[i]); }
  RECO_Jets             = (TDirectoryFile*) RECO->mkdir("RECO_Jets",   "RECO_Jets");
  RECO_Event            = (TDirectoryFile*) RECO->mkdir("RECO_Event", "RECO_Event");

  for(unsigned int h = 0; h<photoncontainers.size(); ++h) {
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Mult",          "Photon Multiplicity Distribution", n_mult, n1_mult, n2_mult));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Pt",            "Photon Pt Distribution", n_et, n1_et, n2_et));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Et",            "Photon Et Distribution", n_et, n1_et, n2_et));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Eta",           "Photon Eta Distribution", n_eta, n1_eta, n2_eta));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Phi",           "Photon Phi Distribution", n_phi, n1_phi, n2_phi));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("SolTrkIso",     "Photon Solid Trk Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("HolTrkIso",     "Photon Hollow Trk Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("EcalIso",       "Photon Ecal Iso Deposit (DR = 0.4) Distribution", n_iso, n1_iso, n2_iso));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("HcalIso",       "Photon Hcal Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Conversions",   "Has Conversions", n_ps, n1_ps, n2_ps));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("PoverE",        "p of e^{+}e^{-} pair divided by SuperCluster E :: ConversionMethod", n_pe, n1_pe, n2_pe));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("PairPOverSCE",  "p of e^{+}e^{-} pair divided by SuperCluster E :: ManualMethod", n_pe, n1_pe, n2_pe));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("nConv",         "Number of Conversions", n_cv, n1_cv, n2_cv));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("nTrack",        "Number of Tracks of Conversion", n_cv, n1_cv, n2_cv));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("V_Conversions", "Has Valid Conversions", n_ps, n1_ps, n2_ps));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("V_PoverE",      "Valid p of e^{+}e^{-} pair divided by SuperCluster E :: ConversionMethod", n_pe, n1_pe, n2_pe));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("V_PairPOverSCE","Valid p of e^{+}e^{-} pair divided by SuperCluster E :: ManualMethod", n_pe, n1_pe, n2_pe));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("V_nConv",       "Number of Valid Conversions", n_cv, n1_cv, n2_cv));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("V_nTrack",      "Number of Tracks of Valid Conversion", n_cv, n1_cv, n2_cv));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("PixelSeeds",    "Has Electron Pixel Seeds", n_ps, n1_ps, n2_ps));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("HE",            "Hadronic over Electromagnetic Energy Fraction", n_he, n1_he, n2_he));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("SigmaIEtaIEta", "#sigma_{i #eta i #eta} ShowerShape Variable", n_s, n1_s, n2_s));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("SigmaIPhiIPhi", "#sigma_{i #phi i #phi} ShowerShape Variable", n_s, n1_s, n2_s));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("nTrackHollow",  "number of Tracks in (DR = 0.4) Cone", n_trk, n1_trk, n2_trk));
    // Analysis
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("anaBookKeeping",  "| 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk));
    // Photon Bump 
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("PFMET",              "PF MET Distribution", n_et, n1_et, n2_et));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("DR_Ph_MET",          "#Delta R Photon PF MET Distribution", n_rd, n1_rd, n2_rd));  
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("DP_Ph_MET",          "#Delta Phi Photon PF MET Distribution", n_pd, n1_pd, n2_pd));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("E",                  "Photon E Distribution", n_ht, n1_ht, n2_ht));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_eMax",        "Barrel eMax distribution", n_et, n1_et, n2_et));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_seedTime",    "Barrel Seed Time distribution", n_sT, n1_sT, n2_sT));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_RecoFlag",    "Barrel Reco Flag distribution", n_RF, n1_RF, n2_RF));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_seedSeverity","Barrel Seed Severity distribution", n_sS, n1_sS, n2_sS));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_SwissCross",  "Barrel Swiss Cross (E1/E4) distribution", n_SC, n1_SC, n2_SC));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_E2E9",        "Barrel E2/E9 distribution", n_SC, n1_SC, n2_SC));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_SigmaIEtaIEta", "Barrel #sigma_{i #eta i #eta} ShowerShape Variable", n_s, n1_s, n2_s));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Barrel_SigmaIPhiIPhi", "Barrel #sigma_{i #phi i #phi} ShowerShape Variable", n_s, n1_s, n2_s));

    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_eMax",        "Endcap eMax distribution ", n_et, n1_et, n2_et));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_seedTime",    "Endcap Seed Time distribution ", n_sT, n1_sT, n2_sT));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_RecoFlag",    "Endcap Reco Flag distribution ", n_RF, n1_RF, n2_RF));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_seedSeverity","Endcap Seed Severity distribution ", n_sS, n1_sS, n2_sS));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_SwissCross",  "Endcap Swiss Cross (E1/E4) distribution ", n_SC, n1_SC, n2_SC));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_E2E9",        "Endcap E2/E9 distribution ", n_SC, n1_SC, n2_SC));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_SigmaIEtaIEta", "Endcap #sigma_{i #eta i #eta} ShowerShape Variable", n_s, n1_s, n2_s));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("Endcap_SigmaIPhiIPhi", "Endcap #sigma_{i #phi i #phi} ShowerShape Variable", n_s, n1_s, n2_s));
    // nice binlables for RecoFlag and SeedSeverity
    for(int i=0; i<n_RF; ++i) {
      photoncontainers[h].get(PhotonContainNameArray[h],"Barrel_RecoFlag")->GetXaxis()->SetBinLabel(i+1, RecoFlagString[i].c_str()); 
      photoncontainers[h].get(PhotonContainNameArray[h],"Endcap_RecoFlag")->GetXaxis()->SetBinLabel(i+1, RecoFlagString[i].c_str());
    }
    for(int i=0; i<n_sS; ++i) {
      photoncontainers[h].get(PhotonContainNameArray[h],"Barrel_seedSeverity")->GetXaxis()->SetBinLabel(i+1, severityString[i].c_str()); 
      photoncontainers[h].get(PhotonContainNameArray[h],"Endcap_seedSeverity")->GetXaxis()->SetBinLabel(i+1, severityString[i].c_str());
    }
  }

  for(unsigned int h = 0; h<histcontainers.size(); ++h) {
    histcontainers[h].book("Jets",   new TH1F("JetsMult", "Jets Multiplicity Distribution (with pt > 50)", n_mult, n1_mult, n2_mult));
    histcontainers[h].book("Jets",   new TH1F("JetsPtUn", "Unweighted Jets Pt Distribution", n_et, n1_et, n2_et));
    histcontainers[h].book("Jets",   new TH1F("JetsPt",   "Jets Pt Distribution",  n_et, n1_et, n2_et));
    histcontainers[h].book("Jets",   new TH1F("JetsEta",  "Jets Eta Distribution", n_eta, n1_eta, n2_eta));
    histcontainers[h].book("Jets",   new TH1F("JetsPhi",  "Jets Phi Distribution", n_phi, n1_phi, n2_phi));
    histcontainers[h].book("Jets",   new TH1F("HT",       "HT Distribution",  n_ht, n1_ht, n2_ht));
    histcontainers[h].book("Jets",   new TH1F("MHT",      "MHT Distribution",  n_et, n1_et, n2_et));
    histcontainers[h].book("Jets",   new TH1F("MHT_J1",   "#Delta Phi (MHT - 1st Jet) Distribution", n_phi, n1_phi, n2_phi));
    histcontainers[h].book("Jets",   new TH1F("MHT_J2",   "#Delta Phi (MHT - 2nd Jet) Distribution", n_phi, n1_phi, n2_phi));
    histcontainers[h].book("Jets",   new TH1F("MHT_J3",   "#Delta Phi (MHT - 3rd Jet) Distribution", n_phi, n1_phi, n2_phi));

    // event information
    histcontainers[h].book("Event",  new TH1F("nConvEvt",  "Number of Conversions in Event", n_mult, n1_mult, n2_mult));
    histcontainers[h].book("Event",  new TH1F("nVertEvt",  "Number of Vertices in Event",    n_mult, n1_mult, n2_mult));
    histcontainers[h].book("Event",  new TH1F("PATMatch",  "Photon is PATMatched?",          n_ps, n1_ps, n2_ps));
    // histcontainers[h].book("Event",  new TH1F("Distance",  "Distance between Hard Photon and Photon in #eta, #varphi:", n_dist, n1_dist, n2_dist));
    // histcontainers[h].book("Event",  new TH1F("DR_Before", "Distance between #gamma and Jet before removing #gamma out the jetcollection by #eta, #varphi: #Delta R (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));
    // histcontainers[h].book("Event",  new TH1F("DP_Before", "Distance between #gamma and Jet before removing #gamma out the jetcollection by p_{T}: #Delta p_{T} / p_{T} (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));
    // histcontainers[h].book("Event",  new TH1F("DR_After",  "Distance between #gamma and Jet in #eta, #varphi: #Delta R (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));
    // histcontainers[h].book("Event",  new TH1F("DP_After",  "Distance between #gamma and Jet in p_{T}: #Delta p_{T} / p_{T} (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));

    histcontainers[0].book("Event", new TH2F("PtJet_PtPhoton_All",      "Pt Closesth Jet vs Pt Photon for All Events", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_PtPhoton_NF",       "Pt Closesth Jet vs Pt Photon for Non Found Jets", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_PtPhoton_F",        "Pt Closesth Jet vs Pt Photon for Found Jets", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_PtPhoton_DeltaR01", "Pt Closesth Jet vs Pt Photon for Jets in 0.1 cone of Photons", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_PtPhoton_DeltaR04", "Pt Closesth Jet vs Pt Photon for Jets in 0.1-0.4 cone of Photons", 100,0,500,100,0,500));

    histcontainers[0].book("Event", new TH2F("PtJet_GenPtPhoton_All",      "Pt Closesth Jet vs Gen Pt Photon for All Events", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_GenPtPhoton_NF",       "Pt Closesth Jet vs Gen Pt Photon for Non Found Jets", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_GenPtPhoton_F",        "Pt Closesth Jet vs Gen Pt Photon for Found Jets", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_GenPtPhoton_DeltaR01", "Pt Closesth Jet vs Gen Pt Photon for Jets in 0.1 cone of Photons", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("PtJet_GenPtPhoton_DeltaR04", "Pt Closesth Jet vs Gen Pt Photon for Jets in 0.1-0.4 cone of Photons", 100,0,500,100,0,500));

    histcontainers[0].book("Event", new TH2F("TotJetPt_PtPhoton_DeltaR01", "Total Jet Pt vs Pt Photon for Jets in 0.1 cone of Photons", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH2F("TotJetPt_PtPhoton_DeltaR04", "Total Jet Pt vs Pt Photon for Jets in 0.1-0.4 cone of Photons", 100,0,500,100,0,500));
    histcontainers[0].book("Event", new TH1F("N_Jets_in_DeltaR01",          "Number of Jets in 0.1 cone of Photons", n_mult, n1_mult, n2_mult));
    histcontainers[0].book("Event", new TH1F("N_Jets_in_DeltaR04",          "Number of Jets in 0.1-0.4 cone of Photons",n_mult, n1_mult, n2_mult));

    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_All",      "DR Closesth Jet - Photon for All Events", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_NF",       "DR Closesth Jet - Photon for Non Found Jets", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_F",        "DR Closesth Jet - Photon for Found Jets", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_DeltaR01", "DR Closesth Jet - Photon for Jets in 0.1 cone of Photons", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_DeltaR04", "DR Closesth Jet - Photon for Jets in 0.1-0.4 cone of Photons", n_rdd, n1_rdd, n2_rdd));

    histcontainers[0].book("Event", new TH1F("DR_Parton_Photon_All",      "DR Closesth Parton - Photon for All Events", n_pd, n1_pd, n2_pd));
    /*
    histcontainers[0].book("Event", new TH1F("DR_Parton_Photon_NF",       "DR Closesth Parton - Photon for Non Found Jets", n_dist, n1_dist, n2_dist));
    histcontainers[0].book("Event", new TH1F("DR_Parton_Photon_F",        "DR Closesth Parton - Photon for Found Jets", n_dist, n1_dist, n2_dist));
    histcontainers[0].book("Event", new TH1F("DR_Parton_Photon_DeltaR01", "DR Closesth Parton - Photon for Jets in 0.1 cone of Photons", n_dist, n1_dist, n2_dist));
    histcontainers[0].book("Event", new TH1F("DR_Parton_Photon_DeltaR04", "DR Closesth Parton - Photon for Jets in 0.1-0.4 cone of Photons", n_dist, n1_dist, n2_dist));
    */
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_Direct",    "DR Mother - Direct Photon", n_pd, n1_pd, n2_pd));
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_Fragment",  "DR Mother - Fragment Photon", n_pd, n1_pd, n2_pd));
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_ISR",       "DR Mother - ISR Photon", n_pd, n1_pd, n2_pd));
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_FSR",       "DR Mother - FSR Photon", n_pd, n1_pd, n2_pd));
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_Decay",     "DR Mother - Decayed Meson", n_pd, n1_pd, n2_pd));
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_Ewk",       "DR Mother - Ewk Rad Photon", n_pd, n1_pd, n2_pd));
    histcontainers[0].book("Event", new TH1F("DR_Mother_Photon_Mistag",    "DR Mother - Mistagged Electron", n_pd, n1_pd, n2_pd));

    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_Direct",    "DR Closesth Jet - Direct Photon", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_Fragment",  "DR Closesth Jet - Fragment Photon", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_ISR",       "DR Closesth Jet - ISR Photon", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_FSR",       "DR Closesth Jet - FSR Photon", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_Decay",     "DR Closesth Jet - Decayed Meson", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_Ewk",       "DR Closesth Jet - Ewk Rad Photon", n_rdd, n1_rdd, n2_rdd));
    histcontainers[0].book("Event", new TH1F("DR_Jet_Photon_Mistag",    "DR Closesth Jet - Mistagged Electron", n_rdd, n1_rdd, n2_rdd));



    histcontainers[h].book("Timing", new TH1F("tSel_Real", "Real Time SelectCollections", n_t1, n1_t1, n2_t1));
    histcontainers[h].book("Timing", new TH1F("tIso_Real", "Real Time Isolate",           n_t1, n1_t1, n2_t1));
    histcontainers[h].book("Timing", new TH1F("tAna_Real", "Real Time Analyze",           n_t1, n1_t1, n2_t1));
    histcontainers[h].book("Timing", new TH1F("tSel_CPU",  "CPU Time SelectCollections",  n_t2, n1_t2, n2_t2));
    histcontainers[h].book("Timing", new TH1F("tIso_CPU",  "CPU Time Isolate",            n_t2, n1_t2, n2_t2));
    histcontainers[h].book("Timing", new TH1F("tAna_CPU",  "CPU Time Analyze",            n_t2, n1_t2, n2_t2));

    // RA2 Selection Plots
    JetMult_All= new TH1F("0_JetMult_All", "PT50ETA25 Jet Multiplicity", n_mult, n1_mult, n2_mult);
    PT_All     = new TH1F("0_PT_All" ,     "Boson PT Distribution",  n_et, n1_et, n2_et);
    HT_All     = new TH1F("0_HT_All" ,     "HT Distribution",  n_ht, n1_ht, n2_ht);
    MHT_All    = new TH1F("0_MHT_All",     "MHT Distribution", n_et, n1_et, n2_et);
    MHT_J1_All = new TH1F("0_MHT_J1_All" , "#Delta Phi (MHT - 1st Jet) Distribution",  n_pd, n1_pd, n2_pd);
    MHT_J2_All = new TH1F("0_MHT_J2_All" , "#Delta Phi (MHT - 2nd Jet) Distribution",  n_pd, n1_pd, n2_pd);
    MHT_J3_All = new TH1F("0_MHT_J3_All" , "#Delta Phi (MHT - 3rd Jet) Distribution",  n_pd, n1_pd, n2_pd);
    Book_All   = new TH1F("0_Book_All"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_All  = new TH1F("0_Book1_All"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_All  = new TH1F("0_Book2_All"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_All = new TH1F("0_PH_JCl_All",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_All = new TH1F("0_PH_RJC_All",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_All = new TH1F("0_PH_MHT_All" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    JetMult_AJC= new TH1F("1_JetMult_AJC", "PT50ETA25 Jet Multiplicity After Jet Cuts", n_mult, n1_mult, n2_mult);
    PT_AJC     = new TH1F("1_PT_AJC" ,     "Boson PT Distribution After Jet Cuts",  n_et, n1_et, n2_et);
    HT_AJC     = new TH1F("1_HT_AJC" ,     "HT Distribution After Jet Cuts", n_ht, n1_ht, n2_ht);
    MHT_AJC    = new TH1F("1_MHT_AJC",     "MHT Distribution After Jet Cuts", n_et, n1_et, n2_et);
    MHT_J1_AJC = new TH1F("1_MHT_J1_AJC" , "#Delta Phi (MHT - 1st Jet) Distribution After Jet Cuts", n_pd, n1_pd, n2_pd);
    MHT_J2_AJC = new TH1F("1_MHT_J2_AJC" , "#Delta Phi (MHT - 2nd Jet) Distribution After Jet Cuts", n_pd, n1_pd, n2_pd);
    MHT_J3_AJC = new TH1F("1_MHT_J3_AJC" , "#Delta Phi (MHT - 3rd Jet) Distribution After Jet Cuts", n_pd, n1_pd, n2_pd);
    Book_AJC   = new TH1F("1_Book_AJC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_AJC  = new TH1F("1_Book1_AJC"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_AJC  = new TH1F("1_Book2_AJC"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_AJC = new TH1F("1_PH_JCl_AJC",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_AJC = new TH1F("1_PH_RJC_AJC",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_AJC = new TH1F("1_PH_MHT_AJC" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    JetMult_AHC= new TH1F("2_JetMult_AHC", "PT50ETA25 Jet Multiplicity After HT Cut", n_mult, n1_mult, n2_mult);
    PT_AHC     = new TH1F("2_PT_AHC" ,     "Boson PT Distribution After HT Cut",  n_et, n1_et, n2_et);
    HT_AHC     = new TH1F("2_HT_AHC" ,     "HT Distribution After HT Cut", n_ht, n1_ht, n2_ht);
    MHT_AHC    = new TH1F("2_MHT_AHC",     "MHT Distribution After HT Cut", n_et, n1_et, n2_et);
    MHT_J1_AHC = new TH1F("2_MHT_J1_AHC" , "#Delta Phi (MHT - 1st Jet) Distribution After HT Cut", n_pd, n1_pd, n2_pd);
    MHT_J2_AHC = new TH1F("2_MHT_J2_AHC" , "#Delta Phi (MHT - 2nd Jet) Distribution After HT Cut", n_pd, n1_pd, n2_pd);
    MHT_J3_AHC = new TH1F("2_MHT_J3_AHC" , "#Delta Phi (MHT - 3rd Jet) Distribution After HT Cut", n_pd, n1_pd, n2_pd);
    Book_AHC   = new TH1F("2_Book_AHC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_AHC  = new TH1F("2_Book1_AHC"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_AHC  = new TH1F("2_Book2_AHC"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_AHC = new TH1F("2_PH_JCl_AHC",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_AHC = new TH1F("2_PH_RJC_AHC",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_AHC = new TH1F("2_PH_MHT_AHC" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    JetMult_AAC= new TH1F("3_JetMult_AAC", "PT50ETA25 Jet Multiplicity After Angular Cuts", n_mult, n1_mult, n2_mult);
    PT_AAC     = new TH1F("3_PT_AAC" ,     "Boson PT Distribution After Angular Cuts",  n_et, n1_et, n2_et);
    HT_AAC     = new TH1F("3_HT_AAC" ,     "HT Distribution After Angular Cuts", n_ht, n1_ht, n2_ht);
    MHT_AAC    = new TH1F("3_MHT_AAC",     "MHT Distribution After Angular Cuts", n_et, n1_et, n2_et);
    MHT_J1_AAC = new TH1F("3_MHT_J1_AAC" , "#Delta Phi (MHT - 1st Jet) Distribution After Angular Cuts", n_pd, n1_pd, n2_pd);
    MHT_J2_AAC = new TH1F("3_MHT_J2_AAC" , "#Delta Phi (MHT - 2nd Jet) Distribution After Angular Cuts", n_pd, n1_pd, n2_pd);
    MHT_J3_AAC = new TH1F("3_MHT_J3_AAC" , "#Delta Phi (MHT - 3rd Jet) Distribution After Angular Cuts", n_pd, n1_pd, n2_pd);
    Book_AAC   = new TH1F("3_Book_AAC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_AAC  = new TH1F("3_Book1_AAC"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_AAC  = new TH1F("3_Book2_AAC"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_AAC = new TH1F("3_PH_JCl_AAC",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_AAC = new TH1F("3_PH_RJC_AAC",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_AAC = new TH1F("3_PH_MHT_AAC" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    JetMult_AMC= new TH1F("4_JetMult_AMC", "PT50ETA25 Jet Multiplicity After MHT Cut: Baseline Selection", n_mult, n1_mult, n2_mult);
    PT_AMC     = new TH1F("4_PT_AMC" ,     "Boson PT Distribution after MHT Cut: Baseline Selection",  n_et, n1_et, n2_et);
    HT_AMC     = new TH1F("4_HT_AMC" ,     "HT Distribution After MHT Cut: Baseline Selection", n_ht, n1_ht, n2_ht);
    MHT_AMC    = new TH1F("4_MHT_AMC",     "MHT Distribution After MHT Cut: Baseline Selection", n_et, n1_et, n2_et);
    MHT_J1_AMC = new TH1F("4_MHT_J1_AMC" , "#Delta Phi (MHT - 1st Jet) Distribution After MHT Cut: Baseline Selection", n_pd, n1_pd, n2_pd);
    MHT_J2_AMC = new TH1F("4_MHT_J2_AMC" , "#Delta Phi (MHT - 2nd Jet) Distribution After MHT Cut: Baseline Selection", n_pd, n1_pd, n2_pd);
    MHT_J3_AMC = new TH1F("4_MHT_J3_AMC" , "#Delta Phi (MHT - 3rd Jet) Distribution After MHT Cut: Baseline Selection", n_pd, n1_pd, n2_pd);
    Book_AMC   = new TH1F("4_Book_AMC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_AMC  = new TH1F("4_Book1_AMC"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_AMC  = new TH1F("4_Book2_AMC"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_AMC = new TH1F("4_PH_JCl_AMC",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_AMC = new TH1F("4_PH_RJC_AMC",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_AMC = new TH1F("4_PH_MHT_AMC" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    JetMult_PTC= new TH1F("5_JetMult_PTC", "PT50ETA25 Jet Multiplicity Baseline Selection After Boson Pt Cut", n_mult, n1_mult, n2_mult);
    PT_PTC     = new TH1F("5_PT_PTC" ,     "Boson PT Distribution Baseline Selection After Boson Pt Cut",  n_et, n1_et, n2_et);
    HT_PTC     = new TH1F("5_HT_PTC" ,     "HT Distribution Baseline Selection After Boson Pt Cut",  n_ht, n1_ht, n2_ht);
    MHT_PTC    = new TH1F("5_MHT_PTC",     "MHT Distribution Baseline Selection After Boson Pt Cut", n_et, n1_et, n2_et);
    MHT_J1_PTC = new TH1F("5_MHT_J1_PTC" , "#Delta Phi (MHT - 1st Jet) Distribution Baseline Selection After Boson Pt Cut",  n_pd, n1_pd, n2_pd);
    MHT_J2_PTC = new TH1F("5_MHT_J2_PTC" , "#Delta Phi (MHT - 2nd Jet) Distribution Baseline Selection After Boson Pt Cut",  n_pd, n1_pd, n2_pd);
    MHT_J3_PTC = new TH1F("5_MHT_J3_PTC" , "#Delta Phi (MHT - 3rd Jet) Distribution Baseline Selection After Boson Pt Cut",  n_pd, n1_pd, n2_pd);
    Book_PTC   = new TH1F("5_Book_PTC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_PTC  = new TH1F("5_Book1_PTC"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_PTC  = new TH1F("5_Book2_PTC"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_PTC = new TH1F("5_PH_JCl_PTC",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_PTC = new TH1F("5_PH_RJC_PTC",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_PTC = new TH1F("5_PH_MHT_PTC" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);    

    JetMult_HHS= new TH1F("6_JetMult_HHS", "PT50ETA25 Jet Multiplicity: High HT Search Selection", n_mult, n1_mult, n2_mult);
    PT_HHS     = new TH1F("6_PT_HHS" ,     "Boson PT Distribution: High HT Selection",  n_et, n1_et, n2_et);
    HT_HHS     = new TH1F("6_HT_HHS" ,     "HT Distribution: High HT Selection", n_ht, n1_ht, n2_ht);
    MHT_HHS    = new TH1F("6_MHT_HHS",     "MHT Distribution: High HT Selection", n_et, n1_et, n2_et);
    MHT_J1_HHS = new TH1F("6_MHT_J1_HHS" , "#Delta Phi (MHT - 1st Jet) Distribution: High HT Selection", n_pd, n1_pd, n2_pd);
    MHT_J2_HHS = new TH1F("6_MHT_J2_HHS" , "#Delta Phi (MHT - 2nd Jet) Distribution: High HT Selection", n_pd, n1_pd, n2_pd);
    MHT_J3_HHS = new TH1F("6_MHT_J3_HHS" , "#Delta Phi (MHT - 3rd Jet) Distribution: High HT Selection", n_pd, n1_pd, n2_pd);
    Book_HHS   = new TH1F("6_Book_HHS"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_HHS  = new TH1F("6_Book1_HHS"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_HHS  = new TH1F("6_Book2_HHS"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_HHS = new TH1F("6_PH_JCl_HHS",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_HHS = new TH1F("6_PH_RJC_HHS",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_HHS = new TH1F("6_PH_MHT_HHS" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    JetMult_HMS= new TH1F("7_JetMult_HMS", "PT50ETA25 Jet Multiplicity: High MHT Selection", n_mult, n1_mult, n2_mult);
    PT_HMS     = new TH1F("7_PT_HMS" ,     "Boson PT Distribution: High MHT Selection",  n_et, n1_et, n2_et);
    HT_HMS     = new TH1F("7_HT_HMS" ,     "HT Distribution: High MHT Selection", n_ht, n1_ht, n2_ht);
    MHT_HMS    = new TH1F("7_MHT_HMS",     "MHT Distribution: High MHT Selection", n_et, n1_et, n2_et);
    MHT_J1_HMS = new TH1F("7_MHT_J1_HMS" , "#Delta Phi (MHT - 1st Jet) Distribution: High MHT Selection", n_pd, n1_pd, n2_pd);
    MHT_J2_HMS = new TH1F("7_MHT_J2_HMS" , "#Delta Phi (MHT - 2nd Jet) Distribution: High MHT Selection", n_pd, n1_pd, n2_pd);
    MHT_J3_HMS = new TH1F("7_MHT_J3_HMS" , "#Delta Phi (MHT - 3rd Jet) Distribution: High MHT Selection", n_pd, n1_pd, n2_pd);
    Book_HMS   = new TH1F("7_Book_HMS"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec | ewk | ISR/FSR | ISR | FSR |", n_bk, n1_bk, n2_bk);
    Book1_HMS  = new TH1F("7_Book1_HMS"   ,"Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE", n_b1, n1_b1, n2_b1);
    Book2_HMS  = new TH1F("7_Book2_HMS"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    PH_JCl_HMS = new TH1F("7_PH_JCl_HMS",  "#Delta Phi (Photon - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
    PH_RJC_HMS = new TH1F("7_PH_RJC_HMS",  "#Delta R (Photon - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
    PH_MHT_HMS = new TH1F("7_PH_MHT_HMS" , "#Delta Phi (Photon - MHT) Distribution",  n_pd, n1_pd, n2_pd);

    Pt_All_BaS = new TH1F("Pt_All_BaS", "All Photon Pt Distribution for Baseline Selection", n_et, n1_et, n2_et);
    Pt_All_HHS = new TH1F("Pt_All_HHS", "All Photon Pt Distribution for High HT  Selection", n_et, n1_et, n2_et);
    Pt_All_HMS = new TH1F("Pt_All_HMS", "All Photon Pt Distribution for High MHT Selection", n_et, n1_et, n2_et);

    Pt_R05_BaS = new TH1F("Pt_R05_BaS", "R05 Photon Pt Distribution for Baseline Selection", n_et, n1_et, n2_et);
    Pt_R05_HHS = new TH1F("Pt_R05_HHS", "R05 Photon Pt Distribution for High HT  Selection", n_et, n1_et, n2_et);
    Pt_R05_HMS = new TH1F("Pt_R05_HMS", "R05 Photon Pt Distribution for High MHT Selection", n_et, n1_et, n2_et);

    Pt_R09_BaS = new TH1F("Pt_R09_BaS", "R09 Photon Pt Distribution for Baseline Selection", n_et, n1_et, n2_et);
    Pt_R09_HHS = new TH1F("Pt_R09_HHS", "R09 Photon Pt Distribution for High HT  Selection", n_et, n1_et, n2_et);
    Pt_R09_HMS = new TH1F("Pt_R09_HMS", "R09 Photon Pt Distribution for High MHT Selection", n_et, n1_et, n2_et);


    // 7 counters + 4 more offline checks + 2 more for Barrel and Endcap Final Checks
    for(unsigned int i=0; i<13; ++i) { Count.push_back(0); }
    EventCounters = new TH1F("EventCounters", "Events: All Processed | post HLT | preFilter | postStdClean | postPFClean | postPhoton | postPFJets| ISO #gamma |2J != #gamma | 3J != #gamma | RA2 Selected", 11, 0.5, 11.5);
    // 5 RA2 Selection Criteria, so 6 counters (+ 1 - 1) + additional counters for influence Boson PT Cut, High HT Search Region and High MHT Search Region
    for(unsigned int i=0; i<9; ++i) { RA2Count.push_back(0); }
    RA2SelectionHisto = new TH1F("RA2SelectionHisto", "Events: Processed | Passing PAT | with ISO #gamma and 3 Jets | with HT > 300 | with Angular Cuts | with MHT > 150 | Boson PT Cut | High HT | High MHT", 9, 0.5, 9.5);
    // Additional Prediction
    for(unsigned int i=0; i<4; ++i) { AddCount.push_back(0); }
    AdditionalPrediction = new TH1F("AdditionalPrediction", "ISO #gamma + 2 Jets | 2 Jets + MHT > 150 | ISO #gamma + 3 Jets | 3 Jets + MHT > 150", 4, 0.5, 4.5);
    Book2_2J_150 = new TH1F("Book2_2J_150"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);
    Book2_3J_150 = new TH1F("Book2_3J_150"   ,"Bookkeeping: #eta #epsilon [0.0 | 0.9 | 1.4442 || 1.566 | 2.1 | 2.5] x E_{T} #epsilon [ 100 | 120 | 200 | 300 | +300]", n_b2, n1_b2, n2_b2);

    // Error Flags Histogram
    ErrorFlags = new TH1F("ErrorFlags", "No Matched GEN #gamma | no GEN #gamma in event | No Jet == #gamma | No Jets in event", 4, 0.5, 4.5);
    // Debug Stream
    print_tss = 0, print_dss = 0, print_rss = 0, print_lss = 0, print_fss = 0;
    timer.Stop();
    tss<<" PhotonAnalyzer :: Constructor :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
  }


}


PhotonAnalyzer::~PhotonAnalyzer()
{
  std::cout<<"Amount of Rejected Events: "<<IdIsoLepton<<std::endl;

  timer.Start();
  outputfile->cd();
  EventCounters->Write();
  RA2SelectionHisto->Write();
  AdditionalPrediction->Write();
  Book2_2J_150->Write();
  Book2_3J_150->Write();
  ErrorFlags->Write();
  histcontainers[0].write(RECO);
  for(int h=0; h<n_PhotonDir; ++h) {photoncontainers[h].write(RECO);}

  outputfile->cd();
  RA2SEL->cd();
  JetMult_All->Write();
  PT_All->Write();
  HT_All->Write();
  MHT_All->Write();
  MHT_J1_All->Write();
  MHT_J2_All->Write();
  MHT_J3_All->Write();
  Book_All->Write();
  Book1_All->Write();
  Book2_All->Write();
  PH_JCl_All->Write();
  PH_RJC_All->Write();
  PH_MHT_All->Write();

  JetMult_AJC->Write();
  PT_AJC->Write();
  HT_AJC->Write();
  MHT_AJC->Write();
  MHT_J1_AJC->Write();
  MHT_J2_AJC->Write();
  MHT_J3_AJC->Write();
  Book_AJC->Write();
  Book1_AJC->Write();
  Book2_AJC->Write();
  PH_JCl_AJC->Write();
  PH_RJC_AJC->Write();
  PH_MHT_AJC->Write();

  JetMult_AHC->Write();
  PT_AHC->Write();
  HT_AHC->Write();
  MHT_AHC->Write();
  MHT_J1_AHC->Write();
  MHT_J2_AHC->Write();
  MHT_J3_AHC->Write();
  Book_AHC->Write();
  Book1_AHC->Write();
  Book2_AHC->Write();
  PH_JCl_AHC->Write();
  PH_RJC_AHC->Write();
  PH_MHT_AHC->Write();

  JetMult_AAC->Write();
  PT_AAC->Write();
  HT_AAC->Write();
  MHT_AAC->Write();
  MHT_J1_AAC->Write();
  MHT_J2_AAC->Write();
  MHT_J3_AAC->Write();
  Book_AAC->Write();
  Book1_AAC->Write();
  Book2_AAC->Write();
  PH_JCl_AAC->Write();
  PH_RJC_AAC->Write();
  PH_MHT_AAC->Write();

  JetMult_AMC->Write();
  PT_AMC->Write();
  HT_AMC->Write();
  MHT_AMC->Write();
  MHT_J1_AMC->Write();
  MHT_J2_AMC->Write();
  MHT_J3_AMC->Write();
  Book_AMC->Write();
  Book1_AMC->Write();
  Book2_AMC->Write();
  PH_JCl_AMC->Write();
  PH_RJC_AMC->Write();
  PH_MHT_AMC->Write();

  JetMult_PTC->Write();
  PT_PTC->Write();
  HT_PTC->Write();
  MHT_PTC->Write();
  MHT_J1_PTC->Write();
  MHT_J2_PTC->Write();
  MHT_J3_PTC->Write();
  Book_PTC->Write();
  Book1_PTC->Write();
  Book2_PTC->Write();
  PH_JCl_PTC->Write();
  PH_RJC_PTC->Write();
  PH_MHT_PTC->Write();

  JetMult_HHS->Write();
  PT_HHS->Write();
  HT_HHS->Write();
  MHT_HHS->Write();
  MHT_J1_HHS->Write();
  MHT_J2_HHS->Write();
  MHT_J3_HHS->Write();
  Book_HHS->Write();
  Book1_HHS->Write();
  Book2_HHS->Write();
  PH_JCl_HHS->Write();
  PH_RJC_HHS->Write();
  PH_MHT_HHS->Write();

  JetMult_HMS->Write();
  PT_HMS->Write();
  HT_HMS->Write();
  MHT_HMS->Write();
  MHT_J1_HMS->Write();
  MHT_J2_HMS->Write();
  MHT_J3_HMS->Write();
  Book_HMS->Write();
  Book1_HMS->Write();
  Book2_HMS->Write();
  PH_JCl_HMS->Write();
  PH_RJC_HMS->Write();
  PH_MHT_HMS->Write();

  ExtraPlots->cd();
  Pt_All_BaS->Write();
  Pt_All_HHS->Write();
  Pt_All_HMS->Write();
  Pt_R05_BaS->Write();
  Pt_R05_HHS->Write();
  Pt_R05_HMS->Write();
  Pt_R09_BaS->Write();
  Pt_R09_HHS->Write();
  Pt_R09_HMS->Write();


  timer.Stop();
  tss<<" PhotonAnalyzer :: Destructor :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
  tss<<"closing "<<RootFileName.c_str()<<std::endl;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // EVENT INFORMATION                                                                                                                                                             
  int evNum = (iEvent.id()).event();
  int rnNum = (iEvent.id()).run();

  // if(rnNum == 146430 && evNum == 7882279) {
  dss.str("");
  dss<<"Analyzing Run nr "<<rnNum<<" Event nr "<<evNum<<std::endl; // Debug[1] = 1; Debug[2] = 1; Debug[3] = 1; 

  // REJECT EVENTS IN CASE ID & ISO LEPTONS
  edm::Handle< std::vector<pat::Electron> >  patElectrons;
  edm::Handle< std::vector<pat::Muon> >  patMuons;
  iEvent.getByLabel("patElectronsPFIDIso", patElectrons);
  iEvent.getByLabel("patMuonsPFIDIso", patMuons);
  // std::cout<<"Size of patElectrons: "<<patElectrons->size()<<std::endl;
  // std::cout<<"Size of patMuons: "<<patMuons->size()<<std::endl;
  if (patElectrons->size() != 0 || patMuons->size() != 0) { 
    ++IdIsoLepton; 
    std::cout<<"Event Rejected !!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Electrons: "<<patElectrons->size()<<" PFIDIso Muons: "<<patMuons->size()<<std::endl;
  }
  if (patElectrons->size() == 0 && patMuons->size() == 0) {

  
  globaltimer.Start();
  print_dss = 0, print_rss = 0, print_lss = 0, print_fss = 0;
  p = 0;
  double evWeight = 1;
  PhotonFound = 0;

  // GEN INFORMATION
  if(Data==0) {
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByLabel("generator",genEventInfo);
    evWeight = genEventInfo->weight();
  }
  // EVENT INFORMATION
  edm::Handle< std::vector<reco::Vertex> > Vertices;
  iEvent.getByLabel("offlinePrimaryVertices", Vertices);
  int nVertices = Vertices->size();

  const pat::Photon * Photon = 0;
  std::vector<const pat::Jet*> Jets;
  int Identity = 0;
  Selection(iEvent, iSetup, Data, Photon, Jets, evWeight, nVertices, Identity);
  if(PhotonFound) {
    RA2Selection(iEvent, iSetup, Photon, Jets, evWeight, Identity); 
  }
  globaltimer.Stop();
  tss<<" PhotonAnalyzer :: Analyze :: Real Time; "<<globaltimer.RealTime()<<" CPU Time: "<<globaltimer.CpuTime()<<std::endl;
  Fill(histcontainers[0].get("Timing","tAna_Real"), n_t1, n1_t1, n2_t1, globaltimer.RealTime());
  Fill(histcontainers[0].get("Timing","tAna_CPU"),  n_t2, n1_t2, n2_t2, globaltimer.CpuTime());

  //debug
  if( (PasBaseline && inDR0104) || (PasBaseline && FailEGM) ) {
    std::cout<<" "<<std::endl<<"Matching Debug Stream:    "<<std::endl<<"--------------------------"<<std::endl<<dss.str()<<std::endl;
    std::cout<<" "<<std::endl<<"RA2Selection Debug Stream:"<<std::endl<<"--------------------------"<<std::endl<<rss.str()<<std::endl;
  }

  
  if(Debug[1] || print_dss) std::cout<<" "<<std::endl<<"Matching Debug Stream:    "<<std::endl<<"--------------------------"<<std::endl<<dss.str()<<std::endl; else dss.str("");
  if(Debug[2] || print_rss) std::cout<<" "<<std::endl<<"RA2Selection Debug Stream:"<<std::endl<<"--------------------------"<<std::endl<<rss.str()<<std::endl; else rss.str("");
  if(Debug[3] || print_lss) std::cout<<" "<<std::endl<<"Luminosity Debug Stream:  "<<std::endl<<"--------------------------"<<std::endl<<lss.str()<<std::endl; else lss.str("");
  if(            print_fss) std::cout<<" "<<std::endl<<"Histo Fill Debug Stream:  "<<std::endl<<"--------------------------"<<std::endl<<fss.str()<<std::endl; else fss.str("");
  }
}

void
PhotonAnalyzer::Selection(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool& Data, const pat::Photon * & photon, std::vector<const pat::Jet*>& Jets, double evW, int nVertices, int& Identity)
{
  timer.Start();
  Identity = 0;
  PhotonFound = 0;
  PasBaseline = 0;
  inDR0104 = 0;
  FailEGM = 0;
  //
  // Handles for Data
  //
  // edm::Handle< std::vector<reco::PFJet> >       recoJets;
  edm::Handle< std::vector<pat::Jet> >          recoJets;
  edm::Handle< std::vector<reco::PFCandidate> > pfPhotons;
  edm::Handle< std::vector<pat::Photon> >       patPhotons;
  edm::Handle< std::vector<reco::Photon> >      recoPhotons;
  bool usePfPhotons = 0, usePatPhotons = 0, useRecoPhotons = 0; 

  if      (RecoPhotons.compare("particleFlow")==0 || RecoPhotons.compare("pfAllPhotons")==0       || RecoPhotons.compare("pfAllPhotonsPF")==0)  
    { iEvent.getByLabel(RecoPhotons, pfPhotons);   usePfPhotons = 1;   }
  else if (RecoPhotons.compare("patPhotons")==0   || RecoPhotons.compare("selectedPatPhotons")==0 || RecoPhotons.compare("cleanPatPhotons")==0) 
    { iEvent.getByLabel(RecoPhotons, patPhotons);  usePatPhotons = 1;  }
  else
    { iEvent.getByLabel(RecoPhotons, recoPhotons); useRecoPhotons = 1; }
  iEvent.getByLabel(RecoJets, recoJets);

  // We use PAT Photons!
  int count_good_photons = 0;
  p = 0; b = 0;
  for(unsigned int i=0; i<patPhotons->size(); ++i) {
     p1 = &((*patPhotons)[i]);
     // Cuts on H/E, Et and Eta
     if (p1->hadronicOverEm()<photon_he && p1->hasPixelSeed()==false  && p1->et() > photon_ptcut && ((p1->isEE() && fabs(p1->eta()) < 2.5) || p1->isEB())) {
       // Cut on sigma_IEtaIEta
       double p_sIEIE = p1->sigmaIetaIeta();
       double p_sIPIP = 0;
       // if ((p1->isEB() && p_sIEIE < sigma_barrel) || (p1->isEE() && p_sIEIE < sigma_endcap)) {
       double feta = fabs(p1->superCluster()->eta());
       double peta = fabs(p1->eta());
       if((peta < 1.379 && p_sIEIE < sigma_barrel) || (peta > 1.579 && p_sIEIE < sigma_endcap)) {
	 // Apparently isEB() and isEE() do not work properly?
	 // if((fabs(p1->eta())>1.379 && p1->isEB()) || (fabs(p1->eta())<1.579 && p1->isEB()) ) {
	 // std::cout<<"Photon: Eta = "<<p1->eta()<<" isEB() = "<<p1->isEB()<<" isEE() = "<<p1->isEE()<<" isEBEEGap() = "<<p1->isEBEEGap()<<" isEBGap() = "<<p1->isEBGap()<<"isEEGap() = "<<p1->isEEGap()<<std::endl;
	 // std::cout<<"Photon: Eta = "<<p1->eta()<<" SuperCluster Eta = "<<p1->superCluster()->eta()<<" SuperCluster Position Eta = "<<p1->superCluster()->position().eta()<<std::endl;
	 // }
	 // We only work with Isolated Photons as defined by Egamma POG:
	 double photon_pt = p1->pt(), photon_et = p1->et();
	 float photon_holTrkIso = p1->trkSumPtHollowConeDR04();
	 float photon_ecalIso = p1->ecalRecHitSumEtConeDR04();
	 float photon_hcalIso = p1->hcalTowerSumEtConeDR04();
	 // Traditional Isolation trk < 2.2 ecal < 4.2 hcal < 2.0
	 // if(photon_holTrkIso < (egm_trk + egm_trk_prop*photon_pt) && photon_ecalIso < (egm_ecal + egm_ecal_prop*photon_et) && photon_hcalIso < (egm_hcal + egm_hcal_prop*photon_et)) {
	 // New Isolation trk + ecal + hcal < 5.0
	 if((photon_holTrkIso + photon_ecalIso + photon_hcalIso < 5.0) {
	   // Some Final Spike Cleaning
	   // ECAL Lazy Tools: sigma_IEtaIEta, sigma_IEtaIPhi, sigma_IPhiIPhi and eMax
	   if(SpikeCleaning) {
	     // Ecal RecHits
	     edm::Handle<EcalRecHitCollection> EBReducedRecHits;
	     iEvent.getByLabel("reducedEcalRecHitsEB", EBReducedRecHits);
	     edm::Handle<EcalRecHitCollection> EEReducedRecHits;
	     iEvent.getByLabel("reducedEcalRecHitsEE", EEReducedRecHits);
	     // get the channel status from the DB
	     edm::ESHandle<EcalChannelStatus> chStatus;
	     iSetup.get<EcalChannelStatusRcd>().get(chStatus);
	     // ECAL Lazy Tool
	     EcalClusterLazyTools lazyTool(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));
	     
	     const reco::CaloClusterPtr  seed = p1->superCluster()->seed();
	     std::vector<float> viCov = lazyTool.localCovariances(*seed);
	     float p_sIPIP = sqrt(viCov[2]);
	     float eMax = lazyTool.eMax(*seed), e2nd = lazyTool.e2nd(*seed), e3x3 = lazyTool.e3x3(*seed);
	     float eLeft = lazyTool.eLeft(*seed), eRight = lazyTool.eRight(*seed), eTop = lazyTool.eTop(*seed), eBottom = lazyTool.eBottom(*seed);
	     float e2overe9 = (e3x3)? 0: (eMax + e2nd)/e3x3;
	     // ECAL Severity Level: severity and swissCross
	     DetId id     = lazyTool.getMaximum(*seed).first;
	     float seedTime  = -999.;
	     int   seedSeverity = -1;
	     const EcalRecHitCollection & rechits = ( p1->isEB() ? *EBReducedRecHits : *EEReducedRecHits);
	     EcalRecHitCollection::const_iterator it = rechits.find( id );
	     if( it != rechits.end() ) {
	       seedSeverity = EcalSeverityLevelAlgo::severityLevel( id, rechits, *chStatus );
	       seedTime = it->time();
	     }
	     if((p1->isEB() && e2overe9 < sc_spike && seedSeverity != 3) || p1->isEE()) {
	       ++count_good_photons;
	       if (count_good_photons == 1) p = &(*p1);
	       else if (p1->et() > p->et()) p = &(*p1);
	       else {}
	     }
	   }
	   else {
	     if((p1->isEB() && p_sIEIE > sigma_spike) || p1->isEE()) {
	       ++count_good_photons;
	       if (count_good_photons == 1) p = &(*p1);
	       else if (p1->et() > p->et()) p = &(*p1);
	       else {}
	     }
	   }                      
	 }
       }
     }
  }
  if (p != 0) {
    Identity = 0;
    PhotonFound = 1;
    photon = p;
    double photon_pt, photon_eta, photon_phi;
    photon_pt = p->pt(); photon_eta = p->eta(); photon_phi = p->phi();
    double p_sIEIE = p->sigmaIetaIeta();
    // const reco::CaloClusterPtr  seed = p->superCluster()->seed();
    // std::vector<float> viCov = lazyTool.localCovariances(*seed);
    // float p_sIPIP = sqrt(viCov[2]);
    double p_sIPIP = 0;
    double p_holTrkIso = p->trkSumPtHollowConeDR04(), p_ecalIso = p->ecalRecHitSumEtConeDR04(), p_hcalIso = p->hcalTowerSumEtConeDR04();
    double p_pt = p->pt(), p_et = p->et(), p_eta = p->eta(), p_phi = p->phi();
    double p_maxholTrkIso = egm_trk + egm_trk_prop*p_pt, p_maxecalIso = egm_ecal + egm_ecal_prop*p_et, p_maxhcalIso = egm_hcal + egm_hcal_prop*p_et;
    // If MC --> who's your mother?
    //           1) who is the related genparticle?
    //           2) who is the mother of the genparticle?
    g = 0; b = 0;
    if(!Data) {
      edm::Handle<reco::GenParticleCollection>      genParticles;
      iEvent.getByLabel(GenParticles, genParticles);
      // Print Status 3 Gen Particle Collection
      /*
      for(unsigned int i=0; i<genParticles->size(); ++i) {
	g1 = &((*genParticles)[i]);
	if(g1->status() != 3) continue;
	std::cout<<"Gen Particle: | index = "<<std::setw(5)<<i<<" | id = "<<std::setw(5)<<g1->pdgId()<<" | st = "<<std::setw(5)<<g1->status()<<" | pt = ";
	std::cout<<std::setw(12)<<g1->pt()<<" GeV/c | et = "<<std::setw(12)<<g1->et()<<" GeV | eta = "<<std::setw(12)<<g1->eta()<<" | phi = "<<std::setw(12)<<g1->phi()<<std::endl;
      }
      */
      const reco::GenParticle & gen = (reco::GenParticle &) * (p->genPhoton());
      g = &gen;
      double p_mot = 0.0; // distance photon - mother 
      if (g != 0) {
	Fill(histcontainers[0].get("Event","PATMatch"), n_ps, n1_ps, n2_ps, 1);
	// Fill(histcontainers[0].get("Event","Distance"), n_dist, n1_dist, n2_dist, reco::deltaR(g->eta(),g->phi(),p->eta(),p->phi()));
	// Matching Debug Stream
	dss<<"GEN Photon found: pt = "<<std::setw(12)<<g->pt();
	dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	dss<<"PAT Photon found: pt = "<<std::setw(12)<<p->pt();
	dss<<" GeV/c | et = "<<std::setw(12)<<p->et()<<" GeV | eta = "<<std::setw(12)<<p->eta()<<" | phi = "<<std::setw(12)<<p->phi()<<std::endl; 
	// If Electron
	if(gen.pdgId() == 11 || gen.pdgId() == -11)     { 
	  Identity = 2; 
	  double p_e = reco::deltaR(g->eta(),g->phi(),p->eta(),p->phi());
	  Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Mistag"), n_pd, n1_pd, n2_pd,p_e);
	}
	// If Photon
	if(gen.pdgId() == 22) {
	  // mother?
	  const reco::GenParticle & mot = (reco::GenParticle &) *(g->mother());
	  // const reco::GenParticle * motptr = (reco::GenParticle *) g->mother();
	  dss<<"GenParticle | id = "<<std::setw(5)<<g->pdgId()<<" | st = "<<std::setw(5)<<g->status()<<" | pt = "<<std::setw(12)<<g->pt();
	  dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	  if(&mot != 0) {
	    dss<<"MotherPart  | id = "<<std::setw(5)<<mot.pdgId()<<" | st = "<<std::setw(5)<<mot.status()<<" | pt = "<<std::setw(12)<<mot.pt();
	    dss<<" GeV/c | et = "<<std::setw(12)<<mot.et()<<" GeV | eta = "<<std::setw(12)<<mot.eta()<<" | phi = "<<std::setw(12)<<mot.phi()<<std::endl;
	    // std::cout<<"Mother Pointer: | index = "<<std::setw(5)<<g->motherRef().key()<<" | id = "<<std::setw(5)<<motptr->pdgId()<<" | st = "<<std::setw(5)<<motptr->status()<<" | pt = ";
	    //std::cout<<std::setw(12)<<motptr->pt()<<" GeV/c | et = "<<std::setw(12)<<motptr->et()<<" GeV | eta = "<<std::setw(12)<<motptr->eta()<<" | phi = "<<std::setw(12)<<motptr->phi()<<std::endl;
	    p_mot = reco::deltaR(mot.eta(),mot.phi(),p->eta(),p->phi());
	    // Distance to closesth parton
	    double r_g_mot = 0.0, minr_g_mot = 0.0;
	    for(unsigned int i=0; i<genParticles->size(); ++i) {
	      g1 = &((*genParticles)[i]);
	      if (g1->status() != 3) continue;
	      r_g_mot = reco::deltaR(mot.eta(), mot.phi(), g1->eta(), g1->phi());
	      if(i==0) minr_g_mot = r_g_mot;
	      if(r_g_mot < minr_g_mot) minr_g_mot = r_g_mot;
	    }
	    Fill(histcontainers[0].get("Event", "DR_Parton_Photon_All"), n_pd, n1_pd, n2_pd, minr_g_mot);


	    // Pythia 6.4 manual p67
	    // ISTHEP(IHEP): status code for entry IHEP, with the following meanings:
	    // = 0 : null entry.
	    // = 1 : an existing entry, which has not decayed or fragmented. This is the main class of entries, which represents final  given by the generator.
	    // = 2 : an entry which has decayed or fragmented and is therefore not appearing in the final state, but is retained for event history information.
	    // = 3 : a documentation line, defined separately from the event history. This could include the two incoming reacting particles, etc
	    
	    // PDG chap 34 Monte Carlo Particle Numbering Scheme
	    // Quarks: 1-6          (duscbd)
	    // Leptons: 11-16       (e- nu_e u-, nu_u, t-, nu_t)
	    // Gauge Bosons: 21-24  (gluon, gamma, Z0, W+)
	    // mesons:  pi0 = 111, pi+ = 211, eta = 221
	    // baryons: p = 2212
	    
	    // My Definition for IdIsoLepton
	    // 0 = IdIsoLepton = N/A
	    // 1 = Direct Photon
	    // 2 = Mistagged Elektron
	    // 3 = Fragmentation Photon
	    // 4 = Decay Photons
	    // 5 = Ewk
	    // 6 = ISR
	    // 7 = FSR
	    // Lets take ISR and FSR together so far

	    if      (mot.pdgId() == 22) { 
	      Identity = 1; b = &mot; 
	      Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Direct"), n_pd, n1_pd, n2_pd, p_mot); DedicatedPlots(iEvent, iSetup,b);
	    }
	    else if (g->motherRef().key() == 2 || g->motherRef().key() == 3)          { Identity = 6; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_ISR"), n_pd, n1_pd, n2_pd, p_mot);}
	    else if (abs(mot.pdgId()) < 25 && mot.pdgId() != 22 && mot.status() == 3) { Identity = 7; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_FSR"), n_pd, n1_pd, n2_pd, p_mot);}
	    else if (abs(mot.pdgId()) < 10 || mot.pdgId() == 21)                      { Identity = 3; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Fragment"), n_pd, n1_pd, n2_pd, p_mot);}
	    else if (mot.status() == 2 && abs(mot.pdgId()) > 25)                      { Identity = 4; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Decay"), n_pd, n1_pd, n2_pd, p_mot);}
	    else if (abs(mot.pdgId()) >10 && mot.pdgId() != 21 && mot.pdgId() != 22)  { Identity = 5; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Ewk"), n_pd, n1_pd, n2_pd, p_mot);}
	    else {}
	  }
	}
	dss<<"PAT MATCH directPhoton  = "<<(Identity==1)<<" | electroPhoton = "<<(Identity==2)<<" | fragPhoton = "<<(Identity==3)<<" | decayPhoton = "<<(Identity==4);
	dss<<" | ewkPhoton = "<<(Identity==5)<<" | ISR/FSR Photon = "<<(Identity==6)<<std::endl;
      }
      else { 
	ErrorFlags->Fill(1);
	Fill(histcontainers[0].get("Event","PATMatch"), n_ps, n1_ps, n2_ps, 0);
	dss<<"No GEN Photon found by the PAT that matches the PAT Photon, trying manually now"<<std::endl; 
	dss<<"... Not yet implemented"<<std::endl;
	dss<<"... Assumed to be a non Prompt Photon"<<std::endl;
	dss<<"... MANUAL ..."<<std::endl;
	bool secondtry = false;
	for(unsigned int i=0; i<genParticles->size(); ++i) {
	  g1 = &((*genParticles)[i]);
	  if (g1->pdgId() == 22) {
	    // mother?
	    if(reco::deltaR(g1->eta(),g1->phi(),p->eta(),p->phi()) < 0.5 && g1->pt()> 10) {
	      g = *(&g1);
	      dss<<"Gen Photon Candidate | id = "<<std::setw(5)<<g->pdgId()<<" | st = "<<std::setw(5)<<g->status()<<" | pt = "<<std::setw(12)<<g->pt();
	      dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	      const reco::GenParticle & mot = (reco::GenParticle &) *(g->mother());
	      // const reco::GenParticle * motptr = (reco::GenParticle *) g->mother();
	      if(&mot != 0) {
		dss<<" Photon Match Found! "<<std::endl;
		dss<<"MotherPart  | id = "<<std::setw(5)<<mot.pdgId()<<" | st = "<<std::setw(5)<<mot.status()<<" | pt = "<<std::setw(12)<<mot.pt();
		dss<<" GeV/c | et = "<<std::setw(12)<<mot.et()<<" GeV | eta = "<<std::setw(12)<<mot.eta()<<" | phi = "<<std::setw(12)<<mot.phi()<<std::endl;
		// std::cout<<"Mother Pointer: | index = "<<std::setw(5)<<g->motherRef().key()<<" | id = "<<std::setw(5)<<motptr->pdgId()<<" | st = "<<std::setw(5)<<motptr->status()<<" | pt = ";
		// std::cout<<std::setw(12)<<motptr->pt()<<" GeV/c | et = "<<std::setw(12)<<motptr->et()<<" GeV | eta = "<<std::setw(12)<<motptr->eta()<<" | phi = "<<std::setw(12)<<motptr->phi()<<std::endl;
		// classiffy
		if      (mot.pdgId() == 22) { Identity = 1; b = &mot; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Direct"), n_pd, n1_pd, n2_pd, p_mot); DedicatedPlots(iEvent, iSetup, b);}
		else if (g->motherRef().key() == 2 || g->motherRef().key() == 3)          { Identity = 6; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_ISR"), n_pd, n1_pd, n2_pd, p_mot);}
		else if (abs(mot.pdgId()) < 25 && mot.pdgId() != 22 && mot.status() == 3) { Identity = 7; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_FSR"), n_pd, n1_pd, n2_pd, p_mot);}
		else if (abs(mot.pdgId()) < 10 || mot.pdgId() == 21)                      { Identity = 3; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Fragment"), n_pd, n1_pd, n2_pd, p_mot);}
		else if (mot.status() == 2 && abs(mot.pdgId()) > 25)                      { Identity = 4; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Decay"), n_pd, n1_pd, n2_pd, p_mot);}
		else if (abs(mot.pdgId()) >10 && mot.pdgId() != 21 && mot.pdgId() != 22)  { Identity = 5; Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Ewk"), n_pd, n1_pd, n2_pd, p_mot);}
		else {}
		secondtry = true;
	      }
	    }
	  }
	  if(g1->pdgId() == 11 || g1->pdgId() == -11) {
	    if(reco::deltaR(g1->eta(),g1->phi(),p->eta(),p->phi()) < 0.5 && g1->pt()> 70) {
	      g = *(&g1);
	      dss<<"Gen Electron Candidate | id = "<<std::setw(5)<<g->pdgId()<<" | st = "<<std::setw(5)<<g->status()<<" | pt = "<<std::setw(12)<<g->pt();
	      dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	      dss<<" Electron Match Found! "<<std::endl;
	      Identity = 2; // Elektro Photon
	      double p_e = reco::deltaR(g->eta(),g->phi(),p->eta(),p->phi());
	      Fill(histcontainers[0].get("Event", "DR_Mother_Photon_Mistag"), n_pd, n1_pd, n2_pd,p_e);
	      secondtry = true;
	    }
	  }
	}
	// if(decayPhoton==1 && fragPhoton==1) {fragPhoton=0;}
	if(secondtry == false) { 
	  dss<<"No Idea"<<std::endl; ErrorFlags->Fill(2); 
	  int count_cands = 0;
	  // Loop over all genparticles with pt > 10 and Delta R < 0.4
	  for(unsigned int i=0; i<genParticles->size(); ++i) {
	    g1 = &((*genParticles)[i]);
	    if(g1->pt()<10) continue;
	    double dr = reco::deltaR(g1->eta(),g1->phi(),p->eta(),p->phi());
	    if(dr > 1.0) continue;
	    dss<<"Close Gen Particle | id = "<<std::setw(5)<<g1->pdgId()<<" | st = "<<std::setw(5)<<g1->status()<<" | dr = "<<std::setw(9)<<dr<<" | pt = "<<std::setw(12)<<g1->pt();
	    dss<<" GeV/c | et = "<<std::setw(12)<<g1->et()<<" GeV | eta = "<<std::setw(12)<<g1->eta()<<" | phi = "<<std::setw(12)<<g1->phi()<<std::endl;
	    if(g1->pdgId()==22) {
	      const reco::GenParticle & mot1 = (reco::GenParticle &) *(g1->mother());
	      if(&mot1!=0) {
		dss<<"MotherPart  | id = "<<std::setw(5)<<mot1.pdgId()<<" | st = "<<std::setw(5)<<mot1.status()<<" | pt = "<<std::setw(12)<<mot1.pt();
		dss<<" GeV/c | et = "<<std::setw(12)<<mot1.et()<<" GeV | eta = "<<std::setw(12)<<mot1.eta()<<" | phi = "<<std::setw(12)<<mot1.phi()<<std::endl;
	      }
	    }
	    if(g1->status()==3) {
	      reco::GenParticle::daughters d = g1->daughterRefVector();
	      std::cout<<"Status 3 particle has "<<d.size()<<" daughter particles"<<std::endl;
              for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
                dss<<"DaughterPart  | id = "<<std::setw(5)<<(**it_d).pdgId()<<" | st = "<<std::setw(5)<<(**it_d).status()<<" | pt = "<<std::setw(12)<<(**it_d).pt();
                dss<<" GeV/c | et = "<<std::setw(12)<<(**it_d).et()<<" GeV | eta = "<<std::setw(12)<<(**it_d).eta()<<" | phi = "<<std::setw(12)<<(**it_d).phi()<<std::endl;
	      }
	    }
	    if(g1->status()==1 && dr < 0.1) {
	      if(g1->pdgId() > 100) Identity = 4; dss<<"Put Identity = 4 (Decay) "<<std::endl;
	    }
	    ++count_cands;
	  }
	  std::cout<<count_cands<<" candidates found"<<std::endl;
	  if(count_cands == 1) {
	    for(unsigned int i=0; i<genParticles->size(); ++i) {
	      g1 = &((*genParticles)[i]);
	      double dr = reco::deltaR(g1->eta(),g1->phi(),p->eta(),p->phi());
	      dss<<"Gen Particle | id = "<<std::setw(5)<<g1->pdgId()<<" | st = "<<std::setw(5)<<g1->status()<<" | dr = "<<dr<<" | pt = "<<std::setw(12)<<g1->pt();
	      dss<<" GeV/c | et = "<<std::setw(12)<<g1->et()<<" GeV | eta = "<<std::setw(12)<<g1->eta()<<" | phi = "<<std::setw(12)<<g1->phi()<<std::endl;
	      Identity = 4; // set strange event to decay (need to investigate this @ sim level)
	    }
	  }
	}
	dss<<"PAT MATCH directPhoton  = "<<(Identity==1)<<" | electroPhoton = "<<(Identity==2)<<" | fragPhoton = "<<(Identity==3)<<" | decayPhoton = "<<(Identity==4);
        dss<<" | ewkPhoton = "<<(Identity==5)<<" | ISR/FSR Photon = "<<(Identity==6)<<std::endl;
	print_dss=1;
      }
    }
    // Count before we throw away events because we can't find them 
    ++Count[7];
    int PhotonJetIndex = -1, ClosesthJetIndex = -1, AllRClosesthJetIndex = -1; 
    dss<<"PAT Photon we work with: | pt = "<<std::setw(12)<<p->pt()<<" GeV/c | et = "<<std::setw(12)<<p->et()<<" GeV | eta = "<<std::setw(12)<<p->eta()<<" | phi = "<<std::setw(12)<<p->phi()<<std::endl;
    // dss<<"Isolation values for this photon:      | trk = "<<std::setw(12)<<p_holTrkIso<<" GeV/c | ecal = "<<std::setw(12)<<p_ecalIso<<" GeV | hcal = "<<std::setw(12)<<p_hcalIso<<" GeV"<<std::endl;
    // dss<<"Isolation max values for this photon : | trk = "<<std::setw(12)<<p_maxholTrkIso<<" GeV/c | ecal = "<<std::setw(12)<<p_maxecalIso<<" GeV | hcal = "<<std::setw(12)<<p_maxhcalIso<<" GeV"<<std::endl;
    // dss<<"         Would it have been Isolated?: | Iso? = "<<(p_holTrkIso < p_maxholTrkIso && p_ecalIso < p_maxecalIso  && p_hcalIso < p_maxhcalIso)<<std::endl;
    // Make a Good Jet Selection (removing the Photon) 

    double jetpt_0104 = 0.0, jetpt_01 = 0.0;
    int n_jets_in0104 = 0, n_jets_in01 = 0; 
    for(unsigned int i=0; i<recoJets->size(); ++i) {
      r = &((*recoJets)[i]);
      double jet_pt, jet_eta, jet_phi;
      jet_pt = r->pt();    jet_eta = r->eta();    jet_phi = r->phi();
      double r_pi = reco::deltaR(photon_eta, photon_phi,jet_eta, jet_phi);
      if(i==0) {min_allr_pi = r_pi; AllRClosesthJetIndex = 0;}
      else if (min_allr_pi > r_pi) {min_allr_pi = r_pi; AllRClosesthJetIndex = i;}
      if(r->pt() > 30) {
	// dss<<"PF Jet; pt = "<<std::setw(12)<<r->pt()<<" GeV/c | et = "<<std::setw(12)<<r->et()<<" GeV | eta = "<<std::setw(12)<<r->eta()<<" | phi = "<<std::setw(12)<<r->phi();
	// dss<<"| dR = "<<std::setw(12)<<r_pi<<" | dPt/Pt = "<<std::setw(12)<<pt_pi<<std::endl;
	// Kinematic & Geometric Cut
	if( r_pi < 0.40 && r_pi >= 0.1) {
	  /*if (g!=0)*/ histcontainers[0].get("Event", "PtJet_PtPhoton_DeltaR04")->Fill(p->pt(),r->pt());
	  if (g!=0)     histcontainers[0].get("Event", "PtJet_GenPtPhoton_DeltaR04")->Fill(g->pt(),r->pt());
	  Fill(histcontainers[0].get("Event", "DR_Jet_Photon_DeltaR04"), n_rdd, n1_rdd, n2_rdd, r_pi, evW);
	  jetpt_0104 += r->pt();
	  ++n_jets_in0104;
	  inDR0104 = true;
	  std::cout<<"Event in DR0104 Cone: photon pt = "<<p->pt()<<" pfjet pt = "<<r->pt()<<std::endl;
	}
	if(r_pi < 0.10 ) {
	  PhotonJetIndex = i; 
	  dss<<"Jet == Photon Found"<<std::endl;
	  /*if (g!=0)*/ histcontainers[0].get("Event", "PtJet_PtPhoton_DeltaR01")->Fill(p->pt(),r->pt());
	  if (g!=0)     histcontainers[0].get("Event", "PtJet_GenPtPhoton_DeltaR01")->Fill(g->pt(),r->pt());
	  Fill(histcontainers[0].get("Event", "DR_Jet_Photon_DeltaR01"), n_rdd, n1_rdd, n2_rdd, r_pi, evW);
	  jetpt_01 += r->pt();
	  /*if (g!=0)*/ histcontainers[0].get("Event", "PtJet_PtPhoton_All")->Fill(p->pt(),r->pt()); // All: found or not found
	  histcontainers[0].get("Event", "PtJet_PtPhoton_F")->Fill(p->pt(),r->pt());   // Found
	  Fill(histcontainers[0].get("Event", "DR_Jet_Photon_All"), n_rdd, n1_rdd, n2_rdd,  r_pi, evW); // All: Found or not found
	  Fill(histcontainers[0].get("Event", "DR_Jet_Photon_F"), n_rdd, n1_rdd, n2_rdd, r_pi, evW); // Found
	  if (g!= 0) histcontainers[0].get("Event", "PtJet_GenPtPhoton_All")->Fill(g->pt(),r->pt()); // All: found or not found
	  if (g!= 0) histcontainers[0].get("Event", "PtJet_GenPtPhoton_F")->Fill(g->pt(),r->pt());   // Found
	  ++n_jets_in01;
	}
	if(r_pi < 0.10 ) { /*nothing*/ }
	else {
	  Jets.push_back(r);
	  // Delta R and Delta P before Photon Jet Cleaning
	  if(i==0) {min_r_pi = r_pi; ClosesthJetIndex = 0;}
	  else if (min_r_pi > r_pi) {min_r_pi = r_pi; ClosesthJetIndex = i;}
	  else {}
	}
      }
    }
    if(ClosesthJetIndex != -1 && ((unsigned) ClosesthJetIndex) < recoJets->size()) {
      r = &((*recoJets)[ClosesthJetIndex]);
      std::cout<<"Distance Photon - Closesth Jet: "<<min_allr_pi<<"  photon pt = "<<p->pt()<<" pfjet pt = "<<r->pt()<<std::endl;
    }
    if(jetpt_0104 > 0.0) histcontainers[0].get("Event", "TotJetPt_PtPhoton_DeltaR04")->Fill(p->pt(),jetpt_0104);
    if(jetpt_01 > 0.0)   histcontainers[0].get("Event", "TotJetPt_PtPhoton_DeltaR01")->Fill(p->pt(),jetpt_01);
    Fill(histcontainers[0].get("Event", "N_Jets_in_DeltaR04"), n_mult, n1_mult, n2_mult, n_jets_in0104);
    Fill(histcontainers[0].get("Event", "N_Jets_in_DeltaR01"), n_mult, n1_mult, n2_mult, n_jets_in01);
    if(PhotonJetIndex == -1) {  // !!! This is a problem !!!
                                // Actually I should throw these events away ... 
                                // Setting p = 0 here would not help, some histograms are already filled
                                // maybe make a separate Directory for this
      dss<<"Jet == Photon *Not* Found"<<std::endl; print_dss = 1; 
      ErrorFlags->Fill(3);
      if(ClosesthJetIndex != -1) { 
	r = &((*recoJets)[ClosesthJetIndex]);
	dss<<"Closesth Jet:  pt = "<<std::setw(12)<<r->pt()<<" GeV/c | et = "<<std::setw(12)<<r->et()<<" GeV | eta = "<<std::setw(12)<<r->eta()<<" | phi = "<<std::setw(12)<<r->phi()<<std::endl; 
	/* if (g!= 0)*/ histcontainers[0].get("Event", "PtJet_PtPhoton_All")->Fill(p->pt(),r->pt()); // All: found or not found
	histcontainers[0].get("Event", "PtJet_PtPhoton_NF")->Fill(p->pt(),r->pt());  // Not found
	Fill(histcontainers[0].get("Event", "DR_Jet_Photon_All"), n_rdd, n1_rdd, n2_rdd,  min_r_pi, evW); // All: Found or not found
	Fill(histcontainers[0].get("Event", "DR_Jet_Photon_NF"), n_rdd, n1_rdd, n2_rdd, min_r_pi, evW); // Not Found
	if (g!= 0) histcontainers[0].get("Event", "PtJet_GenPtPhoton_All")->Fill(g->pt(),r->pt()); // All: found or not found
	if (g!= 0) histcontainers[0].get("Event", "PtJet_GenPtPhoton_NF")->Fill(g->pt(),r->pt());  // Not found
      }
      else { dss<<"No Jets Found!"<<std::endl; ErrorFlags->Fill(4);} 
    }
    // Fill(histcontainers[0].get("Event", "DR_Before"), n_dist, n1_dist, n2_dist, min_r_pi, evW);                                 
    min_r_pi = 0;
    for(unsigned int i=0; i<Jets.size(); ++i) {
      r = Jets[i];
      // in case you want all jets of all events ... 
      // histcontainers[0].get("Event", "PtJet_PtPhoton_All")->Fill(p->pt(),r->pt()); // All: found or not found
      // histcontainers[0].get("Event", "PtJet_PtPhoton_F")->Fill(p->pt(),r->pt());   // Found                                  
      double jet_pt, jet_eta, jet_phi;
      jet_pt = r->pt();    jet_eta = r->eta();    jet_phi = r->phi();
      double r_pi = reco::deltaR(photon_eta, photon_phi, jet_eta,jet_phi);
      // Delta R and Delta P before Photon Jet Cleaning
      if(i==0) {min_r_pi = r_pi;}
      else if (min_r_pi > r_pi) {min_r_pi = r_pi;}
      else {}
    }
    // Fill(histcontainers[0].get("Event", "DR_After"), n_dist, n1_dist, n2_dist, min_r_pi, evW);                             
    Fill(histcontainers[0].get("Event", "nVertEvt"), n_mult, n1_mult, n2_mult, nVertices);

    // Plots
    FillHistograms(iEvent, iSetup, 0, photon, Jets, evW); dss<<"Photon + 2J  Filled"<<std::endl;   
    if(Identity==1) { FillHistograms(iEvent, iSetup, 1, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_Direct"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"Direct Photon"<<std::endl;}
    if(Identity==2) { FillHistograms(iEvent, iSetup, 2, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_Mistag"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"Mistagged Electron"<<std::endl;}
    if(Identity==3) { FillHistograms(iEvent, iSetup, 3, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_Fragment"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"Fragmentation Photon"<<std::endl;}
    if(Identity==4) { FillHistograms(iEvent, iSetup, 4, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_Decay"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"Decay Photon"<<std::endl;}
    if(Identity==5) { FillHistograms(iEvent, iSetup, 3, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_Ewk"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"Ewk Photon"<<std::endl;}
    if(Identity==6) { FillHistograms(iEvent, iSetup, 3, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_ISR"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"ISR Photon"<<std::endl;}
    if(Identity==7) { FillHistograms(iEvent, iSetup, 3, photon, Jets, evW); Fill(histcontainers[0].get("Event", "DR_Jet_Photon_FSR"), n_rdd, n1_rdd, n2_rdd, min_allr_pi); dss<<"FSR Photon"<<std::endl;}
    if(p->isEB())                   { FillHistograms(iEvent, iSetup, 5, photon, Jets, evW);  ++Count[11]; }  // 06_Egamma_Barrel
    if(p->isEE())                   { FillHistograms(iEvent, iSetup, 6, photon, Jets, evW);  ++Count[12]; }  // 07_Egamma_Endcap
    if(p->pt() <= 130)              { FillHistograms(iEvent, iSetup, 10,photon, Jets, evW);}
    if(p->pt() <= 260)              { FillHistograms(iEvent, iSetup, 11,photon, Jets, evW);}
    if(p->pt() >  260)              { FillHistograms(iEvent, iSetup, 12,photon, Jets, evW);}

    timer.Stop();
    tss<<" PhotonAnalyzer :: Selection :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
    histcontainers[0].get("Timing","tSel_Real")->Fill(timer.RealTime());
    histcontainers[0].get("Timing","tSel_CPU")->Fill(timer.CpuTime());
    std::cout<<"Selection: Identity = "<<Identity<<std::endl;
  }
}

// ------------ Promptch Photons in Case of Monte Carlo -----

// ------------ Make Table ------------                                                                                                                                                                  
void 
PhotonAnalyzer::RA2Selection(const edm::Event& iEvent, const edm::EventSetup& iSetup, const pat::Photon * & photon, std::vector<const pat::Jet*>& Jets, double evW, int Identity)
{
  std::cout<<"RA2Selection: Identity = "<<Identity<<std::endl;
  // Variables
  double ppt = photon->pt();
  double peta = fabs(photon->eta());
  bool bar = photon->isEB();
  bool end = photon->isEE();
  int n_jets = 0, n_ra2  = 0;
  double ht  = 0.0;
  reco::MET::LorentzVector mht(0,0,0,0);
  double min_phi_pi = 0.0, min_r_pi = 0.0;
  int ClosesthJetIndex = -1;
  rss<<"Jets with pt > 50 and eta < 2.5):"<<std::endl;
  for(unsigned int i=0; i<Jets.size(); ++i) {
    r = Jets[i];
    if(Jets[i]->pt() > 50 && fabs(Jets[i]->eta()) < 2.50) {
      ++n_jets;
      ht  += Jets[i]->pt();
    }
    if(Jets[i]->pt() > 30) { 
      mht  -= Jets[i]->p4();
      double phi_pi = reco::deltaPhi(photon->phi(),r->phi());
      double r_pi = reco::deltaR(photon->phi(), photon->eta(), r->phi(), r->eta());
      if(i==0) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = 0;}
      else if (min_phi_pi > phi_pi) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = i;}
      else {}
    }
  }
  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  double mht_value = MHT.pt();
  double mht_j1 = 0, mht_j2 = 0, mht_j3 = 0;
  if(n_jets >= 3) {
    r1 = Jets[0]; r2 = Jets[1]; r3 = Jets[2];
    mht_j1 = fabs(reco::deltaPhi(r1->phi(),MHT.phi()));
    mht_j2 = fabs(reco::deltaPhi(r2->phi(),MHT.phi()));
    mht_j3 = fabs(reco::deltaPhi(r3->phi(),MHT.phi()));
  }
  if(n_jets == 2) {
    r1 = Jets[0]; r2 = Jets[1];
    mht_j1 = fabs(reco::deltaPhi(r1->phi(),MHT.phi()));
    mht_j2 = fabs(reco::deltaPhi(r2->phi(),MHT.phi()));
  }
  double ph_jcl = fabs(min_phi_pi);
  double ph_rjc = fabs(min_r_pi);
  double ph_mht = fabs(reco::deltaPhi(photon->phi(), MHT.phi()));
  rss<<" Amount of Jets passing pt and eta criteria: "<<n_jets<<" | ht of event = "<<ht<< " GeV/c | mht of event = "<<mht_value<<" GeV/c"<<std::endl;
  // if (n_jets > 0) { ++Count[7];}
  if (n_jets > 1) { ++Count[8];}
  if (n_jets > 2) { ++Count[9];}

  // Cuts and Plots
  // ++RA2Count[2]; wordt niet geteld
  if(ppt > 100) {
    Fill(JetMult_All, n_mult, n1_mult, n2_mult, n_jets);
    Fill(PT_All,      n_et, n1_et, n2_et, ppt);
    Fill(HT_All,      n_ht, n1_ht, n2_ht, ht);
    Fill(MHT_All,     n_et, n1_et, n2_et, mht_value);
    Fill(MHT_J1_All,  n_pd, n1_pd, n2_pd, mht_j1);
    Fill(MHT_J2_All,  n_pd, n1_pd, n2_pd, mht_j2);
    Fill(MHT_J3_All,  n_pd, n1_pd, n2_pd, mht_j3);
    RA2Book (Book_All, bar, end, ppt, Identity);
    RA2Book1(Book1_All, bar, end, ppt);
    RA2Book2(Book2_All, peta, ppt);

    Fill(PH_JCl_All,  n_pd, n1_pd, n2_pd, ph_jcl);
    Fill(PH_RJC_All,  n_rd, n1_rd, n2_rd, ph_rjc);
    Fill(PH_MHT_All,  n_pd, n1_pd, n2_pd, ph_mht);
    // N Jets > 3
    if(n_jets >= AmountOfJets) {
      ++RA2Count[2];
      Fill(JetMult_AJC, n_mult, n1_mult, n2_mult, n_jets);
      Fill(PT_AJC,      n_et, n1_et, n2_et, ppt);
      Fill(HT_AJC,      n_ht, n1_ht, n2_ht, ht);
      Fill(MHT_AJC,     n_et, n1_et, n2_et, mht_value);
      Fill(MHT_J1_AJC,  n_pd, n1_pd, n2_pd, mht_j1);
      Fill(MHT_J2_AJC,  n_pd, n1_pd, n2_pd, mht_j2);
      Fill(MHT_J3_AJC,  n_pd, n1_pd, n2_pd, mht_j3);
      RA2Book (Book_AJC, bar, end, ppt, Identity);
      RA2Book1(Book1_AJC, bar, end, ppt);
      RA2Book2(Book2_AJC, peta, ppt);
      Fill(PH_JCl_AJC,  n_pd, n1_pd, n2_pd, ph_jcl);
      Fill(PH_RJC_AJC,  n_rd, n1_rd, n2_rd, ph_rjc);
      Fill(PH_MHT_AJC,  n_pd, n1_pd, n2_pd, ph_mht);
      // HT > 300 GeV/c
      if (ht > 300) {
	++RA2Count[3];
	Fill(JetMult_AHC, n_mult, n1_mult, n2_mult, n_jets);
	Fill(PT_AHC,      n_et, n1_et, n2_et, ppt);
	Fill(HT_AHC,      n_ht, n1_ht, n2_ht, ht);
	Fill(MHT_AHC,     n_et, n1_et, n2_et, mht_value);
	Fill(MHT_J1_AHC,  n_pd, n1_pd, n2_pd, mht_j1);
	Fill(MHT_J2_AHC,  n_pd, n1_pd, n2_pd, mht_j2);
	Fill(MHT_J3_AHC,  n_pd, n1_pd, n2_pd, mht_j3);
	RA2Book (Book_AHC, bar, end, ppt, Identity);
	RA2Book1(Book1_AHC, bar, end, ppt);
	RA2Book2(Book2_AHC, peta, ppt);
	Fill(PH_JCl_AHC,  n_pd, n1_pd, n2_pd, ph_jcl);
	Fill(PH_RJC_AHC,  n_rd, n1_rd, n2_rd, ph_rjc);
	Fill(PH_MHT_AHC,  n_pd, n1_pd, n2_pd, ph_mht);
	// MHT Angular Cuts
	if (mht_j1 > 0.5 && mht_j2 > 0.5 && (AmountOfJets == 2 || (AmountOfJets == 3 && mht_j3 > 0.3))) {
	  ++RA2Count[4];
	  rss<<"Event passed Angular Cuts"<<std::endl;
	  Fill(JetMult_AAC, n_mult, n1_mult, n2_mult, n_jets);
	  Fill(PT_AAC,      n_et, n1_et, n2_et, ppt);
	  Fill(HT_AAC,      n_ht, n1_ht, n2_ht, ht);
	  Fill(MHT_AAC,     n_et, n1_et, n2_et, mht_value);
	  Fill(MHT_J1_AAC,  n_pd, n1_pd, n2_pd, mht_j1);
	  Fill(MHT_J2_AAC,  n_pd, n1_pd, n2_pd, mht_j2);
	  Fill(MHT_J3_AAC,  n_pd, n1_pd, n2_pd, mht_j3);
	  RA2Book (Book_AAC, bar, end, ppt, Identity);
	  RA2Book1(Book1_AAC, bar, end, ppt);
	  RA2Book2(Book2_AAC, peta, ppt);
	  Fill(PH_JCl_AAC,  n_pd, n1_pd, n2_pd, ph_jcl);
	  Fill(PH_RJC_AAC,  n_rd, n1_rd, n2_rd, ph_rjc);
	  Fill(PH_MHT_AAC,  n_pd, n1_pd, n2_pd, ph_mht);
	  // MHT > 150 GeV/c
	  if(mht_value > 150) {
	    ++RA2Count[5]; ++Count[10];
	    n_ra2 = 1;
	    Fill(JetMult_AMC, n_mult, n1_mult, n2_mult, n_jets);
	    Fill(PT_AMC,      n_et, n1_et, n2_et, ppt);
	    Fill(HT_AMC,      n_ht, n1_ht, n2_ht, ht);
	    Fill(MHT_AMC,     n_et, n1_et, n2_et, mht_value);
	    Fill(MHT_J1_AMC,  n_pd, n1_pd, n2_pd, mht_j1);
	    Fill(MHT_J2_AMC,  n_pd, n1_pd, n2_pd, mht_j2);
	    Fill(MHT_J3_AMC,  n_pd, n1_pd, n2_pd, mht_j3);
	    RA2Book (Book_AMC, bar, end, ppt, Identity);
	    RA2Book1(Book1_AMC, bar, end, ppt);
	    RA2Book2(Book2_AMC, peta, ppt);
	    Fill(PH_JCl_AMC,  n_pd, n1_pd, n2_pd, ph_jcl);
	    Fill(PH_RJC_AMC,  n_rd, n1_rd, n2_rd, ph_rjc);
	    Fill(PH_MHT_AMC,  n_pd, n1_pd, n2_pd, ph_mht);
	    Fill(histcontainers[0].get("Jets","HT"),  n_ht, n1_ht, n2_ht, ht, evW);
	    Fill(histcontainers[0].get("Jets","MHT"), n_et, n1_et, n2_et, mht_value, evW);
	    Fill(histcontainers[0].get("Jets","MHT_J1"), n_pd, n1_pd, n2_pd, mht_j1, evW);
	    Fill(histcontainers[0].get("Jets","MHT_J2"), n_pd, n1_pd, n2_pd, mht_j2, evW);
	    Fill(histcontainers[0].get("Jets","MHT_J3"), n_pd, n1_pd, n2_pd, mht_j3, evW);
	    PasBaseline = true;
	    std::cout<<"This event passed the baseline: photonpt = "<<ppt<<" ht = "<<ht<<" mht = "<<mht_value<<std::endl;
	    float photon_holTrkIso = photon->trkSumPtHollowConeDR04();
	    float photon_ecalIso = photon->ecalRecHitSumEtConeDR04();
	    float photon_hcalIso = photon->hcalTowerSumEtConeDR04();
	    if(photon_holTrkIso < egm_trk && photon_ecalIso < egm_ecal  && photon_hcalIso < egm_hcal) { // ok, this event is also in EGM Selection
	    }
	    else { 
	      FailEGM = true; std::cout<<"This event failed EGM ID: photonpt = "<<ppt<<" ht = "<<ht<<" mht = "<<mht_value<<std::endl; 
	      std::cout<<"Iso_TRK = "<<photon_holTrkIso<<" Iso_ECAL = "<<photon_ecalIso<<" Iso_HCAL = "<<photon_hcalIso<<" check: "<<photon->hcalTowerSumEtConeDR04()<<std::cout;
	    }
	    // Influence of Boson PT Cut (in case we lower pt cut of photon to 70 GeV or want to study the influence of Z boson pt contribution
	    if(ppt > 100) {
	      ++RA2Count[6];
	      Fill(JetMult_PTC, n_mult, n1_mult, n2_mult, n_jets);
	      Fill(PT_PTC,      n_et, n1_et, n2_et, ppt);
	      Fill(HT_PTC,      n_ht, n1_ht, n2_ht, ht);
	      Fill(MHT_PTC,     n_et, n1_et, n2_et, mht_value);
	      Fill(MHT_J1_PTC,  n_pd, n1_pd, n2_pd, mht_j1);
	      Fill(MHT_J2_PTC,  n_pd, n1_pd, n2_pd, mht_j2);
	      Fill(MHT_J3_PTC,  n_pd, n1_pd, n2_pd, mht_j3);
	      RA2Book (Book_PTC, bar, end, ppt, Identity);
	      RA2Book1(Book1_PTC, bar, end, ppt);
	      RA2Book2(Book2_PTC, peta, ppt);
	      Fill(PH_JCl_PTC,  n_pd, n1_pd, n2_pd, ph_jcl);
	      Fill(PH_RJC_PTC,  n_rd, n1_rd, n2_rd, ph_rjc);
	      Fill(PH_MHT_PTC,  n_pd, n1_pd, n2_pd, ph_mht);
	    }
	    // Hight HT Selection Region
	    if(ht > 500) {
	      ++RA2Count[7];
	      Fill(JetMult_HHS, n_mult, n1_mult, n2_mult, n_jets);
	      Fill(PT_HHS,      n_et, n1_et, n2_et, ppt);
	      Fill(HT_HHS,      n_ht, n1_ht, n2_ht, ht);
	      Fill(MHT_HHS,     n_et, n1_et, n2_et, mht_value);
	      Fill(MHT_J1_HHS,  n_pd, n1_pd, n2_pd, mht_j1);
	      Fill(MHT_J2_HHS,  n_pd, n1_pd, n2_pd, mht_j2);
	      Fill(MHT_J3_HHS,  n_pd, n1_pd, n2_pd, mht_j3);
	      RA2Book (Book_HHS, bar, end, ppt, Identity);
	      RA2Book1(Book1_HHS, bar, end, ppt);
	      RA2Book2(Book2_HHS, peta, ppt);
	      Fill(PH_JCl_HHS,  n_pd, n1_pd, n2_pd, ph_jcl);
	      Fill(PH_RJC_HHS,  n_rd, n1_rd, n2_rd, ph_rjc);
	      Fill(PH_MHT_HHS,  n_pd, n1_pd, n2_pd, ph_mht);
	    }
	    // High MHT Selection Region
	    if(mht_value > 250) {
	      ++RA2Count[8];
	      Fill(JetMult_HMS, n_mult, n1_mult, n2_mult, n_jets);
	      Fill(PT_HMS,      n_et, n1_et, n2_et, ppt);
	      Fill(HT_HMS,      n_ht, n1_ht, n2_ht, ht);
	      Fill(MHT_HMS,     n_et, n1_et, n2_et, mht_value);
	      Fill(MHT_J1_HMS,  n_pd, n1_pd, n2_pd, mht_j1);
	      Fill(MHT_J2_HMS,  n_pd, n1_pd, n2_pd, mht_j2);
	      Fill(MHT_J3_HMS,  n_pd, n1_pd, n2_pd, mht_j3);
	      RA2Book (Book_HMS, bar, end, ppt, Identity);
	      RA2Book1(Book1_HMS, bar, end, ppt);
	      RA2Book2(Book2_HMS, peta, ppt);
	      Fill(PH_JCl_HMS,  n_pd, n1_pd, n2_pd, ph_jcl);
	      Fill(PH_RJC_HMS,  n_rd, n1_rd, n2_rd, ph_rjc);
	      Fill(PH_MHT_HMS,  n_pd, n1_pd, n2_pd, ph_mht);
	    }
	  }
	}
      }
    }
    // Additional Prediction
    // 2 Jets + Boson
    ++AddCount[0];
    if(mht_value > 150) {
      ++AddCount[1];
      RA2Book2(Book2_2J_150, peta, ppt);
    }
    if(n_jets >= 3) {
      ++AddCount[2];
      if(mht_value > 150) {
	++AddCount[3];
	RA2Book2(Book2_3J_150, peta, ppt);
      }
    }
  }

  // Plots
  if(n_jets >= 2) FillHistograms(iEvent, iSetup, 7, photon, Jets, evW);  // 08_2Jets_Egamma
  if(n_jets >= 3) FillHistograms(iEvent, iSetup, 8, photon, Jets, evW); // 09_3Jets_Egamma
  if(n_ra2  == 1) FillHistograms(iEvent, iSetup, 9, photon, Jets, evW);  // 10_RA2Sel_Egamma

}



// ------------ Fill Photon And Jet Information ----------
void 
PhotonAnalyzer::FillHistograms(const edm::Event& iEvent, const edm::EventSetup& iSetup, int h, const pat::Photon * & p, std::vector<const pat::Jet*>& Jets, double evW)
{
  timer.Start();
 

  for(unsigned int i=0; i<Jets.size(); ++i) {
    r = Jets[i];
    // r->();
  }
  // p->();
  // -----------------
  double pt = p->pt(), et = p->et(), eta = p->eta(), phi = p->phi();
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"Pt"), n_et, n1_et, n2_et, pt,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"Et"), n_et, n1_et, n2_et, et,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"Eta"), n_eta, n1_eta, n2_eta, eta,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"Phi"), n_phi, n1_phi, n2_phi, phi,evW);
  // -----------------
  double holTrkIso = p->trkSumPtHollowConeDR04(), solTrkIso = p->trkSumPtSolidConeDR04(), ecalIso = p->ecalRecHitSumEtConeDR04(), hcalIso = p->hcalTowerSumEtConeDR04();
  double p_sIPIP = 0, p_sIEIE = p->sigmaIetaIeta(), p_HE = p->hadronicOverEm();
  bool p_hPS = p->hasPixelSeed(), p_hCT = p->hasConversionTracks();
  int p_nTH = p->nTrkHollowConeDR04();
  double p_Energy = p->energy();
  // Try Photon Conversion ... does not work for GenParticles ...
  double TPSCE, PTSC, dPhi, dCotThe, chi2;
  bool validvertex=0, p_hVCT=0;
  int count_conv=0, count_all_conv=0, ntracks=0;
  if(p_hCT) {
    reco::ConversionRefVector conversions = p->conversions();
    count_conv = 0;
    for(unsigned int j=0; j<conversions.size(); ++j) {
      if (j>0) { fss<<"This Photon has more than one conversion"<<std::endl; print_fss = 1;}
      // do some printout of conversion to fss
      // not yet implemented
      TPSCE = pow(conversions[j]->EoverP(),-1);
      PTSC  = conversions[j]->pairMomentum().mag()/conversions[j]->caloCluster()[0]->energy();
      dPhi=conversions[j]->dPhiTracksAtVtx();
      dCotThe=conversions[j]->pairCotThetaSeparation();
      chi2=conversions[j]->conversionVertex().chi2();
      validvertex = conversions[j]->conversionVertex().isValid();
      ntracks = conversions[j]->nTracks();
      if(validvertex && chi2 > chi2_cut && fabs(dPhi) < dphi_cut && fabs(dCotThe) < dcottheta_cut) p_hVCT = 1;
      ++count_conv;
    }
    count_all_conv += count_conv;
  }
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"HolTrkIso"), n_iso, n1_iso, n2_iso, holTrkIso,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"SolTrkIso"), n_iso, n1_iso, n2_iso, solTrkIso,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"EcalIso"), n_iso, n1_iso, n2_iso, ecalIso,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"HcalIso"), n_iso, n1_iso, n2_iso, hcalIso,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"Conversions"), n_ps, n1_ps, n2_ps, p_hCT,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"V_Conversions"), n_ps, n1_ps, n2_ps, p_hVCT,evW);
  if(p_hCT) {
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"PairPOverSCE"), n_pe, n1_pe, n2_pe, TPSCE, evW);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"PoverE"), n_pe, n1_pe, n2_pe, PTSC, evW);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"nConv"), n_cv, n1_cv, n2_cv, count_conv,evW);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"nTrack"), n_cv, n1_cv, n2_cv, ntracks,evW);
  }
  if(p_hVCT) {
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"V_PairPOverSCE"), n_pe, n1_pe, n2_pe, TPSCE,evW);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"V_PoverE"), n_pe, n1_pe, n2_pe, PTSC,evW);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"V_nConv"), n_cv, n1_cv, n2_cv, count_conv,evW);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h],"V_nTrack"), n_cv, n1_cv, n2_cv, ntracks,evW);
  }
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"PixelSeeds"), n_ps, n1_ps, n2_ps, p_hPS,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"HE"), n_he, n1_he, n2_he, p_HE,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"SigmaIEtaIEta"), n_s, n1_s, n2_s, p_sIEIE,evW);
  // Fill(photoncontainers[h].get(PhotonContainNameArray[h],"SigmaIPhiIPhi"), n_s, n1_s, n2_s, p_sIPIP,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"nTrackHollow"), n_trk, n1_trk, n2_trk, p_nTH,evW);
  // Analysis | Bin1: 100 < EB < 120 | 120 < EB | 100 < EE < 120 | Bin4: 120 < EE |
  if(p->isEB()) {
    if(p->pt()<120) { 
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "anaBookKeeping"), n_bk, n1_bk, n2_bk,1);
    }
    else {
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "anaBookKeeping"), n_bk, n1_bk, n2_bk,2);
    }
  }
  else {
    if(p->pt()<120) { 
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "anaBookKeeping"), n_bk, n1_bk, n2_bk,3);
    }
    else {
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "anaBookKeeping"), n_bk, n1_bk, n2_bk,4);
    }
  }
  // Photon Bump
  if(SpikeCleaning) {
    // Not really related but only available in newest processing
    // and anyway only important related to the additional spike
    double dr_p_met = 0.0, dphi_p_met = 0.0, pfmet = 0.0;
    // patMETs_patMETsPF__PAT
    const pat::MET * met;
    edm::Handle< std::vector<pat::MET> >  patMET;
    iEvent.getByLabel("patMETsPF", patMET);
    // std::cout<<"Size of patMET vector: "<<patMET->size()<<std::endl;
    for(unsigned int i=0; i<patMET->size(); ++i) {
      if (i>0) continue; // just take first element -- is highest 
	met = &((*patMET)[i]);
    }
    pfmet = met->pt();
    dr_p_met = reco::deltaR(met->eta(),met->phi(),p->eta(),p->phi());
    dphi_p_met = reco::deltaPhi(met->phi(),p->phi());
    Fill(photoncontainers[h].get(PhotonContainNameArray[h], "PFMET"), n_et, n1_et, n2_et, pfmet);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h], "DR_Ph_MET"), n_rd, n1_rd, n2_rd, dr_p_met);  // 0 .. 4.5 in steps of 0.01
    Fill(photoncontainers[h].get(PhotonContainNameArray[h], "DP_Ph_MET"), n_pd, n1_pd, n2_pd, dphi_p_met);
    Fill(photoncontainers[h].get(PhotonContainNameArray[h], "E"), n_ht, n1_ht, n2_ht, p_Energy);
    // Ecal RecHits                                                                                                                                                                              
    edm::Handle<EcalRecHitCollection> EBReducedRecHits;
    iEvent.getByLabel("reducedEcalRecHitsEB", EBReducedRecHits);
    edm::Handle<EcalRecHitCollection> EEReducedRecHits;
    iEvent.getByLabel("reducedEcalRecHitsEE", EEReducedRecHits);
    // get the channel status from the DB                                                                                                                                                        
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);
    // ECAL Lazy Tool                                                                                                                                                                            
    EcalClusterLazyTools lazyTool(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));

    // ECAL Lazy Tools: sigma_IEtaIEta, sigma_IEtaIPhi, sigma_IPhiIPhi and eMax
    const reco::CaloClusterPtr  seed = p->superCluster()->seed();
    std::vector<float> vCov = lazyTool.covariances(*seed);
    float covPhiPhi = vCov[2];
    float covEtaPhi = vCov[1];
    float covEtaEta = vCov[0];
    std::vector<float> viCov = lazyTool.localCovariances(*seed);
    float sigmaIphiIphi = sqrt(viCov[2]);
    float sigmaIetaIphi = sqrt(viCov[1]);
    float sigmaIetaIeta = sqrt(viCov[0]);
    float eMax  = lazyTool.eMax(*seed),  e2nd = lazyTool.e2nd(*seed),     e3x3 = lazyTool.e3x3(*seed);
    float eLeft = lazyTool.eLeft(*seed), eRight = lazyTool.eRight(*seed), eTop = lazyTool.eTop(*seed), eBottom = lazyTool.eBottom(*seed);
    float e2overe5 = ((eLeft + eRight + eTop + eBottom)==0)? 0: (eMax + e2nd)/( eMax + eLeft + eRight + eTop + eBottom);
    float e2overe8 = ((e3x3-eMax)==0)? 0: (eMax + e2nd)/(e3x3-eMax);
    float e2overe9 = (e3x3==0)? 0: (eMax + e2nd)/e3x3;

    // ECAL Severity Level: severity and swissCross
    DetId id     = lazyTool.getMaximum(*seed).first;
    float seedTime  = -999., seedOutOfTimeChi2 = -999., seedChi2 = -999., seedSwissCross = -999.;
    int   seedRecoFlag = -1, seedSeverity = -1;
    const EcalRecHitCollection & rechits = ( p->isEB() ? *EBReducedRecHits : *EEReducedRecHits);
    EcalRecHitCollection::const_iterator it = rechits.find( id );
    if( it != rechits.end() ) {
      seedTime = it->time();
      seedOutOfTimeChi2 = it->outOfTimeChi2();
      seedChi2 = it->chi2();
      seedRecoFlag = it->recoFlag();
      seedSeverity = EcalSeverityLevelAlgo::severityLevel( id, rechits, *chStatus );
      seedSwissCross = EcalSeverityLevelAlgo::swissCross( id, rechits);
    }
    if(p->isEB()) {
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_eMax"), n_et, n1_et, n2_et, eMax);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_seedTime"), n_sT, n1_sT, n2_sT, seedTime);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_RecoFlag"), n_RF, n1_RF, n2_RF, seedRecoFlag);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_seedSeverity"), n_sS, n1_sS, n2_sS, seedSeverity);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_SwissCross"), n_SC, n1_SC, n2_SC, seedSwissCross );
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_E2E9"), n_SC, n1_SC, n2_SC, e2overe9);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_SigmaIEtaIEta"), n_s, n1_s, n2_s, sigmaIetaIeta);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Barrel_SigmaIPhiIPhi"), n_s, n1_s, n2_s, sigmaIphiIphi);
    }
    else if(p->isEE()) {
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_eMax"), n_et, n1_et, n2_et, eMax);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_seedTime"), n_sT, n1_sT, n2_sT, seedTime);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_RecoFlag"), n_RF, n1_RF, n2_RF, seedRecoFlag);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_seedSeverity"), n_sS, n1_sS, n2_sS, seedSeverity);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_SwissCross"), n_SC, n1_SC, n2_SC, seedSwissCross);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_E2E9"), n_SC, n1_SC, n2_SC, e2overe9);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_SigmaIEtaIEta"), n_s, n1_s, n2_s, sigmaIetaIeta);
      Fill(photoncontainers[h].get(PhotonContainNameArray[h], "Endcap_SigmaIPhiIPhi"), n_s, n1_s, n2_s, sigmaIphiIphi);
    }
    else {}
    p_sIPIP = sigmaIphiIphi;
  }
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"SigmaIPhiIPhi"), n_s, n1_s, n2_s, p_sIPIP,evW);
  if (h==0) {
    // Jets
    double jet_pt = 0.0, jet_eta = 0.0, jet_phi = 0.0;
    for (unsigned int i=0; i<Jets.size(); ++i) {
      r = Jets[i];
      jet_pt = r->pt(); jet_eta = r->eta(); jet_phi = r->phi();
      Fill(histcontainers[h].get("Jets","JetsPt"), n_et, n1_et, n2_et, jet_pt,evW); 
      Fill(histcontainers[h].get("Jets","JetsPtUn"), n_et, n1_et, n2_et, jet_pt);
      Fill(histcontainers[h].get("Jets","JetsEta"), n_eta, n1_eta, n2_eta, jet_eta,evW);
      Fill(histcontainers[h].get("Jets","JetsPhi"), n_phi, n1_phi, n2_phi, jet_phi,evW);
    }
    Fill(histcontainers[h].get("Jets","JetsMult"), n_mult, n1_mult, n2_mult, Jets.size(),evW); 
    // Event
    Fill(histcontainers[h].get("Event", "nConvEvt"), n_mult, n1_mult, n2_mult, count_all_conv);
  }
  
  timer.Stop();
  tss<<" PhotonAnalyzer :: FillPhotonAndJetHistograms :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
}

void
PhotonAnalyzer::DedicatedPlots(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::GenParticle*& b)
{
  // Require 3 GenJets with pt > 50 and |eta| < 2.5                                                                                                                              
  int goodGenJets_AK5PT50ETA25 = 0, goodGenJets_AK5PT30ETA50 = 0;
  double ht = 0;
  reco::MET::LorentzVector mht(0,0,0,0);
  double min_phi_pi = 0.0, min_r_pi = 0.0;
  // GEN JETS
  edm::Handle< std::vector<reco::GenJet> >      genJets;
  iEvent.getByLabel(GenJets, genJets);
  for(unsigned int i=0; i<genJets->size(); ++i) {
    gj = &((*genJets)[i]);
    if (gj->pt() > 50 && fabs(gj->eta()) < 2.5) {++goodGenJets_AK5PT50ETA25; ht  += gj->pt(); }
    if (gj->pt() > 30) {
      ++goodGenJets_AK5PT30ETA50; mht -= gj->p4();
      double phi_pi = reco::deltaPhi(b->phi(),gj->phi());
      double r_pi = reco::deltaR(b->phi(), b->eta(), gj->phi(), gj->eta());
      if(i==0) {min_phi_pi = phi_pi; min_r_pi = r_pi;}
      else if (min_phi_pi > phi_pi) {min_phi_pi = phi_pi; min_r_pi = r_pi;}
      else {}
    }
  }
  int goodGenJets = goodGenJets_AK5PT50ETA25;
  double b_jcl = fabs(min_phi_pi);
  double b_rjc = fabs(min_r_pi);
  // RA2 Selection Criteria                                                                                                                                                      
  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  double mht_value = MHT.pt();
  double mht_j1 = 0, mht_j2 = 0, mht_j3 = 0;
  if(goodGenJets_AK5PT50ETA25 >= 3) {
    gj1 = &((*genJets)[0]); gj2 = &((*genJets)[1]); gj3 = &((*genJets)[2]);
    mht_j1 = fabs(reco::deltaPhi(gj1->phi(),MHT.phi()));
    mht_j2 = fabs(reco::deltaPhi(gj2->phi(),MHT.phi()));
    mht_j3 = fabs(reco::deltaPhi(gj3->phi(),MHT.phi()));
  }
  if(goodGenJets > 2) {
    if(ht > 300) {
      if ( mht_j1 > 0.5 && mht_j2 > 0.5 && mht_j3 > 0.3) {
	if(mht_value > 150) {
	  Fill(Pt_All_BaS, n_et, n1_et, n2_et, b->pt());
	  if(b_rjc > 0.5) Fill(Pt_R05_BaS, n_et, n1_et, n2_et, b->pt());
	  if(b_rjc > 0.9) Fill(Pt_R09_BaS, n_et, n1_et, n2_et, b->pt());
	  if(ht > 500) {
	    Fill(Pt_All_HHS, n_et, n1_et, n2_et, b->pt());
	    if(b_rjc > 0.5) Fill(Pt_R05_HHS, n_et, n1_et, n2_et, b->pt());
	    if(b_rjc > 0.9) Fill(Pt_R09_HHS, n_et, n1_et, n2_et, b->pt());
	  }
	  if(mht_value > 250) {
	    Fill(Pt_R05_HMS, n_et, n1_et, n2_et, b->pt());
	    if(b_rjc > 0.5) Fill(Pt_All_HMS, n_et, n1_et, n2_et, b->pt());
	    if(b_rjc > 0.9) Fill(Pt_R09_HMS, n_et, n1_et, n2_et, b->pt());
	  }
	}
      }
    }
  }
}

void
PhotonAnalyzer::RA2Book(TH1F * Histo, bool bar, bool end, double ppt, int Identity)
{
  if(bar && ppt <  120) Fill(Histo, n_bk, n1_bk, n2_bk, 1);
  if(bar && ppt >= 120) Fill(Histo, n_bk, n1_bk, n2_bk, 2);
  if(end && ppt <  120) Fill(Histo, n_bk, n1_bk, n2_bk, 3);
  if(end && ppt >= 120) Fill(Histo, n_bk, n1_bk, n2_bk, 4);
  // std::cout<<"RA2Book for Histo = "<<Histo->GetName()<<" Identity = "<<Identity<<std::endl;
  if(Identity <6)     { 
    Fill(Histo, n_bk, n1_bk, n2_bk,Identity+5); 
    // std::cout<<"RA2Book for Histo = "<<Histo->GetName()<<" Fill("<<Histo<<", "<<n_bk<<", "<<n1_bk<<", "<<n2_bk<<", Identity+5 = "<<Identity+5<<std::endl; 
  }
  else {
    Fill(Histo, n_bk, n1_bk, n2_bk,11);           // ISR + FSR                                                                                                                                               
    Fill(Histo, n_bk, n1_bk, n2_bk,Identity+6);   // ISR or FSR
    // std::cout<<"RA2Book for Histo = "<<Histo->GetName()<<" Fill("<<Histo<<", "<<n_bk<<", "<<n1_bk<<", "<<n2_bk<<", Identity+6 = "<<Identity+6<<std::endl;
  }
}

void
PhotonAnalyzer::RA2Book1(TH1F * Histo, bool bar, bool end, double ppt)
{
  if(bar && ppt <  120) Fill(Histo, n_b1, n1_b1, n2_b1, 1);
  if(bar && ppt >= 120) Fill(Histo, n_b1, n1_b1, n2_b1, 2);
  if(end && ppt <  120) Fill(Histo, n_b1, n1_b1, n2_b1, 3);
  if(end && ppt >= 120) Fill(Histo, n_b1, n1_b1, n2_b1, 4);
}

void
PhotonAnalyzer::RA2Book2(TH1F * Histo, double peta, double ppt)
{
  if(peta > 0 && peta < 0.9) {
    if      (ppt > 100 && ppt < 120) Fill(Histo, n_b2, n1_b2, n2_b2, 1);
    else if (ppt > 120 && ppt < 200) Fill(Histo, n_b2, n1_b2, n2_b2, 2);
    else if (ppt > 200 && ppt < 300) Fill(Histo, n_b2, n1_b2, n2_b2, 3);
    else if (ppt > 300)              Fill(Histo, n_b2, n1_b2, n2_b2, 4);
    else {}
  }
  else if (peta > 0.9 && peta < 1.4442) {
    if      (ppt > 100 && ppt < 120) Fill(Histo, n_b2, n1_b2, n2_b2, 5);
    else if (ppt > 120 && ppt < 200) Fill(Histo, n_b2, n1_b2, n2_b2, 6);
    else if (ppt > 200 && ppt < 300) Fill(Histo, n_b2, n1_b2, n2_b2, 7);
    else if (ppt > 300)              Fill(Histo, n_b2, n1_b2, n2_b2, 8);
    else {}
  }
  else if (peta > 1.566 && peta < 2.1) {
    if      (ppt > 100 && ppt < 120) Fill(Histo, n_b2, n1_b2, n2_b2, 9);
    else if (ppt > 120 && ppt < 200) Fill(Histo, n_b2, n1_b2, n2_b2, 10);
    else if (ppt > 200 && ppt < 300) Fill(Histo, n_b2, n1_b2, n2_b2, 11);
    else if (ppt > 300)              Fill(Histo, n_b2, n1_b2, n2_b2, 12);
    else {}
  }
  else if (peta > 2.1 && peta < 2.5) {
    if      (ppt > 100 && ppt < 120) Fill(Histo, n_b2, n1_b2, n2_b2, 13);
    else if (ppt > 120 && ppt < 200) Fill(Histo, n_b2, n1_b2, n2_b2, 14);
    else if (ppt > 200 && ppt < 300) Fill(Histo, n_b2, n1_b2, n2_b2, 15);
    else if (ppt > 300)              Fill(Histo, n_b2, n1_b2, n2_b2, 16);
    else {}
  }
  else {
    // std::cout<<"Not filled in RA2Book2 for Histo = "<<Histo->GetName()<<" photon eta = "<<peta<<" pt = "<<ppt<<std::endl;
  }
}




void 
PhotonAnalyzer::endLuminosityBlock(const edm::LuminosityBlock& iLuminosityBlock, const edm::EventSetup& iSetup) 
{ 

  //  EVENT SELECTION BY PAT
  edm::Handle<edm::MergeableCounter> allEventsCounter, postHLTCounter, prefilterCounter, postStdCleaningCounter, postPFCleaningCounter, postPhotonCounter, postPFJetsCounter;
  iLuminosityBlock.getByLabel("allEventsCounter", allEventsCounter);
  iLuminosityBlock.getByLabel("postHLTCounter",postHLTCounter);
  iLuminosityBlock.getByLabel("prefilterCounter",prefilterCounter);
  iLuminosityBlock.getByLabel("postStdCleaningCounter",postStdCleaningCounter);
  iLuminosityBlock.getByLabel("postPFCleaningCounter",postPFCleaningCounter);
  iLuminosityBlock.getByLabel("postPhotonCounter",postPhotonCounter);
  iLuminosityBlock.getByLabel("postPFJetsCounter",postPFJetsCounter);
                              

  if(allEventsCounter.isValid()) {
    Count[0] += allEventsCounter->value; RA2Count[0] += allEventsCounter->value;
    lss<<"In endLuminosityBlock; allEventsCounter :: adding "<<std::setw(10)<<allEventsCounter->value<<" to total events"<<std::endl; 
  }
  if(postHLTCounter.isValid()) {
    Count[1] += postHLTCounter->value; 
    lss<<"In endLuminosityBlock; postHLTCounter :: adding "<<std::setw(10)<<postHLTCounter->value<<" to total events"<<std::endl; 
  }
  if(prefilterCounter.isValid())       Count[2] += prefilterCounter->value;
  if(postStdCleaningCounter.isValid()) Count[3] += postStdCleaningCounter->value;
  if(postPFCleaningCounter.isValid())  Count[4] += postPFCleaningCounter->value;
  if(postPhotonCounter.isValid()) {
    Count[5] += postPhotonCounter->value;
    lss<<"In endLuminosityBlock; postPhotonCounter :: adding "<<std::setw(10)<<postPhotonCounter->value<<" to total events"<<std::endl; 
  }
  if(postPFJetsCounter.isValid()) {
    Count[6] += postPFJetsCounter->value; RA2Count[1] += postPFJetsCounter->value;
    lss<<"In endLuminosityBlock; postPFJetsCounter :: adding "<<std::setw(10)<<postPFJetsCounter->value<<" to total events"<<std::endl; 
  }

}


void
PhotonAnalyzer::Fill(TH1 * Histo, int n, double n1, double n2, double value) {
  double binwidth = (n2-n1)/n;
  double n1_ = n1+binwidth/2;
  double n2_ = n2-binwidth/2;
  if(value > n1 && value < n2) Histo->Fill(value);
  else if (value < n1) Histo->Fill(n1_);
  else if (value > n2) Histo->Fill(n2_);
  else {}
}
void
PhotonAnalyzer::Fill(TH1 * Histo, int n, double n1, double n2, double value, double evWeight) {
  double binwidth = (n2-n1)/n;
  double n1_ = n1+binwidth/2;
  double n2_ = n2-binwidth/2;
  if(value > n1 && value < n2) Histo->Fill(value, evWeight);
  else if (value < n1) Histo->Fill(n1_, evWeight);
  else if (value > n2) Histo->Fill(n2_, evWeight);
  else {}
}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonAnalyzer::endJob() 
{ 

  if(Debug[0] || print_tss) std::cout<<" "<<std::endl<<"Timing Debug Stream:      "<<std::endl<<"--------------------------"<<std::endl<<tss.str()<<std::endl; else tss.str("");

  std::cout<<"In endJob; allEventsCounter         :: Total amount of Events : "<<std::setw(10)<<Count[0]<<std::endl; 
  std::cout<<"In endJob; postHLTCounter           :: Total amount of Events : "<<std::setw(10)<<Count[1]<<std::endl; 
  std::cout<<"In endJob; postStdCleaningCounter   :: Total amount of Events : "<<std::setw(10)<<Count[3]<<std::endl; 
  std::cout<<"In endJob; postPFCleaningCounter    :: Total amount of Events : "<<std::setw(10)<<Count[4]<<std::endl; 
  std::cout<<"In endJob; postPhotonCounter        :: Total amount of Events : "<<std::setw(10)<<Count[5]<<std::endl; 
  std::cout<<"In endJob; postPFJetCounter         :: Total amount of Events : "<<std::setw(10)<<Count[6]<<std::endl;   
  std::cout<<"In endJob; Offline Cuts             :: Total amount of Events : "<<std::setw(10)<<Count[7]<<std::endl;   
  std::cout<<"In endJob; Photon + 2 Jets          :: Total amount of Events : "<<std::setw(10)<<Count[8]<<std::endl;   
  std::cout<<"In endJob; Photon + 3 Jets          :: Total amount of Events : "<<std::setw(10)<<Count[9]<<std::endl;   

  std::cout<<"In endJob; Barrel + Endcap          :: Total amount of Events : "<<std::setw(10)<<Count[7]<<std::endl;   
  std::cout<<"In endJob; Barrel                   :: Total amount of Events : "<<std::setw(10)<<Count[11]<<std::endl;   
  std::cout<<"In endJob; Endcap                   :: Total amount of Events : "<<std::setw(10)<<Count[12]<<std::endl;   
  std::cout<<""<<std::endl;
  std::cout<<"Amount of Rejected Events: "<<IdIsoLepton<<std::endl;

  for(unsigned int i=0; i<Count.size(); ++i)    { EventCounters->SetBinContent(i+1, Count[i]); }
  for(unsigned int i=0; i<RA2Count.size(); ++i) { RA2SelectionHisto->SetBinContent(i+1, RA2Count[i]); }
  for(unsigned int i=0; i<AddCount.size(); ++i) { AdditionalPrediction->SetBinContent(i+1, AddCount[i]); }
  
  std::string ECBinLabels [] = {"Processed Events","post HLT","preFilter","postStdClean","postPFClean","postPhoton","postPFJets","ISO #gamma > 100 GeV","2J != #gamma","3J != #gamma","RA2 Selected"};
  std::string RSBinLabels [] = {"Processed Events","After PAT Selection", "ISO #gamma > 100 + 3J", "with HT > 300 GeV/c", "with Angular Cuts", "with MHT > 150 GeV/c", "Boson PT cut", "High HT", "High MHT"};
  std::string ACBinLabels [] = {"ISO #gamma + 2 Jets", "2 Jets + MHT > 150", "ISO #gamma + 3 Jets", "3 Jets + MHT > 150"};
  std::string EFBinLabels [] = {"No Matched GEN #gamma", "no GEN #gamma in event", "No Jet == #gamma", "No Jets in event"};
  
  for(unsigned int i=0; i<(Count.size()-2); ++i) { EventCounters->GetXaxis()->SetBinLabel(i+1, ECBinLabels[i].c_str()); }
  for(unsigned int i=0; i<RA2Count.size(); ++i)  { RA2SelectionHisto->GetXaxis()->SetBinLabel(i+1, RSBinLabels[i].c_str()); }
  for(unsigned int i=0; i<AddCount.size(); ++i) { AdditionalPrediction->GetXaxis()->SetBinLabel(i+1, ACBinLabels[i].c_str());}
  for(int i=0; i<4; ++i)  { ErrorFlags->GetXaxis()->SetBinLabel(i+1, EFBinLabels[i].c_str()); }
  

}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonAnalyzer);
