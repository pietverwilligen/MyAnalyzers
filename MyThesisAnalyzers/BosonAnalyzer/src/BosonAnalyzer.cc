// -*- C++ -*-
//
// Package:    BosonAnalyzer
// Class:      BosonAnalyzer
// 
/**\class BosonAnalyzer BosonAnalyzer.cc MyAnalyzers/BosonAnalyzer/src/BosonAnalyzer.cc

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
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include <DataFormats/EgammaReco/interface/ElectronSeed.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/Photon.h>
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/METReco/interface/MET.h>
#include "MyAnalyzers/BosonAnalyzer/interface/SHistContainer.h"
#include "MyAnalyzers/BosonAnalyzer/interface/SGraphContainer.h"


#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

//
// class declaration
//

class BosonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit BosonAnalyzer(const edm::ParameterSet&);
      ~BosonAnalyzer();


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);

  // could be optimized by using pointers to vectors instead of passing the vector as a whole

  void Selection     (const edm::Event&, const edm::EventSetup&, bool&, const pat::Photon*&, std::vector<const pat::Jet*>&, double, int);
  void FillHistograms(int, const pat::Photon*&, std::vector<const pat::Jet*>&, double);
  void Fill          (TH1*, int, double, double, double);
  void Fill          (TH1*, int, double, double, double, double);

  void ZBosonSelection(const edm::Event&, const edm::EventSetup&, const reco::GenParticle*&, std::vector<const pat::Jet*>&);
  void WBosonSelection(const edm::Event&, const edm::EventSetup&, const reco::GenParticle*&, const pat::Muon*&, std::vector<const pat::Jet*>&);
  // Unified Approach
  void RA2Selection  (const pat::Photon*&, const reco::GenParticle*&, const pat::Muon*&, std::vector<const pat::Jet*>&);


  double angle(double eta1, double eta2, double phi1, double phi2);
  double angle(const reco::LeafCandidate * a, const reco::LeafCandidate * b);
  double angle(const reco::LeafCandidate * a, const pat::Jet * b);

  // ----------member data ---------------------------

  // WRITE
  TFile * outputfile;
  TDirectoryFile *RECO, *RA2SEL; 
  TDirectoryFile *RECO_Photon, *RECO_Event, *RECO_Jets;
  int n_PhotonDir;
  TDirectoryFile *RECO_Egamma, *RECO_Egamma_1PV, *RECO_Prompt_Egamma, *RECO_NonPrompt_Egamma, *RECO_Direct_Egamma, *RECO_Fragment_Egamma, *RECO_Barrel_Egamma, *RECO_Endcap_Egamma, *RECO_2Jets_Egamma, *RECO_3Jets_Egamma;
  std::vector<SHistContainer> histcontainers, photoncontainers;
  TH1F * EventCounters, * RA2SelectionHisto, * AdditionalPrediction, *ErrorFlags;
  TH1F *JetMult_All, *PT_All, *HT_All, *MHT_All, *MHT_J1_All, *MHT_J2_All, *MHT_J3_All, *Book_All, *JetMult_PTC, *PT_PTC, *HT_PTC, *MHT_PTC, *MHT_J1_PTC, *MHT_J2_PTC, *MHT_J3_PTC, *Book_PTC;
  TH1F *JetMult_AJC, *PT_AJC, *HT_AJC, *MHT_AJC, *MHT_J1_AJC, *MHT_J2_AJC, *MHT_J3_AJC, *Book_AJC, *JetMult_AHC, *PT_AHC, *HT_AHC, *MHT_AHC, *MHT_J1_AHC, *MHT_J2_AHC, *MHT_J3_AHC, *Book_AHC;
  TH1F *JetMult_AAC, *PT_AAC, *HT_AAC, *MHT_AAC, *MHT_J1_AAC, *MHT_J2_AAC, *MHT_J3_AAC, *Book_AAC, *JetMult_AMC, *PT_AMC, *HT_AMC, *MHT_AMC, *MHT_J1_AMC, *MHT_J2_AMC, *MHT_J3_AMC, *Book_AMC;
  TH1F *JetMult_HHS, *PT_HHS, *HT_HHS, *MHT_HHS, *MHT_J1_HHS, *MHT_J2_HHS, *MHT_J3_HHS, *Book_HHS, *JetMult_HMS, *PT_HMS, *HT_HMS, *MHT_HMS, *MHT_J1_HMS, *MHT_J2_HMS, *MHT_J3_HMS, *Book_HMS;
  TH1F *PH_JCl_All, *PH_JCl_PTC, *PH_JCl_AJC, *PH_JCl_AHC, *PH_JCl_AAC, *PH_JCl_AMC, *PH_JCl_HHS, *PH_JCl_HMS;
  TH1F *PH_RJC_All, *PH_RJC_PTC, *PH_RJC_AJC, *PH_RJC_AHC, *PH_RJC_AAC, *PH_RJC_AMC, *PH_RJC_HHS, *PH_RJC_HMS;
  TH1F *PH_MHT_All, *PH_MHT_PTC, *PH_MHT_AJC, *PH_MHT_AHC, *PH_MHT_AAC, *PH_MHT_AMC, *PH_MHT_HHS, *PH_MHT_HMS;
  TH1F *Wstat;
  // READ
  std::vector<int> Debug;
  int evNum, rnNum;
  int Boson, AmountOfJets;
  bool Data;
  double PtCut;
  std::string RootFileName, GenParticles, GenJets, RecoPhotons, RecoJets;
  // INTERNAL USE
  const reco::GenParticle     *g;
  const pat::Photon           *p, *p1, *p2;
  const pat::Jet              *r, *r1, *r2, *r3;
  // const reco::Jet              *r, *r1, *r2, *r3;
  //
  TStopwatch timer, globaltimer, isolatetimer;
  std::vector<int> Count, RA2Count, AddCount;
  bool photonFound, promptPhoton, directPhoton, fragPhoton, decayPhoton;
  // DEBUG STREAM
  bool print_tss, print_dss, print_rss, print_lss, print_fss;
  std::stringstream tss, dss, rss, lss, fss;  // stringstreams: tss = Timing Debug Stream | dss = Matching Debug Stream | rss = RA2Selection Debug Stream | lss = Luminosity Debug Stream | fss = Fill Debug Stream
};

//
// constants, enums and typedefs
//
double pi = 3.1415926535;

// Photon ID cuts:
double photon_ptcut = 100.0;
double photon_etacut = 2.50;
double egm_trk  = 2.0; double egm_trk_prop  = 0.001;
double egm_ecal = 4.2; double egm_ecal_prop = 0.003;
double egm_hcal = 2.2; double egm_hcal_prop = 0.001;
double photon_he = 0.05;
double sigma_barrel = 0.01;
double sigma_endcap = 0.03;
double sigma_spike = 0.001; // sigma_IEtaIEta and sigma_IPhiIPhi spike cut: p_sIEIE > 0.001 p_sIPIP > 0.001 (manual)
double sc_spike = 0.95;     // Swiss Cross and E2 over E9 spike cut:        E1/E4 < 0.95 (38X rereco: ok) E2/E9 < 0.95 (manual)
double st_spike = -3.5;     // Seed Time spike cut (ns):                    seedTime > -3.5 ns (manual)
// ECAL Geom def
double EBEE_bord = 1.479;
double EB_accept = 1.4442;  // use method: isEB()  --> implemented for supercluster Eta
double EE_accept = 1.566;   // use method: isEE()

// Conversion cuts:                                                                                                                                                                                              
double dphi_cut = 0.2;
double dcottheta_cut = 0.3;
double chi2_cut = 0.0005;

// TH1F Settings
int n_eta  = 30;  double n1_eta  = -3.0,  n2_eta  = 3.0;
int n_phi  = 36;  double n1_phi  = -3.14, n2_phi  = 3.14;
int n_et   = 200; double n1_et   = 0.0,   n2_et   = 1000;
int n_ht   = 500; double n1_ht   = 0.0,   n2_ht   = 2500;
int n_iso  = 32;  double n1_iso  = -1.0,  n2_iso  = 15;
int n_mult = 11;  double n1_mult = -0.5,  n2_mult = 10.5;
int n_pe   = 15;  double n1_pe   = 0.0,   n2_pe   = 3.0;
int n_he   = 20;  double n1_he   = 0.00,  n2_he   = 0.10;
int n_sb   = 30;  double n1_sb   = 0.00,  n2_sb   = 0.03;
int n_se   = 30;  double n1_se   = 0.015, n2_se   = 0.030;
int n_s2e  = 45;  double n1_s2e  = 0.015, n2_s2e  = 0.060;
int n_ps   = 2;   double n1_ps   = -0.5,  n2_ps   = 1.5;
int n_cv   = 3;   double n1_cv   = -0.5,  n2_cv   = 3.5;
int n_trk  = 31;  double n1_trk  = 0.00,  n2_trk  = 30.0;

int n_pd   = 31;  double n1_pd   = 0.00,  n2_pd   = 3.14;
int n_rdd  = 50;  double n1_rdd  = 0.00,  n2_rdd  = 0.50;
int n_rd   = 50;  double n1_rd   = 0.00,  n2_rd   = 5.00;

int n_dist = 150; double n1_dist = 0.00,  n2_dist = 1.50;

int n_t1   = 100; double n1_t1   = 0.00,  n2_t1   = 0.10;
int n_t2   = 10;  double n1_t2   = 0.00,  n2_t2   = 0.10;

int n_sT   = 100; double n1_sT   = -25,   n2_sT   = 25;
int n_RF   = 17;  double n1_RF   = -0.5,  n2_RF   = 15.5;
int n_sS   = 6;   double n1_sS   = -0.5,  n2_sS   = 5.5;
int n_SC   = 60;  double n1_SC   = 0.0,   n2_SC   = 1.2;

int n_bk   = 13;   double n1_bk   = 0.5,   n2_bk   = 13.5;

const char * PhotonContainNameArray [] = { "RECO_Egamma", "RECO_Egamma_1PV", "RECO_Prompt_Egamma", "RECO_NonPrompt_Egamma", 
					   "RECO_Direct_Egamma", "RECO_Fragment_Egamma", "RECO_Barrel_Egamma","RECO_Endcap_Egamma", "RECO_2Jets_Egamma", "RECO_3Jets_Egamma"};


//
// static data member definitions
//

//
// constructors and destructor
//
BosonAnalyzer::BosonAnalyzer(const edm::ParameterSet& iConfig)
{
  timer.Start();

  // read parameters from config file
  Debug         = iConfig.getParameter< std::vector<int> >("Debug");
  Boson         = iConfig.getParameter<int>("Boson");
  Data          = iConfig.getParameter<bool>("Data");
  PtCut         = iConfig.getParameter<double>("PtCut");
  AmountOfJets  = iConfig.getParameter<int>("AmountOfJets");
  RootFileName  = iConfig.getParameter<std::string>("RootFileName");
  GenParticles  = iConfig.getParameter<std::string>("GenParticles");
  GenJets       = iConfig.getParameter<std::string>("GenJets");
  RecoPhotons   = iConfig.getParameter<std::string>("RecoPhotons");
  RecoJets      = iConfig.getParameter<std::string>("RecoJets");

   //now do what ever initialization is needed
  outputfile = new TFile(RootFileName.c_str(), "RECREATE" );
  RECO                  = (TDirectoryFile*) outputfile->mkdir("RECO","RECO");
  RA2SEL                = (TDirectoryFile*) outputfile->mkdir("RA2SEL","RA2SEL");

  SHistContainer RECO_hist; histcontainers.push_back(RECO_hist);
  RECO_Photon           = (TDirectoryFile*) RECO->mkdir("RECO_Photon", "RECO_Photon");
  n_PhotonDir           = 10;    // Keep up to date !!!
  RECO_Egamma           = (TDirectoryFile*) RECO_Photon->mkdir("01_Egamma",           "01_Egamma");
  RECO_Egamma_1PV       = (TDirectoryFile*) RECO_Photon->mkdir("02_Egamma_1PV",       "02_Egamma_1PV");
  RECO_Prompt_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("03_Prompt_Egamma",    "03_Prompt_Egamma");
  RECO_NonPrompt_Egamma = (TDirectoryFile*) RECO_Photon->mkdir("04_NonPrompt_Egamma", "04_NonPrompt_Egamma"); 
  RECO_Direct_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("05_Direct_Egamma",    "05_Direct_Egamma");
  RECO_Fragment_Egamma  = (TDirectoryFile*) RECO_Photon->mkdir("06_Fragment_Egamma",  "06_Fragment_Egamma");
  RECO_Barrel_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("07_Barrel_Egamma",    "07_Barrel_Egamma");
  RECO_Endcap_Egamma    = (TDirectoryFile*) RECO_Photon->mkdir("08_Endcap_Egamma",    "08_Endcap_Egamma");    
  RECO_2Jets_Egamma     = (TDirectoryFile*) RECO_Photon->mkdir("09_2Jets_Egamma",     "09_2Jets_Egamma");
  RECO_3Jets_Egamma     = (TDirectoryFile*) RECO_Photon->mkdir("10_RA2Sel_Egamma",    "10_RA2Sel_Egamma");
  SHistContainer RECO_Egamma_c, RECO_Egamma_1PV_c, RECO_Prompt_Egamma_c, RECO_NonPrompt_Egamma_c, 
                 RECO_Direct_Egamma_c, RECO_Fragment_Egamma_c, RECO_Barrel_Egamma_c, RECO_Endcap_Egamma_c, RECO_2Jets_Egamma_c, RECO_3Jets_Egamma_c;
  SHistContainer PhotonContainerArray [] = { RECO_Egamma_c, RECO_Egamma_1PV_c, RECO_Prompt_Egamma_c, RECO_NonPrompt_Egamma_c, 
					     RECO_Direct_Egamma_c, RECO_Fragment_Egamma_c, RECO_Barrel_Egamma_c,RECO_Endcap_Egamma_c, RECO_2Jets_Egamma_c, RECO_3Jets_Egamma_c };
  // const char * PhotonContainNameArray [] = { "RECO_Egamma", "RECO_Egamma_1PV", "RECO_Prompt_Egamma", "RECO_NonPrompt_Egamma", "RECO_Barrel_Egamma","RECO_Endcap_Egamma", "RECO_2Jets_Egamma", "RECO_3Jets_Egamma"};
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
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("SigmaIEtaIEta", "#sigma_{i #eta i #eta} ShowerShape Variable", n_se,  n1_se, n2_se));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("SigmaIPhiIPhi", "#sigma_{i #phi i #phi} ShowerShape Variable", n_s2e, n1_s2e, n2_s2e));
    photoncontainers[h].book(PhotonContainNameArray[h], new TH1F("nTrackHollow",  "number of Tracks in (DR = 0.4) Cone", n_trk, n1_trk, n2_trk));
  }

  for(unsigned int h = 0; h<histcontainers.size(); ++h) {
    histcontainers[h].book("Jets",   new TH1F("JetsMult", "Jets Multiplicity Distribution (with pt > 50)", n_mult, n1_mult, n2_mult));
    histcontainers[h].book("Jets",   new TH1F("JetsPtUn", "Unweighted Jets Pt Distribution", n_et, n1_et, n2_et));
    histcontainers[h].book("Jets",   new TH1F("JetsPt",   "Jets Pt Distribution",  n_et, n1_et, n2_et));
    histcontainers[h].book("Jets",   new TH1F("JetsEta",  "Jets Eta Distribution", n_eta, n1_eta, n2_eta));
    histcontainers[h].book("Jets",   new TH1F("JetsPhi",  "Jets Phi Distribution", n_phi, n1_phi, n2_phi));
    histcontainers[h].book("Jets",   new TH1F("HT",       "HT Distribution",  n_ht, n1_ht, n2_ht));
    histcontainers[h].book("Jets",   new TH1F("MHT",      "MHT Distribution",  n_et, n1_et, n2_et));
    histcontainers[h].book("Jets",   new TH1F("MHT_J1",   "#Delta Phi (MHT - 1st Jet) Distribution", n_pd, n1_pd, n2_pd));
    histcontainers[h].book("Jets",   new TH1F("MHT_J2",   "#Delta Phi (MHT - 2nd Jet) Distribution", n_pd, n1_pd, n2_pd));
    histcontainers[h].book("Jets",   new TH1F("MHT_J3",   "#Delta Phi (MHT - 3rd Jet) Distribution", n_pd, n1_pd, n2_pd));

    // event information
    histcontainers[h].book("Event",  new TH1F("nConvEvt",  "Number of Conversions in Event", n_mult, n1_mult, n2_mult));
    histcontainers[h].book("Event",  new TH1F("nVertEvt",  "Number of Vertices in Event",    n_mult, n1_mult, n2_mult));
    histcontainers[h].book("Event",  new TH1F("PATMatch",  "Photon is PATMatched?",          n_ps, n1_ps, n2_ps));
    histcontainers[h].book("Event",  new TH1F("Distance",  "Distance between Hard Photon and Photon in #eta, #varphi:", n_dist, n1_dist, n2_dist));
    histcontainers[h].book("Event",  new TH1F("DR_Before", "Distance between #gamma and Jet before removing #gamma out the jetcollection by #eta, #varphi: #Delta R (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));
    histcontainers[h].book("Event",  new TH1F("DP_Before", "Distance between #gamma and Jet before removing #gamma out the jetcollection by p_{T}: #Delta p_{T} / p_{T} (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));
    histcontainers[h].book("Event",  new TH1F("DR_After",  "Distance between #gamma and Jet in #eta, #varphi: #Delta R (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));
    histcontainers[h].book("Event",  new TH1F("DP_After",  "Distance between #gamma and Jet in p_{T}: #Delta p_{T} / p_{T} (#gamma, Closesth Jet)", n_dist, n1_dist, n2_dist));


    histcontainers[h].book("Timing", new TH1F("tSel_Real", "Real Time SelectCollections", n_t1, n1_t1, n2_t1));
    histcontainers[h].book("Timing", new TH1F("tIso_Real", "Real Time Isolate",           n_t1, n1_t1, n2_t1));
    histcontainers[h].book("Timing", new TH1F("tAna_Real", "Real Time Analyze",           n_t1, n1_t1, n2_t1));
    histcontainers[h].book("Timing", new TH1F("tSel_CPU",  "CPU Time SelectCollections",  n_t2, n1_t2, n2_t2));
    histcontainers[h].book("Timing", new TH1F("tIso_CPU",  "CPU Time Isolate",            n_t2, n1_t2, n2_t2));
    histcontainers[h].book("Timing", new TH1F("tAna_CPU",  "CPU Time Analyze",            n_t2, n1_t2, n2_t2));
  } 

  // RA2 Selection Plots
  JetMult_All= new TH1F("0_JetMult_All", "PT50ETA25 Jet Multiplicity", n_mult, n1_mult, n2_mult);
  PT_All     = new TH1F("0_PT_All" ,     "Boson PT Distribution",  n_et, n1_et, n2_et);
  HT_All     = new TH1F("0_HT_All" ,     "HT Distribution",  n_ht, n1_ht, n2_ht);
  MHT_All    = new TH1F("0_MHT_All",     "MHT Distribution", n_et, n1_et, n2_et);
  MHT_J1_All = new TH1F("0_MHT_J1_All" , "#Delta Phi (MHT - 1st Jet) Distribution",  n_pd, n1_pd, n2_pd);
  MHT_J2_All = new TH1F("0_MHT_J2_All" , "#Delta Phi (MHT - 2nd Jet) Distribution",  n_pd, n1_pd, n2_pd);
  MHT_J3_All = new TH1F("0_MHT_J3_All" , "#Delta Phi (MHT - 3rd Jet) Distribution",  n_pd, n1_pd, n2_pd);
  Book_All   = new TH1F("0_Book_All"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_All = new TH1F("0_PH_JCl_All",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_All = new TH1F("0_PH_RJC_All",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_All = new TH1F("0_PH_MHT_All" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);

  JetMult_AJC= new TH1F("1_JetMult_AJC", "PT50ETA25 Jet Multiplicity After Jet Cuts", n_mult, n1_mult, n2_mult);
  PT_AJC     = new TH1F("1_PT_AJC" ,     "Boson PT Distribution After Jet Cuts",  n_et, n1_et, n2_et);
  HT_AJC     = new TH1F("1_HT_AJC" ,     "HT Distribution After Jet Cuts", n_ht, n1_ht, n2_ht);
  MHT_AJC    = new TH1F("1_MHT_AJC",     "MHT Distribution After Jet Cuts", n_et, n1_et, n2_et);
  MHT_J1_AJC = new TH1F("1_MHT_J1_AJC" , "#Delta Phi (MHT - 1st Jet) Distribution After Jet Cuts", n_pd, n1_pd, n2_pd);
  MHT_J2_AJC = new TH1F("1_MHT_J2_AJC" , "#Delta Phi (MHT - 2nd Jet) Distribution After Jet Cuts", n_pd, n1_pd, n2_pd);
  MHT_J3_AJC = new TH1F("1_MHT_J3_AJC" , "#Delta Phi (MHT - 3rd Jet) Distribution After Jet Cuts", n_pd, n1_pd, n2_pd);
  Book_AJC   = new TH1F("1_Book_AJC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_AJC = new TH1F("1_PH_JCl_AJC",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_AJC = new TH1F("1_PH_RJC_AJC",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_AJC = new TH1F("1_PH_MHT_AJC" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);
  
  JetMult_AHC= new TH1F("2_JetMult_AHC", "PT50ETA25 Jet Multiplicity After HT Cut", n_mult, n1_mult, n2_mult);
  PT_AHC     = new TH1F("2_PT_AHC" ,     "Boson PT Distribution After HT Cut",  n_et, n1_et, n2_et);
  HT_AHC     = new TH1F("2_HT_AHC" ,     "HT Distribution After HT Cut", n_ht, n1_ht, n2_ht);
  MHT_AHC    = new TH1F("2_MHT_AHC",     "MHT Distribution After HT Cut", n_et, n1_et, n2_et);
  MHT_J1_AHC = new TH1F("2_MHT_J1_AHC" , "#Delta Phi (MHT - 1st Jet) Distribution After HT Cut", n_pd, n1_pd, n2_pd);
  MHT_J2_AHC = new TH1F("2_MHT_J2_AHC" , "#Delta Phi (MHT - 2nd Jet) Distribution After HT Cut", n_pd, n1_pd, n2_pd);
  MHT_J3_AHC = new TH1F("2_MHT_J3_AHC" , "#Delta Phi (MHT - 3rd Jet) Distribution After HT Cut", n_pd, n1_pd, n2_pd);
  Book_AHC   = new TH1F("2_Book_AHC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_AHC = new TH1F("2_PH_JCl_AHC",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_AHC = new TH1F("2_PH_RJC_AHC",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_AHC = new TH1F("2_PH_MHT_AHC" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);

  JetMult_AAC= new TH1F("3_JetMult_AAC", "PT50ETA25 Jet Multiplicity After Angular Cuts", n_mult, n1_mult, n2_mult);
  PT_AAC     = new TH1F("3_PT_AAC" ,     "Boson PT Distribution After Angular Cuts",  n_et, n1_et, n2_et);
  HT_AAC     = new TH1F("3_HT_AAC" ,     "HT Distribution After Angular Cuts", n_ht, n1_ht, n2_ht);
  MHT_AAC    = new TH1F("3_MHT_AAC",     "MHT Distribution After Angular Cuts", n_et, n1_et, n2_et);
  MHT_J1_AAC = new TH1F("3_MHT_J1_AAC" , "#Delta Phi (MHT - 1st Jet) Distribution After Angular Cuts", n_pd, n1_pd, n2_pd);
  MHT_J2_AAC = new TH1F("3_MHT_J2_AAC" , "#Delta Phi (MHT - 2nd Jet) Distribution After Angular Cuts", n_pd, n1_pd, n2_pd);
  MHT_J3_AAC = new TH1F("3_MHT_J3_AAC" , "#Delta Phi (MHT - 3rd Jet) Distribution After Angular Cuts", n_pd, n1_pd, n2_pd);
  Book_AAC   = new TH1F("3_Book_AAC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_AAC = new TH1F("3_PH_JCl_AAC",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_AAC = new TH1F("3_PH_RJC_AAC",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_AAC = new TH1F("3_PH_MHT_AAC" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);
  
  JetMult_AMC= new TH1F("4_JetMult_AMC", "PT50ETA25 Jet Multiplicity After MHT Cut: Baseline Selection", n_mult, n1_mult, n2_mult);
  PT_AMC     = new TH1F("4_PT_AMC" ,     "Boson PT Distribution after MHT Cut: Baseline Selection",  n_et, n1_et, n2_et);
  HT_AMC     = new TH1F("4_HT_AMC" ,     "HT Distribution After MHT Cut: Baseline Selection", n_ht, n1_ht, n2_ht);
  MHT_AMC    = new TH1F("4_MHT_AMC",     "MHT Distribution After MHT Cut: Baseline Selection", n_et, n1_et, n2_et);
  MHT_J1_AMC = new TH1F("4_MHT_J1_AMC" , "#Delta Phi (MHT - 1st Jet) Distribution After MHT Cut: Baseline Selection", n_pd, n1_pd, n2_pd);
  MHT_J2_AMC = new TH1F("4_MHT_J2_AMC" , "#Delta Phi (MHT - 2nd Jet) Distribution After MHT Cut: Baseline Selection", n_pd, n1_pd, n2_pd);
  MHT_J3_AMC = new TH1F("4_MHT_J3_AMC" , "#Delta Phi (MHT - 3rd Jet) Distribution After MHT Cut: Baseline Selection", n_pd, n1_pd, n2_pd);
  Book_AMC   = new TH1F("4_Book_AMC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_AMC = new TH1F("4_PH_JCl_AMC",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_AMC = new TH1F("4_PH_RJC_AMC",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_AMC = new TH1F("4_PH_MHT_AMC" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);

  JetMult_PTC= new TH1F("5_JetMult_PTC", "PT50ETA25 Jet Multiplicity After Boson Pt Cut", n_mult, n1_mult, n2_mult);
  PT_PTC     = new TH1F("5_PT_PTC" ,     "Boson PT Distribution After Boson Pt Cut",  n_et, n1_et, n2_et);
  HT_PTC     = new TH1F("5_HT_PTC" ,     "HT Distribution After Boson Pt Cut",  n_ht, n1_ht, n2_ht);
  MHT_PTC    = new TH1F("5_MHT_PTC",     "MHT Distribution After Boson Pt Cut", n_et, n1_et, n2_et);
  MHT_J1_PTC = new TH1F("5_MHT_J1_PTC" , "#Delta Phi (MHT - 1st Jet) Distribution After Boson Pt Cut",  n_pd, n1_pd, n2_pd);
  MHT_J2_PTC = new TH1F("5_MHT_J2_PTC" , "#Delta Phi (MHT - 2nd Jet) Distribution After Boson Pt Cut",  n_pd, n1_pd, n2_pd);
  MHT_J3_PTC = new TH1F("5_MHT_J3_PTC" , "#Delta Phi (MHT - 3rd Jet) Distribution After Boson Pt Cut",  n_pd, n1_pd, n2_pd);
  Book_PTC   = new TH1F("5_Book_PTC"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_PTC = new TH1F("5_PH_JCl_PTC",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_PTC = new TH1F("5_PH_RJC_PTC",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_PTC = new TH1F("5_PH_MHT_PTC" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);

  JetMult_HHS= new TH1F("6_JetMult_HHS", "PT50ETA25 Jet Multiplicity: High HT Search Selection", n_mult, n1_mult, n2_mult);
  PT_HHS     = new TH1F("6_PT_HHS" ,     "Boson PT Distribution: High HT Selection",  n_et, n1_et, n2_et);
  HT_HHS     = new TH1F("6_HT_HHS" ,     "HT Distribution: High HT Selection", n_ht, n1_ht, n2_ht);
  MHT_HHS    = new TH1F("6_MHT_HHS",     "MHT Distribution: High HT Selection", n_et, n1_et, n2_et);
  MHT_J1_HHS = new TH1F("6_MHT_J1_HHS" , "#Delta Phi (MHT - 1st Jet) Distribution: High HT Selection", n_pd, n1_pd, n2_pd);
  MHT_J2_HHS = new TH1F("6_MHT_J2_HHS" , "#Delta Phi (MHT - 2nd Jet) Distribution: High HT Selection", n_pd, n1_pd, n2_pd);
  MHT_J3_HHS = new TH1F("6_MHT_J3_HHS" , "#Delta Phi (MHT - 3rd Jet) Distribution: High HT Selection", n_pd, n1_pd, n2_pd);
  Book_HHS   = new TH1F("6_Book_HHS"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);
  PH_JCl_HHS = new TH1F("6_PH_JCl_HHS",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_HHS = new TH1F("6_PH_RJC_HHS",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_HHS = new TH1F("6_PH_MHT_HHS" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);

  JetMult_HMS= new TH1F("7_JetMult_HMS", "PT50ETA25 Jet Multiplicity: High MHT Selection", n_mult, n1_mult, n2_mult);
  PT_HMS     = new TH1F("7_PT_HMS" ,     "Boson PT Distribution: High MHT Selection",  n_et, n1_et, n2_et);
  HT_HMS     = new TH1F("7_HT_HMS" ,     "HT Distribution: High MHT Selection", n_ht, n1_ht, n2_ht);
  MHT_HMS    = new TH1F("7_MHT_HMS",     "MHT Distribution: High MHT Selection", n_et, n1_et, n2_et);
  MHT_J1_HMS = new TH1F("7_MHT_J1_HMS" , "#Delta Phi (MHT - 1st Jet) Distribution: High MHT Selection", n_pd, n1_pd, n2_pd);
  MHT_J2_HMS = new TH1F("7_MHT_J2_HMS" , "#Delta Phi (MHT - 2nd Jet) Distribution: High MHT Selection", n_pd, n1_pd, n2_pd);
  MHT_J3_HMS = new TH1F("7_MHT_J3_HMS" , "#Delta Phi (MHT - 3rd Jet) Distribution: High MHT Selection", n_pd, n1_pd, n2_pd);
  Book_HMS   = new TH1F("7_Book_HMS"   , "Bookkeeping: | 100 < EB < 120 | 120 < EB | 100 < EE < 120 | 120 < EE || ? | dir | ele | frag | sec |", n_bk, n1_bk, n2_bk);  
  PH_JCl_HMS = new TH1F("7_PH_JCl_HMS",  "#Delta Phi (Boson - closesth Jet) Distribution",  n_pd, n1_pd, n2_pd);
  PH_RJC_HMS = new TH1F("7_PH_RJC_HMS",  "#Delta R (Boson - closesth Jet) Distribution",  n_rd, n1_rd, n2_rd);
  PH_MHT_HMS = new TH1F("7_PH_MHT_HMS" , "#Delta Phi (Boson - MHT) Distribution",  n_pd, n1_pd, n2_pd);
  
  // 7 counters + 4 more offline checks + 2 more for Barrel and Endcap Final Checks
  for(unsigned int i=0; i<13; ++i) { Count.push_back(0); }
  EventCounters = new TH1F("EventCounters", "Events: All Processed | post HLT | preFilter | postStdClean | postPFClean | postPhoton | postPFJets| ISO #gamma |2J != #gamma | 3J != #gamma | RA2 Selected", 11, 0.5, 11.5);
  // 5 RA2 Selection Criteria, so 6 counters (+ 1 - 1) + additional counters for influence Boson PT Cut, High HT Search Region and High MHT Search Region
  for(unsigned int i=0; i<9; ++i) { RA2Count.push_back(0); }
  RA2SelectionHisto = new TH1F("RA2SelectionHisto", "Events: Processed | Passing PAT | with ISO #gamma and 3 Jets | with HT > 300 | with Angular Cuts | with MHT > 150 | Boson PT Cut | High HT | High MHT", 9, 0.5, 9.5);
  // Additional Prediction
  for(unsigned int i=0; i<4; ++i) { AddCount.push_back(0); }
  AdditionalPrediction = new TH1F("AdditionalPrediction", "ISO #gamma + 2 Jets | 2 Jets + MHT > 150 | ISO #gamma + 3 Jets | 3 Jets + MHT > 150", 4, 0.5, 4.5);
  // W statistics
  Wstat = new TH1F("Wstatistics", "Decay Channels: e #nu | #mu #nu | #tau #nu | hadrons | ?", 5, 0.5, 5.5);


  // Error Flags Histogram
  ErrorFlags = new TH1F("ErrorFlags", "No Matched GEN #gamma | no GEN #gamma in event | No Jet == #gamma | No Jets in event", 4, 0.5, 4.5);
  // Debug Stream
  print_tss = 0, print_dss = 0, print_rss = 0, print_lss = 0, print_fss = 0;
  timer.Stop();
  tss<<" BosonAnalyzer :: Constructor :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
  photonFound = 0, promptPhoton = 0, directPhoton = 0, fragPhoton = 0, decayPhoton = 0; 

  evNum = 0; rnNum = 0;
}


BosonAnalyzer::~BosonAnalyzer()
{
  timer.Start();
  outputfile->cd();
  EventCounters->Write();
  RA2SelectionHisto->Write();
  AdditionalPrediction->Write();
  ErrorFlags->Write();
  Wstat->Write();
  histcontainers[0].write(RECO);
  if(Boson==22) {for(int h=0; h<n_PhotonDir; ++h) {photoncontainers[h].write(RECO);}}

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
  PH_JCl_HMS->Write();
  PH_RJC_HMS->Write();
  PH_MHT_HMS->Write();

  timer.Stop();
  tss<<" BosonAnalyzer :: Destructor :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
  tss<<"closing "<<RootFileName.c_str()<<std::endl;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BosonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // EVENT INFORMATION                                                                                                                                                             
  evNum = (iEvent.id()).event();
  rnNum = (iEvent.id()).run();

  // if(rnNum == 146430 && evNum == 7882279) {
  //   std::cout<<"Analyzing Run nr "<<rnNum<<" Event nr "<<evNum<<std::endl; Debug[1] = 1; Debug[2] = 1; Debug[3] = 1; 

  globaltimer.Start();
  print_dss = 0, print_rss = 0, print_lss = 0, print_fss = 0;
  p = 0;
  double evWeight = 1;

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
  const reco::GenParticle * boson = 0;
  const pat::Muon * muon = 0;
  std::vector<const pat::Jet*> Jets;

  if(Boson==22) {
    Selection(iEvent, iSetup, Data, Photon, Jets, evWeight, nVertices);
    if(Photon != 0) RA2Selection(Photon, boson, muon, Jets);
  }
  else if(Boson==23) {
    ZBosonSelection(iEvent, iSetup, boson, Jets);
    if(boson != 0) RA2Selection(Photon, boson, muon, Jets);
  }
  else if(Boson==24) {
    WBosonSelection(iEvent, iSetup, boson, muon, Jets);
    if(boson != 0 || muon != 0) RA2Selection(Photon, boson, muon, Jets);
  }

  globaltimer.Stop();
  tss<<" BosonAnalyzer :: Analyze :: Real Time; "<<globaltimer.RealTime()<<" CPU Time: "<<globaltimer.CpuTime()<<std::endl;
  Fill(histcontainers[0].get("Timing","tAna_Real"), n_t1, n1_t1, n2_t1, globaltimer.RealTime());
  Fill(histcontainers[0].get("Timing","tAna_CPU"),  n_t2, n1_t2, n2_t2, globaltimer.CpuTime());
  
  if(Debug[1] || print_dss) std::cout<<" "<<std::endl<<"Matching Debug Stream:    "<<std::endl<<"--------------------------"<<std::endl<<dss.str()<<std::endl; else dss.str("");
  if(Debug[2] || print_rss) std::cout<<" "<<std::endl<<"RA2Selection Debug Stream:"<<std::endl<<"--------------------------"<<std::endl<<rss.str()<<std::endl; else rss.str("");
  if(Debug[3] || print_lss) std::cout<<" "<<std::endl<<"Luminosity Debug Stream:  "<<std::endl<<"--------------------------"<<std::endl<<lss.str()<<std::endl; else lss.str("");
  if(            print_fss) std::cout<<" "<<std::endl<<"Histo Fill Debug Stream:  "<<std::endl<<"--------------------------"<<std::endl<<fss.str()<<std::endl; else fss.str("");
}//}

void
BosonAnalyzer::Selection(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool& Data, const pat::Photon * & photon, std::vector<const pat::Jet*>& Jets, double evW, int nVertices)
{
  timer.Start();
  photonFound = 0, promptPhoton = 0, directPhoton = 0, fragPhoton = 0, decayPhoton = 0;

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

  // We use PAT Photons!
  for(unsigned int i=0; i<patPhotons->size(); ++i) {
     p1 = &((*patPhotons)[i]);
     int count_good_photons = 0;
     // Cuts on H/E, Et and Eta
     if (p1->hadronicOverEm()<photon_he && p1->hasPixelSeed()==false  && p1->et() > photon_ptcut && ((p1->isEE() && fabs(p1->eta()) < 2.5) || p1->isEB())) {
       // Cut on sigma_IEtaIEta
       double p_sIEIE = p1->sigmaIetaIeta();
       if ((p1->isEB() && p_sIEIE < sigma_barrel) || (p1->isEE() && p_sIEIE < sigma_endcap)) {
	 // We only work with Isolated Photons as defined by Egamma POG:
	 double photon_pt = p1->pt(), photon_et = p1->et();
	 float photon_holTrkIso = p1->trkSumPtHollowConeDR04();
	 float photon_ecalIso = p1->ecalRecHitSumEtConeDR04();
	 float photon_hcalIso = p1->hcalTowerSumEtConeDR04();
	 if(photon_holTrkIso < (egm_trk + egm_trk_prop*photon_pt) && photon_ecalIso < (egm_ecal + egm_ecal_prop*photon_et) && photon_hcalIso < (egm_hcal + egm_hcal_prop*photon_et)) {
	   // Some Final Spike Cleaning
	   // ECAL Lazy Tools: sigma_IEtaIEta, sigma_IEtaIPhi, sigma_IPhiIPhi and eMax
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
       }
     }
  }
  if (p != 0) {
    photonFound = 1;
    photon = p;
    double p_sIEIE = p->sigmaIetaIeta();
    const reco::CaloClusterPtr  seed = p->superCluster()->seed();
    std::vector<float> viCov = lazyTool.localCovariances(*seed);
    float p_sIPIP = sqrt(viCov[2]);
    double p_holTrkIso = p->trkSumPtHollowConeDR04(), p_ecalIso = p->ecalRecHitSumEtConeDR04(), p_hcalIso = p->hcalTowerSumEtConeDR04();
    double p_pt = p->pt(), p_et = p->et(), p_eta = p->eta(), p_phi = p->phi();
    double p_maxholTrkIso = egm_trk + egm_trk_prop*p_pt, p_maxecalIso = egm_ecal + egm_ecal_prop*p_et, p_maxhcalIso = egm_hcal + egm_hcal_prop*p_et;
    // If MC --> who's your mother?
    //           1) who is the related genparticle?
    //           2) who is the mother of the genparticle?
    if(!Data) {
      const reco::GenParticle & gen = (reco::GenParticle &) * (p->genPhoton());
      g = &gen; 
      if (g != 0) {
	Fill(histcontainers[0].get("Event","PATMatch"), n_ps, n1_ps, n2_ps, 1);
	Fill(histcontainers[0].get("Event","Distance"), n_dist, n1_dist, n2_dist, angle(g,p));
	// Matching Debug Stream
	dss<<"GEN Photon found: ";
	dss<<" | pt = "<<std::setw(12)<<g->pt()<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	dss<<"PAT Photon found: ";
	dss<<" | pt = "<<std::setw(12)<<p->pt()<<" GeV/c | et = "<<std::setw(12)<<p->et()<<" GeV | eta = "<<std::setw(12)<<p->eta()<<" | phi = "<<std::setw(12)<<p->phi()<<std::endl; 
	// mother?
	const reco::GenParticle & mot = (reco::GenParticle &) *(g->mother());
	dss<<"GenParticle | id = "<<std::setw(5)<<g->pdgId()<<" | st = "<<std::setw(5)<<g->status()<<" | pt = "<<std::setw(12)<<g->pt();
	dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	if(&mot != 0) {
	  dss<<"MotherPart  | id = "<<std::setw(5)<<mot.pdgId()<<" | st = "<<std::setw(5)<<mot.status()<<" | pt = "<<std::setw(12)<<mot.pt();
	  dss<<" GeV/c | et = "<<std::setw(12)<<mot.et()<<" GeV | eta = "<<std::setw(12)<<mot.eta()<<" | phi = "<<std::setw(12)<<mot.phi()<<std::endl;
	  
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
	  
	  if( (mot.status() == 2 || mot.status() == 3) && mot.pdgId() == 22) { promptPhoton = 1; directPhoton = 1;} 
	  else if (mot.pdgId() < 25 && mot.pdgId() != 22)                    { promptPhoton = 1; fragPhoton = 1; }  
	  else if (mot.status() == 2 && mot.pdgId() > 25)                    { promptPhoton = 0; decayPhoton = 1; }  
	  else {}
	}
      }
      else { 
	ErrorFlags->Fill(1);
	Fill(histcontainers[0].get("Event","PATMatch"), n_ps, n1_ps, n2_ps, 0);
	dss<<"No GEN Photon found by the PAT that matches the PAT Photon, trying manually now"<<std::endl; 
	dss<<"... Not yet implemented"<<std::endl;
	dss<<"... Assumed to be a non Prompt Photon"<<std::endl;
	dss<<"... MANUAL ..."<<std::endl;
	int man_counter = 0;
	edm::Handle<reco::GenParticleCollection>      genParticles;
	iEvent.getByLabel(GenParticles, genParticles);                                                                                                                                 
	for(unsigned int i=0; i<genParticles->size(); ++i) {
	  g = &((*genParticles)[i]);
	  if (g->pdgId() == 22) {
	    ++man_counter;
	    dss<<"Gen Photon Candidate | id = "<<std::setw(5)<<g->pdgId()<<" | st = "<<std::setw(5)<<g->status()<<" | pt = "<<std::setw(12)<<g->pt();
	    dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	    // mother?
	    const reco::GenParticle & mot = (reco::GenParticle &) *(g->mother());
	    dss<<"GenParticle | id = "<<std::setw(5)<<g->pdgId()<<" | st = "<<std::setw(5)<<g->status()<<" | pt = "<<std::setw(12)<<g->pt();
	    dss<<" GeV/c | et = "<<std::setw(12)<<g->et()<<" GeV | eta = "<<std::setw(12)<<g->eta()<<" | phi = "<<std::setw(12)<<g->phi()<<std::endl;
	    if(&mot != 0) {
	      dss<<"MotherPart  | id = "<<std::setw(5)<<mot.pdgId()<<" | st = "<<std::setw(5)<<mot.status()<<" | pt = "<<std::setw(12)<<mot.pt();
	      dss<<" GeV/c | et = "<<std::setw(12)<<mot.et()<<" GeV | eta = "<<std::setw(12)<<mot.eta()<<" | phi = "<<std::setw(12)<<mot.phi()<<std::endl;
	    }
	  }
	}
	if(man_counter==0) ErrorFlags->Fill(2);
	print_dss=1;
      }
    }
    /* Misceleanous
    // PF Photons
    if(usePfPhotons) {
    for(unsigned int i=0; i<recoParticles->size(); ++i) {
    p = &((*recoParticles)[i]); // equivalent: recoParticles+i
    int pid = p->pdgId(); double ppt = p->pt(); double p_eta = fabs(p->eta());
    if (pid == 22 && ppt >= photon_ptcut && p_eta < photon_etacut && (p_eta < EB_accept || p_eta > EE_accept)) Photons.push_back(p); // Since There is no Reference to original reco::photon, not much photon ID can be applied ...
    }
    }
    // PAT Photons
    else if (usePatPhotons) {
    }
    // RECO Photons
    else if (useRecoPhotons){
    for(unsigned int i=0; i<photons->size(); ++i) {
      p = &((*photons)[i]);
      const reco::Photon * ph = &((*photons)[i]);
      double p_sIEIE = ph->sigmaIetaIeta();
      double p_eta = fabs(ph->eta());
      if (ph->hadronicOverEm()<photon_he && ph->hasPixelSeed()==false && ((ph->isEB() && p_sIEIE < sigma_barrel) || (ph->isEE() && p_sIEIE < sigma_endcap)) && ph->et() > photon_ptcut && 
      p_eta < photon_etacut && (p_eta < EB_accept || p_eta > EE_accept)) Photons.push_back(p);
      }  
      }
    */
    int PhotonJetIndex = -1, ClosesthJetIndex = -1; 
    double min_r_pi = 0, min_pt_pi = 0;
    dss<<"PAT Photon we work with: | pt = "<<std::setw(12)<<p->pt()<<" GeV/c | et = "<<std::setw(12)<<p->et()<<" GeV | eta = "<<std::setw(12)<<p->eta()<<" | phi = "<<std::setw(12)<<p->phi()<<std::endl;
    dss<<"Isolation values for this photon:      | trk = "<<std::setw(12)<<p_holTrkIso<<" GeV/c | ecal = "<<std::setw(12)<<p_ecalIso<<" GeV | hcal = "<<std::setw(12)<<p_hcalIso<<" GeV"<<std::endl;
    dss<<"Isolation max values for this photon : | trk = "<<std::setw(12)<<p_maxholTrkIso<<" GeV/c | ecal = "<<std::setw(12)<<p_maxecalIso<<" GeV | hcal = "<<std::setw(12)<<p_maxhcalIso<<" GeV"<<std::endl;
    dss<<"         Would it have been Isolated?: | Iso? = "<<(p_holTrkIso < p_maxholTrkIso && p_ecalIso < p_maxecalIso  && p_hcalIso < p_maxhcalIso)<<std::endl;


    // Make a Good Jet Selection (removing the Photon) 
    double photon_pt, photon_eta, photon_phi;
    photon_pt = p->pt(); photon_eta = p->eta(); photon_phi = p->phi();
    for(unsigned int i=0; i<recoJets->size(); ++i) {
      r = &((*recoJets)[i]);
      double jet_pt, jet_eta, jet_phi;
      jet_pt = r->pt();    jet_eta = r->eta();    jet_phi = r->phi();
      double r_pi = angle(photon_eta, jet_eta, photon_phi, jet_phi);
      double pt_pi = (photon_pt-jet_pt)/photon_pt;
      if(r->pt() > 25) {
	dss<<"PF Jet; pt = "<<std::setw(12)<<r->pt()<<" GeV/c | et = "<<std::setw(12)<<r->et()<<" GeV | eta = "<<std::setw(12)<<r->eta()<<" | phi = "<<std::setw(12)<<r->phi();
	dss<<"| dR = "<<std::setw(12)<<r_pi<<" | dPt/Pt = "<<std::setw(12)<<pt_pi<<std::endl;
      }
      // Kinematic & Geometric Cut
      if(/*pt_pi < 0.25 &&*/ r_pi < 0.10 ) {
	PhotonJetIndex = i; 
	dss<<"Jet == Photon Found"<<std::endl;
      }
      else {
	Jets.push_back(r);
	// Delta R and Delta P before Photon Jet Cleaning
	if(i==0) {min_r_pi = r_pi; ClosesthJetIndex = 0;}
	else if (min_r_pi > r_pi) {min_r_pi = r_pi; ClosesthJetIndex = i;}
	else {}
	if(i==0) {min_pt_pi = pt_pi;}
	else if (min_pt_pi > pt_pi) {min_pt_pi = pt_pi;}
	else {}
      }
    }
    if(PhotonJetIndex == -1) {  // !!! This is a problem !!!
      dss<<"Jet == Photon *Not* Found"<<std::endl; print_dss = 1; 
      ErrorFlags->Fill(3);
      r = &((*recoJets)[ClosesthJetIndex]);
      if(ClosesthJetIndex != -1) { 
	dss<<"Closesth Jet:  pt = "<<std::setw(12)<<r->pt()<<" GeV/c | et = "<<std::setw(12)<<r->et()<<" GeV | eta = "<<std::setw(12)<<r->eta()<<" | phi = "<<std::setw(12)<<r->phi()<<std::endl; 
      }
      else { dss<<"No Jets Found!"<<std::endl; ErrorFlags->Fill(4);} 
    }
    Fill(histcontainers[0].get("Event", "DR_Before"), n_dist, n1_dist, n2_dist, min_r_pi, evW);                             
    Fill(histcontainers[0].get("Event", "DP_Before"), n_dist, n1_dist, n2_dist, min_pt_pi, evW);                             
    
    min_r_pi = 0, min_pt_pi = 0;
    for(unsigned int i=0; i<Jets.size(); ++i) {
      r = Jets[i];
      double jet_pt, jet_eta, jet_phi;
      jet_pt = r->pt();    jet_eta = r->eta();    jet_phi = r->phi();
      double r_pi = angle(photon_eta, jet_eta, photon_phi, jet_phi);
      double pt_pi = (photon_pt-jet_pt)/photon_pt;
      // Delta R and Delta P before Photon Jet Cleaning
      if(i==0) {min_r_pi = r_pi;}
      else if (min_r_pi > r_pi) {min_r_pi = r_pi;}
      else {}
      if(i==0) {min_pt_pi = pt_pi;}
      else if (min_pt_pi > pt_pi) {min_pt_pi = pt_pi;}
      else {}
    }
    Fill(histcontainers[0].get("Event", "DR_After"), n_dist, n1_dist, n2_dist, min_r_pi, evW);                             
    Fill(histcontainers[0].get("Event", "DP_After"), n_dist, n1_dist, n2_dist, min_pt_pi, evW);                             
    Fill(histcontainers[0].get("Event", "nVertEvt"), n_mult, n1_mult, n2_mult, nVertices);

    // Plots
    FillHistograms(0, photon, Jets, evW);                                     // 01_Egamma
    if(nVertices == 1) FillHistograms(1, photon, Jets, evW);                  // 02_Egamma_1PV 
    if(promptPhoton)   FillHistograms(2, photon, Jets, evW);                  // 03_Prompt_Egamma
    else               FillHistograms(3, photon, Jets, evW);                  // 04_NonPrompt_Egamma
    if(directPhoton)   FillHistograms(4, photon, Jets, evW);                  // 05_Direct_Egamma
    if(fragPhoton)     FillHistograms(5, photon, Jets, evW);                  // 06_Fragment_Egamma
    if(p->isEB())    { FillHistograms(6, photon, Jets, evW);  ++Count[11]; }  // 07_Egamma_Barrel
    if(p->isEE())    { FillHistograms(7, photon, Jets, evW);  ++Count[12]; }  // 08_Egamma_Endcap

    timer.Stop();
    tss<<" BosonAnalyzer :: Selection :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
    histcontainers[0].get("Timing","tSel_Real")->Fill(timer.RealTime());
    histcontainers[0].get("Timing","tSel_CPU")->Fill(timer.CpuTime());

  }
}

// ------------ Prompt Photons in Case of Monte Carlo -----
void
BosonAnalyzer::ZBosonSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::GenParticle * & zboson, std::vector<const pat::Jet*>& Jets) 
{
  // We assume that every event in the Zinv sample contains a Z boson
  // We do not need the properties, only the result for the RA2 Histograms
  zboson = 0;
  muon = 0;
  if(Data == 0) {
    edm::Handle<reco::GenParticleCollection>      genParticles;
    iEvent.getByLabel(GenParticles, genParticles);
    bool MEZbosonFound = 0;
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      if(MEZbosonFound) continue;
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if (g->pdgId() == 23) {
	std::cout<<"ME Zboson: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = "<<std::setw(9)<<g->phi();
	std::cout<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	zboson = &(*g);
      }
    }
    edm::Handle< std::vector<pat::Jet> >          recoJets;
    iEvent.getByLabel(RecoJets, recoJets);
    for(unsigned int i=0; i<recoJets->size(); ++i) {
      r = &((*recoJets)[i]);
      Jets.push_back(r);
    }
  }
  if (Data == 1) {
    // REJECT EVENTS IN CASE ID & ISO LEPTONS
    // Changed: no electrons allowed, only 1 (!!! strictly 1 !!!) muon allowed
    edm::Handle< std::vector<pat::Electron> >  patElectrons;
    edm::Handle< std::vector<pat::Muon> >  patMuons;
    iEvent.getByLabel("patElectronsPFIDIso", patElectrons);
    iEvent.getByLabel("patMuonsPFIDIso", patMuons);
    if (patElectrons->size() != 0 || patMuons->size() > 2 || patMuons->size() < 2) {
      std::cout<<"Event Rejected !!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Electrons: "<<patElectrons->size()<<" PFIDIso Muons: "<<patMuons->size()<<std::endl;
    }
    else if (patMuons->size() == 0) {
      std::cout<<"Empty Event, no muon found!!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Muons: "<<patMuons->size()<<std::endl;
    }
    else if (patElectrons->size() != 0 || patMuons->size() == 2) {
      for(unsigned int i=0; i<patMuons->size(); ++i) {
	std::cout<<"Muon Found !!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Electrons: "<<patElectrons->size()<<" PFIDIso Muons: "<<patMuons->size();
	std::cout<<" Muon pt = "<<((*patMuons)[i]).pt()<<" GeV/c | Eta = "<<((*patMuons)[i]).eta()<<" Phi = "<<((*patMuons)[i]).phi();
        if(((*patMuons)[i]).pt() > 20 && abs(((*patMuons)[i]).eta()) < 2.1) {
          muon = &((*patMuons)[i]);
        }
      }
    }
    else {
      std::cout<<"This should not happen"<<std::endl;
    }
  }
  if(muon != 0 || zboson != 0) {
    edm::Handle< std::vector<pat::Jet> >          recoJets;
    iEvent.getByLabel(RecoJets, recoJets);
    for(unsigned int i=0; i<recoJets->size(); ++i) {
      r = &((*recoJets)[i]);
      // Look for the muon/electron and remove the jets in a 0.1 cone around gen muon / electron
      if(Data == 0) {
	reco::GenParticle::daughters d = zboson->daughterRefVector();
        for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
          if(abs((*it_d)->pdgId()) == 11 || abs((*it_d)->pdgId()) == 13 || abs((*it_d)->pdgId()) == 15) {
            if(angle((*it_d)->eta(), r->eta(), (*it_d)->phi(), r->phi()) < 0.10) {
	      std::cout<<"PF Jet Related to W boson decay product Found: PFJet is not taken into account"<<std::endl;
            }
            else {
              Jets.push_back(r);
            }
          }
        }
      }
      if(Data == 1) {
        if(angle(muon->eta(), r->eta(), muon->phi(), r->phi()) < 0.10) {
	  std::cout<<"PF Jet Related to muon from W boson Found: PFJet is not taken into account"<<std::endl;
        }
        else {
          Jets.push_back(r);
        }
      }
    }
  }
}

void
BosonAnalyzer::WBosonSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::GenParticle * & wboson, const pat::Muon * & muon, std::vector<const pat::Jet*>& Jets)
{
  // We assume that every event in the W sample contains a W boson
  // We do not need the properties, only the result for the RA2 Histograms
  edm::Handle<reco::GenParticleCollection>      genParticles;
  iEvent.getByLabel(GenParticles, genParticles);
  bool MEWbosonFound = 0;
  wboson = 0;
  muon = 0;
  if(Data == 0) {
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      if(MEWbosonFound) continue;
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if (abs(g->pdgId()) == 24) {
	std::cout<<"ME Wboson: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = "<<std::setw(9)<<g->phi();
	std::cout<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	// keep all (ele, mu, tau)
	// wboson = &(*g);
	// keep only mu from W
	reco::GenParticle::daughters d = g->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
	  std::cout<<"     "<<"ME Wboson Daughter: id = "<<std::setw(2)<<(*it_d)->pdgId()<<" | st = "<<std::setw(2)<<(*it_d)->status()<<" | eta = "<<std::setw(9)<<(*it_d)->eta();
	  std::cout<<" | phi = "<<std::setw(9)<<(*it_d)->phi()<<" | pt = "<<std::setw(9)<<(*it_d)->pt()<<" GeV/c"<<std::endl;
	  if(abs((*it_d)->pdgId()) == 11) {
	    Wstat->Fill(1);
	  }
	  else if(abs((*it_d)->pdgId()) == 13) {
	    Wstat->Fill(2);
	    if((*it_d)->pt() > 20 && abs((*it_d)->eta()) < 2.1) {
	      wboson = &(*g);
	    }
	  }
	  else if(abs((*it_d)->pdgId()) == 15) {
	    Wstat->Fill(3);
	  }
	  else if(abs((*it_d)->pdgId()) < 10) {
	    Wstat->Fill(4);
	  }
	  else if(abs((*it_d)->pdgId()) != 24) {
	    Wstat->Fill(5);
	  }
	}
	// keep all mu (from W and from tau)
	// not implemented
      }
    }
  }
  if (Data == 1) {
    // REJECT EVENTS IN CASE ID & ISO LEPTONS
    // Changed: no electrons allowed, only 1 (!!! strictly 1 !!!) muon allowed
    edm::Handle< std::vector<pat::Electron> >  patElectrons;
    edm::Handle< std::vector<pat::Muon> >  patMuons;
    iEvent.getByLabel("patElectronsPFIDIso", patElectrons);
    iEvent.getByLabel("patMuonsPFIDIso", patMuons);
    if (patElectrons->size() != 0 || patMuons->size() > 1) {
      std::cout<<"Event Rejected !!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Electrons: "<<patElectrons->size()<<" PFIDIso Muons: "<<patMuons->size()<<std::endl;
    }
    else if (patMuons->size() == 0) {
      std::cout<<"Empty Event, no muon found!!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Muons: "<<patMuons->size()<<std::endl;
    }
    else if (patElectrons->size() != 0 || patMuons->size() == 1) {
      for(unsigned int i=0; i<patMuons->size(); ++i) {
	std::cout<<"Muon Found !!! Run nr "<<rnNum<<" Event nr "<<evNum<<" PFIDIso Electrons: "<<patElectrons->size()<<" PFIDIso Muons: "<<patMuons->size();
	std::cout<<" Muon pt = "<<((*patMuons)[i]).pt()<<" GeV/c | Eta = "<<((*patMuons)[i]).eta()<<" Phi = "<<((*patMuons)[i]).phi(); 
	if(((*patMuons)[i]).pt() > 20 && abs(((*patMuons)[i]).eta()) < 2.1) {
	  muon = &((*patMuons)[i]);
	}
      }
    }
    else {
      std::cout<<"This should not happen"<<std::endl;
    }
  }
  if(muon != 0 || wboson != 0) { 
    edm::Handle< std::vector<pat::Jet> >          recoJets;
    iEvent.getByLabel(RecoJets, recoJets);
    for(unsigned int i=0; i<recoJets->size(); ++i) {
      r = &((*recoJets)[i]);
      // Look for the muon/electron and remove the jets in a 0.1 cone around gen muon / electron
      if(Data == 0) {
	reco::GenParticle::daughters d = wboson->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
	  if(abs((*it_d)->pdgId()) == 11 || abs((*it_d)->pdgId()) == 13 || abs((*it_d)->pdgId()) == 15) {
	    if(angle((*it_d)->eta(), r->eta(), (*it_d)->phi(), r->phi()) < 0.10) {
	      std::cout<<"PF Jet Related to W boson decay product Found: PFJet is not taken into account"<<std::endl;
	    }
	    else {
	      Jets.push_back(r);
	    }
	  }
	}
      }
      if(Data == 1) {
	if(angle(muon->eta(), r->eta(), muon->phi(), r->phi()) < 0.10) {
	  std::cout<<"PF Jet Related to muon from W boson Found: PFJet is not taken into account"<<std::endl;
	}
	else {
	  Jets.push_back(r);
	}
      }
    }
  }
}


// ------------ Make Table ------------                                                                                                                                                                  
void 
BosonAnalyzer::RA2Selection(const pat::Photon * & photon, const reco::GenParticle * & boson, const pat::Muon * & muon, std::vector<const pat::Jet*>& Jets)
{

  // Variables
  double bosonpt = 0.0;
  int n_jets = 0, n_ra2  = 0;
  double ht  = 0.0;
  reco::MET::LorentzVector mht(0,0,0,0);
  double min_phi_pi = 0.0, min_r_pi = 0.0;
  int ClosesthJetIndex = -1;

  std::cout<<"Photon, Boson or Muon: P = "<<photon<<" B = "<<boson<<" M = "<<muon<<std::endl;
  if(photon != 0) { bosonpt = photon->pt(); std::cout<<"Test: photon = "<<photon<<" pt = "<<photon->pt()<<std::endl; } 
  if(boson != 0)  { bosonpt = boson->pt();  std::cout<<"Test: boson = "<<boson<<" pt = "<<boson->pt()<<std::endl; } 
  if(muon != 0)   { bosonpt = muon->pt();   std::cout<<"Test: muon = "<<muon<<" pt = "<<muon->pt()<<std::endl; } 

  rss<<"Jets with pt > 50 and eta < 2.5):"<<std::endl;
  for(unsigned int i=0; i<Jets.size(); ++i) {
    r = Jets[i];
    if(Jets[i]->pt() > 50 && fabs(Jets[i]->eta()) < 2.50) {
      ++n_jets;
      ht  += Jets[i]->pt();
    }
    if(Jets[i]->pt() > 30)  {
      mht  -= Jets[i]->p4();
    }
    if(photon) {
      double phi_pi = reco::deltaPhi(photon->phi(),r->phi());
      double r_pi = reco::deltaR(photon->phi(), photon->eta(), r->phi(), r->eta());
      if(i==0) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = 0;}
      else if (min_phi_pi > phi_pi) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = i;}
    }
    if(boson) {
      double phi_pi = reco::deltaPhi(boson->phi(),r->phi());
      double r_pi = reco::deltaR(boson->phi(), boson->eta(), r->phi(), r->eta());
      if(i==0) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = 0;}
      else if (min_phi_pi > phi_pi) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = i;}
    }
    if(muon) {
      double phi_pi = reco::deltaPhi(muon->phi(),r->phi());
      double r_pi = reco::deltaR(muon->phi(), muon->eta(), r->phi(), r->eta());
      if(i==0) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = 0;}
      else if (min_phi_pi > phi_pi) {min_phi_pi = phi_pi; min_r_pi = r_pi; ClosesthJetIndex = i;}
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
  double ph_mht = 0.0;
  if(photon) ph_mht = fabs(reco::deltaPhi(photon->phi(), MHT.phi()));
  if(boson)  ph_mht = fabs(reco::deltaPhi(boson->phi(), MHT.phi()));
  if(muon)   ph_mht = fabs(reco::deltaPhi(muon->phi(), MHT.phi()));

  rss<<" Amount of Jets passing pt and eta criteria: "<<n_jets<<" | ht of event = "<<ht<< " GeV/c | mht of event = "<<mht_value<<" GeV/c"<<std::endl;
  if (n_jets > 0) { ++Count[7];}
  if (n_jets > 1) { ++Count[8];}
  if (n_jets > 2) { ++Count[9];}

  int Identity = 0;
  if(bosonpt > PtCut) {
    Fill(JetMult_All, n_mult, n1_mult, n2_mult, n_jets);
    Fill(PT_All,      n_et, n1_et, n2_et, bosonpt);
    Fill(HT_All,      n_ht, n1_ht, n2_ht, ht);
    Fill(MHT_All,     n_et, n1_et, n2_et, mht_value);
    Fill(MHT_J1_All,  n_pd, n1_pd, n2_pd, mht_j1);
    Fill(MHT_J2_All,  n_pd, n1_pd, n2_pd, mht_j2);
    Fill(MHT_J3_All,  n_pd, n1_pd, n2_pd, mht_j3);
    Fill(Book_All, n_bk, n1_bk, n2_bk,Identity+5);
    Fill(PH_JCl_All,  n_pd, n1_pd, n2_pd, ph_jcl);
    Fill(PH_RJC_All,  n_rd, n1_rd, n2_rd, ph_rjc);
    Fill(PH_MHT_All,  n_pd, n1_pd, n2_pd, ph_mht);
    // N Jets > 3
    if(n_jets >= AmountOfJets) {
      ++RA2Count[2];
      Fill(JetMult_AJC, n_mult, n1_mult, n2_mult, n_jets);
      Fill(PT_AJC,      n_et, n1_et, n2_et, bosonpt);
      Fill(HT_AJC,      n_ht, n1_ht, n2_ht, ht);
      Fill(MHT_AJC,     n_et, n1_et, n2_et, mht_value);
      Fill(MHT_J1_AJC,  n_pd, n1_pd, n2_pd, mht_j1);
      Fill(MHT_J2_AJC,  n_pd, n1_pd, n2_pd, mht_j2);
      Fill(MHT_J3_AJC,  n_pd, n1_pd, n2_pd, mht_j3);
      Fill(Book_AJC, n_bk, n1_bk, n2_bk,Identity+5);
      Fill(PH_JCl_AJC,  n_pd, n1_pd, n2_pd, ph_jcl);
      Fill(PH_RJC_AJC,  n_rd, n1_rd, n2_rd, ph_rjc);
      Fill(PH_MHT_AJC,  n_pd, n1_pd, n2_pd, ph_mht);
      // HT > 300 GeV/c
      if (ht > 300) {
	++RA2Count[3];
	Fill(JetMult_AHC, n_mult, n1_mult, n2_mult, n_jets);
	Fill(PT_AHC,      n_et, n1_et, n2_et, bosonpt);
	Fill(HT_AHC,      n_ht, n1_ht, n2_ht, ht);
	Fill(MHT_AHC,     n_et, n1_et, n2_et, mht_value);
	Fill(MHT_J1_AHC,  n_pd, n1_pd, n2_pd, mht_j1);
	Fill(MHT_J2_AHC,  n_pd, n1_pd, n2_pd, mht_j2);
	Fill(MHT_J3_AHC,  n_pd, n1_pd, n2_pd, mht_j3);
	Fill(Book_AHC, n_bk, n1_bk, n2_bk,Identity+5);
	Fill(PH_JCl_AHC,  n_pd, n1_pd, n2_pd, ph_jcl);
	Fill(PH_RJC_AHC,  n_rd, n1_rd, n2_rd, ph_rjc);
	Fill(PH_MHT_AHC,  n_pd, n1_pd, n2_pd, ph_mht);
	// MHT Angular Cuts
	if (mht_j1 > 0.5 && mht_j2 > 0.5 && (AmountOfJets == 2 || (AmountOfJets == 3 && mht_j3 > 0.3))) {
	  ++RA2Count[4];
	  rss<<"Event passed Angular Cuts"<<std::endl;
	  Fill(JetMult_AAC, n_mult, n1_mult, n2_mult, n_jets);
	  Fill(PT_AAC,      n_et, n1_et, n2_et, bosonpt);
	  Fill(HT_AAC,      n_ht, n1_ht, n2_ht, ht);
	  Fill(MHT_AAC,     n_et, n1_et, n2_et, mht_value);
	  Fill(MHT_J1_AAC,  n_pd, n1_pd, n2_pd, mht_j1);
	  Fill(MHT_J2_AAC,  n_pd, n1_pd, n2_pd, mht_j2);
	  Fill(MHT_J3_AAC,  n_pd, n1_pd, n2_pd, mht_j3);
	  Fill(Book_AAC, n_bk, n1_bk, n2_bk,Identity+5);
	  Fill(PH_JCl_AAC,  n_pd, n1_pd, n2_pd, ph_jcl);
	  Fill(PH_RJC_AAC,  n_rd, n1_rd, n2_rd, ph_rjc);
	  Fill(PH_MHT_AAC,  n_pd, n1_pd, n2_pd, ph_mht);
	  // MHT > 150 GeV/c
	  if(mht_value > 150) {
	    ++RA2Count[5]; ++Count[10];
	    n_ra2 = 1;
	    Fill(JetMult_AMC, n_mult, n1_mult, n2_mult, n_jets);
	    Fill(PT_AMC,      n_et, n1_et, n2_et, bosonpt);
	    Fill(HT_AMC,      n_ht, n1_ht, n2_ht, ht);
	    Fill(MHT_AMC,     n_et, n1_et, n2_et, mht_value);
	    Fill(MHT_J1_AMC,  n_pd, n1_pd, n2_pd, mht_j1);
	    Fill(MHT_J2_AMC,  n_pd, n1_pd, n2_pd, mht_j2);
	    Fill(MHT_J3_AMC,  n_pd, n1_pd, n2_pd, mht_j3);
	    Fill(Book_AMC, n_bk, n1_bk, n2_bk,Identity+5);
	    Fill(PH_JCl_AMC,  n_pd, n1_pd, n2_pd, ph_jcl);
	    Fill(PH_RJC_AMC,  n_rd, n1_rd, n2_rd, ph_rjc);
	    Fill(PH_MHT_AMC,  n_pd, n1_pd, n2_pd, ph_mht);
	    // Influence of Boson PT Cut (in case we lower pt cut of photon to 70 GeV or want to study the influence of Z boson pt contribution
	    if(bosonpt > 100) {
	      ++RA2Count[6];
	      Fill(JetMult_PTC, n_mult, n1_mult, n2_mult, n_jets);
	      Fill(PT_PTC,      n_et, n1_et, n2_et, bosonpt);
	      Fill(HT_PTC,      n_ht, n1_ht, n2_ht, ht);
	      Fill(MHT_PTC,     n_et, n1_et, n2_et, mht_value);
	      Fill(MHT_J1_PTC,  n_pd, n1_pd, n2_pd, mht_j1);
	      Fill(MHT_J2_PTC,  n_pd, n1_pd, n2_pd, mht_j2);
	      Fill(MHT_J3_PTC,  n_pd, n1_pd, n2_pd, mht_j3);
	      Fill(Book_PTC, n_bk, n1_bk, n2_bk,Identity+5);
	      Fill(PH_JCl_PTC,  n_pd, n1_pd, n2_pd, ph_jcl);
	      Fill(PH_RJC_PTC,  n_rd, n1_rd, n2_rd, ph_rjc);
	      Fill(PH_MHT_PTC,  n_pd, n1_pd, n2_pd, ph_mht);
	    }
	    // Hight HT Selection Region
	    if(ht > 500) {
	      ++RA2Count[7];
	      Fill(JetMult_HHS, n_mult, n1_mult, n2_mult, n_jets);
	      Fill(PT_HHS,      n_et, n1_et, n2_et, bosonpt);
	      Fill(HT_HHS,      n_ht, n1_ht, n2_ht, ht);
	      Fill(MHT_HHS,     n_et, n1_et, n2_et, mht_value);
	      Fill(MHT_J1_HHS,  n_pd, n1_pd, n2_pd, mht_j1);
	      Fill(MHT_J2_HHS,  n_pd, n1_pd, n2_pd, mht_j2);
	      Fill(MHT_J3_HHS,  n_pd, n1_pd, n2_pd, mht_j3);
	      Fill(Book_HHS, n_bk, n1_bk, n2_bk,Identity+5);
	      Fill(PH_JCl_HHS,  n_pd, n1_pd, n2_pd, ph_jcl);
	      Fill(PH_RJC_HHS,  n_rd, n1_rd, n2_rd, ph_rjc);
	      Fill(PH_MHT_HHS,  n_pd, n1_pd, n2_pd, ph_mht);
	    }
	    // High MHT Selection Region
	    if(mht_value > 250) {
	      ++RA2Count[8];
	      Fill(JetMult_HMS, n_mult, n1_mult, n2_mult, n_jets);
	      Fill(PT_HMS,      n_et, n1_et, n2_et, bosonpt);
	      Fill(HT_HMS,      n_ht, n1_ht, n2_ht, ht);
	      Fill(MHT_HMS,     n_et, n1_et, n2_et, mht_value);
	      Fill(MHT_J1_HMS,  n_pd, n1_pd, n2_pd, mht_j1);
	      Fill(MHT_J2_HMS,  n_pd, n1_pd, n2_pd, mht_j2);
	      Fill(MHT_J3_HMS,  n_pd, n1_pd, n2_pd, mht_j3);
	      Fill(Book_HMS, n_bk, n1_bk, n2_bk,Identity+5);
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
    if(mht_value > 150) ++AddCount[1];
    if(n_jets >= 3) {
      ++AddCount[2];
      if(mht_value > 150) ++AddCount[3];
    }
  }
  /*
  if(photon) {
    if(n_jets >= 2) FillHistograms(8, photon, Jets, evW); // 09_2Jets_Egamma
    // if(n_jets  >= 3) FillHistograms(9, photon, Jets, evW); // 10_3Jets_Egamma
    if(n_ra2  == 1) FillHistograms(9, photon, Jets, evW); // 10_3Jets_Egamma ... --> 10_RA2_Egamma
  }
  */
}

// ------------ Fill Photon And Jet Information ----------
void 
BosonAnalyzer::FillHistograms(int h, const pat::Photon * & p, std::vector<const pat::Jet*>& Jets, double evW)
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
  double p_sIEIE = p->sigmaIetaIeta(), p_HE = p->hadronicOverEm();
  bool p_hPS = p->hasPixelSeed(), p_hCT = p->hasConversionTracks();
  int p_nTH = p->nTrkHollowConeDR04();
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
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"SigmaIEtaIEta"), n_se,  n1_se, n2_se, p_sIEIE,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"SigmaIPhiIPhi"), n_s2e, n1_s2e, n2_s2e, p_sIEIE,evW);
  Fill(photoncontainers[h].get(PhotonContainNameArray[h],"nTrackHollow"), n_trk, n1_trk, n2_trk, p_nTH,evW);
  
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
  tss<<" BosonAnalyzer :: FillPhotonAndJetHistograms :: Real Time; "<<timer.RealTime()<<" CPU Time: "<<timer.CpuTime()<<std::endl;
}



void 
BosonAnalyzer::endLuminosityBlock(const edm::LuminosityBlock& iLuminosityBlock, const edm::EventSetup& iSetup) 
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
BosonAnalyzer::Fill(TH1 * Histo, int n, double n1, double n2, double value) {
  double binwidth = (n2-n1)/n;
  double n1_ = n1+binwidth/2;
  double n2_ = n2-binwidth/2;
  if(value > n1 && value < n2) Histo->Fill(value);
  else if (value < n1) Histo->Fill(n1_);
  else if (value > n2) Histo->Fill(n2_);
  else {}
}
void
BosonAnalyzer::Fill(TH1 * Histo, int n, double n1, double n2, double value, double evWeight) {
  double binwidth = (n2-n1)/n;
  double n1_ = n1+binwidth/2;
  double n2_ = n2-binwidth/2;
  if(value > n1 && value < n2) Histo->Fill(value, evWeight);
  else if (value < n1) Histo->Fill(n1_, evWeight);
  else if (value > n2) Histo->Fill(n2_, evWeight);
  else {}
}


// ------------ method for calculating distance between photon en genJet            -----
double
BosonAnalyzer::angle(double eta1, double eta2, double phi1, double phi2) {
  double angle;
  double abs_diff_eta = fabs(eta1-eta2);
  // double abs_diff_phi = fabs(phi1-phi2);
  // if (abs_diff_phi > pi) abs_diff_phi-=pi;
  double abs_diff_phi = fabs(reco::deltaPhi(phi1,phi2));
  angle = sqrt(pow(abs_diff_eta,2) + pow(abs_diff_phi,2));
  return angle;
}

double 
BosonAnalyzer::angle(const reco::LeafCandidate * a, const reco::LeafCandidate * b) {
  double angle;
  double abs_diff_eta = fabs(a->eta()-b->eta());
  double abs_diff_phi = fabs (reco::deltaPhi(a->phi(),b->phi()));
  angle = sqrt(pow(abs_diff_eta,2) + pow(abs_diff_phi,2));
  return angle;
}

double
BosonAnalyzer::angle(const reco::LeafCandidate * a, const pat::Jet * b) {
  double angle;
  double abs_diff_eta = fabs(a->eta()-b->eta());
  double abs_diff_phi = fabs (reco::deltaPhi(a->phi(),b->phi()));
  angle = sqrt(pow(abs_diff_eta,2) + pow(abs_diff_phi,2));
  return angle;
}

// ------------ method called once each job just before starting event loop  ------------
void 
BosonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BosonAnalyzer::endJob() 
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


  for(unsigned int i=0; i<Count.size(); ++i)    { EventCounters->SetBinContent(i+1, Count[i]); }
  for(unsigned int i=0; i<RA2Count.size(); ++i) { RA2SelectionHisto->SetBinContent(i+1, RA2Count[i]); }
  for(unsigned int i=0; i<AddCount.size(); ++i) { AdditionalPrediction->SetBinContent(i+1, AddCount[i]); }

  
  std::string ECBinLabels [] = {"Processed Events","post HLT","preFilter","postStdClean","postPFClean","postPhoton","postPFJets","ISO #gamma > 100 GeV","2J != #gamma","3J != #gamma","RA2 Selected"};
  std::string RSBinLabels [] = {"Processed Events","After PAT Selection", "ISO #gamma > 100 + 3J", "HT > 300 GeV/c", "Angular Cuts", "MHT > 150 GeV/c", "Boson PT Cut", "High HT", "High MHT"};
  std::string ACBinLabels [] = {"ISO #gamma + 2 Jets", "2 Jets + MHT > 150", "ISO #gamma + 3 Jets", "3 Jets + MHT > 150"};
  std::string EFBinLabels [] = {"No Matched GEN #gamma", "no GEN #gamma in event", "No Jet == #gamma", "No Jets in event"};
  
  for(unsigned int i=0; i<(Count.size()-2); ++i) { EventCounters->GetXaxis()->SetBinLabel(i+1, ECBinLabels[i].c_str()); }
  for(unsigned int i=0; i<RA2Count.size(); ++i)  { RA2SelectionHisto->GetXaxis()->SetBinLabel(i+1, RSBinLabels[i].c_str()); }
  for(unsigned int i=0; i<AddCount.size(); ++i)  { AdditionalPrediction->GetXaxis()->SetBinLabel(i+1, ACBinLabels[i].c_str());}
  for(int i=0; i<4; ++i)  { ErrorFlags->GetXaxis()->SetBinLabel(i+1, EFBinLabels[i].c_str()); }
  

}
    
//define this as a plug-in
DEFINE_FWK_MODULE(BosonAnalyzer);
