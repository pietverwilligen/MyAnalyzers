// -*- C++ -*-
//
// Package:    PhotonPreSelectionPlotter
// Class:      PhotonPreSelectionPlotter
// 
/**\class PhotonPreSelectionPlotter PhotonPreSelectionPlotter.cc MyAnalyzers/PhotonPreSelectionPlotter/src/PhotonPreSelectionPlotter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  local user
//         Created:  Thu Nov 25 19:02:51 CET 2010
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

// ROOT include files
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectoryFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"





#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/EgammaCandidates/interface/Photon.h>
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"


//
// class declaration
//

class PhotonPreSelectionPlotter : public edm::EDAnalyzer {
   public:
      explicit PhotonPreSelectionPlotter(const edm::ParameterSet&);
      ~PhotonPreSelectionPlotter();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  void Fill(TH1F *, int, double, double, double);
  
      // ----------member data ---------------------------
  TFile * outputfile;
  TDirectoryFile *PHOTON, *PSPIKE, *PSTRAT;

  TH1F *PB_Et_B, *PB_Eta_B, *PB_TRK, *PB_ECAL, *PB_HCAL, *PB_PS, *PB_HE, *PB_SIEIE, *PB_SIPIP, *PB_Et_A, *PB_Eta_A;
  TH1F *PB_eM_B, *PB_sT_B, *PB_RF_B, *PB_sS_B, *PB_SC_B, *PB_EE_B, *PB_eM_A, *PB_sT_A, *PB_RF_A, *PB_sS_A, *PB_SC_A, *PB_EE_A; 
  TH1F *PE_Et_B, *PE_Eta_B, *PE_TRK, *PE_ECAL, *PE_HCAL, *PE_PS, *PE_HE, *PE_SIEIE, *PE_SIPIP, *PE_Et_A, *PE_Eta_A;
  TH1F *PE_eM_B, *PE_sT_B, *PE_RF_B, *PE_sS_B, *PE_SC_B, *PE_EE_B, *PE_eM_A, *PE_sT_A, *PE_RF_A, *PE_sS_A, *PE_SC_A, *PE_EE_A;
  TH1F *P_Et_B, *P_Eta_B, *P_Et_A, *P_Eta_A;

  TH1F *PB_EE_Strat1,                      *PB_sT_Strat1;
  TH1F *PB_SIEIE_Strat2, *PB_SIPIP_Strat2, *PB_sS_Strat2; 
  TH1F *PB_SIEIE_Strat3, *PB_SIPIP_Strat3, *PB_sT_Strat3;

  TH1F *Selection1, *Selection2;
  const reco::Photon           *p, *p1, *p2;
};

//
// constants, enums and typedefs
//
double pi = 3.1415926535;

// Photon ID cuts:                                                                                                                                                                                         
int    photon_ptcut = 30.0;
int    photon_etacut = 2.50;
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
int n_phi  = 18;  double n1_phi  = -3.14, n2_phi  = 3.14;
int n_phi2 = 360; double n1_phi2 = -3.14, n2_phi2 = 3.14;
int n_et   = 100; double n1_et   = 0.0,   n2_et   = 500;
int n_iso  = 32;  double n1_iso  = -1.0,   n2_iso = 15;
int n_mult = 11;  double n1_mult = -0.5,  n2_mult = 10.5;
int n_pe   = 15;  double n1_pe   = 0.0,   n2_pe   = 3.0;
int n_he   = 20;  double n1_he   = 0.00,  n2_he   = 0.10;  
int n_sb   = 30;  double n1_sb   = 0.00,  n2_sb   = 0.03; 
int n_se   = 30;  double n1_se   = 0.015, n2_se   = 0.03;
int n_s2e  = 45;  double n1_s2e  = 0.015, n2_s2e  = 0.06;
int n_ps   = 2;   double n1_ps   = -0.5,  n2_ps   = 1.5;  

int n_sT   = 100; double n1_sT   = -25,   n2_sT   = 25;
int n_RF   = 17;  double n1_RF   = -0.5,  n2_RF   = 15.5;
int n_sS   = 6;   double n1_sS   = -0.5,  n2_sS   = 5.5;
int n_SC   = 60;  double n1_SC   = 0.0,   n2_SC   = 1.2;
//
// static data member definitions
//

//
// constructors and destructor
//
PhotonPreSelectionPlotter::PhotonPreSelectionPlotter(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  outputfile = new TFile("PhotonPreSelectionPlots.root", "RECREATE" );

  PHOTON = (TDirectoryFile*) outputfile->mkdir("PHOTON","PHOTON");
  PSPIKE = (TDirectoryFile*) outputfile->mkdir("PSPIKE","PSPIKE");
  PSTRAT = (TDirectoryFile*)     PSPIKE->mkdir("PSTRAT","PSTRAT");

  PB_TRK   = new TH1F("PB_TRK",        "N-1 Plot Photon Hollow Trk Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso);
  PB_ECAL  = new TH1F("PB_ECAL",       "N-1 Plot Photon Ecal Iso Deposit (DR = 0.4) Distribution", n_iso, n1_iso, n2_iso);
  PB_HCAL  = new TH1F("PB_HCAL",       "N-1 Plot Photon Hcal Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso);
  PB_PS    = new TH1F("PB_PS",         "N-1 Plot Has Electron Pixel Seeds", n_ps, n1_ps, n2_ps);
  PB_HE    = new TH1F("PB_HE",         "N-1 Plot Hadronic over Electromagnetic Energy Fraction", n_he, n1_he, n2_he);
  PB_SIEIE = new TH1F("PB_SIEIE",      "N Plot #sigma_{i #eta i #eta} ShowerShape Variable", n_sb, n1_sb, n2_sb);
  PB_SIPIP = new TH1F("PB_SIPIP",      "N Plot #sigma_{i #phi i #phi} ShowerShape Variable", n_sb, n1_sb, n2_sb);
  PB_Et_B  = new TH1F("PB_Et_B",       "N Plot Photon Et Distribution before #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PB_Eta_B = new TH1F("PB_Eta_B",      "N Plot Photon Eta Distribution before #sigma_{i #eta i #eta} cut", n_eta, n1_eta, n2_eta);
  PB_Et_A  = new TH1F("PB_Et_A",       "N+1 Plot Photon Et Distribution after #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PB_Eta_A = new TH1F("PB_Eta_A",      "N+1 Plot Photon Eta Distribution after #sigma_{i #eta i #eta} cut", n_eta, n1_eta, n2_eta);

  PE_TRK   = new TH1F("PE_TRK",        "N-1 Plot Photon Hollow Trk Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso);
  PE_ECAL  = new TH1F("PE_ECAL",       "N-1 Plot Photon Ecal Iso Deposit (DR = 0.4) Distribution", n_iso, n1_iso, n2_iso);
  PE_HCAL  = new TH1F("PE_HCAL",       "N-1 Plot Photon Hcal Iso (DR = 0.4) Deposit Distribution", n_iso, n1_iso, n2_iso);
  PE_PS    = new TH1F("PE_PS",         "N-1 Plot Has Electron Pixel Seeds", n_ps, n1_ps, n2_ps);
  PE_HE    = new TH1F("PE_HE",         "N-1 Plot Hadronic over Electromagnetic Energy Fraction", n_he, n1_he, n2_he);
  PE_SIEIE = new TH1F("PE_SIEIE",      "N Plot #sigma_{i #eta i #eta} ShowerShape Variable", n_se, n1_se, n2_se);
  PE_SIPIP = new TH1F("PE_SIPIP",      "N Plot #sigma_{i #phi i #phi} ShowerShape Variable", n_s2e, n1_s2e, n2_s2e);
  PE_Et_B  = new TH1F("PE_Et_B",       "N Plot Photon Et Distribution before #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PE_Eta_B = new TH1F("PE_Eta_B",      "N plot Photon Eta Distribution before #sigma_{i #eta i #eta} cut", n_eta, n1_eta, n2_eta);
  PE_Et_A  = new TH1F("PE_Et_A",       "N+1 Plot Photon Et Distribution after #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PE_Eta_A = new TH1F("PE_Eta_A",      "N+1 Plot Photon Eta Distribution after #sigma_{i #eta i #eta} cut", n_eta, n1_eta, n2_eta);

  P_Et_B  = new TH1F("P_Et_B",       "N Plot Photon Et Distribution before #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  P_Eta_B = new TH1F("P_Eta_B",      "N Plot Photon Eta Distribution before #sigma_{i #eta i #eta} cut", n_eta, n1_eta, n2_eta);
  P_Et_A  = new TH1F("P_Et_A",       "N+1 Plot Photon Et Distribution after #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  P_Eta_A = new TH1F("P_Eta_A",      "N+1 Plot Photon Eta Distribution after #sigma_{i #eta i #eta} cut", n_eta, n1_eta, n2_eta);

  Selection1 = new TH1F("Selection", "[-] | HE | PS | SIEIE | ISO TRK | ISO ECAL | ISO HCAL | No Spike | #gamma > 100 GeV", 9, 0.5, 9.5); 
  Selection2 = new TH1F("Selection", "[-] | HE | PS | ISO TRK | ISO ECAL | ISO HCAL | SIEIE | No Spike | #gamma > 100 GeV", 9, 0.5, 9.5); 

  PB_eM_B = new TH1F("PB_eM_B",       "N Plot Photon eMax distribution before #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PB_sT_B = new TH1F("PB_sT_B",       "N Plot Photon Seed Time distribution before #sigma_{i #eta i #eta} cut", n_sT, n1_sT, n2_sT);
  PB_RF_B = new TH1F("PB_RF_B",       "N Plot Photon Reco Flag distribution before #sigma_{i #eta i #eta} cut", n_RF, n1_RF, n2_RF);
  PB_sS_B = new TH1F("PB_sS_B",       "N Plot Photon Severity distribution before #sigma_{i #eta i #eta} cut", n_sS, n1_sS, n2_sS);
  PB_SC_B = new TH1F("PB_SC_B",       "N Plot Photon Swiss Cross (E1/E4) distribution before #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);
  PB_EE_B = new TH1F("PB_EE_B",       "N Plot Photon E2/E9 distribution before #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);
  PE_eM_B = new TH1F("PE_eM_B",       "N Plot Photon eMax distribution before #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PE_sT_B = new TH1F("PE_sT_B",       "N Plot Photon Seed Time distribution before #sigma_{i #eta i #eta} cut", n_sT, n1_sT, n2_sT);
  PE_RF_B = new TH1F("PE_RF_B",       "N Plot Photon Reco Flag distribution before #sigma_{i #eta i #eta} cut", n_RF, n1_RF, n2_RF);
  PE_sS_B = new TH1F("PE_sS_B",       "N Plot Photon Severity distribution before #sigma_{i #eta i #eta} cut", n_sS, n1_sS, n2_sS);
  PE_SC_B = new TH1F("PE_SC_B",       "N Plot Photon Swiss Cross (E1/E4) distribution before #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);
  PE_EE_B = new TH1F("PE_EE_B",       "N Plot Photon E2/E9 distribution before #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);

  PB_eM_A = new TH1F("PB_eM_A",       "N Plot Photon eMax distribution after #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PB_sT_A = new TH1F("PB_sT_A",       "N Plot Photon Seed Time distribution after #sigma_{i #eta i #eta} cut", n_sT, n1_sT, n2_sT);
  PB_RF_A = new TH1F("PB_RF_A",       "N Plot Photon Reco Flag distribution after #sigma_{i #eta i #eta} cut", n_RF, n1_RF, n2_RF);
  PB_sS_A = new TH1F("PB_sS_A",       "N Plot Photon Severity distribution after #sigma_{i #eta i #eta} cut", n_sS, n1_sS, n2_sS);
  PB_SC_A = new TH1F("PB_SC_A",       "N Plot Photon Swiss Cross (E1/E4) distribution after #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);
  PB_EE_A = new TH1F("PB_EE_A",       "N Plot Photon E2/E9distribution after #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);
  PE_eM_A = new TH1F("PE_eM_A",       "N Plot Photon eMax distribution after #sigma_{i #eta i #eta} cut", n_et, n1_et, n2_et);
  PE_sT_A = new TH1F("PE_sT_A",       "N Plot Photon Seed Time distribution after #sigma_{i #eta i #eta} cut", n_sT, n1_sT, n2_sT);
  PE_RF_A = new TH1F("PE_RF_A",       "N Plot Photon Reco Flag after #sigma_{i #eta i #eta} cut", n_RF, n1_RF, n2_RF);
  PE_sS_A = new TH1F("PE_sS_A",       "N Plot Photon Severity Flag after #sigma_{i #eta i #eta} cut", n_sS, n1_sS, n2_sS);
  PE_SC_A = new TH1F("PE_SC_A",       "N Plot Photon Swiss Cross (E1/E4) distribution after #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);
  PE_EE_A = new TH1F("PE_EE_A",       "N Plot Photon E2/E9 distribution after #sigma_{i #eta i #eta} cut", n_SC, n1_SC, n2_SC);


  // Spike Cut Strategy 1: sigma_IEtaIEta > 0.001 && sigma_IPhiIPhi > 0.001 && severityFlag != 3
  PB_EE_Strat1 = new TH1F("PB_EE_Strat1",       "N + 1 Plot Photon E2/E9distribution after Strat1", n_SC, n1_SC, n2_SC);
  PB_sT_Strat1 = new TH1F("PB_sT_Strat1",       "N + 1 Plot Photon Seed Time distribution after Strat1", n_sT, n1_sT, n2_sT);

  // Spike Cut Strategy 2: E2/E9 < 0.95 && seedTime > -3.5 ns
  PB_SIEIE_Strat2 = new TH1F("PB_SIEIE_Strat2",      "N + 1 Plot #sigma_{i #eta i #eta} after Strat 2", n_sb, n1_sb, n2_sb);
  PB_SIPIP_Strat2 = new TH1F("PB_SIPIP_Strat2",      "N + 1 Plot #sigma_{i #phi i #phi} after Strat 2", n_sb, n1_sb, n2_sb);
  PB_sS_Strat2    = new TH1F("PB_sS_Strat2",         "N + 1 Plot Photon Severity distribution after Strat 2", n_sS, n1_sS, n2_sS);

  // Spike Cut Strategy 3: E2/E9 < 0.95 && severityFlag != 3
  PB_SIEIE_Strat3 = new TH1F("PB_SIEIE_Strat3",      "N + 1 Plot #sigma_{i #eta i #eta} ShowerShape Variable Strat 3", n_sb, n1_sb, n2_sb);
  PB_SIPIP_Strat3 = new TH1F("PB_SIPIP_Strat3",      "N + 1 Plot #sigma_{i #phi i #phi} ShowerShape Variable Strat 3", n_sb, n1_sb, n2_sb);
  PB_sT_Strat3    = new TH1F("PB_sT_Strat3",         "N + 1 Plot Photon Seed Time distribution Strat 3", n_sT, n1_sT, n2_sT);




  std::string RecoFlagString [] = {"kGood", "kPoorReco", "kOutOfTime", "kFaultyHardware", "kNoisy", "kPoorCalib", "kSaturated", "kLeadingEdgeRecovered", "kNeighboursRecovered", "kTowerRecovered", "kFake", "kFakeNeighbours", "kDead", "kKilled", "kTPSaturated", "kL1SpikeFlag", "kUnknown"};
  std::string severityString [] = {"kGood", "kProblematic", "kRecovered", "kTime", "kWeird", "kBad"};

  for(int i=0; i<n_RF; ++i) {
    PB_RF_B->GetXaxis()->SetBinLabel(i+1, RecoFlagString[i].c_str()); PE_RF_B->GetXaxis()->SetBinLabel(i+1, RecoFlagString[i].c_str());
    PB_RF_A->GetXaxis()->SetBinLabel(i+1, RecoFlagString[i].c_str()); PE_RF_A->GetXaxis()->SetBinLabel(i+1, RecoFlagString[i].c_str());
  }
  for(int i=0; i<n_sS; ++i) {
    PB_sS_B->GetXaxis()->SetBinLabel(i+1, severityString[i].c_str()); PE_sS_B->GetXaxis()->SetBinLabel(i+1, severityString[i].c_str());
    PB_sS_A->GetXaxis()->SetBinLabel(i+1, severityString[i].c_str()); PE_sS_A->GetXaxis()->SetBinLabel(i+1, severityString[i].c_str());
  }


}


PhotonPreSelectionPlotter::~PhotonPreSelectionPlotter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();
  PHOTON->cd();
  PB_Et_B->Write();
  PB_Eta_B->Write();
  PB_TRK->Write();
  PB_ECAL->Write();
  PB_HCAL->Write();
  PB_PS->Write();
  PB_HE->Write();
  PB_SIEIE->Write();
  PB_SIPIP->Write();
  PB_Et_A->Write();
  PB_Eta_A->Write();
  PE_Et_B->Write();
  PE_Eta_B->Write();
  PE_TRK->Write();
  PE_ECAL->Write();
  PE_HCAL->Write();
  PE_PS->Write();
  PE_HE->Write();
  PE_SIEIE->Write();
  PE_SIPIP->Write();
  PE_Et_A->Write();
  PE_Eta_A->Write();
  P_Et_B->Write();
  P_Eta_B->Write();
  P_Et_A->Write();
  P_Eta_A->Write();

  Selection1->Write();
  Selection2->Write();

  outputfile->cd();
  PSPIKE->cd();
  PB_eM_B->Write();
  PB_sT_B->Write();
  PB_RF_B->Write();
  PB_sS_B->Write();
  PB_SC_B->Write();
  PB_EE_B->Write();
  PE_eM_B->Write();
  PE_sT_B->Write();
  PE_RF_B->Write();
  PE_sS_B->Write();
  PE_SC_B->Write();
  PE_EE_B->Write();

  PB_eM_A->Write();
  PB_sT_A->Write();
  PB_RF_A->Write();
  PB_sS_A->Write();
  PB_SC_A->Write();
  PB_EE_A->Write();
  PE_eM_A->Write();
  PE_sT_A->Write();
  PE_RF_A->Write();
  PE_sS_A->Write();
  PE_SC_A->Write();
  PE_EE_A->Write();

  PSTRAT->cd();
  PB_EE_Strat1->Write();
  PB_sT_Strat1->Write();
  PB_SIEIE_Strat2->Write();
  PB_SIPIP_Strat2->Write();
  PB_sS_Strat2->Write();
  PB_SIEIE_Strat3->Write();
  PB_SIPIP_Strat3->Write();
  PB_sT_Strat3->Write();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PhotonPreSelectionPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Selection counter: [-]
  Selection1->Fill(1); Selection2->Fill(1);

  // edm::Handle< std::vector<pat::Photon> >       patPhotons;
  // iEvent.getByLabel("patPhotons", patPhotons);
  edm::Handle< std::vector<reco::Photon> >       recoPhotons;
  iEvent.getByLabel("photons", recoPhotons);

  // Ecal RecHits
  edm::Handle<EcalRecHitCollection> EBReducedRecHits;
  iEvent.getByLabel("reducedEcalRecHitsEB", EBReducedRecHits);
  edm::Handle<EcalRecHitCollection> EEReducedRecHits;
  iEvent.getByLabel("reducedEcalRecHitsEE", EEReducedRecHits); 
  // get the channel status from the DB
  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);

  EcalClusterLazyTools lazyTool(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));

  int n_photons = 0;
  p = 0;
 // Selecting the highest Et Photon
  for(unsigned int i=0; i<recoPhotons->size(); ++i) {
    p1 = &((*recoPhotons)[i]);
    // within Kinematic & Geometric Cuts
    if(p1->et() > photon_ptcut && ((p1->isEE() && fabs(p1->eta()) < 2.5) || p1->isEB())) {
      ++n_photons;
      if (n_photons == 1) p = &(*p1);
      else if (p1->et() > p->et()) p = &(*p1);
      else {}
    }
    else {}
  }
  if(p != 0) {
    // Photon Variables
    double p_sIEIE = p->sigmaIetaIeta(), p_HE = p->hadronicOverEm();
    bool p_hPS = p->hasPixelSeed();
    double p_trkIso = p->trkSumPtHollowConeDR04(), p_ecalIso = p->ecalRecHitSumEtConeDR04(), p_hcalIso = p->hcalTowerSumEtConeDR04();
    double p_pt = p->pt(), p_et = p->et(), p_eta = p->eta(), p_phi = p->phi();
    double p_maxTrkIso = egm_trk + egm_trk_prop*p_pt, p_maxEcalIso = egm_ecal + egm_ecal_prop*p_et, p_maxHcalIso = egm_hcal + egm_hcal_prop*p_et;

    // ECAL Lazy Tools: sigma_IEtaIEta, sigma_IEtaIPhi, sigma_IPhiIPhi and eMax
    const reco::CaloClusterPtr  seed = p->superCluster()->seed();
    std::vector<float> vCov = lazyTool.covariances(*seed);
    float covPhiPhi = vCov[2];
    float covEtaPhi = vCov[1];
    float covEtaEta = vCov[0];
    std::vector<float> viCov = lazyTool.localCovariances(*seed);
    float sigmaIphiIphi = sqrt(viCov[2]);
    float sigmaIetaIphi = sqrt(viCov[1]);
    float eMax = lazyTool.eMax(*seed), e2nd = lazyTool.e2nd(*seed), e3x3 = lazyTool.e3x3(*seed);
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

    // From: EcalRecHit.h: seedRecoFlag
    // --------------------------------
    // recHit flags
    /* enum Flags { 
            kGood=0,                   // channel ok, the energy and time measurement are reliable
            kPoorReco,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
            kOutOfTime,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
            kFaultyHardware,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
            kNoisy,                    // the channel is very noisy
            kPoorCalib,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
            kSaturated,                // saturated channel (recovery not tried)
            kLeadingEdgeRecovered,     // saturated channel: energy estimated from the leading edge before saturation
            kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
            kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
            kFake,                     // the signal in the channel is a fake (e.g. a so-called spike)
            kFakeNeighbours,           // the signal in the channel is a fake and it is detected by looking at the neighbours
            kDead,                     // channel is dead and any recovery fails
            kKilled,                   // MC only flag: the channel is killed in the real detector
            kTPSaturated,              // only for flagBits_: the channel is in a region with saturated TP
            kL1SpikeFlag,              // only for flagBits_: the channel is in a region with TP with sFGVB = 0
                                       // pro tempore, this will become obsolete when the online protection against spikes will be activated
                                       //
            kUnknown                   // to easy the interface with functions returning flags
    }*/


    // From: EcalSeverityLevelAlgo.h: seedSeverity
    // -------------------------------------------
    // give the severity level from the EcalRecHit flags + the DB information stored in EcalChannelStatus
    // Levels of severity:
    // - kGood        --> good channel
    // - kProblematic --> problematic (e.g. noisy)
    // - kRecovered   --> recovered (e.g. an originally dead or saturated)
    // - kTime        --> the channel is out of time (e.g. spike)
    // - kWeird       --> weird (e.g. spike)
    // - kBad         --> bad, not suitable to be used in the reconstruction





    // Spike Cleaning : AN-2010-268 for CMSSW_3_6_X
    // Spike <=> (isEB && (seedSeverity!=3 && seedSeverity!=4 ) && (seedRecoFlag != 2 || (eMax>130&&seedTime>0) ) && sigmaIetaIeta > 0.001 && sqrt(covPhiPhi)>0.001 )
    double p_sIPIP = sigmaIphiIphi;

    /*
    std::cout<<"PAT Photon we work with: | pt = "<<std::setw(12)<<p->pt()<<" GeV/c | et = "<<std::setw(12)<<p->et()<<" GeV | eta = "<<std::setw(12)<<p->eta()<<" | phi = "<<std::setw(12)<<p->phi()<<std::endl;
    std::cout<<"Isolation values for this photon:      | trk = "<<std::setw(12)<<p_trkIso<<" GeV/c | ecal = "<<std::setw(12)<<p_ecalIso<<" GeV | hcal = "<<std::setw(12)<<p_hcalIso<<" GeV"<<std::endl;
    std::cout<<"Isolation max values for this photon : | trk = "<<std::setw(12)<<p_maxTrkIso<<" GeV/c | ecal = "<<std::setw(12)<<p_maxEcalIso<<" GeV | hcal = "<<std::setw(12)<<p_maxHcalIso<<" GeV"<<std::endl;
    std::cout<<"         Would it have been Isolated?: | Iso? = "<<(p_trkIso < p_maxTrkIso && p_ecalIso < p_maxEcalIso  && p_hcalIso < p_maxHcalIso)<<std::endl;
    std::cout<<"        Does it has good showershape?: | Pix? = "<<p_hPS<<" | H/E = "<<p_HE<<" | p_sIEIE = "<<p_sIEIE<<" | p_sIPIP = "<<p_sIPIP<<std::endl; 
    */

    // N-1 Cuts without p_sIEIE
    // N - 1 Cut: H/E
    if (!p_hPS && p_trkIso < p_maxTrkIso && p_ecalIso < p_maxEcalIso && p_hcalIso < p_maxHcalIso) { p->isEB()? Fill(PB_HE, n_he, n1_he, n2_he, p_HE) : Fill(PE_HE, n_he, n1_he, n2_he, p_HE); }
    // N - 1 Cut: hPS
    if(p_HE < photon_he && p_trkIso < p_maxTrkIso && p_ecalIso < p_maxEcalIso && p_hcalIso < p_maxHcalIso) { p->isEB()? Fill(PB_PS, n_ps, n1_ps, n2_ps, p_hPS) : Fill(PE_PS, n_ps, n1_ps, n2_ps, p_hPS); }
    // N - 1 Cut: Trk Iso
    if(!p_hPS && p_HE < photon_he && p_ecalIso < p_maxEcalIso && p_hcalIso < p_maxHcalIso) { p->isEB()? Fill(PB_TRK, n_iso, n1_iso, n2_iso, p_trkIso) : Fill(PE_TRK, n_iso, n1_iso, n2_iso, p_trkIso); }
    // N - 1 Cut: ECAL Iso
    if(!p_hPS && p_HE < photon_he && p_trkIso < p_maxTrkIso && p_hcalIso < p_maxHcalIso) { p->isEB()? Fill(PB_ECAL, n_iso, n1_iso, n2_iso, p_ecalIso) : Fill(PE_ECAL, n_iso, n1_iso, n2_iso, p_ecalIso); }
    // N - 1 Cut: HCAL Iso
    if(!p_hPS && p_HE < photon_he && p_trkIso < p_maxTrkIso && p_ecalIso < p_maxEcalIso) { p->isEB()? Fill(PB_HCAL, n_iso, n1_iso, n2_iso, p_hcalIso) : Fill(PE_HCAL, n_iso, n1_iso, n2_iso, p_hcalIso); }
    // N Cuts:
    if(!p_hPS && p_HE < photon_he && p_trkIso < p_maxTrkIso && p_ecalIso < p_maxEcalIso && p_hcalIso < p_maxHcalIso) { 
      p->isEB()? Fill(PB_SIEIE, n_sb, n1_sb, n2_sb, p_sIEIE) : Fill(PE_SIEIE, n_se, n1_se, n2_se, p_sIEIE);
      p->isEB()? Fill(PB_SIPIP, n_sb, n1_sb, n2_sb, p_sIPIP) : Fill(PE_SIPIP, n_s2e, n1_s2e, n2_s2e, p_sIPIP);
      if(seedRecoFlag != -1) {
	p->isEB()? Fill(PB_eM_B, n_et, n1_et, n2_et, eMax)           : Fill(PE_eM_B, n_et, n1_et, n2_et, eMax);
	p->isEB()? Fill(PB_sT_B, n_sT, n1_sT, n2_sT, seedTime)       : Fill(PE_sT_B, n_sT, n1_sT, n2_sT, seedTime);
	p->isEB()? Fill(PB_RF_B, n_RF, n1_RF, n2_RF, seedRecoFlag)   : Fill(PE_RF_B, n_RF, n1_RF, n2_RF, seedRecoFlag);
	p->isEB()? Fill(PB_sS_B, n_sS, n1_sS, n2_sS, seedSeverity)   : Fill(PE_sS_B, n_sS, n1_sS, n2_sS, seedSeverity);
	p->isEB()? Fill(PB_SC_B, n_SC, n1_SC, n2_SC, seedSwissCross) : Fill(PE_SC_B, n_SC, n1_SC, n2_SC, seedSwissCross);
	p->isEB()? Fill(PB_EE_B, n_SC, n1_SC, n2_SC, e2overe9)       : Fill(PE_EE_B, n_SC, n1_SC, n2_SC, e2overe9);
      }
      Fill(P_Et_B, n_et, n1_et, n2_et, p_et);     p->isEB()? Fill(PB_Et_B, n_et, n1_et, n2_et, p_et)     : Fill(PE_Et_B, n_et, n1_et, n2_et, p_et);
      Fill(P_Eta_B, n_eta, n1_eta, n2_eta, p_eta); p->isEB()? Fill(PB_Eta_B, n_eta, n1_eta, n2_eta, p_eta) : Fill(PE_Eta_B, n_eta, n1_eta, n2_eta, p_eta);
    }
    // N + 1 Cuts: p_sIEIE 
    if(!p_hPS && p_HE < photon_he && p_trkIso < p_maxTrkIso && p_ecalIso < p_maxEcalIso && p_hcalIso < p_maxHcalIso && ((p->isEB() && p_sIEIE < sigma_barrel) || (p->isEE() && p_sIEIE < sigma_endcap))) {
      // spike removal with E2/E9 and seedSeverity
      // Recommended most powerful method, but only partially documented (Strategy 3)
      if((p->isEB() && e2overe9 < sc_spike && seedSeverity != 3) || p->isEE()) {
	Fill(P_Et_A, n_et, n1_et, n2_et, p_et);     p->isEB()? Fill(PB_Et_A, n_et, n1_et, n2_et, p_et)     : Fill(PE_Et_A, n_et, n1_et, n2_et, p_et);
	Fill(P_Eta_A, n_eta, n1_eta, n2_eta, p_eta); p->isEB()? Fill(PB_Eta_A, n_eta, n1_eta, n2_eta, p_eta) : Fill(PE_Eta_A, n_eta, n1_eta, n2_eta, p_eta);
	if(seedRecoFlag != -1) {
	  p->isEB()? Fill(PB_eM_A, n_et, n1_et, n2_et, eMax)           : Fill(PE_eM_A, n_et, n1_et, n2_et, eMax);
	  p->isEB()? Fill(PB_sT_A, n_sT, n1_sT, n2_sT, seedTime)       : Fill(PE_sT_A, n_sT, n1_sT, n2_sT, seedTime);
	  p->isEB()? Fill(PB_RF_A, n_RF, n1_RF, n2_RF, seedRecoFlag)   : Fill(PE_RF_A, n_RF, n1_RF, n2_RF, seedRecoFlag);
	  p->isEB()? Fill(PB_sS_A, n_sS, n1_sS, n2_sS, seedSeverity)   : Fill(PE_sS_A, n_sS, n1_sS, n2_sS, seedSeverity);
	  p->isEB()? Fill(PB_SC_A, n_SC, n1_SC, n2_SC, seedSwissCross) : Fill(PE_SC_A, n_SC, n1_SC, n2_SC, seedSwissCross);
	  p->isEB()? Fill(PB_EE_A, n_SC, n1_SC, n2_SC, e2overe9)       : Fill(PE_EE_A, n_SC, n1_SC, n2_SC, e2overe9);
	  // Print out of what is still left
	  if(seedRecoFlag != 0) {
	    std::cout<<"PAT Photon we work with: | eta = "<<std::setw(9)<<p->eta()<<" | et = "<<std::setw(9)<<p->et()<<" GeV | eMax = "<<std::setw(9)<<eMax;
	    std::cout<<" | seedRecoFlag = "<<std::setw(9)<<seedRecoFlag<<std::endl;
	  }
	}
      }
      // Specific Strategy 3 Plots
      if (p->isEB() && e2overe9 < sc_spike && seedSeverity != 3){
	  Fill(PB_SIEIE_Strat3, n_sb, n1_sb, n2_sb, p_sIEIE);
	  Fill(PB_SIPIP_Strat3, n_sb, n1_sb, n2_sb, p_sIPIP);
	  Fill(PB_sT_Strat3,    n_sT, n1_sT, n2_sT, seedTime); 
      }
      // spike removal with p_sIEIE, p_sIPIP and seedSeverity
      // AN-2010/268 method
      if(p->isEB() && p_sIEIE > sigma_spike && p_sIPIP > sigma_spike && seedSeverity != 3) {
	Fill(PB_EE_Strat1, n_SC, n1_SC, n2_SC, e2overe9);
	Fill(PB_sT_Strat1, n_sT, n1_sT, n2_sT, seedTime);
      }
      // spike removal with E2/E9 and seedTime
      // AN-2010/340
      if(p->isEB() && e2overe9 < sc_spike && seedTime > st_spike) {
	Fill(PB_SIEIE_Strat2, n_sb, n1_sb, n2_sb, p_sIEIE);
	Fill(PB_SIPIP_Strat2, n_sb, n1_sb, n2_sb, p_sIPIP);
	Fill(PB_sS_Strat2, n_sS, n1_sS, n2_sS, seedSeverity);
      }
    }
    // Selection 1 counter: [-] | HE | PS | SIEIE | ISO TRK | ISO ECAL | ISO HCAL
    if(p_HE < photon_he) { Selection1->Fill(2);
      if(!p_hPS) { Selection1->Fill(3); 
	if ((p->isEB() && p_sIEIE < sigma_barrel) || (p->isEE() && p_sIEIE < sigma_endcap)) { Selection1->Fill(4); 
	  if(p_trkIso < p_maxTrkIso) { Selection1->Fill(5);
	    if(p_ecalIso < p_maxEcalIso) { Selection1->Fill(6);
	      if(p_hcalIso < p_maxHcalIso) { Selection1->Fill(7);
		if((p->isEB() && e2overe9 < sc_spike && seedSeverity != 3) || p->isEE()) { Selection1->Fill(8);
		  if(p->et() > 100) { Selection1->Fill(9);
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // Selection 2 counter: [-] | HE | PS | ISO TRK | ISO ECAL | ISO HCAL | SIEIE
    if(p_HE < photon_he) { Selection2->Fill(2);
      if(!p_hPS) { Selection2->Fill(3); 
	if(p_trkIso < p_maxTrkIso) { Selection2->Fill(4);
	  if(p_ecalIso < p_maxEcalIso) { Selection2->Fill(5);
	    if(p_hcalIso < p_maxHcalIso) { Selection2->Fill(6);
	      if ((p->isEB() && p_sIEIE < sigma_barrel) || (p->isEE() && p_sIEIE < sigma_endcap)) { Selection2->Fill(7);
		if((p->isEB() && e2overe9 < sc_spike && seedSeverity != 3) || p->isEE()) { Selection2->Fill(8);
		  if(p->et() > 100) { Selection2->Fill(9);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } 
}

void 
PhotonPreSelectionPlotter::Fill(TH1F * Histo, int n, double n1, double n2, double value) {
  double binwidth = (n2-n1)/n;
  double n1_ = n1+binwidth/2;
  double n2_ = n2-binwidth/2;
  if(value > n1 && value < n2) Histo->Fill(value);
  else if (value < n1) Histo->Fill(n1_);
  else if (value > n2) Histo->Fill(n2_);
  else {}
}

// ------------ method called once each job just before starting event loop  ------------
void 
PhotonPreSelectionPlotter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonPreSelectionPlotter::endJob() {

  std::string S1BinLabels [] = {"Cleaned Events", "HE", "PS", "SIEIE", "ISO TRK", "ISO ECAL", "ISO HCAL", "No Spike", "#gamma > 100 GeV"};
  std::string S2BinLabels [] = {"Cleaned Events", "HE", "PS", "ISO TRK", "ISO ECAL", "ISO HCAL", "SIEIE", "No Spike", "#gamma > 100 GeV"};
  for(int i=0; i<9; ++i) { 
    Selection1->GetXaxis()->SetBinLabel(i+1, S1BinLabels[i].c_str());
    Selection2->GetXaxis()->SetBinLabel(i+1, S2BinLabels[i].c_str());
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonPreSelectionPlotter);
