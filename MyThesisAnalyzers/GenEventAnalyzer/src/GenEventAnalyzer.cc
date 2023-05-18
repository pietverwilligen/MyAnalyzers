// -*- C++ -*-
//
// Package:    GenEventAnalyzer
// Class:      GenEventAnalyzer
// 
/**\class GenEventAnalyzer GenEventAnalyzer.cc MyAnalyzers/GenEventAnalyzer/src/GenEventAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen,40 4-B15,+41227671568,
//         Created:  Sat Dec  4 14:27:21 CET 2010
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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <DataFormats/METReco/interface/MET.h>

#include <DataFormats/PatCandidates/interface/Jet.h>

// own include files
#include "MyAnalyzers/PhotonAnalyzer/interface/SHistContainer.h"

//
// class declaration
//

class GenEventAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenEventAnalyzer(const edm::ParameterSet&);
      ~GenEventAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  void BosonFill(int, const reco::GenParticle*);
  void JetsFill(int h, edm::Handle< std::vector<reco::GenJet> > genJets);
  void METFill(int h, reco::MET::LorentzVector mht, double ht);
  void Fill(TH1*, int, double, double, double);

      // ----------member data ---------------------------
  // WRITE
  TFile * outputfile;
  TDirectoryFile *Boson_0J_BA, *Boson_0J_AA, *Boson_0J_BA_100, *Boson_0J_AA_100;
  TDirectoryFile *Boson_1J_BA, *Boson_1J_AA, *Boson_1J_BA_100, *Boson_1J_AA_100;
  TDirectoryFile *Boson_2J_BA, *Boson_2J_AA, *Boson_2J_BA_100, *Boson_2J_AA_100;
  TDirectoryFile *Boson_3J_BA, *Boson_3J_AA, *Boson_3J_BA_100, *Boson_3J_AA_100;
  TDirectoryFile *ExtraPlots;
  std::vector<SHistContainer> histcontainers;
  TH1F *EventCounters_B000, *EventCounters_B100, *NJets_AK5PT50ETA25, *NJets_AK5PT30ETA50, *GenEvQScale, *GenEvPTHat, *GenEvQScale_AfterHTCut, *GenEvPTHat_AfterHTCut;
  TH1F *J1Counter, *J2Counter, *J3Counter;
  TH1F *Pt_All_BaS, *Pt_All_HHS, *Pt_All_HMS,*Pt_R05_BaS, *Pt_R05_HHS, *Pt_R05_HMS,*Pt_R09_BaS, *Pt_R09_HHS, *Pt_R09_HMS;
  // READ
  bool Debug, Zinv, Pat;
  double HTCut;
  std::string RootFileName, GenParticles, GenJets, RecoJets, PatJets;
  // INTERNAL USE
  const reco::GenParticle     *g, *b, *d1, *d2;
  const reco::GenJet          *gj, *gj1, *gj2, *gj3;
  std::vector<int> Count;
};

//
// constants, enums and typedefs
//
double pi = 3.1415926535;

double jet_ptcut = 50.0;
double jet_etacut = 2.5;

double photon_ptcut = 100.0;
double photon_etacut = 2.50;
double EBEE_bord = 1.479;
// kapoenen!
// isEB() and isEE() methods:
double EB_accept = 1.379;
double EE_accept = 1.579;
// if you do it yourself:

// double EB_accept = 1.4442;
// double EE_accept = 1.566;



int n_eta = 50;  double n1_eta = -5.0,  n2_eta  = 5.0;
int n_phi = 18;  double n1_phi = -3.14, n2_phi  = 3.14;
int n_et  = 100; double n1_et  = 0.0,   n2_et   = 500;
int n_ht  = 200; double n1_ht  = 0.0,   n2_ht   = 1000;
int n_m   = 300; double n1_m   = 0.0,   n2_m = 150.0;
int n_n   = 6;   double n1_n   = -0.5,  n2_n = 5.5;

std::string ContainNameArray [] = {"Boson_0J_BA", "Boson_0J_AA", "Boson_0J_BA_100", "Boson_0J_AA_100", "Boson_1J_BA", "Boson_1J_AA", "Boson_1J_BA_100", "Boson_1J_AA_100", "Boson_2J_BA", "Boson_2J_AA", "Boson_2J_BA_100", "Boson_2J_AA_100", "Boson_3J_BA", "Boson_3J_AA", "Boson_3J_BA_100", "Boson_3J_AA_100"};

//
// static data member definitions
//

//
// constructors and destructor
//
GenEventAnalyzer::GenEventAnalyzer(const edm::ParameterSet& iConfig)

{
  
  HTCut = 0.0;
  //now do what ever initialization is needed
  Debug         = iConfig.getParameter<bool>("Debug");
  Zinv          = iConfig.getParameter<bool>("Zinv");
  Pat           = iConfig.getParameter<bool>("Pat");
  HTCut         = iConfig.getParameter<double>("HTCut");
  RootFileName  = iConfig.getParameter<std::string>("RootFileName");
  GenParticles  = iConfig.getParameter<std::string>("GenParticles");
  GenJets       = iConfig.getParameter<std::string>("GenJets");
  RecoJets      = iConfig.getParameter<std::string>("RecoJets");
  PatJets       = iConfig.getParameter<std::string>("PatJets");

  outputfile = new TFile(RootFileName.c_str(), "RECREATE" );

  Boson_0J_BA     = (TDirectoryFile*) outputfile->mkdir("Boson_0J_BA",    "Boson_0J_BA");
  Boson_0J_AA     = (TDirectoryFile*) outputfile->mkdir("Boson_0J_AA",    "Boson_0J_AA");
  Boson_0J_BA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_0J_BA_100","Boson_0J_BA_100");
  Boson_0J_AA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_0J_AA_100","Boson_0J_AA_100");
  Boson_1J_BA     = (TDirectoryFile*) outputfile->mkdir("Boson_1J_BA",    "Boson_1J_BA");
  Boson_1J_AA     = (TDirectoryFile*) outputfile->mkdir("Boson_1J_AA",    "Boson_1J_AA");
  Boson_1J_BA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_1J_BA_100","Boson_1J_BA_100");
  Boson_1J_AA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_1J_AA_100","Boson_1J_AA_100");
  Boson_2J_BA     = (TDirectoryFile*) outputfile->mkdir("Boson_2J_BA",    "Boson_2J_BA");
  Boson_2J_AA     = (TDirectoryFile*) outputfile->mkdir("Boson_2J_AA",    "Boson_2J_AA");
  Boson_2J_BA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_2J_BA_100","Boson_2J_BA_100");
  Boson_2J_AA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_2J_AA_100","Boson_2J_AA_100");
  Boson_3J_BA     = (TDirectoryFile*) outputfile->mkdir("Boson_3J_BA",    "Boson_3J_BA");
  Boson_3J_AA     = (TDirectoryFile*) outputfile->mkdir("Boson_3J_AA",    "Boson_3J_AA");
  Boson_3J_BA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_3J_BA_100","Boson_3J_BA_100");
  Boson_3J_AA_100 = (TDirectoryFile*) outputfile->mkdir("Boson_3J_AA_100","Boson_3J_AA_100");

  ExtraPlots = (TDirectoryFile*) outputfile->mkdir("ExtraPlots","ExtraPlots");

  // # Jets:                          |   0  |   1   |   2   |   3   |
  // -----------------------------------------------------------------
  // Before Acceptance                |   0  |   4   |   8   |  12   |
  // After  Acceptance                |   1  |   5   |   9   |  13   |
  // Before Acceptance & > 100 GeV/c  |   2  |   6   |  10   |  14   |
  // After  Acceptance & > 100 GeV/c  |   3  |   7   |  11   |  15   |

  SHistContainer Boson_0J_BA_c, Boson_0J_AA_c, Boson_0J_BA_100_c, Boson_0J_AA_100_c, Boson_1J_BA_c, Boson_1J_AA_c, Boson_1J_BA_100_c, Boson_1J_AA_100_c, Boson_2J_BA_c, Boson_2J_AA_c, Boson_2J_BA_100_c, Boson_2J_AA_100_c, Boson_3J_BA_c, Boson_3J_AA_c, Boson_3J_BA_100_c, Boson_3J_AA_100_c;
  SHistContainer ContainerArray [] = {Boson_0J_BA_c, Boson_0J_AA_c, Boson_0J_BA_100_c, Boson_0J_AA_100_c, Boson_1J_BA_c, Boson_1J_AA_c, Boson_1J_BA_100_c, Boson_1J_AA_100_c, Boson_2J_BA_c, Boson_2J_AA_c, Boson_2J_BA_100_c, Boson_2J_AA_100_c, Boson_3J_BA_c, Boson_3J_AA_c, Boson_3J_BA_100_c, Boson_3J_AA_100_c};
  for(int i=0; i<16; ++i) { histcontainers.push_back(ContainerArray[i]); }
  
  for(unsigned int h = 0; h<histcontainers.size(); ++h) {
    histcontainers[h].book(ContainNameArray[h], new TH1F("Eta", "Boson Eta",  n_eta, n1_eta, n2_eta));
    histcontainers[h].book(ContainNameArray[h], new TH1F("Phi", "Boson Phi",  n_phi, n1_phi, n2_phi));
    histcontainers[h].book(ContainNameArray[h], new TH1F("Et",  "Boson Et",   n_et,  n1_et,  n2_et));
    histcontainers[h].book(ContainNameArray[h], new TH1F("Pt",  "Boson Pt",   n_et,  n1_et,  n2_et));
    histcontainers[h].book(ContainNameArray[h], new TH1F("Mt",  "Boson Mt",   n_m,   n1_m,   n2_m));
    histcontainers[h].book(ContainNameArray[h], new TH1F("Ms",  "Boson Mass", n_m,   n1_m,   n2_m));
    histcontainers[h].book(ContainNameArray[h], new TH1F("NJ",  "N Jets",     n_n,   n1_n,   n2_n));
    histcontainers[h].book(ContainNameArray[h], new TH1F("HMt", "Hadronic Mt",   n_et,   n1_et,   n2_et));
    histcontainers[h].book(ContainNameArray[h], new TH1F("HMtCheck", "Hadronic Mt By Hand: M^{2} = E^{2}-p_{x}^{2}-p_{y}^{2}", n_et,   n1_et,   n2_et));
    histcontainers[h].book(ContainNameArray[h], new TH1F("HMs", "Hadronic Mass", n_et,   n1_et,   n2_et));
    histcontainers[h].book(ContainNameArray[h], new TH1F("HMsCheck", "Hadronic Mass By Hand: M^{2} = E^{2}-p_{x}^{2}-p_{y}^{2}-p_{z}^{2}", n_et,   n1_et,   n2_et));
    histcontainers[h].book(ContainNameArray[h], new TH1F("HT",  "HT: Transverse Momentum",          n_ht,   n1_ht,   n2_ht));
    histcontainers[h].book(ContainNameArray[h], new TH1F("MHT", "MHT: Missing Transverse Momentum", n_et,   n1_et,   n2_et));

    histcontainers[h].book(ContainNameArray[h], new TH1F("J1_Eta", "Jet 1 Eta",  n_eta, n1_eta, n2_eta));
    histcontainers[h].book(ContainNameArray[h], new TH1F("J1_Phi", "Jet 1 Phi",  n_phi, n1_phi, n2_phi));
    histcontainers[h].book(ContainNameArray[h], new TH1F("J1_Pt",  "Jet 1 Pt",   n_et,  n1_et,  n2_et));

    histcontainers[h].book(ContainNameArray[h], new TH1F("J2_Eta", "Jet 2 Eta",  n_eta, n1_eta, n2_eta));
    histcontainers[h].book(ContainNameArray[h], new TH1F("J2_Phi", "Jet 2 Phi",  n_phi, n1_phi, n2_phi));
    histcontainers[h].book(ContainNameArray[h], new TH1F("J2_Pt",  "Jet 2 Pt",   n_et,  n1_et,  n2_et));

    histcontainers[h].book(ContainNameArray[h], new TH1F("J3_Eta", "Jet 3 Eta",  n_eta, n1_eta, n2_eta));
    histcontainers[h].book(ContainNameArray[h], new TH1F("J3_Phi", "Jet 3 Phi",  n_phi, n1_phi, n2_phi));
    histcontainers[h].book(ContainNameArray[h], new TH1F("J3_Pt",  "Jet 3 Pt",   n_et,  n1_et,  n2_et));

  }
  EventCounters_B000 = new TH1F("EventCounters_B000", "All Events | B + 0J | B + 1J | B + 2J | B + 3J | HT | Ang | MHT | HHS | HMS", 10,0.5,10.5);
  EventCounters_B100 = new TH1F("EventCounters_B100", "All Events | B100 + 0J | B100 + 1J | B100 + 2J | B100 + 3J | HT | Ang | MHT | HHS | HMS", 10,0.5,10.5);
  NJets_AK5PT50ETA25 = new TH1F("NJets_AK5PT50ETA25", "Number Of AK5 Jets with PT > 50 and ETA < 2.5", n_n, n1_n, n2_n);
  NJets_AK5PT30ETA50 = new TH1F("NJets_AK5PT30ETA50", "Number Of AK5 Jets with PT > 30 and ETA < 5.0", n_n, n1_n, n2_n);

  GenEvQScale = new TH1F("GenEvQScale", "qScale of the Event (related to HTBinning at generator level)", n_ht, n1_ht, n2_ht);
  GenEvPTHat  = new TH1F("GenEvPTHat",  "PTHat of the Event (related to HTBinning at generator level)", n_ht, n1_ht, n2_ht);
  GenEvQScale_AfterHTCut = new TH1F("GenEvQScale_AfterHTCut", "qScale of the Event after HT cut (related to HTBinning at generator level)", n_ht, n1_ht, n2_ht);
  GenEvPTHat_AfterHTCut  = new TH1F("GenEvPTHat_AfterHTCut",  "PTHat of the Event after HT cut (related to HTBinning at generator level)", n_ht, n1_ht, n2_ht);

  J1Counter = new TH1F("J1Counter", "All Events | B + 0J | B + 1J | HT | Ang | MHT | HHS | HMS", 8,0.5,8.5);
  J2Counter = new TH1F("J2counter", "All Events | B + 0J | B + 2J | HT | Ang | MHT | HHS | HMS", 8,0.5,8.5);
  J3Counter = new TH1F("J3Counter", "All Events | B + 0J | B + 3J | HT | Ang | MHT | HHS | HMS", 8,0.5,8.5);

  Pt_All_BaS = new TH1F("Pt_All_BaS", "All Photon Pt Distribution for Baseline Selection", n_et, n1_et, n2_et);
  Pt_All_HHS = new TH1F("Pt_All_HHS", "All Photon Pt Distribution for High HT  Selection", n_et, n1_et, n2_et);
  Pt_All_HMS = new TH1F("Pt_All_HMS", "All Photon Pt Distribution for High MHT Selection", n_et, n1_et, n2_et);

  Pt_R05_BaS = new TH1F("Pt_R05_BaS", "R05 Photon Pt Distribution for Baseline Selection", n_et, n1_et, n2_et);
  Pt_R05_HHS = new TH1F("Pt_R05_HHS", "R05 Photon Pt Distribution for High HT  Selection", n_et, n1_et, n2_et);
  Pt_R05_HMS = new TH1F("Pt_R05_HMS", "R05 Photon Pt Distribution for High MHT Selection", n_et, n1_et, n2_et);

  Pt_R09_BaS = new TH1F("Pt_R09_BaS", "R09 Photon Pt Distribution for Baseline Selection", n_et, n1_et, n2_et);
  Pt_R09_HHS = new TH1F("Pt_R09_HHS", "R09 Photon Pt Distribution for High HT  Selection", n_et, n1_et, n2_et);
  Pt_R09_HMS = new TH1F("Pt_R09_HMS", "R09 Photon Pt Distribution for High MHT Selection", n_et, n1_et, n2_et);


  for(int i=0; i<20; ++i) {Count.push_back(0);}
  std::string binlabel [] = {"All Events", "B + 0J", "B + 1J", "B + 2J", "B + 3J", "HT", "Ang", "MHT", "HHS", "HMS"};
  for(int i=1; i<=10; ++i) {EventCounters_B000->GetXaxis()->SetBinLabel(i,(binlabel[i-1]).c_str()); EventCounters_B100->GetXaxis()->SetBinLabel(i,(binlabel[i-1]).c_str());}
}


GenEventAnalyzer::~GenEventAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();
  for(unsigned int h = 0; h<histcontainers.size(); ++h) {histcontainers[h].write(outputfile);}
  outputfile->cd();
  EventCounters_B000->Write();
  EventCounters_B100->Write();
  NJets_AK5PT50ETA25->Write();
  NJets_AK5PT30ETA50->Write();
  GenEvQScale->Write();
  GenEvPTHat->Write();
  GenEvQScale_AfterHTCut->Write();
  GenEvPTHat_AfterHTCut->Write();
  J1Counter->Write();
  J2Counter->Write();
  J3Counter->Write();
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
  // std::cout<<"Total Amount of Z/Y + 3 Jets = "<<Count3JetEvents<<std::endl;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GenEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // EVENT INFORMATION
  int evNum = (iEvent.id()).event();
  int rnNum = (iEvent.id()).run();

  // GEN INFORMATION
  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByLabel("generator",genEventInfo);
  double evWeight = genEventInfo->weight();
  double qScale   = genEventInfo->qScale();
  double pthat    = (genEventInfo->hasBinningValues() ? (genEventInfo->binningValues())[0] : 0.0);
  // std::cout<<"evWeight = "<<evWeight<<" qScale = "<<qScale<<" pthat = "<<pthat<<std::endl;

  GenEvQScale->Fill(qScale);
  GenEvPTHat->Fill(pthat);

  if(qScale < HTCut) return;

  GenEvQScale_AfterHTCut->Fill(qScale);
  GenEvPTHat_AfterHTCut->Fill(pthat);

  // GEN PARTICLES
  edm::Handle<reco::GenParticleCollection>      genParticles;
  iEvent.getByLabel(GenParticles, genParticles);
  // GEN JETS
  edm::Handle< std::vector<reco::GenJet> >      genJets;
  iEvent.getByLabel(GenJets, genJets);
  
  // RESET
  b = 0;
  
  EventCounters_B000->Fill(1); // Analyzed Events
  EventCounters_B100->Fill(1);
  
  J1Counter->Fill(1);
  J2Counter->Fill(1);
  J3Counter->Fill(1);
  
  
  // Find the bosons you need in the GEN collection
  // ----------------------------------------------
  if(Zinv==0) {
    // Find the ME Photon and the GenPhoton linked to it
    bool MEPhotonFound = 0;
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      if(MEPhotonFound) continue;
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if (g->pdgId() == 22) {
	if(Debug) std::cout<<"ME Photon: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = "<<std::setw(9)<<g->phi();
	if(Debug) std::cout<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	reco::GenParticle::daughters d = g->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
	  if((*it_d)->status()==1 && (*it_d)->pdgId()==22) {
	    b = &(*(*it_d));
	    if(Debug) std::cout<<"GenPhoton: id = "<<std::setw(2)<<b->pdgId()<<" | st = "<<std::setw(2)<<b->status()<<" | eta = "<<std::setw(9)<<b->eta()<<" | phi = "<<std::setw(9)<<b->phi();
	    if(Debug) std::cout<<" | pt = "<<std::setw(9)<<b->pt()<<" GeV/c | et = "<<std::setw(9)<<b->et()<<" GeV | Energy = "<<b->energy()<<" GeV ";
	    if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<b->mt()<<" GeV/c^2"<<std::endl;
	  }
	}
      }
    }
  }
  else if (Zinv==1) {
    // Find the ME Zboson
    bool MEZbosonFound = 0;
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      if(MEZbosonFound) continue;
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if (g->pdgId() == 23) {
	if(Debug) std::cout<<"ME Zboson: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = "<<std::setw(9)<<g->phi();
	if(Debug) std::cout<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	b = &(*g);
      }
    }
  }
  else if (Zinv==2) { 
    // Drell Yan: Find the ME Y* Boson
    // And Scale to the invisible Branching Ratios
    // We will ignore the mix with Photon, since it only affects the left tail (left side of the resonance)
    // while we are only interested in the right tail
    
    // 2B Implemented
    // std::cout<<"Zinv > 1 is not implemented"<<std::endl;
    
    // Find the ME Zboson
    bool MEZbosonFound = 0;
    for(unsigned int i=0; i<genParticles->size(); ++i) {
      if(MEZbosonFound) continue;
      g = &((*genParticles)[i]);
      if (g->status() != 3) continue;
      if (g->pdgId() == 23) {
	if(Debug) std::cout<<"ME Zboson: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = "<<std::setw(9)<<g->phi();
	if(Debug) std::cout<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	reco::GenParticle::daughters d = g->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
	  std::cout<<"     "<<"ME Zboson Daughter: id = "<<std::setw(2)<<(*it_d)->pdgId()<<" | st = "<<std::setw(2)<<(*it_d)->status()<<" | eta = "<<std::setw(9)<<(*it_d)->eta();
	  std::cout<<" | phi = "<<std::setw(9)<<(*it_d)->phi()<<" | pt = "<<std::setw(9)<<(*it_d)->pt()<<" GeV/c"<<std::endl;
	}
	b = &(*g);
      }
    }
    // Tag two jets that have Z decay products clustered in
    /* for(unsigned int i=0; i<genJets->size(); ++i) {
       gj = &((*genJets)[i]);
       std::vector<const GenParticle*> jetconstituents = gj->getGenConstituents();
       for(unsigned int j=0; j<jetconstituents.size(); ++j) {
       g = &((*jetconstituents)[i]);
       }
       }*/
  }
  else {std::cout<<"Zinv > 2 is not implemented"<<std::endl;}
  
  
  
  // If boson found, fill histograms
  // -------------------------------
  
  if(b != 0) {
    // Require 3 GenJets with pt > 50 and |eta| < 2.5
    int goodGenJets_AK5PT50ETA25 = 0, goodGenJets_AK5PT30ETA50 = 0;
    double ht = 0;
    reco::MET::LorentzVector mht(0,0,0,0);
    double min_phi_pi = 0.0, min_r_pi = 0.0;
    for(unsigned int i=0; i<genJets->size(); ++i) {
      gj = &((*genJets)[i]);
      if (gj->pt() > jet_ptcut && fabs(gj->eta()) < jet_etacut) {++goodGenJets_AK5PT50ETA25; ht  += gj->pt(); }
      if (gj->pt() > 30)                                        {
	++goodGenJets_AK5PT30ETA50; mht -= gj->p4();
	double phi_pi = reco::deltaPhi(b->phi(),gj->phi());
	double r_pi = reco::deltaR(b->phi(), b->eta(), gj->phi(), gj->eta());
	if(i==0) {min_phi_pi = phi_pi; min_r_pi = r_pi;}
	else if (min_phi_pi > phi_pi) {min_phi_pi = phi_pi; min_r_pi = r_pi;}
	else {}
      }
    }
    Fill(NJets_AK5PT50ETA25, n_n, n1_n, n2_n, goodGenJets_AK5PT50ETA25);
    Fill(NJets_AK5PT30ETA50, n_n, n1_n, n2_n, goodGenJets_AK5PT30ETA50);
    int goodGenJets = goodGenJets_AK5PT50ETA25;
    double b_jcl = fabs(min_phi_pi);
    double b_rjc = fabs(min_r_pi);

    // Before Acceptance
    BosonFill(0, b);  JetsFill(0, genJets); METFill(0, mht, ht); EventCounters_B000->Fill(2); J1Counter->Fill(2); J2Counter->Fill(2); J3Counter->Fill(2);
    if(goodGenJets > 0) { BosonFill(4, b);  JetsFill(4, genJets);  METFill(4, mht, ht); EventCounters_B000->Fill(3); }
    if(goodGenJets > 1) { BosonFill(8, b);  JetsFill(8, genJets);  METFill(8, mht, ht); EventCounters_B000->Fill(4); }
    if(goodGenJets > 2) { BosonFill(12, b); JetsFill(12, genJets); METFill(12, mht, ht); EventCounters_B000->Fill(5); }
    // After Acceptance
    if((fabs(b->eta()) < EB_accept) || (fabs(b->eta()) > EE_accept && fabs(b->eta()) < photon_etacut )) {
      BosonFill(1, b); JetsFill(1, genJets); METFill(1, mht, ht);
      if(goodGenJets > 0) { BosonFill(5, b);  JetsFill(5, genJets);  METFill(5, mht, ht);}
      if(goodGenJets > 1) { BosonFill(9, b);  JetsFill(9, genJets);  METFill(9, mht, ht);}
      if(goodGenJets > 2) { BosonFill(13, b); JetsFill(13, genJets); METFill(13, mht, ht);}
    }
    // Before Acceptance && Boson Pt > 100
    if(b->pt() > 100) {
      BosonFill(2, b); JetsFill(2, genJets); METFill(2, mht, ht);    EventCounters_B100->Fill(2);
      if(goodGenJets > 0) { BosonFill(6, b);  JetsFill(6, genJets);  METFill(6, mht, ht); EventCounters_B100->Fill(3); }
      if(goodGenJets > 1) { BosonFill(10, b); JetsFill(10, genJets); METFill(10, mht, ht); EventCounters_B100->Fill(4); } 
      if(goodGenJets > 2) { BosonFill(14, b); JetsFill(14, genJets); METFill(14, mht, ht); EventCounters_B100->Fill(5); } 
    }
    // After Acceptance  && Boson Pt > 100
    if(((fabs(b->eta()) < EB_accept) || (fabs(b->eta()) > EE_accept && fabs(b->eta()) < photon_etacut )) && b->pt() > 100) {
      BosonFill(3, b); JetsFill(3, genJets); METFill(3, mht, ht);
      if(goodGenJets > 0) { BosonFill(7, b);  JetsFill(7, genJets);  METFill(7, mht, ht);}
      if(goodGenJets > 1) { BosonFill(11, b); JetsFill(11, genJets); METFill(11, mht, ht);}
      if(goodGenJets > 2) { BosonFill(15, b); JetsFill(15, genJets); METFill(15, mht, ht);}
    }

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
	EventCounters_B000->Fill(6); if(b->pt() > 100) EventCounters_B100->Fill(6);
	if ( mht_j1 > 0.5 && mht_j2 > 0.5 && mht_j3 > 0.3) {
	  EventCounters_B000->Fill(7); if(b->pt() > 100) EventCounters_B100->Fill(7);
	  if(mht_value > 150) {
	    EventCounters_B000->Fill(8); if(b->pt() > 100) EventCounters_B100->Fill(8);
	    Fill(Pt_All_BaS, n_et, n1_et, n2_et, b->pt());
	    if(b_rjc > 0.5) Fill(Pt_R05_BaS, n_et, n1_et, n2_et, b->pt());
	    if(b_rjc > 0.9) Fill(Pt_R09_BaS, n_et, n1_et, n2_et, b->pt());
	    if(ht > 500) {
	      EventCounters_B000->Fill(9); if(b->pt() > 100) EventCounters_B100->Fill(9);
	      Fill(Pt_All_HHS, n_et, n1_et, n2_et, b->pt());
	      if(b_rjc > 0.5) Fill(Pt_R05_HHS, n_et, n1_et, n2_et, b->pt());
	      if(b_rjc > 0.9) Fill(Pt_R09_HHS, n_et, n1_et, n2_et, b->pt());
	    }
	    if(mht_value > 250) {
	      EventCounters_B000->Fill(10); if(b->pt() > 100) EventCounters_B100->Fill(10);
	      Fill(Pt_R05_HMS, n_et, n1_et, n2_et, b->pt());
	      if(b_rjc > 0.5) Fill(Pt_All_HMS, n_et, n1_et, n2_et, b->pt());
	      if(b_rjc > 0.9) Fill(Pt_R09_HMS, n_et, n1_et, n2_et, b->pt());
	    }
	  }
	}
      }
    }

    // Additional:
    if(goodGenJets_AK5PT50ETA25 > 2) {
      gj1 = &((*genJets)[0]); gj2 = &((*genJets)[1]); gj3 = &((*genJets)[2]);
      mht_j1 = fabs(reco::deltaPhi(gj1->phi(),MHT.phi()));
      mht_j2 = fabs(reco::deltaPhi(gj2->phi(),MHT.phi()));
      mht_j3 = fabs(reco::deltaPhi(gj3->phi(),MHT.phi()));
    }
    else if(goodGenJets_AK5PT50ETA25 > 1) {
      gj1 = &((*genJets)[0]); gj2 = &((*genJets)[1]); 
      mht_j1 = fabs(reco::deltaPhi(gj1->phi(),MHT.phi()));
      mht_j2 = fabs(reco::deltaPhi(gj2->phi(),MHT.phi()));
    }
    else if(goodGenJets_AK5PT50ETA25 > 0) {
      gj1 = &((*genJets)[0]); 
      mht_j1 = fabs(reco::deltaPhi(gj1->phi(),MHT.phi()));
    }

    if(goodGenJets > 2) {
      J3Counter->Fill(3);
      if(ht > 300) {
        J3Counter->Fill(4); 
        if ( mht_j1 > 0.5 && mht_j2 > 0.5 && mht_j3 > 0.3) {
          J3Counter->Fill(5);
          if(mht_value > 150) {
            J3Counter->Fill(6);
            if(ht > 500) {
              J3Counter->Fill(7);
            }
            if(mht_value > 250) {
              J3Counter->Fill(8); 
            }
          }
        }
      }
    }

    if(goodGenJets > 1) {
      J2Counter->Fill(3);
      if(ht > 300) {
        J2Counter->Fill(4);
        if ( mht_j1 > 0.5 && mht_j2 > 0.5) {
          J2Counter->Fill(5);
          if(mht_value > 150) {
            J2Counter->Fill(6);
            if(ht > 500) {
              J2Counter->Fill(7);
            }
            if(mht_value > 250) {
              J2Counter->Fill(8);
            }
          }
        }
      }
    }

    if(goodGenJets > 0) {
      J1Counter->Fill(3);
      if(ht > 300) {
        J1Counter->Fill(4);
        if ( mht_j1 > 0.5 ) {
          J1Counter->Fill(5);
          if(mht_value > 150) {
            J1Counter->Fill(6);
            if(ht > 500) {
              J1Counter->Fill(7);
            }
            if(mht_value > 250) {
              J1Counter->Fill(8);
            }
          }
        }
      }
    }

  }
  else {
    std::cout<<"No ME Boson found!"<<std::endl;
  }

}

void
GenEventAnalyzer::BosonFill(int h, const reco::GenParticle *g) {
      Fill(histcontainers[h].get(ContainNameArray[h],"Eta"), n_eta, n1_eta, n2_eta, g->eta());
      Fill(histcontainers[h].get(ContainNameArray[h],"Phi"), n_phi, n1_phi, n2_phi, g->phi());
      Fill(histcontainers[h].get(ContainNameArray[h],"Et"),  n_et,  n1_et,  n2_et,  g->et());
      Fill(histcontainers[h].get(ContainNameArray[h],"Pt"),  n_et,  n1_et,  n2_et,  g->pt());
      Fill(histcontainers[h].get(ContainNameArray[h],"Mt"),  n_m,   n1_m,   n2_m,   g->mt());
      Fill(histcontainers[h].get(ContainNameArray[h],"Ms"),  n_m,   n1_m,   n2_m,   g->mass());
}


void
GenEventAnalyzer::JetsFill(int h, edm::Handle< std::vector<reco::GenJet> > genJets) {
  int counter = 0;
  for(unsigned int i=0; i<genJets->size(); ++i) {
    gj = &((*genJets)[i]);
    if (gj->pt() > jet_ptcut && fabs(gj->eta()) < jet_etacut) {
      ++counter;
      if(counter==1) {
	Fill(histcontainers[h].get(ContainNameArray[h],"J1_Eta"), n_eta, n1_eta, n2_eta, gj->eta());
	Fill(histcontainers[h].get(ContainNameArray[h],"J1_Phi"), n_phi, n1_phi, n2_phi, gj->phi());
	Fill(histcontainers[h].get(ContainNameArray[h],"J1_Pt"),  n_et,  n1_et,  n2_et, gj->pt());
      }
      if(counter==2) {
	Fill(histcontainers[h].get(ContainNameArray[h],"J2_Eta"), n_eta, n1_eta, n2_eta, gj->eta());
	Fill(histcontainers[h].get(ContainNameArray[h],"J2_Phi"), n_phi, n1_phi, n2_phi, gj->phi());
	Fill(histcontainers[h].get(ContainNameArray[h],"J2_Pt"),  n_et,  n1_et,  n2_et, gj->pt());
      }
      if(counter==3) {
	Fill(histcontainers[h].get(ContainNameArray[h],"J3_Eta"), n_eta, n1_eta, n2_eta, gj->eta());
	Fill(histcontainers[h].get(ContainNameArray[h],"J3_Phi"), n_phi, n1_phi, n2_phi, gj->phi());
	Fill(histcontainers[h].get(ContainNameArray[h],"J3_Pt"),  n_et,  n1_et,  n2_et, gj->pt());
      }
    }
  }
  Fill(histcontainers[h].get(ContainNameArray[h],"NJ"),  n_n,   n1_n,   n2_n, counter);
}

void 
GenEventAnalyzer::METFill(int h, reco::MET::LorentzVector mht, double ht) {
  // Hadronic Mass by Hand:
  double mt   = sqrt(pow(mht.energy(),2)-pow(mht.px(),2)-pow(mht.py(),2));
  double mass = sqrt(pow(mht.energy(),2)-pow(mht.px(),2)-pow(mht.py(),2)-pow(mht.pz(),2));

  Fill(histcontainers[h].get(ContainNameArray[h],"HMt"),  n_et,   n1_et,   n2_et,   mht.mt());
  Fill(histcontainers[h].get(ContainNameArray[h],"HMtCheck"),  n_et,   n1_et,   n2_et,   mt);
  Fill(histcontainers[h].get(ContainNameArray[h],"HMs"),  n_et,   n1_et,   n2_et,   mht.mass());
  Fill(histcontainers[h].get(ContainNameArray[h],"HMsCheck"),  n_et,   n1_et,   n2_et,   mass);
  Fill(histcontainers[h].get(ContainNameArray[h],"HT"),   n_ht,  n1_ht,  n2_ht,  ht);
  Fill(histcontainers[h].get(ContainNameArray[h],"MHT"),  n_et,  n1_et,  n2_et,  mht.pt());
}

void
GenEventAnalyzer::Fill(TH1 * Histo, int n, double n1, double n2, double value) {
  if(Debug) std::cout<<"Inside Filling Method"<<std::endl;
  double binwidth = (n2-n1)/n;
  double n1_ = n1+binwidth/2;
  double n2_ = n2-binwidth/2;
  if(value > n1 && value < n2) Histo->Fill(value);
  else if (value <= n1) Histo->Fill(n1_);
  else if (value >= n2) Histo->Fill(n2_);
  else {}
}




// ------------ method called once each job just before starting event loop  ------------
void 
GenEventAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenEventAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenEventAnalyzer);
