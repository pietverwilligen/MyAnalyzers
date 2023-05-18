// -*- C++ -*-
//
// Package:    Acceptance
// Class:      Acceptance
// 
/**\class Acceptance Acceptance.cc MyAnalyzers/Acceptance/src/Acceptance.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  local user
//         Created:  Sun Jan 23 11:42:31 CET 2011
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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include <DataFormats/EgammaCandidates/interface/Photon.h>

//
// class declaration
//

class Acceptance : public edm::EDAnalyzer {
   public:
      explicit Acceptance(const edm::ParameterSet&);
      ~Acceptance();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  void Fill(TH1*, int, double, double, double);
      // ----------member data ---------------------------
  TFile * outputfile;
  TH1F *Eta_All_Gen, *Phi_All_Gen, *Pt_All_Gen, *Eta_All_Reco, *Phi_All_Reco, *Pt_All_Reco, *Pt_Res_All;
  TH1F *Eta_3Jets_Gen, *Phi_3Jets_Gen, *Pt_3Jets_Gen, *Eta_3Jets_Reco, *Phi_3Jets_Reco, *Pt_3Jets_Reco, *Pt_Res_3Jets;
  std::string RootFileName;
  const reco::GenParticle     *g, *b;
  const reco::GenJet          *gj;
  const reco::Photon          *p;
};

//
// constants, enums and typedefs
//

int n_eta = 100; double n1_eta = -5.0,  n2_eta = 5.0;
int n_phi = 72;  double n1_phi = -3.14, n2_phi = 3.14;
int n_et  = 100; double n1_et  = 0.0,   n2_et  = 500;
int n_res = 200; double n1_res = -1.0,  n2_res = 1.0;

//
// static data member definitions
//

//
// constructors and destructor
//
Acceptance::Acceptance(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  RootFileName  = iConfig.getParameter<std::string>("RootFileName");
  outputfile = new TFile(RootFileName.c_str(), "RECREATE" );

  Eta_All_Gen  = new TH1F("Eta_All_Gen",  "Boson Eta Gen Photon",  n_eta, n1_eta, n2_eta);
  Phi_All_Gen  = new TH1F("Phi_All_Gen",  "Boson Phi Gen Photon",  n_phi, n1_phi, n2_phi);
  Pt_All_Gen   = new TH1F("Pt_All_Gen",   "Boson Pt Gen Photon",   n_et,  n1_et,  n2_et);
  Eta_All_Reco = new TH1F("Eta_All_Reco", "Boson Eta Reco Photon", n_eta, n1_eta, n2_eta);
  Phi_All_Reco = new TH1F("Phi_All_Reco", "Boson Phi Reco Photon", n_phi, n1_phi, n2_phi);
  Pt_All_Reco  = new TH1F("Pt_All_Reco",  "Boson Pt Reco Photon",  n_et,  n1_et,  n2_et);
  Pt_Res_All   = new TH1F("Pt_Res_All",   "Boson Pt Resolution",   n_res,  n1_res,  n2_res);

  Eta_3Jets_Gen  = new TH1F("Eta_3Jets_Gen",  "Boson Eta Gen Photon",  n_eta, n1_eta, n2_eta);
  Phi_3Jets_Gen  = new TH1F("Phi_3Jets_Gen",  "Boson Phi Gen Photon",  n_phi, n1_phi, n2_phi);
  Pt_3Jets_Gen   = new TH1F("Pt_3Jets_Gen",   "Boson Pt Gen Photon",   n_et,  n1_et,  n2_et);
  Eta_3Jets_Reco = new TH1F("Eta_3Jets_Reco", "Boson Eta Reco Photon", n_eta, n1_eta, n2_eta);
  Phi_3Jets_Reco = new TH1F("Phi_3Jets_Reco", "Boson Phi Reco Photon", n_phi, n1_phi, n2_phi);
  Pt_3Jets_Reco  = new TH1F("Pt_3Jets_Reco",  "Boson Pt Reco Photon",  n_et,  n1_et,  n2_et);
  Pt_Res_3Jets   = new TH1F("Pt_Res_3Jets",   "Boson Pt Resolution",   n_res,  n1_res,  n2_res);

}


Acceptance::~Acceptance()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  Eta_All_Gen->Write();
  Phi_All_Gen->Write();
  Pt_All_Gen->Write();
  Eta_All_Reco->Write();
  Phi_All_Reco->Write();
  Pt_All_Reco->Write();
  Pt_Res_All->Write();

  Eta_3Jets_Gen->Write();
  Phi_3Jets_Gen->Write();
  Pt_3Jets_Gen->Write();
  Eta_3Jets_Reco->Write();
  Phi_3Jets_Reco->Write();
  Pt_3Jets_Reco->Write();
  Pt_Res_3Jets->Write();
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
Acceptance::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenParticleCollection>      genParticles;
  edm::Handle< std::vector<reco::Photon> >      recoPhotons;
  //  edm::Handle< std::vector<reco::GenJet> >      genJets;
  // edm::Handle< std::vector<pat::Jet> >          recoJets;

  iEvent.getByLabel("genParticles", genParticles);                   
  iEvent.getByLabel("photons", recoPhotons); 
  // iEvent.getByLabel("ak5GenJets", genJets);
  // iEvent.getByLabel("patJetsAK5PF", recoJets);


  bool MEPhotonFound = 0;
  bool GENPhotonFound = 0;
  bool RECOPhotonFound = 0;
  
  bool Debug = 1; 

  double genEta =-25.0, genPhi =-25.0, genPt = -25.0;
  double recoEta=-25.0, recoPhi=-25.0, recoPt= -25.0;

  for(unsigned int i=0; i<genParticles->size(); ++i) {
    if(MEPhotonFound) continue;
    g = &((*genParticles)[i]);
    if (g->status() != 3) continue;
    if (g->pdgId() == 22) {
      MEPhotonFound = 1;
      reco::GenParticle::daughters d = g->daughterRefVector();
      for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
	if(GENPhotonFound) continue;
	if((*it_d)->status()==1 && (*it_d)->pdgId()==22) {
	  b = &(*(*it_d));
	  if(b->pt() > 100) {
	    GENPhotonFound = 1;
	    if(Debug) std::cout<<"ME Photon: id = "<<std::setw(2)<<g->pdgId()<<" | st = "<<std::setw(2)<<g->status()<<" | eta = "<<std::setw(9)<<g->eta()<<" | phi = ";
	    if(Debug) std::cout<<std::setw(9)<<g->phi()<<" | pt = "<<std::setw(9)<<g->pt()<<" GeV/c | et = "<<std::setw(9)<<g->et()<<" GeV | Energy = "<<g->energy()<<" GeV";
	    if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<g->mt()<<" GeV/c^2"<<std::endl;
	    if(Debug) std::cout<<"GEN Photon: id = "<<std::setw(2)<<b->pdgId()<<" | st = "<<std::setw(2)<<b->status()<<" | eta = "<<std::setw(9)<<b->eta()<<" | phi = ";
	    if(Debug) std::cout<<std::setw(9)<<b->phi()<<" | pt = "<<std::setw(9)<<b->pt()<<" GeV/c | et = "<<std::setw(9)<<b->et()<<" GeV | Energy = "<<b->energy()<<" GeV ";
	    if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<b->mt()<<" GeV/c^2"<<std::endl;
	    genEta = b->eta();
	    genPhi = b->phi();
	    genPt = b->pt();
	    
	    // Check if This Photon is the photon found in RECO
	    p = 0;
	    for(unsigned int i=0; i<recoPhotons->size(); ++i) {
	      if(RECOPhotonFound) continue;
	      p = &((*recoPhotons)[i]);
	      // if(reco::deltaR(b->eta(),b->phi(),p->eta(),p->phi()) < 0.5 && p->pt() > 100) {
	      if(p->pt() > 100) {
		RECOPhotonFound = 1;
		if(Debug) std::cout<<"RECO Photon: id = "<<std::setw(2)<<p->pdgId()<<" | st = "<<std::setw(2)<<p->status()<<" | eta = "<<std::setw(9)<<p->eta()<<" | phi = ";
		if(Debug) std::cout<<std::setw(9)<<p->phi()<<" | pt = "<<std::setw(9)<<p->pt()<<" GeV/c | et = "<<std::setw(9)<<p->et()<<" GeV | Energy = "<<p->energy()<<" GeV";
		if(Debug) std::cout<<" | transverse mass = "<<std::setw(9)<<p->mt()<<" GeV/c^2"<<std::endl;
		recoEta = p->eta();
		recoPhi = p->phi();
		recoPt = p->pt();
	      }
	    }
	  }
	  // In case you did not save reco::Photon but fortunately saved PhotonCore and the three SuperCluster Collections
	  //     reco::SuperClusterRef superCluster() const {return this->photonCore()->superCluster();}
	  //     reco::SuperClusterRef superCluster() const {return this->photonCore()->superCluster();}
	  // To bad. SuperCluster does only contain energy and not pt, so requires additional computation
	}
      }
      // Require 3 GenJets with pt > 50 and |eta| < 2.5
      /*
      int goodGenJets_AK5PT50ETA25 = 0, goodGenJets_AK5PT30ETA50 = 0;
      for(unsigned int i=0; i<genJets->size(); ++i) {
	gj = &((*genJets)[i]);
	if (gj->pt() > 50 && fabs(gj->eta()) < 2.5) {++goodGenJets_AK5PT50ETA25;}
	if (gj->pt() > 30)                          {++goodGenJets_AK5PT30ETA50;}
      }
      int goodGenJets = goodGenJets_AK5PT50ETA25;
      */
      if(GENPhotonFound) {
	Fill(Eta_All_Gen, n_eta, n1_eta, n2_eta, genEta);
	Fill(Phi_All_Gen, n_phi, n1_phi, n2_phi, genPhi);
	Fill(Pt_All_Gen,  n_et,  n1_et,  n2_et,  genPt);
      }
      if(RECOPhotonFound) {
	Fill(Eta_All_Reco, n_eta, n1_eta, n2_eta, genEta);
	Fill(Phi_All_Reco, n_phi, n1_phi, n2_phi, genPhi);
	Fill(Pt_All_Reco,  n_et,  n1_et,  n2_et,  genPt);
	Fill(Pt_Res_All,  n_res,  n1_res,  n2_res, (genPt-recoPt)/genPt);
      }
      /*
      if(goodGenJets >= 3) {
	if(GENPhotonFound) {
	  Fill(Eta_3Jets_Gen, n_eta, n1_eta, n2_eta, genEta);
	  Fill(Phi_3Jets_Gen, n_phi, n1_phi, n2_phi, genPhi);
	  Fill(Pt_3Jets_Gen,  n_et,  n1_et,  n2_et,  genPt);
	}
	if(RECOPhotonFound) {
	  Fill(Eta_3Jets_Reco, n_eta, n1_eta, n2_eta, genEta);
	  Fill(Phi_3Jets_Reco, n_phi, n1_phi, n2_phi, genPhi);
	  Fill(Pt_3Jets_Reco,  n_et,  n1_et,  n2_et,  genPt);
	  Fill(Pt_Res_3Jets,  n_res,  n1_res,  n2_res, (genPt-recoPt)/genPt);
	}
      }
      */
    }
  }
}

void
Acceptance::Fill(TH1 * Histo, int n, double n1, double n2, double value) {
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
Acceptance::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Acceptance::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(Acceptance);
