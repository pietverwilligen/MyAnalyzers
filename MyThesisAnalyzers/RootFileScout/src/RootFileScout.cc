// -*- C++ -*-
//
// Package:    RootFileScout
// Class:      RootFileScout
// 
/**\class RootFileScout RootFileScout.cc MyAnalyzers/RootFileScout/src/RootFileScout.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  local user
//         Created:  Tue Jul 13 10:23:27 CEST 2010
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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/Photon.h>
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"

#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/EgammaReco/interface/ElectronSeed.h>
//
// class declaration
//

class RootFileScout : public edm::EDAnalyzer {
   public:
      explicit RootFileScout(const edm::ParameterSet&);
      ~RootFileScout();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RootFileScout::RootFileScout(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


RootFileScout::~RootFileScout()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
RootFileScout::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //
  // Event Info
  //
  int evNum = (iEvent.id()).event();
  int rnNum = (iEvent.id()).run();


  //
  // Handles for Data
  //
  edm::Handle< reco::GenParticleCollection >    genParticles;
  // edm::Handle< reco::GenParticleCollection >    prunedGenParticles;
  edm::Handle< reco::GenJetCollection >         genJets;
  edm::Handle< GenEventInfoProduct >            genEventInfo;
  edm::Handle< std::vector<reco::PFCandidate> > recoParticles;
  edm::Handle< std::vector<reco::PFJet> >       recoJets;
  edm::Handle< CaloTowerCollection >            caloTowers;
  edm::Handle< std::vector<reco::Photon> >      photons;
  edm::Handle< std::vector<pat::Photon> >       patPhotons;


  iEvent.getByLabel("genParticles",       genParticles);
  // iEvent.getByLabel("prunedGenParticles", prunedGenParticles);
  iEvent.getByLabel("ak5GenJets",         genJets);
  iEvent.getByLabel("generator",          genEventInfo);
  iEvent.getByLabel("particleFlow",       recoParticles);
  iEvent.getByLabel("ak5PFJets",          recoJets);
  iEvent.getByLabel("towerMaker",         caloTowers);
  iEvent.getByLabel("photons",            photons);
  // iEvent.getByLabel("cleanPatPhotons",    patPhotons);

  double evWeight = genEventInfo->weight();
  double qScale   = genEventInfo->qScale();
  double pthat    = (genEventInfo->hasBinningValues() ? (genEventInfo->binningValues())[0] : 0.0);
  std::cout<<"evWeight = "<<evWeight<<" qScale = "<<qScale<<" pthat = "<<pthat<<std::endl;


  double jets_pt, jets_et, jets_eta, jets_phi, jets_Eem, jets_Ehad;
  double p_id, p_st, p_pt, p_et, p_eta, p_phi, p_sIEIE, p_HE, p_eIso, p_hTIso, p_hIso;
  int p_nTH;
  bool p_iP, p_hPS, p_hCT;
  double p_ecalIso, p_hcalIso, p_trackIso;

  std::cout<<"--------------- Run: "<<std::setw(6)<<rnNum<<" Event: "<<std::setw(6)<<evNum<<" Weight: "<<std::setw(12)<<evWeight<<" ---------------"<<std::endl;
  for (unsigned int i=0; i<genJets->size(); ++i) {
    const reco::GenJet & p = (*genJets)[i];
    jets_pt = p.pt(); jets_et = p.et(); jets_eta = p.eta(); jets_phi = p.phi();
    jets_Eem = p.emEnergy(); jets_Ehad = p.hadEnergy();
    std::cout<<"GenJet "<<std::setw(6)<<i<<" | pt = "<<std::setw(12)<<jets_pt<<" GeV/c | et = "<<std::setw(12)<<jets_et<<" GeV | eta = "<<std::setw(12)<<jets_eta;
    std::cout<<" | phi = "<<std::setw(12)<<jets_phi;
    std::cout<<" | EM energy = "<<std::setw(12)<<jets_Eem<<" | HAD energy = "<<std::setw(12)<<jets_Ehad<<std::endl;
  }
  std::cout<<"-----------------------------------------------------------------------"<<std::endl;
  /*
  for (unsigned int i=0; i<prunedGenParticles->size(); ++i) {
    const reco::GenParticle & p = (*prunedGenParticles)[i];
    p_id = p.pdgId(); p_st = p.status(); p_pt = p.pt(); p_et = p.et(); p_eta = p.eta(); p_phi = p.phi(); 
    if (p_pt>=15 && fabs(p_eta)<=2.5) {
      std::cout<<"prunedGenParticle "<<std::setw(3)<<i<<" | id = "<<std::setw(5)<<p_id<<" | st = "<<std::setw(1)<<p_st<<" | pt = "<<std::setw(12)<<p_pt;
      std::cout<<" GeV/c | et = "<<std::setw(12)<<p_et<<" GeV | eta = "<<std::setw(12)<<p_eta<<" | phi = "<<std::setw(12)<<p_phi<<" |"<<std::endl;
    }
  }
  */
  std::cout<<"-----------------------------------------------------------------------"<<std::endl;
  std::cout<<"------     PDG chap 34 Monte Carlo Particle Numbering Scheme     ------"<<std::endl;
  std::cout<<"      Quarks:        1- 6  (duscbd)"<<std::endl;
  std::cout<<"      Leptons:      11-16  (e- nu_e u-, nu_u, t-, nu_t)"<<std::endl;
  std::cout<<"      Gauge Bosons: 21-24 (gluon, gamma, Z0, W+)"<<std::endl;
  std::cout<<"      Mesons:       pi0 = 111, pi+ = 211, eta = 221"<<std::endl;
  std::cout<<"      Baryons:      p = 2212"<<std::endl;
  std::cout<<"      Diquarks:     (ud)_0 = 2101"<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"------     Particles Involved in Matrix Element Calculation      ------"<<std::endl;
  for (unsigned int i=0; i<genParticles->size(); ++i) {
    const reco::GenParticle & p = (*genParticles)[i];
    if(p.status() != 3) continue;
    p_id = p.pdgId(); p_st = p.status(); p_pt = p.pt(); p_et = p.et(); p_eta = p.eta(); p_phi = p.phi();
    std::cout<<"genParticle "<<std::setw(3)<<i<<" | id = "<<std::setw(5)<<p_id<<" | pt = "<<std::setw(12)<<p_pt;
    std::cout<<" GeV/c | et = "<<std::setw(12)<<p_et<<" GeV | eta = "<<std::setw(12)<<p_eta<<" | phi = "<<std::setw(12)<<p_phi<<" |";
    std::cout<<" mass = "<<std::setw(12)<<p.mass()<<" GeV/c^2 | transverse mass = "<<std::setw(12)<<p.mt()<<" GeV/c^2"<<std::endl;
  }
  std::cout<<"------     Matrix Element Event Reconstruction                   ------"<<std::endl;
  for (unsigned int i=0; i<genParticles->size(); ++i) {
    const reco::GenParticle & p = (*genParticles)[i];
    if(p.status() != 3) continue;
    std::cout<<"genParticle "<<std::setw(3)<<i<</*" | index "<<std::setw(3)<<p.index()<<*/" | id = "<<std::setw(5)<<p.pdgId()<<" | mothers = ";
    reco::GenParticle::mothers m = p.motherRefVector();
    reco::GenParticle::daughters d = p.daughterRefVector();
    int setwidthsize = 0;
    if (m.size()==0) std::cout<<std::setw(10)<<"";
    else if (m.size() == 1) setwidthsize = 8;
    else setwidthsize = 3;
    for (reco::GenParticle::mothers::const_iterator it_m = m.begin(), e = m.end(); it_m != e; ++it_m) {
      std::cout<<std::setw(setwidthsize)<<(**it_m).pdgId()<<", ";
    }
    std::cout<<" | daughters = ";
    for (reco::GenParticle::daughters::const_iterator it_d = d.begin(), e = d.end(); it_d != e; ++it_d) {
      if((**it_d).status()==3) std::cout<<(**it_d).pdgId()<<", ";
    }
    std::cout<<""<<std::endl;
  }
  std::cout<<"------     Status 1 Photon                                       ------"<<std::endl;
  for (unsigned int i=0; i<genParticles->size(); ++i) {
    const reco::GenParticle & p = (*genParticles)[i];
    if(p.status() != 1 || p.pdgId() != 22 || p.pt() < 15) continue;
    p_id = p.pdgId(); p_st = p.status(); p_pt = p.pt(); p_et = p.et(); p_eta = p.eta(); p_phi = p.phi();
    std::cout<<"genPhoton   "<<std::setw(3)<<i<<" | id = "<<std::setw(5)<<p_id<<" | pt = "<<std::setw(12)<<p_pt;
    std::cout<<" GeV/c | et = "<<std::setw(12)<<p_et<<" GeV | eta = "<<std::setw(12)<<p_eta<<" | phi = "<<std::setw(12)<<p_phi<<" |";
    std::cout<<" mass = "<<std::setw(12)<<p.mass()<<" GeV/c^2 | transverse mass = "<<std::setw(12)<<p.mt()<<" GeV/c^2"<<std::endl;
  }
  std::cout<<"-----------------------------------------------------------------------"<<std::endl;
  for (unsigned int i=0; i<photons->size(); ++i) {
    const reco::Photon & p = (*photons)[i];
    p_pt = p.pt(); p_et = p.et(); p_eta = p.eta(); p_phi = p.phi(); p_sIEIE = p.sigmaIetaIeta(); p_HE = p.hadronicOverEm(); 
    p_eIso = p.ecalRecHitSumEtConeDR04(); p_hTIso = p.trkSumPtHollowConeDR04(); p_hIso = p.hcalTowerSumEtConeDR04();
    p_nTH = p.nTrkHollowConeDR04(); p_iP = p.isPhoton(); p_hPS = p.hasPixelSeed(); p_hCT = p.hasConversionTracks();
    if (p_pt>=15 && fabs(p_eta)<=2.5) {
      std::cout<<"Photon "<<std::setw(6)<<i<<" | pt = "<<std::setw(12)<<p_pt<<" GeV/c | et = "<<std::setw(12)<<p_et<<" GeV | eta = "<<std::setw(12)<<p_eta;
      std::cout<<" | phi = "<<std::setw(12)<<p_phi<<" | sigmaIEtaIEta = "<<std::setw(12)<<p_sIEIE<<" | Had/Em = "<<std::setw(12)<<p_HE<<std::endl;
      std::cout<<" | ECAL Iso = "<<std::setw(12)<<p_eIso<<" GeV | TRK Iso = "<<std::setw(12)<<p_hTIso<<" GeV/c| HCAL Iso = "<<std::setw(12)<<p_hIso;
      std::cout<<" GeV | # tracks in cone = "<<std::setw(2)<<p_nTH<<" | Photon? = "<<std::setw(1)<<p_iP<<" | PixelSeed? = "<<std::setw(1)<<p_hPS;
      std::cout<<" | ConversionTks? = "<<std::setw(1)<<p_hCT<<std::endl;
    }
  }
  std::cout<<"-----------------------------------------------------------------------"<<std::endl;
  /*
  for (unsigned int i=0; i<patPhotons->size(); ++i) {
    const pat::Photon & p = (*patPhotons)[i];
    p_pt = p.pt(); p_et = p.et(); p_eta = p.eta(); p_phi = p.phi(); p_sIEIE = p.sigmaIetaIeta(); p_HE = p.hadronicOverEm(); 
    p_eIso = p.ecalRecHitSumEtConeDR04(); p_hTIso = p.trkSumPtHollowConeDR04(); p_hIso = p.hcalTowerSumEtConeDR04();
    p_nTH = p.nTrkHollowConeDR04(); p_iP = p.isPhoton(); p_hPS = p.hasPixelSeed(); p_hCT = p.hasConversionTracks();
    p_ecalIso = p.ecalIso(); p_hcalIso = p.hcalIso(); p_trackIso = p.trackIso();
    const reco::GenParticle & gen = (reco::GenParticle &) * (p.genPhoton());
    if (p_pt>=15 && fabs(p_eta)<=2.5) {
      std::cout<<"PAT Photon "<<std::setw(6)<<i<<" | pt = "<<std::setw(12)<<p_pt<<" GeV/c | et = "<<std::setw(12)<<p_et<<" GeV | eta = "<<std::setw(12)<<p_eta;
      std::cout<<" | phi = "<<std::setw(12)<<p_phi<<" | sigmaIEtaIEta = "<<std::setw(12)<<p_sIEIE<<" | Had/Em = "<<std::setw(12)<<p_HE<<std::endl;
      std::cout<<" | ECAL Iso = "<<std::setw(12)<<p_eIso<<" GeV | TRK Iso = "<<std::setw(12)<<p_hTIso<<" GeV/c| HCAL Iso = "<<std::setw(12)<<p_hIso;
      std::cout<<" GeV | # tracks in cone = "<<std::setw(2)<<p_nTH<<" | Photon? = "<<std::setw(1)<<p_iP<<" | PixelSeed? = "<<std::setw(1)<<p_hPS;
      std::cout<<" | ConversionTks? = "<<std::setw(1)<<p_hCT<<std::endl;
      std::cout<<" | PAT ecalIso = "<<std::setw(12)<<p_ecalIso<<" GeV | PAT trackIso = "<<std::setw(12)<<p_trackIso<<" GeV/c| PAT hcalIso = ";
      std::cout<<std::setw(12)<<p_hcalIso<<std::endl;
      if (&gen !=0) {
	const reco::GenParticle & mot = (reco::GenParticle &) *(gen.mother());
	std::cout<<"GenParticle | id = "<<std::setw(5)<<gen.pdgId()<<" | st = "<<std::setw(5)<<gen.status()<<" | pt = "<<std::setw(12)<<gen.pt();
	std::cout<<" GeV/c | et = "<<std::setw(12)<<gen.et()<<" GeV | eta = "<<std::setw(12)<<gen.eta()<<" | phi = "<<std::setw(12)<<gen.phi()<<std::endl;
	if(&mot != 0) {
	  std::cout<<"MotherPart  | id = "<<std::setw(5)<<mot.pdgId()<<" | st = "<<std::setw(5)<<mot.status()<<" | pt = "<<std::setw(12)<<mot.pt();
	  std::cout<<" GeV/c | et = "<<std::setw(12)<<mot.et()<<" GeV | eta = "<<std::setw(12)<<mot.eta()<<" | phi = "<<std::setw(12)<<mot.phi()<<std::endl;
	}
	std::cout<<"- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
      }
    }
  }
  */
  std::cout<<"-----------------------------------------------------------------------"<<std::endl;


}


// ------------ method called once each job just before starting event loop  ------------
void 
RootFileScout::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RootFileScout::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(RootFileScout);
