// -*- C++ -*-
//
// Package:    MySimHitAnalyzer
// Class:      MySimHitAnalyzer
// 
/**\class MySimHitAnalyzer MySimHitAnalyzer.cc MyAnalyzers/MySimHitAnalyzer/src/MySimHitAnalyzer.cc

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


// root include files
#include <TRandom.h>
#include "TROOT.h"
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
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>

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

class MySimHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MySimHitAnalyzer(const edm::ParameterSet&);
      ~MySimHitAnalyzer();
  edm::ESHandle <RPCGeometry> rpcGeom;
  edm::ESHandle <CSCGeometry> cscGeom;

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
  std::string rootFileName;
  TFile * outputfile;

  std::vector<double> x_p, y_p, z_p, x_n, y_n, z_n, r_r, z_r;

  TGraph  * RE4_Plus_XY, * RE4_Minus_XY, * RE_YZ; 
  TCanvas * Canvas_RE4_Plus_XY, * Canvas_RE4_Minus_XY, * Canvas_RE_YZ, * Canvas_VTX; 

};

//
// constants, enums and typedefs
//

// int n_tof  = 75;  double n1_tof  = 30,    n2_tof = 45;
// int n_v_x  = 40;  double n1_v_x  = -2.0,  n2_v_x = 2.0;
// int n_v_z  = 50;  double n1_v_z  = -25.0, n2_v_z = 25.0;
// int p_v_x  = 50;  double p1_v_x  = -500.0,  p2_v_x = 500.0;
// int p_v_z  = 125;  double p1_v_z  = -1250.0, p2_v_z = 1250.0;


//
// static data member definitions
//

//
// constructors and destructor
//
MySimHitAnalyzer::MySimHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");

  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

}


MySimHitAnalyzer::~MySimHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputfile->cd();

  const int n_p = x_p.size();  double x_ap[n_p]; double y_ap[n_p];  /* double z_ap[n_p]; */
  const int n_n = x_n.size();  double x_an[n_n]; double y_an[n_n];  /* double z_an[n_n]; */
  const int n_r = r_r.size();  double r_ar[n_r]; double z_ar[n_r];
  for(int i=0; i< n_p; ++i) { x_ap[i] = x_p[i]; y_ap[i] = y_p[i]; /* z_ap[i] = z_p[i]; */}
  for(int i=0; i< n_n; ++i) { x_an[i] = x_n[i]; y_an[i] = y_n[i]; /* z_an[i] = z_n[i]; */}
  for(int i=0; i< n_r; ++i) { r_ar[i] = r_r[i]; z_ar[i] = z_r[i]; }

  RE4_Plus_XY  = new TGraph(n_p, x_ap, y_ap);
  RE4_Minus_XY = new TGraph(n_n, x_an, y_an);
  RE_YZ        = new TGraph(n_r, z_ar, r_ar);

  Canvas_RE4_Plus_XY  = new TCanvas("Canvas_RE4_Plus_XY",  "Canvas_RE4_Plus_XY", 800, 600);
  Canvas_RE4_Minus_XY = new TCanvas("Canvas_RE4_Minus_XY", "Canvas_RE4_Minus_XY", 800, 600);
  Canvas_RE_YZ        = new TCanvas("Canvas_RE_YZ",        "Canvas_RE_YZ", 800, 600);


  // XY and RZ Graphs
 Canvas_RE4_Plus_XY->cd();  RE4_Plus_XY->SetMarkerStyle(5);  RE4_Plus_XY->Draw("AP");  RE4_Plus_XY->GetXaxis()->SetTitle("X [cm]");  RE4_Plus_XY->GetYaxis()->SetTitle("Y [cm]");   RE4_Plus_XY->SetTitle("RE+4 SimHits");
 Canvas_RE4_Minus_XY->cd(); RE4_Minus_XY->SetMarkerStyle(5); RE4_Minus_XY->Draw("AP"); RE4_Minus_XY->GetXaxis()->SetTitle("X [cm]"); RE4_Minus_XY->GetYaxis()->SetTitle("Y [cm]");  RE4_Minus_XY->SetTitle("RE-4 SimHits");
 Canvas_RE_YZ->cd();        RE_YZ->SetMarkerStyle(5);        RE_YZ->Draw("AP");        RE_YZ->GetXaxis()->SetTitle("Z [cm]");        RE_YZ->GetYaxis()->SetTitle("R [cm]");         RE_YZ->SetTitle("RE SimHits");

  RE4_Plus_XY->Write();
  RE4_Minus_XY->Write();
  RE_YZ->Write();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MySimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // SimHits
  std::cout << " Getting the SimHits " <<std::endl;
  std::vector<edm::Handle<edm::PSimHitContainer> > theSimHitContainers;
  iEvent.getManyByType(theSimHitContainers);
  std::cout << " The number of SimHit Containers is  " << theSimHitContainers.size() <<std::endl;
  std::vector<PSimHit> theSimHits;
  for (int i = 0; i < int(theSimHitContainers.size()); ++i) {
    theSimHits.insert(theSimHits.end(),theSimHitContainers.at(i)->begin(),theSimHitContainers.at(i)->end());
  }
  // SimTracks
  std::vector<SimTrack> theSimTracks;
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel("g4SimHits",SimTk);
  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  std::cout << "This Event has " <<  theSimTracks.size() << " sim tracks " << std::endl;
  // SimVertices
  std::vector<SimVertex> theSimVertices; 
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel("g4SimHits",SimVtx);
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  std::cout << "This Event has " <<  theSimVertices.size() << " sim vertices " << std::endl;


  // ===========================================================
  //      Search for the right SimTrack and right SimVertex
  // ===========================================================
  int vertex_index = -1;
  for (std::vector<SimTrack>::const_iterator iTrack = theSimTracks.begin(); iTrack != theSimTracks.end(); ++iTrack) {
    SimTrack simtrack = (*iTrack);
    // if(fabs(simtrack.type()) != 13 || simtrack.momentum().vect().perp() < 95.0) continue;
    if(fabs(simtrack.type()) != 13) continue;
    else {
      double mu_pt = simtrack.momentum().pt();
      std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<simtrack.momentum().pt()<<" GeV/c from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      // double Px = simtrack.momentum().x(); double Py = simtrack.momentum().y(); double mu_pt = sqrt(pow(Px,2)+pow(Py,2));
      // std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" pt = "<<mu_pt<<" GeV/c from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      // std::cout<<"Muon SimTrack Found: id = "<<simtrack.type()<<" from Vertex no. "<<simtrack.vertIndex()<<std::endl;
      if(mu_pt > 95) {
	vertex_index = simtrack.vertIndex();
      }
    }
  }
  double vtx_x = -99.9, vtx_y = -99.9, vtx_z = -99.9; 
  int counter = 0;
  for (std::vector<SimVertex>::const_iterator iVertex = theSimVertices.begin(); iVertex != theSimVertices.end(); ++iVertex) {
    SimVertex simvertex = (*iVertex);
    if(vertex_index==counter) {
      vtx_x = simvertex.position().x();
      vtx_y = simvertex.position().y();
      vtx_z = simvertex.position().z();
      // std::cout<<"Muon SimVertex Found: x = "<<vtx_x<<" y = "<<vtx_y<<" z = "<<vtx_z<<std::endl;
    }
    else {
      // std::cout<<"PileUp SimVertex Found: x = "<<simvertex.position().x()<<" y = "<<simvertex.position().y()<<" z = "<<simvertex.position().z()<<std::endl;
    }
    ++counter;
  }
  GlobalPoint VtxGlobalPoint(vtx_x,vtx_y,vtx_z);
  std::cout<<"PileUp SimVertex Found :: "; 
  std::cout<<" Local Coordinates :: x = "<<vtx_x<<" y = "<<vtx_y<<" z = "<<vtx_z<<" || ";
  std::cout<<"Global Coordinates :: x = "<<VtxGlobalPoint.x()<<" y = "<<VtxGlobalPoint.y()<<" z = "<<VtxGlobalPoint.z()<<std::endl;



  // ===============================
  //      Loop over the SimHits
  // ===============================
  for (std::vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); ++iHit) {
    DetId theDetUnitId((*iHit).detUnitId());
    DetId simdetid= DetId((*iHit).detUnitId());

    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){ // Only RPCs
      // std::cout<<"RPC SimHit "<<std::endl;
      RPCDetId rollId(theDetUnitId);
      RPCGeomServ rpcsrv(rollId);
      // std::cout << " Reading the Roll"<<std::endl;
      const RPCRoll* rollasociated = rpcGeom->roll(rollId);
      // std::cout << " Getting the Surface"<<std::endl;
      const BoundPlane & RPCSurface = rollasociated->surface();
      GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal((*iHit).localPosition());
      // GlobalPoint VtxGlobalPoint; if(theSimVertices.size() == 0) {;}  else {theSimVertices[0].position();}
      // VTX_X->Fill(VtxGlobalPoint.x()); VTX_Y->Fill(VtxGlobalPoint.y()); VTX_Z->Fill(VtxGlobalPoint.z());
      std::cout<<"RPC SimHit in "<<std::setw(12)<<(int)rollId<<" a.k.a. "<<std::setw(24)<<rpcsrv.name()<<" details: "<<std::setw(24)<<rollId;
      std::cout<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
      std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta()<<std::endl;
      // std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;


      // Only RE4
      if ((rollId.station()==4) && rollId.region() != 0) {
	// std::cout<<"RPC SimHit in "<<std::setw(12)<<(int)rollId<<" a.k.a. "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
	// std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta();
	// std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
	// Fill Histograms
	if(rollId.region() == 1) {
	  // XY Graph RE+4
          x_p.push_back(RPCGlobalPoint.x()); y_p.push_back(RPCGlobalPoint.y()); z_p.push_back(RPCGlobalPoint.z());
	  // TOF_SimHits_RE4_Plus->Fill((*iHit).timeOfFlight());
	}
	if(rollId.region() == -1) {
	  // XY Graph RE-4
	  x_n.push_back(RPCGlobalPoint.x()); y_n.push_back(RPCGlobalPoint.y()); z_n.push_back(RPCGlobalPoint.z());
	  // TOF_SimHits_RE4_Minus->Fill((*iHit).timeOfFlight());
	}
      }
      // Only RE3
      if ((rollId.station()==3) && rollId.region() != 0) {
	// std::cout<<"RPC SimHit in "<<std::setw(12)<<(int)rollId<<" a.k.a. "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<RPCGlobalPoint.z();
	// std::cout<<" | r = "<<std::setw(12)<<RPCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<RPCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<RPCGlobalPoint.eta();
	// std::cout<<" | global position = "<<RPCGlobalPoint<<std::endl;
        // Fill Histograms
        if(rollId.region() == 1) {
          // TOF_SimHits_RE3_Plus->Fill((*iHit).timeOfFlight());
	}
	if(rollId.region() == -1) {
          // TOF_SimHits_RE3_Minus->Fill((*iHit).timeOfFlight());
	}
      }

      if(rollId.region() != 0) {
	// RZ graph endcap
	r_r.push_back(sqrt(pow(RPCGlobalPoint.x(),2) + pow(RPCGlobalPoint.y(),2))); z_r.push_back(RPCGlobalPoint.z());
      }
    }
    /*
    if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::CSC){ // Only CSCs
      CSCDetId rollId(theDetUnitId);
      // CSCGeomServ rpcsrv(rollId);
      GlobalPoint CSCGlobalPoint = cscGeom->idToDet(rollId)->toGlobal((*iHit).localPosition());
      // Only ME4
      if (rollId.station()==4) {
	// std::cout<<"CSC SimHit in "<<std::setw(24)<<rpcsrv.name()<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
	std::cout<<"CSC SimHit in "<<std::setw(24)<<rollId<<" | time t = "<<std::setw(12)<<(*iHit).timeOfFlight()<<" | z = "<<std::setw(12)<<CSCGlobalPoint.z();
	std::cout<<" | r = "<<std::setw(12)<<CSCGlobalPoint.mag()<<" | phi = "<<std::setw(12)<<CSCGlobalPoint.phi()<<" | eta = "<<std::setw(12)<<CSCGlobalPoint.eta();
	std::cout<<" | global position = "<<CSCGlobalPoint<<std::endl;
      }
    }
    */
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MySimHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MySimHitAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MySimHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
}


// ------------ method called when ending the processing of a run  ------------
void 
MySimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MySimHitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MySimHitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MySimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MySimHitAnalyzer);
