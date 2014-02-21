// -*- C++ -*-
//
// Package:    MyRPCTriggerAnalyzer
// Class:      MyRPCTriggerAnalyzer
// 
/**\class MyRPCTriggerAnalyzer MyRPCTriggerAnalyzer.cc MyAnalyzers/MyRPCTriggerAnalyzer/src/MyRPCTriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet Verwilligen,161 R-006,+41227676292,
//         Created:  Sun Nov 10 12:43:09 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// root include files
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

// RPC Geometry
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "DataFormats/Provenance/interface/Timestamp.h"

#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>

// #include "DQMServices/Core/interface/DQMStore.h"
// #include "DQMServices/Core/interface/MonitorElement.h"

#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>


// L1 Trigger
#include <DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h>
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include <DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h>
//
// class declaration
//

class MyRPCTriggerAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyRPCTriggerAnalyzer(const edm::ParameterSet&);
      ~MyRPCTriggerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  std::vector<int> m_GMTcandidatesBx;
  std::vector<int> m_DTcandidatesBx;
  std::vector<int> m_RPCcandidatesBx;

  bool debug;
  std::string rootFileName;
  TFile * outputfile;
 
  TH1F * RPC_Triggers_ETA_All, * RPC_Triggers_ETA_Q0, * RPC_Triggers_ETA_Q1, * RPC_Triggers_ETA_Q2, * RPC_Triggers_ETA_Q3; 
  TH1F * RPC_Triggers_PHI_All, * RPC_Triggers_PHI_Q0, * RPC_Triggers_PHI_Q1, * RPC_Triggers_PHI_Q2, * RPC_Triggers_PHI_Q3; 

  TH2F * RPC_Triggers_ETA_PHI_All;

  // edm::InputTag m_rpcDigiLabel;
  edm::InputTag m_gtReadoutLabel;
  edm::InputTag m_gmtReadoutLabel;
};

//
// constants, enums and typedefs
//
int n_phi = 144;
int n_eta =  48;
double n_phi_1 =  0.0000; 
double n_phi_2 =  6.2832;
double n_eta_1 = -2.40;
double n_eta_2 =  2.40;
double n_eta_exact = 64;
double n_eta_vec[] = {-2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2, -1.95, -1.9, -1.85, -1.8, -1.75, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45};

//
// static data member definitions
//

//
// constructors and destructor
//
MyRPCTriggerAnalyzer::MyRPCTriggerAnalyzer(const edm::ParameterSet& iConfig)

{
   if(debug) std::cout<<"MyRPCTriggerAnalyzer :: Constructor]"<<std::endl; 
   //now do what ever initialization is needed
   m_gtReadoutLabel     = iConfig.getParameter<edm::InputTag>("GTReadoutRcd");
   m_gmtReadoutLabel    = iConfig.getParameter<edm::InputTag>("GMTReadoutRcd");
   rootFileName         = iConfig.getUntrackedParameter<std::string>("RootFileName");
   debug                = iConfig.getUntrackedParameter<bool>("Debug");

   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

   RPC_Triggers_PHI_All = new TH1F("RPC_Triggers_PHI_All,",  "RPC_Triggers_PHI_All,", n_phi, n_phi_1, n_phi_2);
   RPC_Triggers_PHI_Q0  = new TH1F("RPC_Triggers_PHI_Q0,",   "RPC_Triggers_PHI_Q0,",  n_phi, n_phi_1, n_phi_2);
   RPC_Triggers_PHI_Q1  = new TH1F("RPC_Triggers_PHI_Q1,",   "RPC_Triggers_PHI_Q1,",  n_phi, n_phi_1, n_phi_2);
   RPC_Triggers_PHI_Q2  = new TH1F("RPC_Triggers_PHI_Q2,",   "RPC_Triggers_PHI_Q2,",  n_phi, n_phi_1, n_phi_2);
   RPC_Triggers_PHI_Q3  = new TH1F("RPC_Triggers_PHI_Q3,",   "RPC_Triggers_PHI_Q3,",  n_phi, n_phi_1, n_phi_2);
   
   RPC_Triggers_ETA_All = new TH1F("RPC_Triggers_ETA_All,",  "RPC_Triggers_ETA_All,", n_eta_exact, n_eta_vec);
   RPC_Triggers_ETA_Q0  = new TH1F("RPC_Triggers_ETA_Q0,",   "RPC_Triggers_ETA_Q0,",  n_eta_exact, n_eta_vec);
   RPC_Triggers_ETA_Q1  = new TH1F("RPC_Triggers_ETA_Q1,",   "RPC_Triggers_ETA_Q1,",  n_eta_exact, n_eta_vec);
   RPC_Triggers_ETA_Q2  = new TH1F("RPC_Triggers_ETA_Q2,",   "RPC_Triggers_ETA_Q2,",  n_eta_exact, n_eta_vec);
   RPC_Triggers_ETA_Q3  = new TH1F("RPC_Triggers_ETA_Q3,",   "RPC_Triggers_ETA_Q3,",  n_eta_exact, n_eta_vec);

   RPC_Triggers_ETA_PHI_All = new TH2F("RPC_Triggers_ETA_PHI_All,",  "RPC_Triggers_ETA_PHI_All,", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);
}


MyRPCTriggerAnalyzer::~MyRPCTriggerAnalyzer()
{
   if(debug) std::cout<<"MyRPCTriggerAnalyzer :: Destructor :: begin]"<<std::endl; 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   outputfile->cd();

   RPC_Triggers_PHI_All->Write();
   RPC_Triggers_PHI_Q0->Write();
   RPC_Triggers_PHI_Q1->Write();
   RPC_Triggers_PHI_Q2->Write();
   RPC_Triggers_PHI_Q3->Write();
  
   RPC_Triggers_ETA_All->Write();
   RPC_Triggers_ETA_Q0->Write();
   RPC_Triggers_ETA_Q1->Write();
   RPC_Triggers_ETA_Q2->Write();
   RPC_Triggers_ETA_Q3->Write();
   
   RPC_Triggers_ETA_PHI_All->Write();
   
   outputfile->Close();
   if(debug) std::cout<<"MyRPCTriggerAnalyzer :: Destructor :: end]"<<std::endl; 
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyRPCTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
/* 
using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/

  int evNum = (iEvent.id()).event();
  int rnNum = (iEvent.id()).run();

  edm::Handle<L1MuGMTReadoutCollection> pCollection;
  iEvent.getByLabel(m_gmtReadoutLabel,pCollection);

  if ( ! pCollection.isValid() ) {
    // edm::LogError("discriminateGMT") << "can't find L1MuGMTReadoutCollection with label "<< m_gmtReadoutLabel ;
    std::cout<<"can't find L1MuGMTReadoutCollection with label "<< m_gmtReadoutLabel<<std::endl;
    // return -1; 
  }
  
  // get GMT readout collection
  const L1MuGMTReadoutCollection * gmtRC = pCollection.product();
  
  // get record vector
  std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr;
  std::vector<L1MuGMTReadoutRecord> gmt_records = gmtRC->getRecords();
  
  edm::LogInfo("DiscriminateGMT") << "nRecords: " << gmt_records.size() << '\n';
  
  for( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
    
    int BxInEvent = RRItr->getBxInEvent();
    int BxInEventNew = RRItr->getBxNr();
    
    // RPC barrel muon candidates
    int nrpcB = 0;
    int ndtB  = 0;
    
    std::vector<L1MuRegionalCand> BrlRpcCands = RRItr->getBrlRPCCands();
    std::vector<L1MuRegionalCand> BrlDtCands  = RRItr->getDTBXCands ();
    
    std::vector<L1MuRegionalCand>::const_iterator RCItr;
    
    for( RCItr = BrlRpcCands.begin(); RCItr !=BrlRpcCands.end(); ++RCItr) {
      if ( !(*RCItr).empty() ) {
	m_GMTcandidatesBx.push_back( BxInEventNew );
	nrpcB++;
	if(debug) {
	  std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | "; 
	  std::cout<<"RPC Barrel Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<std::endl;
	}
	// std::cout<<"RPC Barrel Trigger :: "<<RCItr->print()<<std::endl;

	if(!debug) {

	  int quality = RCItr->quality();
	  double eta = RCItr->etaValue();
	  double phi = RCItr->phiValue();

	  if(debug) std::cout<<"Fill All Histos"<<std::endl;
	  RPC_Triggers_ETA_PHI_All->Fill(eta, phi);
	  RPC_Triggers_ETA_All->Fill(eta);
	  RPC_Triggers_PHI_All->Fill(phi);
	  
	  switch (quality) {
	  case 0: {
	    if(debug) std::cout<<"Fill Quality 0 Histos"<<std::endl;
	    RPC_Triggers_ETA_Q0->Fill(eta);
	    RPC_Triggers_PHI_Q0->Fill(phi);
	  } break;
	  case 1: {
	    if(debug) std::cout<<"Fill Quality 1 Histos"<<std::endl;
	    RPC_Triggers_ETA_Q1->Fill(eta);
	    RPC_Triggers_PHI_Q1->Fill(phi);
	  } break;
	  case 2: {
	    if(debug) std::cout<<"Fill Quality 2 Histos"<<std::endl;
	    RPC_Triggers_ETA_Q2->Fill(eta);
	    RPC_Triggers_PHI_Q2->Fill(phi);
	  } break;
	  case 3: {
	    if(debug) std::cout<<"Fill Quality 3 Histos"<<std::endl;
	    RPC_Triggers_ETA_Q3->Fill(eta);
	    RPC_Triggers_PHI_Q3->Fill(phi);
	  } break;
	  default : std::cout<<"Quality = "<<quality<<" > 3"<<std::endl;
	  }
	}

      }
    }
    // std::cout<<"End RPC Stuff"<<std::endl;
    for( RCItr = BrlDtCands.begin(); RCItr !=BrlDtCands.end(); ++RCItr) {
      if ( !(*RCItr).empty() ) {
	m_DTcandidatesBx.push_back( BxInEventNew );
	ndtB++;
      }
    } 
  }
  // if(debug) std::cout<<"End Event"<<std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyRPCTriggerAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRPCTriggerAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyRPCTriggerAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MyRPCTriggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyRPCTriggerAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyRPCTriggerAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRPCTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRPCTriggerAnalyzer);
