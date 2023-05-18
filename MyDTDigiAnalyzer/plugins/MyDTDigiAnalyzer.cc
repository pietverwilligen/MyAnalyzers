// -*- C++ -*-
//
// Package:    MyDTDigiAnalyzer
// Class:      MyDTDigiAnalyzer
// 
/**\class MyDTDigiAnalyzer MyDTDigiAnalyzer.cc MyAnalyzers/MyDTDigiAnalyzer/plugins/MyDTDigiAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri, 19 Jun 2015 14:18:37 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/Utilities/interface/InputTag.h"

// Digis
#include <DataFormats/DTDigi/interface/DTDigi.h>
#include <DataFormats/DTDigi/interface/DTDigiCollection.h>
#include <DataFormats/MuonDetId/interface/DTLayerId.h>
#include <DataFormats/MuonDetId/interface/DTChamberId.h>
#include "CondFormats/DataRecord/interface/DTReadOutMappingRcd.h"

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

//
// class declaration
//

class MyDTDigiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyDTDigiAnalyzer(const edm::ParameterSet&);
      ~MyDTDigiAnalyzer();

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
  // edm::InputTag dtDigiLabel;
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
MyDTDigiAnalyzer::MyDTDigiAnalyzer(const edm::ParameterSet& iConfig)

{
  // now do what ever initialization is needed
  // dtDigiLabel = iConfig.getParameter<InputTag>("dtDigiLabel");
}


MyDTDigiAnalyzer::~MyDTDigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyDTDigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get the digis from the event
  edm::Handle<DTDigiCollection> dtdigis;
  // iEvent.getByLabel(dtDigiLabel, dtdigis);
  iEvent.getByLabel("simMuonDTDigis", dtdigis);

  // LOOP OVER ALL THE DIGIS OF THE EVENT
  DTDigiCollection::DigiRangeIterator dtLayerId_It;
  for (dtLayerId_It=dtdigis->begin(); dtLayerId_It!=dtdigis->end(); ++dtLayerId_It){
    for (DTDigiCollection::const_iterator digiIt = ((*dtLayerId_It).second).first;
	 digiIt!=((*dtLayerId_It).second).second; ++digiIt){
 
      DTChamberId chId = ((*dtLayerId_It).first).chamberId();
      DTLayerId lyId = ((*dtLayerId_It).first);

      if(chId.wheel()==0 && chId.sector()==11) {
	std::cout<<"DT Digi in DetId "<<((*dtLayerId_It).first).rawId()<< " = LayerId " << lyId<< " at wire "<<(*digiIt).wire()<<" with TDC count "<<(*digiIt).countsTDC()<<" and time "<<(*digiIt).time()<<" ns"<<std::endl;  
      }

      // Check the TDC trigger width
      /*
      int tdcTime = (*digiIt).countsTDC();
      double upperLimit = tTrigStMap[(*dtLayerId_It).first.superlayerId().chamberId()]-safeMargin;
      if(doTimeBoxHistos)
	tbHistos[(*dtLayerId_It).first.superlayerId()]->Fill(tdcTime);
      if(tdcTime>upperLimit)
	continue;
      */

    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyDTDigiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyDTDigiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MyDTDigiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyDTDigiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MyDTDigiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MyDTDigiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyDTDigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyDTDigiAnalyzer);
