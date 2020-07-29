// -*- C++ -*-
//
// Package:    temp/StealthTriggerEfficiency
// Class:      StealthTriggerEfficiency
//
/**\class StealthTriggerEfficiency StealthTriggerEfficiency.cc temp/StealthTriggerEfficiency/plugins/StealthTriggerEfficiency.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Tanmay Mudholkar
//         Created:  Thu, 20 Feb 2020 02:41:34 GMT
//
//


#include "../interface/StealthTriggerEfficiency.h"

//
// constructors and destructor
//
StealthTriggerEfficiency::StealthTriggerEfficiency(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");
  edm::Service<TFileService> fileService;
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity");
  edm::LogInfo("StealthTriggerEfficiency") << "Starting StealthTriggerEfficiency with verbosity: " << verbosity_;

  eventInfoTree_ = fileService->make<TTree>("eventInfoTree", "eventInfoTree");
  eventInfoTree_->Branch("pT_leadingPhoton", &pT_leadingPhoton_, "pT_leadingPhoton/F");
  eventInfoTree_->Branch("eta_leadingPhoton", &eta_leadingPhoton_, "eta_leadingPhoton/F");
  eventInfoTree_->Branch("pT_subLeadingPhoton", &pT_subLeadingPhoton_, "pT_subLeadingPhoton/F");
  eventInfoTree_->Branch("eta_subLeadingPhoton", &eta_subLeadingPhoton_, "eta_subLeadingPhoton/F");
  eventInfoTree_->Branch("passesSelection", &passesSelection_, "passesSelection/O");
  eventInfoTree_->Branch("passesTrigger", &passesTrigger_, "passesTrigger/O");

  photonCollection_ = consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonSrc"));
}


StealthTriggerEfficiency::~StealthTriggerEfficiency() {
  // no manual cleanup needed
}

//
// member functions
//

// ------------ method called for each event  ------------
void
StealthTriggerEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonCollection_, photonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("MinimalMiniAODAnalyzer") << "no pat::Photons in event";
    return;
  }

  // set event variables to default values
  pT_leadingPhoton_ = -9.;
  eta_leadingPhoton_ = -9.;
  pT_subLeadingPhoton_ = -9.;
  eta_subLeadingPhoton_ = -9.;
  passesSelection_ = false;
  passesTrigger_ = false;

  // for (edm::View<pat::Photon>::const_iterator edmPhoton = photonHandle->begin(); edmPhoton != photonHandle->end(); ++edmPhoton) {
  //   float photon_pT = edmPhoton->et();
  //   float photon_eta = edmPhoton->eta();
  // }
  eventInfoTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
StealthTriggerEfficiency::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
StealthTriggerEfficiency::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
StealthTriggerEfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(StealthTriggerEfficiency);
