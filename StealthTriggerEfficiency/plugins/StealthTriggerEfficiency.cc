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
  eventInfoTree_->Branch("eventRho", &eventRho_, "eventRho/F");
  eventInfoTree_->Branch("pT_leadingPhoton", &pT_leadingPhoton_, "pT_leadingPhoton/F");
  eventInfoTree_->Branch("eta_leadingPhoton", &eta_leadingPhoton_, "eta_leadingPhoton/F");
  eventInfoTree_->Branch("pT_subLeadingPhoton", &pT_subLeadingPhoton_, "pT_subLeadingPhoton/F");
  eventInfoTree_->Branch("eta_subLeadingPhoton", &eta_subLeadingPhoton_, "eta_subLeadingPhoton/F");
  eventInfoTree_->Branch("passesSelection", &passesSelection_, "passesSelection/O");
  for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
    std::string branchName = "passesTrigger_patternIndex_" + std::to_string(patternIndex);
    passesTrigger_[patternIndex] = false;
    eventInfoTree_->Branch(branchName.c_str(), &(passesTrigger_[patternIndex]), (branchName + "/O").c_str());
  }

  rhoCollection_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
  triggerCollection_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerCollection"));
  photonCollection_ = consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonSrc"));
}


StealthTriggerEfficiency::~StealthTriggerEfficiency() {
  // no manual cleanup needed
}

//
// member functions
//

float
StealthTriggerEfficiency::getEffectiveArea_neutIso(const float& absEta) {
  if (absEta < 1.0) return (EAValues::neutIso_central);
  else if (absEta < 1.479) return (EAValues::neutIso_peripheral);
  return 0.;
}

float
StealthTriggerEfficiency::getEffectiveArea_phoIso(const float& absEta) {
  if (absEta < 1.0) return (EAValues::phoIso_central);
  else if (absEta < 1.479) return (EAValues::phoIso_peripheral);
  return 0.;
}

float
StealthTriggerEfficiency::getEffectiveArea_chIso(const float& absEta) {
  if (absEta < 1.0) return (EAValues::chIso_central);
  else if (absEta < 1.479) return (EAValues::chIso_peripheral);
  return 0.;
}

// ------------ method called for each event  ------------
void
StealthTriggerEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // set event variables to default values
  eventRho_ = -9.;
  pT_leadingPhoton_ = -9.;
  eta_leadingPhoton_ = -9.;
  pT_subLeadingPhoton_ = -9.;
  eta_subLeadingPhoton_ = -9.;
  passesSelection_ = false;
  for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
    passesTrigger_[patternIndex] = false;
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoCollection_, rhoHandle);
  if (!rhoHandle.isValid()) {
    edm::LogWarning("StealthTriggerEfficiency") << "rho unavailable in event";
    return;
  }
  eventRho_ = *(rhoHandle.product());

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonCollection_, photonHandle);
  if (!photonHandle.isValid()) {
    edm::LogWarning("StealthTriggerEfficiency") << "no pat::Photons in event";
    return;
  }

  for (edm::View<pat::Photon>::const_iterator edmPhoton = photonHandle->begin(); edmPhoton != photonHandle->end(); ++edmPhoton) {
    float photon_pT = edmPhoton->et();

    float photon_eta = edmPhoton->eta();
    float photon_abs_eta = std::fabs(photon_eta);
    bool inBarrel = (photon_abs_eta < photonCuts::eta);

    bool passes_electronVeto = edmPhoton->passElectronVeto();

    float photon_hOverE = edmPhoton->hadTowOverEm();
    bool passes_hOverE = (photon_hOverE < photonCuts::hOverE);

    float photon_sigma_ieta_ieta = edmPhoton->full5x5_sigmaIetaIeta();
    bool passes_sigma_ieta_ieta = (photon_sigma_ieta_ieta < photonCuts::sigma_ieta_ieta);

    float photon_neutIso_uncorrected = edmPhoton->userFloat("phoNeutralHadronIsolation");
    float photon_neutIso = std::max(0.0f, photon_neutIso_uncorrected - eventRho_*getEffectiveArea_neutIso(photon_abs_eta));
    bool passes_neutIso = (photon_neutIso < ((photonCuts::neutIso_const) + photon_pT*(photonCuts::neutIso_lin) + photon_pT*photon_pT*(photonCuts::neutIso_quad)));

    float photon_phoIso_uncorrected = edmPhoton->userFloat("phoPhotonIsolation");
    float photon_phoIso = std::max(0.0f, photon_phoIso_uncorrected - eventRho_*getEffectiveArea_phoIso(photon_abs_eta));
    bool passes_phoIso = (photon_phoIso < ((photonCuts::phoIso_const) + photon_pT*(photonCuts::phoIso_lin)));

    float photon_chIso_uncorrected = edmPhoton->userFloat("phoChargedIsolation");
    float photon_chIso = std::max(0.0f, photon_chIso_uncorrected - eventRho_*getEffectiveArea_phoIso(photon_abs_eta));
    float passes_chIso = (photon_chIso < (photonCuts::chIso));

    if (inBarrel &&
        passes_electronVeto &&
        passes_hOverE &&
        passes_sigma_ieta_ieta &&
        passes_neutIso &&
        passes_phoIso &&
        passes_chIso) {
      // check pT
      if (photon_pT >= pT_leadingPhoton_) {
        // current leading photon becomes new subleading photon
        pT_subLeadingPhoton_ = pT_leadingPhoton_;
        eta_subLeadingPhoton_ = eta_leadingPhoton_;

        // this photon becomes new leading photon
        pT_leadingPhoton_ = photon_pT;
        eta_leadingPhoton_ = photon_eta;
      }
      else if (photon_pT >= pT_subLeadingPhoton_) {
        // this photon becomes new subleading photon
        pT_subLeadingPhoton_ = photon_pT;
        eta_subLeadingPhoton_ = photon_eta;
      }
    }
  }
  if (pT_leadingPhoton_ > 0.) assert(pT_leadingPhoton_ >= pT_subLeadingPhoton_);
  passesSelection_ = ((pT_leadingPhoton_ > photonCuts::pTLeading) &&
                      (pT_subLeadingPhoton_ > photonCuts::pTSubleading));

  edm::Handle<edm::TriggerResults> triggerHandle;
  iEvent.getByToken(triggerCollection_, triggerHandle);
  edm::TriggerNames triggerNames = iEvent.triggerNames(*triggerHandle);
  if (verbosity_ > 1) edm::LogInfo("StealthTriggerEfficiency") << "N_triggers: " << triggerHandle->size();
  for (unsigned int trigIndex = 0; trigIndex < triggerHandle->size(); ++trigIndex) {
    if (verbosity_ > 2) edm::LogInfo("StealthTriggerEfficiency") << "name[" << trigIndex << "]:" << triggerNames.triggerName(trigIndex);
  }

  for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
    const std::string& triggerName = (triggerPatterns::patternsToSave)[patternIndex];
    std::vector< std::vector<std::string>::const_iterator > triggerMatches = edm::regexMatch(triggerNames.triggerNames(), triggerName);
    if (verbosity_ > 1) edm::LogInfo("StealthTriggerEfficiency") << "For trigger with name: " << triggerName << ", N_matches: " << triggerMatches.size();
    bool passes_trigger = false;
    if (!(triggerMatches.empty())) {
      for (auto const& trigIndex : triggerMatches) {
        if (triggerHandle->accept(triggerNames.triggerIndex(*trigIndex))) {
          passes_trigger = true;
          break;
        };
      }
    }
    passesTrigger_[patternIndex] = passes_trigger;
  }

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
