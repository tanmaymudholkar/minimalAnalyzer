// -*- C++ -*-
//
// Package:    temp/GenLevelDeltaRAnalyzer
// Class:      GenLevelDeltaRAnalyzer
//
/**\class GenLevelDeltaRAnalyzer GenLevelDeltaRAnalyzer.cc temp/GenLevelDeltaRAnalyzer/plugins/GenLevelDeltaRAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Tanmay Mudholkar
//         Created:  Thu, 20 Feb 2020 02:41:34 GMT
//
//


#include "../interface/GenLevelDeltaRAnalyzer.h"

//
// constructors and destructor
//
GenLevelDeltaRAnalyzer::GenLevelDeltaRAnalyzer(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");
  edm::Service<TFileService> fileService;
  eventProgenitor_ = iConfig.getUntrackedParameter<std::string>("eventProgenitor");
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity");
  edm::LogInfo("GenDeltaRAnalyzer") << "Starting GenLevelDeltaRAnalyzer with verbosity: " << verbosity_ << ", eventProgenitor: " << eventProgenitor_;

  eventInfoTree_ = fileService->make<TTree>("eventInfoTree", "eventInfoTree");
  eventInfoTree_->Branch("nStealthPhotons", &nStealthPhotons_, "nStealthPhotons/I");
  eventInfoTree_->Branch("nKinematicStealthPhotons", &nKinematicStealthPhotons_, "nKinematicStealthPhotons/I");
  eventInfoTree_->Branch("eventProgenitorMass", &eventProgenitorMass_, "eventProgenitorMass/F");
  eventInfoTree_->Branch("neutralinoMass", &neutralinoMass_, "neutralinoMass/F");
  eventInfoTree_->Branch("nGenJets", &nGenJets_, "nGenJets/I");

  deltaRTree_ = fileService->make<TTree>("deltaRTree", "deltaRTree");
  deltaRTree_->Branch("evt_eventProgenitorMass", &evt_eventProgenitorMass_, "evt_eventProgenitorMass/F");
  deltaRTree_->Branch("evt_neutralinoMass", &evt_neutralinoMass_, "evt_neutralinoMass/F");
  deltaRTree_->Branch("deltaR_closestGenJet", &deltaR_closestGenJet_, "deltaR_closestGenJet/F");
  deltaRTree_->Branch("deltaR_secondClosestGenJet", &deltaR_secondClosestGenJet_, "deltaR_secondClosestGenJet/F");
  deltaRTree_->Branch("photonMom_pdgId", &photonMom_pdgId_, "photonMom_pdgId/I");
  deltaRTree_->Branch("photonPT", &photonPT_, "photonPT/F");
  deltaRTree_->Branch("closestGenJet_PT", &closestGenJet_PT_, "closestGenJet_PT/F");
  deltaRTree_->Branch("closestGenJet_fraction_EM", &closestGenJet_fraction_EM_, "closestGenJet_fraction_EM/F");
  deltaRTree_->Branch("closestGenJet_fraction_hadronic", &closestGenJet_fraction_hadronic_, "closestGenJet_fraction_hadronic/F");
  deltaRTree_->Branch("closestGenJet_fraction_invisible", &closestGenJet_fraction_invisible_, "closestGenJet_fraction_invisible/F");
  deltaRTree_->Branch("closestGenJet_fraction_aux", &closestGenJet_fraction_aux_, "closestGenJet_fraction_aux/F");
  deltaRTree_->Branch("closestGenJet_totalFraction", &closestGenJet_totalFraction_, "closestGenJet_totalFraction/F");
  deltaRTree_->Branch("secondClosestGenJet_PT", &secondClosestGenJet_PT_, "secondClosestGenJet_PT/F");
  deltaRTree_->Branch("secondClosestGenJet_fraction_EM", &secondClosestGenJet_fraction_EM_, "secondClosestGenJet_fraction_EM/F");
  deltaRTree_->Branch("secondClosestGenJet_fraction_hadronic", &secondClosestGenJet_fraction_hadronic_, "secondClosestGenJet_fraction_hadronic/F");
  deltaRTree_->Branch("secondClosestGenJet_fraction_invisible", &secondClosestGenJet_fraction_invisible_, "secondClosestGenJet_fraction_invisible/F");
  deltaRTree_->Branch("secondClosestGenJet_fraction_aux", &secondClosestGenJet_fraction_aux_, "secondClosestGenJet_fraction_aux/F");
  deltaRTree_->Branch("secondClosestGenJet_totalFraction", &secondClosestGenJet_totalFraction_, "secondClosestGenJet_totalFraction/F");

  prunedGenParticlesCollection_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
  genJetsCollection_ = consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsSrc"));
}


GenLevelDeltaRAnalyzer::~GenLevelDeltaRAnalyzer()
{
  // outputFile_->WriteTObject(eventInfoTree_);
  // outputFile_->WriteTObject(deltaRTree_);
  // outputFile_->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenLevelDeltaRAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<reco::GenParticle> > prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesCollection_, prunedGenParticlesHandle);
  int nPrunedParticles = (*(prunedGenParticlesHandle.product())).size();
  if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Found " << nPrunedParticles << " prunedParticles.";

  edm::Handle<edm::View<reco::GenJet> > genJetsHandle;
  iEvent.getByToken(genJetsCollection_, genJetsHandle);
  nGenJets_ = (*(genJetsHandle.product())).size();
  if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Found " << nGenJets_ << " jets.";

  // First pass: count stealth photons in EB with PT > 25, and set neutralino and event progenitor masses.
  std::vector<int> kinematicStealthPhotonIndices;
  nStealthPhotons_ = 0;
  nKinematicStealthPhotons_ = 0;
  eventProgenitorMass_ = -1.0;
  neutralinoMass_ = -1.0;
  for (int prunedParticleIndex = 0; prunedParticleIndex < nPrunedParticles; ++prunedParticleIndex) {
    const reco::GenParticle& prunedParticle = (*(prunedGenParticlesHandle.product())).at(prunedParticleIndex);
    if (verbosity_ >= 4) edm::LogInfo("GenDeltaRAnalyzer") << "Found pruned particle at prunedParticleIndex = " << prunedParticleIndex << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi() << ", pdgId: " << prunedParticle.pdgId();
    if (PIDUtils::isPhotonPID(prunedParticle.pdgId())) {
      int prunedParticleMomID = ((prunedParticle.mother(0) == nullptr) ? 0 : prunedParticle.mother(0)->pdgId());
      if (verbosity_ >= 3) edm::LogInfo("GenDeltaRAnalyzer") << "Found photon at prunedParticleIndex = " << prunedParticleIndex << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi() << ", pT = " << prunedParticle.pt() << ", mom ID: " << prunedParticleMomID;
      if (PIDUtils::isNeutralinoPID(prunedParticleMomID)) {
	++nStealthPhotons_;
	assert(prunedParticle.mother(0) != nullptr);
	if (neutralinoMass_ < 0.) neutralinoMass_ = prunedParticle.mother(0)->mass();
	if ((prunedParticle.pt() >= 25.) && ((std::fabs(prunedParticle.eta())) < 1.442)) {
	  ++nKinematicStealthPhotons_;
	  kinematicStealthPhotonIndices.push_back(prunedParticleIndex);
	}
      }
    }
    else if (eventProgenitorMass_ < 0.) {
      if (verbosity_ >= 4) edm::LogInfo("GenDeltaRAnalyzer") << "eventProgenitorMass unset. isSquarkPID: " << ((PIDUtils::isSquarkPID(prunedParticle.pdgId()))? "yes": "no") << "; isGluinoPID: " << ((PIDUtils::isGluinoPID(prunedParticle.pdgId()))? "yes": "no");
      if (((eventProgenitor_ == "squark") && (PIDUtils::isSquarkPID(prunedParticle.pdgId()))) ||
	  ((eventProgenitor_ == "gluino") && (PIDUtils::isGluinoPID(prunedParticle.pdgId())))) {
	if (verbosity_ >= 4) edm::LogInfo("GenDeltaRAnalyzer") << "Setting eventProgenitorMass: " << prunedParticle.mass();
	eventProgenitorMass_ = prunedParticle.mass();
      }
    }
  }
  if (neutralinoMass_ > 0.) assert(eventProgenitorMass_ > 0.);
  eventInfoTree_->Fill();
  assert(nStealthPhotons_ <= 2);

  if (!(nKinematicStealthPhotons_ == 2)) return;
  assert(kinematicStealthPhotonIndices.size() == static_cast<unsigned int>(2));

  evt_eventProgenitorMass_ = eventProgenitorMass_;
  evt_neutralinoMass_ = neutralinoMass_;
  for (const int& photonIndex: kinematicStealthPhotonIndices) {
    const reco::GenParticle& prunedStealthPhoton = (*(prunedGenParticlesHandle.product())).at(photonIndex);

    angularVariablesStruct photonEtaPhi = angularVariablesStruct(prunedStealthPhoton.eta(), prunedStealthPhoton.phi());
    photonMom_pdgId_ = ((prunedStealthPhoton.mother(0) == nullptr) ? 0 : prunedStealthPhoton.mother(0)->pdgId()); // just a sanity check...
    photonPT_ = prunedStealthPhoton.pt();
    deltaR_closestGenJet_ = -0.1;
    deltaR_secondClosestGenJet_ = -0.1;
    int index_closestGenJet = -1;
    int index_secondClosestGenJet = -1;
    for (int genJetIndex = 0; genJetIndex < nGenJets_; ++genJetIndex) {
      const reco::GenJet& genJet = (*(genJetsHandle.product())).at(genJetIndex);
      angularVariablesStruct genJetEtaPhi = angularVariablesStruct(genJet.eta(), genJet.phi());
      float deltaR = photonEtaPhi.get_deltaR(genJetEtaPhi);
      if (verbosity_ >= 3) edm::LogInfo("GenDeltaRAnalyzer") << "Found genJet at genJetIndex = " << genJetIndex << ", eta = " << genJet.eta() << ", phi = " << genJet.phi() << ", deltaR = " << deltaR << ", PT = " << genJet.pt() << ", EM energy = " << genJet.emEnergy() << ", Hadronic energy = " << genJet.hadEnergy();
      if (deltaR_closestGenJet_ < 0.) { // closest deltaR is unset
	assert(deltaR_secondClosestGenJet_ < 0.); // check that second-closest deltaR is unset as well
	deltaR_closestGenJet_ = deltaR; // set closest deltaR to this deltaR
	index_closestGenJet = genJetIndex;
      }
      else if (deltaR_secondClosestGenJet_ < 0.) { // second-closest deltaR is unset
	if (deltaR < deltaR_closestGenJet_) { // if this deltaR is closer than previously set closest deltaR, move previously set closest deltaR to second-closest deltaR
	  deltaR_secondClosestGenJet_ = deltaR_closestGenJet_;
	  index_secondClosestGenJet = index_closestGenJet;
	  deltaR_closestGenJet_ = deltaR;
	  index_closestGenJet = genJetIndex;
	}
	else { // second-closest deltaR is unset but this deltaR is not as close as previously set deltaR
	  deltaR_secondClosestGenJet_ = deltaR;
	  index_secondClosestGenJet = genJetIndex;
	}
      }
      else if (deltaR < deltaR_closestGenJet_) { // both closest and second-closest are set, and this deltaR is closer than previously set closest deltaR
	deltaR_secondClosestGenJet_ = deltaR_closestGenJet_;
	index_secondClosestGenJet = index_closestGenJet;
	deltaR_closestGenJet_ = deltaR;
	index_closestGenJet = genJetIndex;
      }
      else if (deltaR < deltaR_secondClosestGenJet_) { // // both closest and second-closest are set, and this deltaR is closer than previously set second-closest deltaR but not as close as previously set closest deltaR
	deltaR_secondClosestGenJet_ = deltaR;
	index_secondClosestGenJet = genJetIndex;
      }
    } // ends loop over genJets
    closestGenJet_PT_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).pt();
    float closestGenJet_energy = ((*(genJetsHandle.product())).at(index_closestGenJet)).energy();
    closestGenJet_fraction_EM_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).emEnergy()/closestGenJet_energy;
    closestGenJet_fraction_hadronic_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).hadEnergy()/closestGenJet_energy;
    closestGenJet_fraction_invisible_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).invisibleEnergy()/closestGenJet_energy;
    closestGenJet_fraction_aux_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).auxiliaryEnergy()/closestGenJet_energy;
    closestGenJet_totalFraction_ = closestGenJet_fraction_EM_ + closestGenJet_fraction_hadronic_ + closestGenJet_fraction_invisible_ + closestGenJet_fraction_aux_;

    secondClosestGenJet_PT_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).pt();
    float secondClosestGenJet_energy = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).energy();
    secondClosestGenJet_fraction_EM_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).emEnergy()/secondClosestGenJet_energy;
    secondClosestGenJet_fraction_hadronic_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).hadEnergy()/secondClosestGenJet_energy;
    secondClosestGenJet_fraction_invisible_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).invisibleEnergy()/secondClosestGenJet_energy;
    secondClosestGenJet_fraction_aux_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).auxiliaryEnergy()/secondClosestGenJet_energy;
    secondClosestGenJet_totalFraction_ = secondClosestGenJet_fraction_EM_ + secondClosestGenJet_fraction_hadronic_ + secondClosestGenJet_fraction_invisible_ + secondClosestGenJet_fraction_aux_;
    
    if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Closest deltaR = " << deltaR_closestGenJet_ << " at genJetIndex = " << index_closestGenJet;
    if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Second closest deltaR = " << deltaR_secondClosestGenJet_ << " at genJetIndex = " << index_secondClosestGenJet;
    deltaRTree_->Fill();
  } // ends loop over all gen particles
}


// ------------ method called once each job just before starting event loop  ------------
void
GenLevelDeltaRAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenLevelDeltaRAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenLevelDeltaRAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(GenLevelDeltaRAnalyzer);
