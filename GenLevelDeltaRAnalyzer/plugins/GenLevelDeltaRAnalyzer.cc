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
  // eventProgenitor_ = iConfig.getUntrackedParameter<std::string>("eventProgenitor");
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity");
  edm::LogInfo("GenDeltaRAnalyzer") << "Starting GenLevelDeltaRAnalyzer with verbosity: " << verbosity_; // << ", eventProgenitor: " << eventProgenitor_;

  event_start_end_banner_ = "";
  for (int event_start_end_banner_counter = 0; event_start_end_banner_counter <= 200; ++event_start_end_banner_counter) {
    event_start_end_banner_ += "=";
  }

  eventInfoTree_ = fileService->make<TTree>("eventInfoTree", "eventInfoTree");
  // eventInfoTree_->Branch("nStealthPhotons", &nStealthPhotons_, "nStealthPhotons/I");
  // eventInfoTree_->Branch("nKinematicStealthPhotons", &nKinematicStealthPhotons_, "nKinematicStealthPhotons/I");
  // eventInfoTree_->Branch("nEnergeticStealthPhotons", &nEnergeticStealthPhotons_, "nEnergeticStealthPhotons/I");
  // eventInfoTree_->Branch("eventProgenitorMass", &eventProgenitorMass_, "eventProgenitorMass/F");
  // eventInfoTree_->Branch("neutralinoMass", &neutralinoMass_, "neutralinoMass/F");
  // eventInfoTree_->Branch("eta_progenitor", &eta_progenitor_);
  // eventInfoTree_->Branch("eta_neutralino_progenitor_child", &eta_neutralino_progenitor_child_);
  // eventInfoTree_->Branch("eta_neutralino_photon_mother", &eta_neutralino_photon_mother_);
  // eventInfoTree_->Branch("eta_photon_leading", &eta_photon_leading_, "eta_photon_leading/F");
  // eventInfoTree_->Branch("eta_photon_subleading", &eta_photon_subleading_, "eta_photon_subleading/F");
  eventInfoTree_->Branch("deltaR_photonPair", &deltaR_photonPair_, "deltaR_photonPair/F");
  eventInfoTree_->Branch("nKinematicGenJets", &nKinematicGenJets_, "nKinematicGenJets/I");
  eventInfoTree_->Branch("nKinematicGenPhotons", &nKinematicGenPhotons_, "nKinematicGenPhotons/I");
  eventInfoTree_->Branch("pt_kinematicGenPhoton", &pt_kinematicGenPhoton_);
  eventInfoTree_->Branch("eta_kinematicGenPhoton", &eta_kinematicGenPhoton_);
  eventInfoTree_->Branch("phi_kinematicGenPhoton", &phi_kinematicGenPhoton_);
  eventInfoTree_->Branch("deltaR_kinematicGenPhoton_closestGenJet", &deltaR_kinematicGenPhoton_closestGenJet_);
  eventInfoTree_->Branch("pt_kinematicGenPhoton_closestGenJet", &pt_kinematicGenPhoton_closestGenJet_);
  eventInfoTree_->Branch("eta_kinematicGenPhoton_closestGenJet", &eta_kinematicGenPhoton_closestGenJet_);
  eventInfoTree_->Branch("phi_kinematicGenPhoton_closestGenJet", &phi_kinematicGenPhoton_closestGenJet_);
  eventInfoTree_->Branch("em_fraction_kinematicGenPhoton_closestGenJet", &em_fraction_kinematicGenPhoton_closestGenJet_);
  eventInfoTree_->Branch("deltaR_kinematicGenPhoton_secondClosestGenJet", &deltaR_kinematicGenPhoton_secondClosestGenJet_);
  eventInfoTree_->Branch("pt_kinematicGenPhoton_secondClosestGenJet", &pt_kinematicGenPhoton_secondClosestGenJet_);
  eventInfoTree_->Branch("eta_kinematicGenPhoton_secondClosestGenJet", &eta_kinematicGenPhoton_secondClosestGenJet_);
  eventInfoTree_->Branch("phi_kinematicGenPhoton_secondClosestGenJet", &phi_kinematicGenPhoton_secondClosestGenJet_);
  eventInfoTree_->Branch("em_fraction_kinematicGenPhoton_secondClosestGenJet", &em_fraction_kinematicGenPhoton_secondClosestGenJet_);
  eventInfoTree_->Branch("nTotalFinalStatePhotons", &nTotalFinalStatePhotons_, "nTotalFinalStatePhotons/I");
  eventInfoTree_->Branch("parentage_finalStatePhoton", &parentage_finalStatePhoton_);
  eventInfoTree_->Branch("nMatchedFinalStatePhotons", &nMatchedFinalStatePhotons_, "nMatchedFinalStatePhotons/I");
  eventInfoTree_->Branch("parentage_leadingPhoton", &parentage_leadingPhoton_, "parentage_leadingPhoton/I");
  eventInfoTree_->Branch("parentage_subLeadingPhoton", &parentage_subLeadingPhoton_, "parentage_subLeadingPhoton/I");

  // deltaRTree_ = fileService->make<TTree>("deltaRTree", "deltaRTree");
  // deltaRTree_->Branch("evt_eventProgenitorMass", &evt_eventProgenitorMass_, "evt_eventProgenitorMass/F");
  // deltaRTree_->Branch("evt_neutralinoMass", &evt_neutralinoMass_, "evt_neutralinoMass/F");
  // deltaRTree_->Branch("deltaR_closestGenJet", &deltaR_closestGenJet_, "deltaR_closestGenJet/F");
  // deltaRTree_->Branch("deltaR_secondClosestGenJet", &deltaR_secondClosestGenJet_, "deltaR_secondClosestGenJet/F");
  // deltaRTree_->Branch("photonMom_pdgId", &photonMom_pdgId_, "photonMom_pdgId/I");
  // deltaRTree_->Branch("photonPT", &photonPT_, "photonPT/F");
  // deltaRTree_->Branch("closestGenJet_PT", &closestGenJet_PT_, "closestGenJet_PT/F");
  // deltaRTree_->Branch("closestGenJet_fraction_EM", &closestGenJet_fraction_EM_, "closestGenJet_fraction_EM/F");
  // deltaRTree_->Branch("closestGenJet_fraction_hadronic", &closestGenJet_fraction_hadronic_, "closestGenJet_fraction_hadronic/F");
  // deltaRTree_->Branch("closestGenJet_fraction_invisible", &closestGenJet_fraction_invisible_, "closestGenJet_fraction_invisible/F");
  // deltaRTree_->Branch("closestGenJet_fraction_aux", &closestGenJet_fraction_aux_, "closestGenJet_fraction_aux/F");
  // deltaRTree_->Branch("closestGenJet_totalFraction", &closestGenJet_totalFraction_, "closestGenJet_totalFraction/F");
  // deltaRTree_->Branch("secondClosestGenJet_PT", &secondClosestGenJet_PT_, "secondClosestGenJet_PT/F");
  // deltaRTree_->Branch("secondClosestGenJet_fraction_EM", &secondClosestGenJet_fraction_EM_, "secondClosestGenJet_fraction_EM/F");
  // deltaRTree_->Branch("secondClosestGenJet_fraction_hadronic", &secondClosestGenJet_fraction_hadronic_, "secondClosestGenJet_fraction_hadronic/F");
  // deltaRTree_->Branch("secondClosestGenJet_fraction_invisible", &secondClosestGenJet_fraction_invisible_, "secondClosestGenJet_fraction_invisible/F");
  // deltaRTree_->Branch("secondClosestGenJet_fraction_aux", &secondClosestGenJet_fraction_aux_, "secondClosestGenJet_fraction_aux/F");
  // deltaRTree_->Branch("secondClosestGenJet_totalFraction", &secondClosestGenJet_totalFraction_, "secondClosestGenJet_totalFraction/F");

  lheEventInfoCollection_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEInfoSrc"));
  prunedGenParticlesCollection_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
  genJetsCollection_ = consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsSrc"));

  selection_map_is_available_ = iConfig.getUntrackedParameter<bool>("selection_map_is_available");
  selection_map_source_ = "";
  if (selection_map_is_available_) {
    selection_map_source_ = iConfig.getUntrackedParameter<std::string>("selection_map_source");
    assert(selection_map_source_ != "/dev/null");
    std::ifstream selection_map_source_handle(selection_map_source_.c_str());
    assert(selection_map_source_handle.is_open());
    int run_identifier, lumi_identifier, event_identifier;
    double pt_leadingPhoton, eta_leadingPhoton, phi_leadingPhoton, pt_subLeadingPhoton, eta_subLeadingPhoton, phi_subLeadingPhoton;
    while (selection_map_source_handle >> run_identifier >> lumi_identifier >> event_identifier
	   >> pt_leadingPhoton >> eta_leadingPhoton >> phi_leadingPhoton
	   >> pt_subLeadingPhoton >> eta_subLeadingPhoton >> phi_subLeadingPhoton) {
      selection_map_leading_photon_[run_identifier][lumi_identifier][event_identifier] = PTEtaPhiStruct(pt_leadingPhoton, eta_leadingPhoton, phi_leadingPhoton);
      selection_map_subLeading_photon_[run_identifier][lumi_identifier][event_identifier] = PTEtaPhiStruct(pt_subLeadingPhoton, eta_subLeadingPhoton, phi_subLeadingPhoton);
    }
    selection_map_source_handle.close();
  }
}


GenLevelDeltaRAnalyzer::~GenLevelDeltaRAnalyzer()
{
  // outputFile_->WriteTObject(eventInfoTree_);
  // outputFile_->WriteTObject(deltaRTree_);
  // outputFile_->Close();
}

std::string
GenLevelDeltaRAnalyzer::daughter_pdgIds(const reco::GenParticle& particle) {
  std::string out = "";
  for (int daughter_index = 0; daughter_index < static_cast<int>(particle.numberOfDaughters()); ++daughter_index) {
    out += (std::to_string((particle.daughter(daughter_index))->pdgId()) + "  ");
  }
  return out;
}

std::string
GenLevelDeltaRAnalyzer::mother_pdgIds(const reco::GenParticle& particle) {
  std::string out = "";
  for (int mother_index = 0; mother_index < static_cast<int>(particle.numberOfMothers()); ++mother_index) {
    out += (std::to_string((particle.mother(mother_index))->pdgId()) + "  ");
  }
  return out;
}

int
GenLevelDeltaRAnalyzer::get_photon_parentage(const reco::Candidate * final_state_photon, const int & nLHEPhotons, const std::vector<float> & LHEPhotonPTs, const std::vector<angularVariablesStruct> & LHEPhotonAngles) {
  // find progenitor photon
  const reco::Candidate * current_photon = final_state_photon;
  while (true) {
    if (current_photon->numberOfMothers() == 1) {
      if (PIDUtils::isPhotonPID(current_photon->mother(0)->pdgId())) {
	current_photon = current_photon->mother(0);
	continue;
      }
      else {
	// current_photon is the progenitor
	break;
      }
    }
    else {
      // current_photon is the progenitor
      break;
    }
  }
  if (verbosity_ >= 3) {
    edm::LogInfo("GenDeltaRAnalyzer") << "Found progenitor of final state photon at eta = " << current_photon->eta() << ", phi = " << current_photon->phi();
  }

  // check if progenitor matches an LHE photon
  if (nLHEPhotons > 0) {
    angularVariablesStruct current_photon_eta_phi = angularVariablesStruct(current_photon->eta(), current_photon->phi());
    for (int LHEPhotonIndex = 0; LHEPhotonIndex < nLHEPhotons; ++LHEPhotonIndex) {
      float ptratio = (LHEPhotonPTs.at(LHEPhotonIndex))/(current_photon->pt());
      if (std::fabs(ptratio - 1.0) < 0.1) { // pt-matching: the two pts must differ by no more than 10%
	float deltaR_wrt_lhe_photon = current_photon_eta_phi.get_deltaR(LHEPhotonAngles.at(LHEPhotonIndex));
	if (deltaR_wrt_lhe_photon < 0.05) { // deltaR-matching: the deltaR must be less than 0.05
	  if (verbosity_ >= 3) {
	    edm::LogInfo("GenDeltaRAnalyzer") << "Found match with LHE photon.";
	  }
	  return 1;
	}
      }
    }
  }

  // check if progenitor's mother is proton with status = 4
  if (current_photon->numberOfMothers() == 1) {
    if ((PIDUtils::isProtonPID(current_photon->mother(0)->pdgId())) && ((current_photon->mother(0)->status()) == 4)) {
      if (verbosity_ >= 3) {
	edm::LogInfo("GenDeltaRAnalyzer") << "Found match with status 4 proton.";
      }
      return 2;
    }
  }

  // check if progenitor's mother is a g or q
  if (current_photon->numberOfMothers() >= 1) {
    if (PIDUtils::isJetCandidatePID(current_photon->mother(0)->pdgId())) {
      if (verbosity_ >= 3) {
	edm::LogInfo("GenDeltaRAnalyzer") << "Found match with jet.";
      }
      return 3;
    }
  }
  return 0;
}

//
// member functions
//

// ------------ method called for each event  ------------
void
GenLevelDeltaRAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (verbosity_ >= 2) {
    edm::LogInfo("GenDeltaRAnalyzer") << event_start_end_banner_;
    edm::LogInfo("GenDeltaRAnalyzer") << "Begin processing event ID: " << iEvent.id().event() << ", run: " << iEvent.id().run() << ", lumi: " << iEvent.luminosityBlock();
  }

  nLHEParticles_ = 0;
  genHT_ = 0;
  id_LHEParticle_.clear();
  px_LHEParticle_.clear();
  py_LHEParticle_.clear();
  pz_LHEParticle_.clear();
  energy_LHEParticle_.clear();
  pt_LHEParticle_.clear(); /* maybe unreliable? */
  eta_LHEParticle_.clear(); /* maybe unreliable? */
  phi_LHEParticle_.clear(); /* maybe unreliable? */
  int nLHEPhotons = 0;
  std::vector<float> LHEPhotonPTs;
  LHEPhotonPTs.clear();
  std::vector<angularVariablesStruct> LHEPhotonAngles;
  LHEPhotonAngles.clear();

  edm::Handle<LHEEventProduct> LHEEventInfoHandle;
  iEvent.getByToken(lheEventInfoCollection_, LHEEventInfoHandle);
  if (LHEEventInfoHandle.isValid()) {
    const lhef::HEPEUP& lheEventInfo = LHEEventInfoHandle->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEventInfo.PUP;
    size_t nLHEParticlesFromHandle = lheParticles.size();
    if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Found " << nLHEParticlesFromHandle << " particles (of all status values) from LHE handle.";
    for ( size_t lhe_particle_index = 0; lhe_particle_index < nLHEParticlesFromHandle; ++lhe_particle_index ) {
      int status = lheEventInfo.ISTUP[lhe_particle_index];
      if (status == 1) {
        ++nLHEParticles_;
        int lhe_id = lheEventInfo.IDUP[lhe_particle_index];
        id_LHEParticle_.push_back(lhe_id);

        // found from LHE product
        float lhe_px = lheParticles[lhe_particle_index][0];
        float lhe_py = lheParticles[lhe_particle_index][1];
        float lhe_pz = lheParticles[lhe_particle_index][2];
        float lhe_e  = lheParticles[lhe_particle_index][3];
        px_LHEParticle_.push_back(lhe_px);
        py_LHEParticle_.push_back(lhe_py);
        pz_LHEParticle_.push_back(lhe_pz);
        energy_LHEParticle_.push_back(lhe_e);

        // calculated explicitly
        float lhe_pt = std::sqrt(std::pow(lhe_px, 2) + std::pow(lhe_py, 2));
        bool onPlusZSide = (lhe_pz >= 0.);
        float lhe_theta_abs = std::atan(std::fabs(lhe_pt/lhe_pz));
        float lhe_eta_abs = -1.0*std::log(std::tan(0.5*lhe_theta_abs));
        float lhe_eta = (onPlusZSide ? lhe_eta_abs : (-1.0*lhe_eta_abs));
        float lhe_phi = 0.;
        float tmp_atan = std::atan(std::fabs(lhe_py/lhe_px));
        if (lhe_px >= 0.) {
          if (lhe_py >= 0.) lhe_phi = tmp_atan;
          else lhe_phi = -1.0*tmp_atan;
        }
        else {
          if (lhe_py >= 0.) lhe_phi = (constants::PI - tmp_atan);
          else lhe_phi = (-1.0*(constants::PI - tmp_atan));
        }
        genHT_ += lhe_pt;
        pt_LHEParticle_.push_back(lhe_pt);
        eta_LHEParticle_.push_back(lhe_eta);
        phi_LHEParticle_.push_back(lhe_phi);

        if (verbosity_ >= 3) {
          edm::LogInfo("GenDeltaRAnalyzer") << "Found LHE particle: " << PIDUtils::getParticleString(lhe_id) << " at index: " << lhe_particle_index << ", eta = " << lhe_eta << ", phi = " << lhe_phi << ", pt = " << lhe_pt << ", px = " << lhe_px << ", py = " << lhe_py << ", pz = " << lhe_pz << ", energy = " << lhe_e;
        }

        if (PIDUtils::isPhotonPID(lhe_id)) {
          ++nLHEPhotons;
          LHEPhotonPTs.push_back(lhe_pt);
          LHEPhotonAngles.push_back(angularVariablesStruct(lhe_eta, lhe_phi));
        }
      }
    }
    if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "number of LHE particles with status 1 = " << nLHEParticles_ << "; genHT = " << genHT_;
  }
  else {
    edm::LogInfo("GenDeltaRAnalyzer") << "WARNING: LHE info unavailable.";
  }

  edm::Handle<edm::View<reco::GenParticle> > prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesCollection_, prunedGenParticlesHandle);
  int nPrunedParticles = (*(prunedGenParticlesHandle.product())).size();
  if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Found " << nPrunedParticles << " prunedParticles.";

  edm::Handle<edm::View<reco::GenJet> > genJetsHandle;
  iEvent.getByToken(genJetsCollection_, genJetsHandle);
  int nGenJets = (*(genJetsHandle.product())).size();
  if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Found " << nGenJets << " gen jets.";

  deltaR_photonPair_ = DEFAULT_DELTAR;
  nKinematicGenJets_ = 0;
  nKinematicGenPhotons_ = 0;
  pt_kinematicGenPhoton_.clear();
  eta_kinematicGenPhoton_.clear();
  phi_kinematicGenPhoton_.clear();
  deltaR_kinematicGenPhoton_closestGenJet_.clear();
  pt_kinematicGenPhoton_closestGenJet_.clear();
  eta_kinematicGenPhoton_closestGenJet_.clear();
  phi_kinematicGenPhoton_closestGenJet_.clear();
  em_fraction_kinematicGenPhoton_closestGenJet_.clear();
  deltaR_kinematicGenPhoton_secondClosestGenJet_.clear();
  pt_kinematicGenPhoton_secondClosestGenJet_.clear();
  eta_kinematicGenPhoton_secondClosestGenJet_.clear();
  phi_kinematicGenPhoton_secondClosestGenJet_.clear();
  em_fraction_kinematicGenPhoton_secondClosestGenJet_.clear();
  nTotalFinalStatePhotons_ = 0;
  parentage_finalStatePhoton_.clear();
  nMatchedFinalStatePhotons_ = 0;
  parentage_leadingPhoton_ = 0;
  parentage_subLeadingPhoton_ = 0;

  // first get list and angles of gen jets
  std::vector<angularVariablesStruct> genJetAngles;
  std::vector<float> genJetPTs;
  std::vector<float> genJetEMFractions;
  for (int genJetIndex = 0; genJetIndex < nGenJets; ++genJetIndex) {
    const reco::GenJet& genJet = (*(genJetsHandle.product())).at(genJetIndex);
    if ((verbosity_ >= 4) || ((verbosity_ >= 3) && (genJet.pt() >= 30.))) edm::LogInfo("GenDeltaRAnalyzer") << "Found gen jet at genJetIndex = " << genJetIndex << ", eta = " << genJet.eta() << ", phi = " << genJet.phi() << ", pt = " << genJet.pt() << ", EM energy = " << genJet.emEnergy() << ", total energy = " << genJet.energy();
    if (genJet.pt() > 30.) {
      ++nKinematicGenJets_;
      angularVariablesStruct genJetEtaPhi = angularVariablesStruct(genJet.eta(), genJet.phi());
      genJetAngles.push_back(genJetEtaPhi);
      genJetPTs.push_back(genJet.pt());
      genJetEMFractions.push_back((genJet.emEnergy())/(genJet.energy()));
    }
  }
  assert(static_cast<int>(genJetAngles.size()) == nKinematicGenJets_);
  assert(static_cast<int>(genJetPTs.size()) == nKinematicGenJets_);
  assert(static_cast<int>(genJetEMFractions.size()) == nKinematicGenJets_);

  PTEtaPhiStruct ptEtaPhi_leadingPhoton, ptEtaPhi_subLeadingPhoton;
  angularVariablesStruct angle_leadingPhoton, angle_subLeadingPhoton;

  if (selection_map_is_available_) {
    ptEtaPhi_leadingPhoton = ((selection_map_leading_photon_.at((iEvent.id().run()))).at((iEvent.luminosityBlock()))).at((iEvent.id().event()));
    ptEtaPhi_subLeadingPhoton = ((selection_map_subLeading_photon_.at((iEvent.id().run()))).at((iEvent.luminosityBlock()))).at((iEvent.id().event()));
    angle_leadingPhoton = angularVariablesStruct(ptEtaPhi_leadingPhoton.eta, ptEtaPhi_leadingPhoton.phi);
    angle_subLeadingPhoton = angularVariablesStruct(ptEtaPhi_subLeadingPhoton.eta, ptEtaPhi_subLeadingPhoton.phi);
    if (verbosity_ >= 3) {
      edm::LogInfo("GenDeltaRAnalyzer") << "Read information from selection map:";
      edm::LogInfo("GenDeltaRAnalyzer") << "Leading photon details:" << ptEtaPhi_leadingPhoton;
      edm::LogInfo("GenDeltaRAnalyzer") << "Subleading photon details:" << ptEtaPhi_subLeadingPhoton;
    }
  }

  // now get list of kinematic stealth photons, using the list of gen jets to get closest and second closest gen jets
  std::vector<int> kinematicGenPhotonIndices;
  for (int prunedParticleIndex = 0; prunedParticleIndex < nPrunedParticles; ++prunedParticleIndex) {
    const reco::GenParticle& prunedParticle = (*(prunedGenParticlesHandle.product())).at(prunedParticleIndex);
    std::string mothers_description = "(";
    if (prunedParticle.numberOfMothers() > 0) {
      for (size_t mother_index = 0; mother_index < prunedParticle.numberOfMothers(); ++mother_index) {
        const reco::Candidate * prunedParticleMother = prunedParticle.mother(mother_index);
        mothers_description += (PIDUtils::getParticleString(prunedParticleMother->pdgId()) + ": (eta: " + std::to_string(prunedParticleMother->eta()) + ", phi: " + std::to_string(prunedParticleMother->phi()) + ", pt: " + std::to_string(prunedParticleMother->pt()) + ", status: " + std::to_string(prunedParticleMother->status())) + ", vertex: " + ("<" + std::to_string(prunedParticleMother->vx()) + ", " + std::to_string(prunedParticleMother->vy()) + ", " + std::to_string(prunedParticleMother->vz()) + ">");
        int mother_nmothers = prunedParticleMother->numberOfMothers();
        mothers_description += (", mother_nmothers: " + std::to_string(mother_nmothers));
        int mother_ndaughters = prunedParticleMother->numberOfDaughters();
        mothers_description += (", mother_ndaughters: " + std::to_string(mother_ndaughters));
        mothers_description += "), ";
      }
    }
    mothers_description += ")";
    std::string daughters_description = "(";
    if (prunedParticle.numberOfDaughters() > 0) {
      for (size_t daughter_index = 0; daughter_index < prunedParticle.numberOfDaughters(); ++daughter_index) {
        // daughters_description += (PIDUtils::getParticleString((prunedParticle.daughter(daughter_index))->pdgId()) + ": (" + std::to_string((prunedParticle.daughter(daughter_index))->eta()) + ", " + std::to_string((prunedParticle.daughter(daughter_index))->phi()) + "), ");
        daughters_description += (PIDUtils::getParticleString((prunedParticle.daughter(daughter_index))->pdgId()) + ": (eta: " + std::to_string((prunedParticle.daughter(daughter_index))->eta()) + ", phi: " + std::to_string((prunedParticle.daughter(daughter_index))->phi()) + ", pt: " + std::to_string((prunedParticle.daughter(daughter_index))->pt()) + ", status: " + std::to_string((prunedParticle.daughter(daughter_index))->status()) + "), ");
      }
    }
    daughters_description += ")";
    if ((verbosity_ >= 4) || ((verbosity_ >= 3) && (prunedParticle.pt() >= 25.))) edm::LogInfo("GenDeltaRAnalyzer") << "Found pruned particle: " << PIDUtils::getParticleString(prunedParticle.pdgId()) << " at prunedParticleIndex = " << prunedParticleIndex << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi() << ", pt: " << prunedParticle.pt() << ", isPromptFinalState: " << ((prunedParticle.isPromptFinalState()) ? "true" : "false") << ", vertex: " << "<" << prunedParticle.vx() << ", " << prunedParticle.vy() << ", " << prunedParticle.vz() << ">" << ", mothers: " << mothers_description << ", daughters: " << daughters_description;
    
    if (PIDUtils::isPhotonPID(prunedParticle.pdgId())) {
      if ((prunedParticle.pt() >= 25.) &&
          (std::fabs(prunedParticle.eta()) <= 1.442) &&
          (prunedParticle.isPromptFinalState())) { // prompt final state photon passing kinematic cuts
        kinematicGenPhotonIndices.push_back(prunedParticleIndex);
        ++nKinematicGenPhotons_;
        pt_kinematicGenPhoton_.push_back(prunedParticle.pt());
        eta_kinematicGenPhoton_.push_back(prunedParticle.eta());
        phi_kinematicGenPhoton_.push_back(prunedParticle.phi());
        angularVariablesStruct photonEtaPhi = angularVariablesStruct(prunedParticle.eta(), prunedParticle.phi());
        std::pair<indexAndDeltaRStruct, indexAndDeltaRStruct> indices_and_deltaRs_closest_and_second_closest_gen_jet = photonEtaPhi.get_indices_and_deltaRs_closest_and_second_closest_from_list_of_angles(genJetAngles);
        indexAndDeltaRStruct index_deltaR_closest_gen_jet = indices_and_deltaRs_closest_and_second_closest_gen_jet.first;
        int index_closestGenJet = index_deltaR_closest_gen_jet.index;
        if (index_closestGenJet >= 0) {
          float deltaR_closestGenJet = index_deltaR_closest_gen_jet.deltaR;
          deltaR_kinematicGenPhoton_closestGenJet_.push_back(deltaR_closestGenJet);
          pt_kinematicGenPhoton_closestGenJet_.push_back(genJetPTs.at(index_closestGenJet));
          eta_kinematicGenPhoton_closestGenJet_.push_back((genJetAngles.at(index_closestGenJet)).eta);
          phi_kinematicGenPhoton_closestGenJet_.push_back((genJetAngles.at(index_closestGenJet)).phi);
          em_fraction_kinematicGenPhoton_closestGenJet_.push_back(genJetEMFractions.at(index_closestGenJet));
        }
        indexAndDeltaRStruct index_deltaR_second_closest_gen_jet = indices_and_deltaRs_closest_and_second_closest_gen_jet.second;
        int index_secondClosestGenJet = index_deltaR_second_closest_gen_jet.index;
        if (index_secondClosestGenJet > 0) {
          float deltaR_secondClosestGenJet = index_deltaR_second_closest_gen_jet.deltaR;
          deltaR_kinematicGenPhoton_secondClosestGenJet_.push_back(deltaR_secondClosestGenJet);
          pt_kinematicGenPhoton_secondClosestGenJet_.push_back(genJetPTs.at(index_secondClosestGenJet));
          eta_kinematicGenPhoton_secondClosestGenJet_.push_back((genJetAngles.at(index_secondClosestGenJet)).eta);
          phi_kinematicGenPhoton_secondClosestGenJet_.push_back((genJetAngles.at(index_secondClosestGenJet)).phi);
          em_fraction_kinematicGenPhoton_secondClosestGenJet_.push_back(genJetEMFractions.at(index_secondClosestGenJet));
        }
        if (verbosity_ >= 3) {
          edm::LogInfo("GenDeltaRAnalyzer") << "Found kinematic gen photon particle at prunedParticleIndex = " << prunedParticleIndex;
          edm::LogInfo("GenDeltaRAnalyzer") << "Index of closest gen jet: " << index_deltaR_closest_gen_jet.index << ", with deltaR: " << index_deltaR_closest_gen_jet.deltaR;
          edm::LogInfo("GenDeltaRAnalyzer") << "Index of second-closest gen jet: " << index_deltaR_second_closest_gen_jet.index << ", with deltaR: " << index_deltaR_second_closest_gen_jet.deltaR;
        }
      } // ends loop over prompt final state photons
      if ((prunedParticle.pt() >= 25.) &&
	  (std::fabs(prunedParticle.eta()) <= 1.442) &&
	  (prunedParticle.status() == 1)) { // final state photon (not necessarily prompt) passing kinematic cuts
	++nTotalFinalStatePhotons_;
	int parentage_photon = get_photon_parentage(&((*(prunedGenParticlesHandle.product())).at(prunedParticleIndex)), nLHEPhotons, LHEPhotonPTs, LHEPhotonAngles);
	parentage_finalStatePhoton_.push_back(parentage_photon);
	if (selection_map_is_available_) {
	  angularVariablesStruct final_photon_angle(prunedParticle.eta(), prunedParticle.phi());
          float deltaR_wrt_leading_photon = final_photon_angle.get_deltaR(angle_leadingPhoton);
          bool matches_leading = (deltaR_wrt_leading_photon < 0.05);
          if (verbosity_ >= 3) {
            edm::LogInfo("GenDeltaRAnalyzer") << "Found match with leading photon at pt = " << prunedParticle.pt() << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi();
          }
          bool matches_subLeading = false;
          if (!(matches_leading)) {
            float deltaR_wrt_subLeading_photon = final_photon_angle.get_deltaR(angle_subLeadingPhoton);
            matches_subLeading = (deltaR_wrt_subLeading_photon < 0.05);
            if (verbosity_ >= 3) {
              edm::LogInfo("GenDeltaRAnalyzer") << "Found match with subleading photon at pt = " << prunedParticle.pt() << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi();
            }
          }
          if (matches_leading || matches_subLeading) {
	    if (verbosity_ >= 3) {
              edm::LogInfo("GenDeltaRAnalyzer") << "Either leading or subleading photon matches the current gen photon.";
            }
            ++nMatchedFinalStatePhotons_;
            int * parentage_ptr = nullptr;
            if (matches_leading) parentage_ptr = &parentage_leadingPhoton_;
            else if (matches_subLeading) parentage_ptr = &parentage_subLeadingPhoton_;
            assert(parentage_ptr != nullptr);
	    *parentage_ptr = parentage_photon;
          } // ends condition entering parentage finding algorithm if leading or subleading is matched
        } // ends condition entering loop over final state photons only if selection map is available
      } // ends loop over final state photons
    } // ends loop over all gen photons
  } // ends loop over all gen particles
  if (verbosity_ >= 2) {
    edm::LogInfo("GenDeltaRAnalyzer") << "parentage_leadingPhoton: " << parentage_leadingPhoton_;
    edm::LogInfo("GenDeltaRAnalyzer") << "parentage_subLeadingPhoton: " << parentage_subLeadingPhoton_;
  }

  eventInfoTree_->Fill();
  if (verbosity_ >= 2) {
    edm::LogInfo("GenDeltaRAnalyzer") << "Finish processing event ID: " << iEvent.id().event() << ", run: " << iEvent.id().run() << ", lumi: " << iEvent.luminosityBlock();
    edm::LogInfo("GenDeltaRAnalyzer") << event_start_end_banner_;
  }

  // // First pass: count stealth photons in EB with PT > 25, and set neutralino and event progenitor masses.
  // std::vector<int> kinematicStealthPhotonIndices;
  // std::vector<int> energeticStealthPhotonIndices;
  // std::vector<int> indices_progenitor;
  // std::vector<int> indices_neutralino_progenitor_child;
  // std::vector<int> indices_neutralino_photon_mother;
  // nStealthPhotons_ = 0;
  // nKinematicStealthPhotons_ = 0;
  // nEnergeticStealthPhotons_ = 0;
  // eventProgenitorMass_ = -1.0;
  // neutralinoMass_ = -1.0;
  // for (int prunedParticleIndex = 0; prunedParticleIndex < nPrunedParticles; ++prunedParticleIndex) {
  //   const reco::GenParticle& prunedParticle = (*(prunedGenParticlesHandle.product())).at(prunedParticleIndex);
  //   if (PIDUtils::isNeutralinoPID(prunedParticle.pdgId())) {
  //     // check if photon is among daughters
  //     if (prunedParticle.numberOfDaughters() >= 1) {
  //       bool break_flag = false;
  //       for (int daughter_index = 0; daughter_index < static_cast<int>(prunedParticle.numberOfDaughters()); ++daughter_index) {
  //         if (break_flag) break;
  //         int daughter_pdgid = (prunedParticle.daughter(daughter_index))->pdgId();
  //         if (PIDUtils::isPhotonPID(daughter_pdgid)) {
  //           // this neutralino is the mom of a final state photon
  //        indices_neutralino_photon_mother.push_back(prunedParticleIndex);
  //        break_flag = true;
  //         }
  //       }
  //     }
  //     // check if progenitor is among mothers
  //     if (prunedParticle.numberOfMothers() >= 1) {
  //       bool break_flag = false;
  //       for (int mother_index = 0; mother_index < static_cast<int>(prunedParticle.numberOfMothers()); ++mother_index) {
  //         if (break_flag) break;
  //         int mother_pdgid = (prunedParticle.mother(mother_index))->pdgId();
  //         if (((eventProgenitor_ == "squark") && (PIDUtils::isSquarkPID(mother_pdgid))) ||
  //             ((eventProgenitor_ == "gluino") && (PIDUtils::isGluinoPID(mother_pdgid)))) {
  //        // this neutralino is the daughter of a progenitor
  //           indices_neutralino_progenitor_child.push_back(prunedParticleIndex);
  //           break_flag = true;
  //         }
  //       }
  //     }
  //   }
  //   else if (((eventProgenitor_ == "squark") && (PIDUtils::isSquarkPID(prunedParticle.pdgId()))) ||
  //            ((eventProgenitor_ == "gluino") && (PIDUtils::isGluinoPID(prunedParticle.pdgId())))) {
  //     // check if neutralino is among daughters
  //     if (prunedParticle.numberOfDaughters() >= 1) {
  //       bool break_flag = false;
  //       for (int daughter_index = 0; daughter_index < static_cast<int>(prunedParticle.numberOfDaughters()); ++daughter_index) {
  //         if (break_flag) break;
  //         int daughter_pdgid = (prunedParticle.daughter(daughter_index))->pdgId();
  //         if (PIDUtils::isNeutralinoPID(daughter_pdgid)) {
  //        // this progenitor is a mother of a neutralino
  //           indices_progenitor.push_back(prunedParticleIndex);
  //        break_flag = true;
  //         }
  //       }
  //     }
  //   }

  //   if (verbosity_ >= 4) edm::LogInfo("GenDeltaRAnalyzer") << "Found pruned particle at prunedParticleIndex = " << prunedParticleIndex << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi() << ", pdgId: " << prunedParticle.pdgId();
  //   if (PIDUtils::isPhotonPID(prunedParticle.pdgId())) {
  //     int prunedParticleMomID = ((prunedParticle.mother(0) == nullptr) ? 0 : prunedParticle.mother(0)->pdgId());
  //     if (verbosity_ >= 3) edm::LogInfo("GenDeltaRAnalyzer") << "Found photon at prunedParticleIndex = " << prunedParticleIndex << ", eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi() << ", pT = " << prunedParticle.pt() << ", mom ID: " << prunedParticleMomID;
  //     if (PIDUtils::isNeutralinoPID(prunedParticleMomID)) {
  //       ++nStealthPhotons_;
  //       assert(prunedParticle.mother(0) != nullptr);
  //       if (neutralinoMass_ < 0.) neutralinoMass_ = prunedParticle.mother(0)->mass();
  //       if ((prunedParticle.pt() >= 25.) && ((std::fabs(prunedParticle.eta())) < 1.442)) {
  //         ++nKinematicStealthPhotons_;
  //         kinematicStealthPhotonIndices.push_back(prunedParticleIndex);
  //       }
  //       if (prunedParticle.pt() >= 25.) {
  //         ++nEnergeticStealthPhotons_;
  //         energeticStealthPhotonIndices.push_back(prunedParticleIndex);
  //       }
  //     }
  //   }
  //   else if (eventProgenitorMass_ < 0.) {
  //     if (verbosity_ >= 4) edm::LogInfo("GenDeltaRAnalyzer") << "eventProgenitorMass unset. isSquarkPID: " << ((PIDUtils::isSquarkPID(prunedParticle.pdgId()))? "yes": "no") << "; isGluinoPID: " << ((PIDUtils::isGluinoPID(prunedParticle.pdgId()))? "yes": "no");
  //     if (((eventProgenitor_ == "squark") && (PIDUtils::isSquarkPID(prunedParticle.pdgId()))) ||
  //         ((eventProgenitor_ == "gluino") && (PIDUtils::isGluinoPID(prunedParticle.pdgId())))) {
  //       if (verbosity_ >= 4) edm::LogInfo("GenDeltaRAnalyzer") << "Setting eventProgenitorMass: " << prunedParticle.mass();
  //       eventProgenitorMass_ = prunedParticle.mass();
  //     }
  //   }
  // }
  // if (neutralinoMass_ > 0.) assert(eventProgenitorMass_ > 0.);
  // deltaR_photonPair_ = DEFAULT_DELTAR;
  // if (nKinematicStealthPhotons_ == 2) {
  //   int index_photon1 = kinematicStealthPhotonIndices.at(0);
  //   int index_photon2 = kinematicStealthPhotonIndices.at(1);
  //   const reco::GenParticle& photon1 = (*(prunedGenParticlesHandle.product())).at(index_photon1);
  //   const reco::GenParticle& photon2 = (*(prunedGenParticlesHandle.product())).at(index_photon2);
  //   angularVariablesStruct photon1_EtaPhi = angularVariablesStruct(photon1.eta(), photon1.phi());
  //   angularVariablesStruct photon2_EtaPhi = angularVariablesStruct(photon2.eta(), photon2.phi());
  //   deltaR_photonPair_ = photon1_EtaPhi.get_deltaR(photon2_EtaPhi);
  // }

  // eta_photon_leading_ = -99.0;
  // eta_photon_subleading_ = -99.0;
  // eta_progenitor_.clear();
  // eta_neutralino_progenitor_child_.clear();
  // eta_neutralino_photon_mother_.clear();
  // if (nEnergeticStealthPhotons_ == 2) {
  //   int index_photon1 = energeticStealthPhotonIndices.at(0);
  //   int index_photon2 = energeticStealthPhotonIndices.at(1);
  //   const reco::GenParticle& photon1 = (*(prunedGenParticlesHandle.product())).at(index_photon1);
  //   const reco::GenParticle& photon2 = (*(prunedGenParticlesHandle.product())).at(index_photon2);
  //   float ET_photon1 = photon1.pt();
  //   float ET_photon2 = photon2.pt();
  //   if (ET_photon1 > ET_photon2) {
  //     eta_photon_leading_ = photon1.eta();
  //     eta_photon_subleading_ = photon2.eta();
  //   }
  //   else {
  //     eta_photon_leading_ = photon2.eta();
  //     eta_photon_subleading_ = photon1.eta();
  //   }
  //   assert(indices_progenitor.size() == 2);
  //   for (const int & progenitor_index : indices_progenitor) {
  //     const reco::GenParticle& progenitor_handle = (*(prunedGenParticlesHandle.product())).at(progenitor_index);
  //     eta_progenitor_.push_back(progenitor_handle.eta());
  //   }
  //   assert(indices_neutralino_progenitor_child.size() == 2);
  //   for (const int & neutralino_progenitor_child_index : indices_neutralino_progenitor_child) {
  //     const reco::GenParticle& neutralino_progenitor_child_handle = (*(prunedGenParticlesHandle.product())).at(neutralino_progenitor_child_index);
  //     eta_neutralino_progenitor_child_.push_back(neutralino_progenitor_child_handle.eta());
  //   }
  //   assert(indices_neutralino_photon_mother.size() == 2);
  //   for (const int & neutralino_photon_mother_index : indices_neutralino_photon_mother) {
  //     const reco::GenParticle& neutralino_photon_mother_handle = (*(prunedGenParticlesHandle.product())).at(neutralino_photon_mother_index);
  //     eta_neutralino_photon_mother_.push_back(neutralino_photon_mother_handle.eta());
  //   }
  // }
  // eventInfoTree_->Fill();
  // assert(nStealthPhotons_ <= 2);

  // if (!(nKinematicStealthPhotons_ == 2)) return;
  // assert(kinematicStealthPhotonIndices.size() == static_cast<unsigned int>(2));

  // evt_eventProgenitorMass_ = eventProgenitorMass_;
  // evt_neutralinoMass_ = neutralinoMass_;
  // for (const int& photonIndex: kinematicStealthPhotonIndices) {
  //   const reco::GenParticle& prunedStealthPhoton = (*(prunedGenParticlesHandle.product())).at(photonIndex);

  //   angularVariablesStruct photonEtaPhi = angularVariablesStruct(prunedStealthPhoton.eta(), prunedStealthPhoton.phi());
  //   photonMom_pdgId_ = ((prunedStealthPhoton.mother(0) == nullptr) ? 0 : prunedStealthPhoton.mother(0)->pdgId()); // just a sanity check...
  //   photonPT_ = prunedStealthPhoton.pt();
  //   deltaR_closestGenJet_ = DEFAULT_DELTAR;
  //   deltaR_secondClosestGenJet_ = DEFAULT_DELTAR;
  //   int index_closestGenJet = -1;
  //   int index_secondClosestGenJet = -1;
  //   for (int genJetIndex = 0; genJetIndex < nGenJets_; ++genJetIndex) {
  //     const reco::GenJet& genJet = (*(genJetsHandle.product())).at(genJetIndex);
  //     angularVariablesStruct genJetEtaPhi = angularVariablesStruct(genJet.eta(), genJet.phi());
  //     float deltaR = photonEtaPhi.get_deltaR(genJetEtaPhi);
  //     if (verbosity_ >= 3) edm::LogInfo("GenDeltaRAnalyzer") << "Found genJet at genJetIndex = " << genJetIndex << ", eta = " << genJet.eta() << ", phi = " << genJet.phi() << ", deltaR = " << deltaR << ", PT = " << genJet.pt() << ", EM energy = " << genJet.emEnergy() << ", Hadronic energy = " << genJet.hadEnergy();
  //     if (deltaR_closestGenJet_ < 0.) { // closest deltaR is unset
  //       assert(deltaR_secondClosestGenJet_ < 0.); // check that second-closest deltaR is unset as well
  //       deltaR_closestGenJet_ = deltaR; // set closest deltaR to this deltaR
  //       index_closestGenJet = genJetIndex;
  //     }
  //     else if (deltaR_secondClosestGenJet_ < 0.) { // second-closest deltaR is unset
  //       if (deltaR < deltaR_closestGenJet_) { // if this deltaR is closer than previously set closest deltaR, move previously set closest deltaR to second-closest deltaR
  //         deltaR_secondClosestGenJet_ = deltaR_closestGenJet_;
  //         index_secondClosestGenJet = index_closestGenJet;
  //         deltaR_closestGenJet_ = deltaR;
  //         index_closestGenJet = genJetIndex;
  //       }
  //       else { // second-closest deltaR is unset but this deltaR is not as close as previously set deltaR
  //         deltaR_secondClosestGenJet_ = deltaR;
  //         index_secondClosestGenJet = genJetIndex;
  //       }
  //     }
  //     else if (deltaR < deltaR_closestGenJet_) { // both closest and second-closest are set, and this deltaR is closer than previously set closest deltaR
  //       deltaR_secondClosestGenJet_ = deltaR_closestGenJet_;
  //       index_secondClosestGenJet = index_closestGenJet;
  //       deltaR_closestGenJet_ = deltaR;
  //       index_closestGenJet = genJetIndex;
  //     }
  //     else if (deltaR < deltaR_secondClosestGenJet_) { // // both closest and second-closest are set, and this deltaR is closer than previously set second-closest deltaR but not as close as previously set closest deltaR
  //       deltaR_secondClosestGenJet_ = deltaR;
  //       index_secondClosestGenJet = genJetIndex;
  //     }
  //   } // ends loop over genJets
  //   closestGenJet_PT_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).pt();
  //   float closestGenJet_energy = ((*(genJetsHandle.product())).at(index_closestGenJet)).energy();
  //   closestGenJet_fraction_EM_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).emEnergy()/closestGenJet_energy;
  //   closestGenJet_fraction_hadronic_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).hadEnergy()/closestGenJet_energy;
  //   closestGenJet_fraction_invisible_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).invisibleEnergy()/closestGenJet_energy;
  //   closestGenJet_fraction_aux_ = ((*(genJetsHandle.product())).at(index_closestGenJet)).auxiliaryEnergy()/closestGenJet_energy;
  //   closestGenJet_totalFraction_ = closestGenJet_fraction_EM_ + closestGenJet_fraction_hadronic_ + closestGenJet_fraction_invisible_ + closestGenJet_fraction_aux_;

  //   secondClosestGenJet_PT_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).pt();
  //   float secondClosestGenJet_energy = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).energy();
  //   secondClosestGenJet_fraction_EM_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).emEnergy()/secondClosestGenJet_energy;
  //   secondClosestGenJet_fraction_hadronic_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).hadEnergy()/secondClosestGenJet_energy;
  //   secondClosestGenJet_fraction_invisible_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).invisibleEnergy()/secondClosestGenJet_energy;
  //   secondClosestGenJet_fraction_aux_ = ((*(genJetsHandle.product())).at(index_secondClosestGenJet)).auxiliaryEnergy()/secondClosestGenJet_energy;
  //   secondClosestGenJet_totalFraction_ = secondClosestGenJet_fraction_EM_ + secondClosestGenJet_fraction_hadronic_ + secondClosestGenJet_fraction_invisible_ + secondClosestGenJet_fraction_aux_;
    
  //   if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Closest deltaR = " << deltaR_closestGenJet_ << " at genJetIndex = " << index_closestGenJet;
  //   if (verbosity_ >= 2) edm::LogInfo("GenDeltaRAnalyzer") << "Second closest deltaR = " << deltaR_secondClosestGenJet_ << " at genJetIndex = " << index_secondClosestGenJet;
  //   deltaRTree_->Fill();
  // } // ends loop over all gen particles
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
