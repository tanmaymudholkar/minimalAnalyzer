// system include files
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "dataStructures.h"
#include "constants.h"
#include "pdgUtils.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#define DEFAULT_DELTAR -0.1

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class GenLevelDeltaRAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit GenLevelDeltaRAnalyzer(const edm::ParameterSet&);
  ~GenLevelDeltaRAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  std::string daughter_pdgIds(const reco::GenParticle& particle);
  std::string mother_pdgIds(const reco::GenParticle& particle);
  int get_photon_parentage(const reco::Candidate *, const int &, const std::vector<float> &, const std::vector<angularVariablesStruct> &);

  // ----------member data ---------------------------
  /* std::string eventProgenitor_; */
  std::string event_start_end_banner_;
  int verbosity_;

  TTree* eventInfoTree_;
  /* int nStealthPhotons_; */
  /* int nKinematicStealthPhotons_; */
  /* int nEnergeticStealthPhotons_; */
  /* float eventProgenitorMass_; */
  /* float neutralinoMass_; */
  /* std::vector<float> eta_progenitor_; */
  /* std::vector<float> eta_neutralino_progenitor_child_; */
  /* std::vector<float> eta_neutralino_photon_mother_; */
  /* float eta_photon_leading_; */
  /* float eta_photon_subleading_; */
  int nLHEParticles_;
  float genHT_;
  std::vector<int> id_LHEParticle_;
  std::vector<float> px_LHEParticle_;
  std::vector<float> py_LHEParticle_;
  std::vector<float> pz_LHEParticle_;
  std::vector<float> energy_LHEParticle_;
  std::vector<float> pt_LHEParticle_; /* maybe unreliable? */
  std::vector<float> eta_LHEParticle_; /* maybe unreliable? */
  std::vector<float> phi_LHEParticle_; /* maybe unreliable? */

  float deltaR_photonPair_;
  int nKinematicGenJets_;
  int nKinematicGenPhotons_;
  std::vector<float> pt_kinematicGenPhoton_;
  std::vector<float> eta_kinematicGenPhoton_;
  std::vector<float> phi_kinematicGenPhoton_;
  std::vector<float> deltaR_kinematicGenPhoton_closestGenJet_;
  std::vector<float> pt_kinematicGenPhoton_closestGenJet_;
  std::vector<float> eta_kinematicGenPhoton_closestGenJet_;
  std::vector<float> phi_kinematicGenPhoton_closestGenJet_;
  std::vector<float> em_fraction_kinematicGenPhoton_closestGenJet_;
  std::vector<float> deltaR_kinematicGenPhoton_secondClosestGenJet_;
  std::vector<float> pt_kinematicGenPhoton_secondClosestGenJet_;
  std::vector<float> eta_kinematicGenPhoton_secondClosestGenJet_;
  std::vector<float> phi_kinematicGenPhoton_secondClosestGenJet_;
  std::vector<float> em_fraction_kinematicGenPhoton_secondClosestGenJet_;

  bool selection_map_is_available_;
  // convention for selection map:
  // (run, lumi, event) details are ((selection_map_<leading|subLeading>.at(run)).at(lumi)).at(event)
  std::string selection_map_source_;
  std::map<int, std::map<int, std::map<int, PTEtaPhiStruct> > > selection_map_leading_photon_;
  std::map<int, std::map<int, std::map<int, PTEtaPhiStruct> > > selection_map_subLeading_photon_;
  // Definition of parentage:
  // Loop over all final state photons with pt > 25 GeV
  // Go up the chain and find earliest photon that doesn't have a single mother photon
  // Call that one the progenitor
  // 0: unknown
  // 1: progenitor matched (by pt and eta/phi) to an LHE photon
  // 2: progenitor's mother is proton with status = 4
  // 3: progenitor's mother is a quark
  int nTotalFinalStatePhotons_;
  std::vector<int> parentage_finalStatePhoton_;
  int nMatchedFinalStatePhotons_;
  int parentage_leadingPhoton_;
  int parentage_subLeadingPhoton_;

  /* TTree* deltaRTree_; */
  /* float evt_eventProgenitorMass_; */
  /* float evt_neutralinoMass_; */
  /* float deltaR_closestGenJet_; */
  /* float deltaR_secondClosestGenJet_; */
  /* int photonMom_pdgId_; */
  /* float photonPT_; */
  /* float closestGenJet_PT_; */
  /* float closestGenJet_fraction_EM_; */
  /* float closestGenJet_fraction_hadronic_; */
  /* float closestGenJet_fraction_invisible_; */
  /* float closestGenJet_fraction_aux_; */
  /* float closestGenJet_totalFraction_; // for debugging only */
  /* float secondClosestGenJet_PT_; */
  /* float secondClosestGenJet_fraction_EM_; */
  /* float secondClosestGenJet_fraction_hadronic_; */
  /* float secondClosestGenJet_fraction_invisible_; */
  /* float secondClosestGenJet_fraction_aux_; */
  /* float secondClosestGenJet_totalFraction_; // for debugging only */

  edm::EDGetTokenT<LHEEventProduct> lheEventInfoCollection_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenParticlesCollection_;
  edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsCollection_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
