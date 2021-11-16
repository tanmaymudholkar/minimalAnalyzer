// system include files
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenParticlesCollection_;
  edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsCollection_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
