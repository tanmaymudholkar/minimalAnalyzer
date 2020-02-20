// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "dataStructures.h"
#include "constants.h"
#include "pdgUtils.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class GenLevelDeltaRAnalyzer : public edm::one::EDAnalyzer<>  {
public:
  explicit GenLevelDeltaRAnalyzer(const edm::ParameterSet&);
  ~GenLevelDeltaRAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  std::string outputPath_;
  TFile *outputFile_;
  int verbosity_;

  // TTree* pairs_packedGenParticles_;
  TTree* pairs_prunedGenParticles_;
  // Information to store in these trees
  // float packed_deltaR_;
  // int packed_photonMom_pdgId_;
  // float packed_photonPT_;
  // int packed_otherParticle_pdgId_;
  // int packed_otherParticleMom_pdgId_;
  // float packed_otherParticlePT_;
  float pruned_deltaR_;
  int pruned_photonMom_pdgId_;
  float pruned_photonPT_;
  int pruned_otherParticle_pdgId_;
  float pruned_otherParticlePT_;
  int pruned_otherParticleMom_pdgId_;

  // edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenParticlesCollection_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenParticlesCollection_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
