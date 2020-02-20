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
  outputPath_ = iConfig.getUntrackedParameter<std::string>("outputPath");
  outputFile_ = TFile::Open(outputPath_.c_str(), "RECREATE");
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity");

  // pairs_packedGenParticles_ = new TTree("pairs_packedGenParticles_", "pairs_packedGenParticles_");
  // pairs_packedGenParticles_->Branch("packed_deltaR", &packed_deltaR_, "packed_deltaR/F");
  // pairs_packedGenParticles_->Branch("packed_photonMom_pdgId", &packed_photonMom_pdgId_, "packed_photonMom_pdgId/I");
  // pairs_packedGenParticles_->Branch("packed_photonPT", &packed_photonPT_, "packed_photonPT/F");
  // pairs_packedGenParticles_->Branch("packed_otherParticle_pdgId", &packed_otherParticle_pdgId_, "packed_otherParticle_pdgId/I");
  // pairs_packedGenParticles_->Branch("packed_otherParticleMom_pdgId", &packed_otherParticleMom_pdgId_, "packed_otherParticleMom_pdgId/I");
  // pairs_packedGenParticles_->Branch("packed_otherParticlePT", &packed_otherParticlePT_, "packed_otherParticlePT/F");
  pairs_prunedGenParticles_ = new TTree("pairs_prunedGenParticles_", "pairs_prunedGenParticles_");
  pairs_prunedGenParticles_->Branch("pruned_deltaR", &pruned_deltaR_, "pruned_deltaR/F");
  pairs_prunedGenParticles_->Branch("pruned_photonMom_pdgId", &pruned_photonMom_pdgId_, "pruned_photonMom_pdgId/I");
  pairs_prunedGenParticles_->Branch("pruned_photonPT", &pruned_photonPT_, "pruned_photonPT/F");
  pairs_prunedGenParticles_->Branch("pruned_otherParticle_pdgId", &pruned_otherParticle_pdgId_, "pruned_otherParticle_pdgId/I");
  pairs_prunedGenParticles_->Branch("pruned_otherParticleMom_pdgId", &pruned_otherParticleMom_pdgId_, "pruned_otherParticleMom_pdgId/I");
  pairs_prunedGenParticles_->Branch("pruned_otherParticlePT", &pruned_otherParticlePT_, "pruned_otherParticlePT/F");

  // packedGenParticlesCollection_ = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGenParticlesSrc"));
  prunedGenParticlesCollection_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
}


GenLevelDeltaRAnalyzer::~GenLevelDeltaRAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // outputFile_->WriteTObject(pairs_packedGenParticles_);
  outputFile_->WriteTObject(pairs_prunedGenParticles_);
  outputFile_->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenLevelDeltaRAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // edm::Handle<edm::View<pat::PackedGenParticle> > packedGenParticlesHandle;
  // iEvent.getByToken(packedGenParticlesCollection_, packedGenParticlesHandle);
  // int nPackedParticles = (*(packedGenParticlesHandle.product())).size();
  // // for (const pat::PackedGenParticle& packedParticle: *(packedGenParticlesHandle.product())) {
  // for (int packedParticleIndex = 0; packedParticleIndex < (nPackedParticles-1); ++packedParticleIndex) {
  //   const pat::PackedGenParticle& packedParticle = (*(packedGenParticlesHandle.product()))[packedParticleIndex];
  //   if (verbosity_ >= 2) std::cout << "Found packed particle at eta = " << packedParticle.eta() << ", phi = " << packedParticle.phi() << ", pdgId: " << packedParticle.pdgId() <<  std::endl;
  //   if (PIDUtils::isPhotonPID(packedParticle.pdgId())) {
  //     angularVariablesStruct photonEtaPhi = angularVariablesStruct(packedParticle.eta(), packedParticle.phi());
  //     packed_photonMom_pdgId_ = ((packedParticle.mother(0) == nullptr) ? 0 : packedParticle.mother(0)->pdgId());
  //     packed_photonPT_ = packedParticle.pt();
  //     for (int otherPackedParticleIndex = (1+packedParticleIndex); otherPackedParticleIndex < nPackedParticles; ++otherPackedParticleIndex) {
  // 	const pat::PackedGenParticle& otherPackedParticle = (*(packedGenParticlesHandle.product()))[otherPackedParticleIndex];
  // 	if (PIDUtils::isInterestingPID(otherPackedParticle.pdgId())) {
  // 	  angularVariablesStruct otherParticleEtaPhi = angularVariablesStruct(otherPackedParticle.eta(), otherPackedParticle.phi());
  // 	  packed_deltaR_ = photonEtaPhi.get_deltaR(otherParticleEtaPhi);
  // 	  packed_otherParticle_pdgId_ = otherPackedParticle.pdgId();
  // 	  packed_otherParticleMom_pdgId_ = ((otherPackedParticle.mother(0) == nullptr) ? 0 : otherPackedParticle.mother(0)->pdgId());
  // 	  packed_otherParticlePT_ = otherPackedParticle.pt();
  // 	  pairs_packedGenParticles_->Fill();
  // 	}
  //     }
  //   }
  // }

  edm::Handle<edm::View<reco::GenParticle> > prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesCollection_, prunedGenParticlesHandle);
  int nPrunedParticles = (*(prunedGenParticlesHandle.product())).size();
  // for (const reco::GenParticle& prunedParticle: *(prunedGenParticlesHandle.product())) {
  for (int prunedParticleIndex = 0; prunedParticleIndex < (nPrunedParticles-1); ++prunedParticleIndex) {
    const reco::GenParticle& prunedParticle = (*(prunedGenParticlesHandle.product()))[prunedParticleIndex];
    if (verbosity_ >= 2) std::cout << "Found pruned particle at eta = " << prunedParticle.eta() << ", phi = " << prunedParticle.phi() << ", pdgId: " << prunedParticle.pdgId() << std::endl;
    if (PIDUtils::isPhotonPID(prunedParticle.pdgId())) {
      angularVariablesStruct photonEtaPhi = angularVariablesStruct(prunedParticle.eta(), prunedParticle.phi());
      pruned_photonMom_pdgId_ = ((prunedParticle.mother(0) == nullptr) ? 0 : prunedParticle.mother(0)->pdgId());
      pruned_photonPT_ = prunedParticle.pt();
      for (int otherPrunedParticleIndex = (1+prunedParticleIndex); otherPrunedParticleIndex < nPrunedParticles; ++otherPrunedParticleIndex) {
	const reco::GenParticle& otherPrunedParticle = (*(prunedGenParticlesHandle.product()))[otherPrunedParticleIndex];
	if (PIDUtils::isInterestingPID(otherPrunedParticle.pdgId())) {
	  angularVariablesStruct otherParticleEtaPhi = angularVariablesStruct(otherPrunedParticle.eta(), otherPrunedParticle.phi());
	  pruned_deltaR_ = photonEtaPhi.get_deltaR(otherParticleEtaPhi);
	  pruned_otherParticle_pdgId_ = otherPrunedParticle.pdgId();
	  pruned_otherParticleMom_pdgId_ = ((otherPrunedParticle.mother(0) == nullptr) ? 0 : otherPrunedParticle.mother(0)->pdgId());
	  pruned_otherParticlePT_ = otherPrunedParticle.pt();
	  pairs_prunedGenParticles_->Fill();
	}
      }
    }
  }

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
