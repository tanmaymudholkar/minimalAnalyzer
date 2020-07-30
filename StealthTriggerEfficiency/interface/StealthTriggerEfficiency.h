// system include files
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <array>

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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "constants.h"
#include "triggers.h"

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

class StealthTriggerEfficiency : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit StealthTriggerEfficiency(const edm::ParameterSet&);
  ~StealthTriggerEfficiency();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  float getEffectiveArea_neutIso(const float&);
  float getEffectiveArea_phoIso(const float&);
  float getEffectiveArea_chIso(const float&);

  // ----------member data ---------------------------
  int verbosity_;

  TTree* eventInfoTree_;
  float eventRho_;
  float pT_leadingPhoton_;
  float eta_leadingPhoton_;
  float phi_leadingPhoton_;
  float energy_leadingPhoton_;
  float pT_subLeadingPhoton_;
  float eta_subLeadingPhoton_;
  float phi_subLeadingPhoton_;
  float energy_subLeadingPhoton_;
  float invariantMass_;
  bool passesSelection_;
  std::array<bool, (triggerPatterns::patternsToSave).size()> passesTrigger_;

  edm::EDGetTokenT<double> rhoCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerCollection_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollection_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
