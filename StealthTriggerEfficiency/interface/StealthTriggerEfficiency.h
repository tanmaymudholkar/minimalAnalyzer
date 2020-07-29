// system include files
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>

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

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "constants.h"

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

  // ----------member data ---------------------------
  int verbosity_;

  TTree* eventInfoTree_;
  float pT_leadingPhoton_;
  float eta_leadingPhoton_;
  float pT_subLeadingPhoton_;
  float eta_subLeadingPhoton_;
  bool passesSelection_;
  bool passesTrigger_;

  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollection_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
