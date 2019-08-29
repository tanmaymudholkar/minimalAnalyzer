// -*- C++ -*-
//
// Package:    temp/MinimalMiniAODAnalyzer
// Class:      MinimalMiniAODAnalyzer
// 
/**\class MinimalMiniAODAnalyzer MinimalMiniAODAnalyzer.cc temp/MinimalMiniAODAnalyzer/plugins/MinimalMiniAODAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Tanmay Mudholkar
//         Created:  Mon, 15 Jul 2019 13:25:52 GMT
//
//


// system include files
#include <memory>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <utility>
#include <bitset>
#include <string>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MinimalMiniAODAnalyzer : public edm::EDAnalyzer  {
public:
  explicit MinimalMiniAODAnalyzer(const edm::ParameterSet&);
  ~MinimalMiniAODAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  edm::EDGetTokenT<double> rhoLabel_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollection_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesCollection_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  std::string outputPath_;
  TFile *outputFile_;
  int verbosity_;
  std::string filterType_;
  int desired_parent_pid_ = -1;
  std::vector<std::string> photonIDCriteria_ = {"hOverE", "sigmaIEtaIEta", "chIso", "neutIso", "phoIso"};
  std::vector<std::vector<std::string> > stepByStepSequences_ = {
    {"chIso", "hOverE", "sigmaIEtaIEta", "neutIso", "phoIso"},
    {"hOverE", "chIso", "sigmaIEtaIEta", "neutIso", "phoIso"},
    {"hOverE", "sigmaIEtaIEta", "chIso", "neutIso", "phoIso"},
    {"hOverE", "sigmaIEtaIEta", "neutIso", "chIso", "phoIso"},
    {"hOverE", "sigmaIEtaIEta", "neutIso", "phoIso", "chIso"},
    {"chIso", "phoIso", "hOverE", "sigmaIEtaIEta", "neutIso"},
    {"phoIso", "chIso", "hOverE", "sigmaIEtaIEta", "neutIso"}
  };
  std::map<std::string, int> nHistBins_ = {
    {"hOverE", 500},
    {"sigmaIEtaIEta", 500},
    {"chIso", 1000},
    {"neutIso", 500},
    {"phoIso", 500}
  };
  std::map<std::string, float> lowerHistLimits_ = {
    {"hOverE", 0.},
    {"sigmaIEtaIEta", 0.},
    {"chIso", 0.},
    {"neutIso", 0.},
    {"phoIso", 0.}
  };
  std::map<std::string, float> upperHistLimits_ = {
    {"hOverE", 0.2},
    {"sigmaIEtaIEta", 0.025},
    {"chIso", 60.},
    {"neutIso", 10.},
    {"phoIso", 10.}
  };

  std::map<std::string, int> nHistBins_HLT_ = {
    {"R9_leading", 500},
    {"hOverE_leading", 500},
    {"sigmaIEtaIEta_leading", 500},
    {"clusIso_leading", 1000},
    {"trkIso_leading", 500},
    {"R9_subLeading", 500},
    {"hOverE_subLeading", 500},
    {"sigmaIEtaIEta_subLeading", 500},
    {"clusIso_subLeading", 1000},
    {"trkIso_subLeading", 500}
  };
  std::map<std::string, float> lowerHistLimits_HLT_ = {
    {"R9_leading", 0.},
    {"hOverE_leading", 0.},
    {"sigmaIEtaIEta_leading", 0.},
    {"clusIso_leading", 0.},
    {"trkIso_leading", 0.},
    {"R9_subLeading", 0.},
    {"hOverE_subLeading", 0.},
    {"sigmaIEtaIEta_subLeading", 0.},
    {"clusIso_subLeading", 0.},
    {"trkIso_subLeading", 0.}
  };
  std::map<std::string, float> upperHistLimits_HLT_ = {
    {"R9_leading", 1.},
    {"hOverE_leading", 0.2},
    {"sigmaIEtaIEta_leading", 0.025},
    {"clusIso_leading", 50.},
    {"trkIso_leading", 50.},
    {"R9_subLeading", 1.},
    {"hOverE_subLeading", 0.2},
    {"sigmaIEtaIEta_subLeading", 0.025},
    {"clusIso_subLeading", 50.},
    {"trkIso_subLeading", 50.}
  };

  std::vector<std::string> diphotonHLTCriteria_ = {"R9_leading", "hOverE_leading", "sigmaIEtaIEta_leading", "clusIso_leading", "trkIso_leading", "R9_subLeading", "hOverE_subLeading", "sigmaIEtaIEta_subLeading", "clusIso_subLeading", "trkIso_subLeading"};
  std::vector<double> efficiencyBinEdges_;
  std::map<std::string, TH1F*> h_global1D_TruthMatched_;
  std::map<std::string, TH1F*> h_global1D_HLT_TruthMatched_;
  std::map<std::string, TEfficiency*> h_global1D_ETEfficiencies_TruthMatched_;
  std::map<unsigned int, std::map<unsigned int, TH1F*> > h_stepByStep_TruthMatched_; // Convention: h_stepByStep_TruthMatched_[sequenceIndex][stepIndex]
  std::map<std::string, TH1F*> h_NMinus1_;
  std::map<std::string, TH1F*> h_NMinus1_TruthMatched_;
  std::map<std::string, TH1F*> h_NMinus1_HLT_TruthMatched_;
  std::map<std::string, TEfficiency*> h_NMinus1_ETEfficiencies_TruthMatched_;
  TEfficiency* h_HLTEfficiency_TruthMatched_;
  std::map<std::string, std::map<std::string, TH2F*> > h_global2D_TruthMatched_;
  std::map<std::string, std::map<std::string, TH2F*> > h_NMinus2_TruthMatched_;
  TH1F* h_nMediumPhotons_;
  TH1F* h_nLoosePhotons_;
  TH1F* h_nMediumPhotons_TruthMatched_;
  TH1F* h_nLoosePhotons_TruthMatched_;
  TH1F* h_selectionRegion_;
  TH1F* h_photonType_;
  TH1F* h_nPU_;
  TH1F* h_rho_;
  TH1F* h_chIso_raw_TruthMatched_;
  TH1F* h_neutIso_raw_TruthMatched_;
  TH1F* h_phoIso_raw_TruthMatched_;
  TH1F* h_phoET_TruthMatched_;
  TH1F* h_phoET_passingID_TruthMatched_;
  TEfficiency* h_overallETEfficiency_TruthMatched_;
  TH2F* h_mediumFakeCriteria_;
  TH2F* h_mediumFakeCriteria_TruthMatched_;
  enum class PFTypeForEA{chargedHadron=0, neutralHadron, photon, nPFTypesForEA};
  std::map<PFTypeForEA, std::string> PFTypeNames_ = {
    {PFTypeForEA::chargedHadron, std::string("chargedHadron")},
    {PFTypeForEA::neutralHadron, std::string("neutralHadron")},
    {PFTypeForEA::photon, std::string("photon")}
  };
  std::map<PFTypeForEA, float> EAValues_barrelCenter_ = {
    {PFTypeForEA::chargedHadron, 0.0112},
    {PFTypeForEA::neutralHadron, 0.0668},
    {PFTypeForEA::photon, 0.1113}
  };
  std::map<PFTypeForEA, float> EAValues_barrelEdge_ = {
    {PFTypeForEA::chargedHadron, 0.0108},
    {PFTypeForEA::neutralHadron, 0.1054},
    {PFTypeForEA::photon, 0.0953}
  };
  float getRhoCorrectedIsolation(const float& absEta, const PFTypeForEA& type, const float& isolation, const double& rho);
  bool passesTruthBasedSelection(const edm::Event& iEvent, std::vector<std::pair<float, float> >& truthPhotonsEtaPhi);
  float get_deltaR(const float& source_eta, const float& source_phi, const float& target_eta, const float& target_phi);
  float getMinDeltaR(const float& eta, const float& phi, std::vector<std::pair<float, float> >& targetEtaPhiList);
  void fillGlobal1DAndNMinus1Histograms(const std::map<std::string, bool>& photonIDBits, const std::map<std::string, float>& photonProperties, const bool& isTruthMatched, const float& ET);
  void fillGlobal1DAndNMinus1HLTHistograms(const std::map<std::string, bool>& eventHLTBits, const std::map<std::string, float>& eventHLTProperties, const float& ET_leading, const float& ET_subLeading);
  void fillGlobal2DAndNMinus2Histograms(const std::map<std::string, bool>& photonIDBits, const std::map<std::string, float>& photonProperties);
  void fillStepByStepHistograms(const std::map<std::string, bool>& photonIDBits, const std::map<std::string, float>& photonProperties);
  template<typename valueType> void checkMapKeysAgainstVector(const std::map<std::string, valueType>& mapToCheck, const std::vector<std::string>& allowedValues);

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MinimalMiniAODAnalyzer::MinimalMiniAODAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  checkMapKeysAgainstVector(nHistBins_, photonIDCriteria_);
  checkMapKeysAgainstVector(lowerHistLimits_, photonIDCriteria_);
  checkMapKeysAgainstVector(upperHistLimits_, photonIDCriteria_);

  for (int i = 0; i <= 12; ++i) {
    efficiencyBinEdges_.push_back(1.0*(10+i*20)); // Bins up to 250 GeV are binned in 20 GeV; first edge at 10 GeV, last edge at 250 GeV
  }
  for (int i = 0; i <= 3; ++i) {
    efficiencyBinEdges_.push_back(1.0*(300+i*50)); // Bins up to 450 GeV are binned in 50 GeV; first edge at 300 GeV, last edge at 450 GeV
  }
  efficiencyBinEdges_.push_back(600.);
  efficiencyBinEdges_.push_back(1000.);

  for (const auto& criterion: photonIDCriteria_) {
    h_global1D_TruthMatched_[criterion] = new TH1F((criterion + "_global_TruthMatched").c_str(), (criterion + "_global_TruthMatched;" + criterion).c_str(), nHistBins_.at(criterion), lowerHistLimits_.at(criterion), upperHistLimits_.at(criterion));
    h_global1D_TruthMatched_[criterion]->StatOverflows(kTRUE);
    h_global1D_ETEfficiencies_TruthMatched_[criterion] = new TEfficiency((criterion + "_ETEfficiency_global_TruthMatched").c_str(), (criterion + "_ETEfficiency_global_TruthMatched;photon ET;#epsilon").c_str(), (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]));
    h_NMinus1_[criterion] = new TH1F((criterion + "_NMinus1").c_str(), (criterion + "_NMinus1;" + criterion).c_str(), nHistBins_.at(criterion), lowerHistLimits_.at(criterion), upperHistLimits_.at(criterion));
    h_NMinus1_[criterion]->StatOverflows(kTRUE);
    h_NMinus1_TruthMatched_[criterion] = new TH1F((criterion + "_NMinus1_TruthMatched").c_str(), (criterion + "_NMinus1_TruthMatched;" + criterion).c_str(), nHistBins_.at(criterion), lowerHistLimits_.at(criterion), upperHistLimits_.at(criterion));
    h_NMinus1_TruthMatched_[criterion]->StatOverflows(kTRUE);
    h_NMinus1_ETEfficiencies_TruthMatched_[criterion] = new TEfficiency((criterion + "_ETEfficiency_NMinus1_TruthMatched").c_str(), (criterion + "_ETEfficiency_NMinus1_TruthMatched;photon ET;#epsilon").c_str(), (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]));
  }

  for (const auto& criterion: diphotonHLTCriteria_) {
    h_global1D_HLT_TruthMatched_[criterion] = new TH1F((criterion + "_global_HLT_TruthMatched_").c_str(), (criterion + "_global_TruthMatched_;" + criterion).c_str(), nHistBins_HLT_.at(criterion), lowerHistLimits_HLT_.at(criterion), upperHistLimits_HLT_.at(criterion));
    h_global1D_HLT_TruthMatched_[criterion]->StatOverflows(kTRUE);
    h_NMinus1_HLT_TruthMatched_[criterion] = new TH1F((criterion + "_NMinus1_HLT_TruthMatched_").c_str(), (criterion + "_NMinus1_TruthMatched_;" + criterion).c_str(), nHistBins_HLT_.at(criterion), lowerHistLimits_HLT_.at(criterion), upperHistLimits_HLT_.at(criterion));
    h_NMinus1_HLT_TruthMatched_[criterion]->StatOverflows(kTRUE);
  }
  h_HLTEfficiency_TruthMatched_ = new TEfficiency("HLTEfficiency_TruthMatched", "HLTEfficiency_TruthMatched;ET_leading;ET_subLeading", (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]), (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]));

  for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria_.size()); ++criterion1Index) {
    std::string& criterion1 = photonIDCriteria_.at(criterion1Index);
    for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria_.size(); ++criterion2Index) {
      std::string& criterion2 = photonIDCriteria_.at(criterion2Index);
      std::string prefix2D = criterion1 + "_" + criterion2;
      h_NMinus2_TruthMatched_[criterion1][criterion2] = new TH2F((prefix2D + "_NMinus2_TruthMatched").c_str(), (prefix2D + "_NMinus2_TruthMatched").c_str(), nHistBins_.at(criterion1), lowerHistLimits_.at(criterion1), upperHistLimits_.at(criterion1), nHistBins_.at(criterion2), lowerHistLimits_.at(criterion2), upperHistLimits_.at(criterion2));
      h_NMinus2_TruthMatched_[criterion1][criterion2]->StatOverflows(kTRUE);
      h_NMinus2_TruthMatched_[criterion1][criterion2]->GetXaxis()->SetTitle(criterion1.c_str());
      h_NMinus2_TruthMatched_[criterion1][criterion2]->GetYaxis()->SetTitle(criterion2.c_str());
      h_global2D_TruthMatched_[criterion1][criterion2] = new TH2F((prefix2D + "_global2D_TruthMatched").c_str(), (prefix2D + "_global2D_TruthMatched").c_str(), nHistBins_.at(criterion1), lowerHistLimits_.at(criterion1), upperHistLimits_.at(criterion1), nHistBins_.at(criterion2), lowerHistLimits_.at(criterion2), upperHistLimits_.at(criterion2));
      h_global2D_TruthMatched_[criterion1][criterion2]->StatOverflows(kTRUE);
      h_global2D_TruthMatched_[criterion1][criterion2]->GetXaxis()->SetTitle(criterion1.c_str());
      h_global2D_TruthMatched_[criterion1][criterion2]->GetYaxis()->SetTitle(criterion2.c_str());
    }
  }

  for (unsigned int sequenceIndex = 0; sequenceIndex < stepByStepSequences_.size(); ++sequenceIndex) {
    const std::vector<std::string>& sequence = stepByStepSequences_.at(sequenceIndex);
    if (!(std::is_permutation(sequence.begin(), sequence.end(), photonIDCriteria_.begin()))) {
      std::cout << "ERROR: sequence is not a permutation of photonIDCriteria." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::string stepByStepPrefix = "sequence_";
    for (unsigned int stepIndex = 0; stepIndex < sequence.size(); ++stepIndex) {
      stepByStepPrefix += (sequence.at(stepIndex) + "_");
    }
    for (unsigned int stepIndex = 0; stepIndex < sequence.size(); ++stepIndex) {
      unsigned int stepNumber = stepIndex + 1;
      const std::string& criterion = sequence.at(stepIndex);
      h_stepByStep_TruthMatched_[sequenceIndex][stepIndex] = new TH1F((stepByStepPrefix + "step" + std::to_string(stepNumber)).c_str(), (stepByStepPrefix + "step" + std::to_string(stepNumber)).c_str(), nHistBins_.at(criterion), lowerHistLimits_.at(criterion), upperHistLimits_.at(criterion));
      h_stepByStep_TruthMatched_[sequenceIndex][stepIndex]->StatOverflows(kTRUE);
      h_stepByStep_TruthMatched_[sequenceIndex][stepIndex]->GetXaxis()->SetTitle(criterion.c_str());
    }
  }

  h_nMediumPhotons_ = new TH1F("nMediumPhotons", "nMediumPhotons;nMediumPhotons", 10, -0.5, 9.5);
  h_nMediumPhotons_->StatOverflows(kTRUE);
  h_nLoosePhotons_ = new TH1F("nLoosePhotons", "nLoosePhotons;nLoosePhotons", 10, -0.5, 9.5);
  h_nLoosePhotons_->StatOverflows(kTRUE);
  h_nMediumPhotons_TruthMatched_ = new TH1F("nMediumPhotons_TruthMatched", "nMediumPhotons_TruthMatched;nMediumPhotons_TruthMatched", 10, -0.5, 9.5);
  h_nMediumPhotons_TruthMatched_->StatOverflows(kTRUE);
  h_nLoosePhotons_TruthMatched_ = new TH1F("nLoosePhotons_TruthMatched", "nLoosePhotons_TruthMatched;nLoosePhotons_TruthMatched", 10, -0.5, 9.5);
  h_nLoosePhotons_TruthMatched_->StatOverflows(kTRUE);
  h_selectionRegion_ = new TH1F("selectionRegion", "selectionRegion;region", 4, -0.5, 3.5);
  h_selectionRegion_->StatOverflows(kTRUE);
  h_photonType_ = new TH1F("photonType", "photonType", 5, 0.5, 5.5);
  h_photonType_->StatOverflows(kTRUE);
  h_photonType_->GetXaxis()->SetBinLabel(h_photonType_->GetXaxis()->FindFixBin(1.0), "medium");
  h_photonType_->GetXaxis()->SetBinLabel(h_photonType_->GetXaxis()->FindFixBin(2.0), "medium, truth-matched");
  h_photonType_->GetXaxis()->SetBinLabel(h_photonType_->GetXaxis()->FindFixBin(3.0), "loose");
  h_photonType_->GetXaxis()->SetBinLabel(h_photonType_->GetXaxis()->FindFixBin(4.0), "loose, truth-matched");
  h_photonType_->GetXaxis()->SetBinLabel(h_photonType_->GetXaxis()->FindFixBin(5.0), "all, truth-matched");
  h_nPU_ = new TH1F("nPU", "nPU;nPU", 200, -0.5, 199.5);
  h_nPU_->StatOverflows(kTRUE);
  h_rho_ = new TH1F("rho", "rho;rho", 500, 0., 100.);
  h_rho_->StatOverflows(kTRUE);
  h_chIso_raw_TruthMatched_ = new TH1F("chIso_raw_TruthMatched", "chIso_raw_TruthMatched;chIso_raw", nHistBins_.at("chIso"), lowerHistLimits_.at("chIso"), upperHistLimits_.at("chIso"));
  h_neutIso_raw_TruthMatched_ = new TH1F("neutIso_raw_TruthMatched", "neutIso_raw_TruthMatched;neutIso_raw", nHistBins_.at("neutIso"), lowerHistLimits_.at("neutIso"), upperHistLimits_.at("neutIso"));
  h_phoIso_raw_TruthMatched_ = new TH1F("phoIso_raw_TruthMatched", "phoIso_raw_TruthMatched;phoIso_raw", nHistBins_.at("phoIso"), lowerHistLimits_.at("phoIso"), upperHistLimits_.at("phoIso"));
  h_phoET_TruthMatched_ = new TH1F("phoET_TruthMatched", "phoET_TruthMatched;photon ET;nEvents/bin", (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]));
  h_phoET_passingID_TruthMatched_ = new TH1F("phoET_passingID_TruthMatched", "phoET_passingID_TruthMatched;photon ET;nEvents/bin", (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]));
  h_overallETEfficiency_TruthMatched_ = new TEfficiency("overallETEfficiency_TruthMatched", "overallETEfficiency_TruthMatched;photon ET;#epsilon", (efficiencyBinEdges_.size()-1), &(efficiencyBinEdges_[0]));
  h_mediumFakeCriteria_ = new TH2F("mediumFakeCriteria", "ID criteria: (N-2) plot;sigmaIEtaIEta;chIso", nHistBins_.at("sigmaIEtaIEta"), lowerHistLimits_.at("sigmaIEtaIEta"), upperHistLimits_.at("sigmaIEtaIEta"), nHistBins_.at("chIso"), lowerHistLimits_.at("chIso"), upperHistLimits_.at("chIso"));
  h_mediumFakeCriteria_->StatOverflows(kTRUE);
  h_mediumFakeCriteria_TruthMatched_ = new TH2F("mediumFakeCriteria_TruthMatched", "ID criteria(truth-matched): (N-2) plot;sigmaIEtaIEta;chIso", nHistBins_.at("sigmaIEtaIEta"), lowerHistLimits_.at("sigmaIEtaIEta"), upperHistLimits_.at("sigmaIEtaIEta"), nHistBins_.at("chIso"), lowerHistLimits_.at("chIso"), upperHistLimits_.at("chIso"));
  h_mediumFakeCriteria_TruthMatched_->StatOverflows(kTRUE);
  outputPath_ = iConfig.getUntrackedParameter<std::string>("outputPath");
  outputFile_ = new TFile(outputPath_.c_str(), "RECREATE");
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity");
  filterType_ = iConfig.getUntrackedParameter<std::string>("filterType");
  if (filterType_ == "stealth") desired_parent_pid_ = 1000022;
  else if (filterType_ == "hgg") desired_parent_pid_ = 25;
  rhoLabel_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"));
  photonCollection_ = consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonSrc"));
  genParticlesCollection_ = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"));
  pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
}


MinimalMiniAODAnalyzer::~MinimalMiniAODAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  for (const auto& criterion: photonIDCriteria_) {
    outputFile_->WriteTObject(h_NMinus1_[criterion]);
    outputFile_->WriteTObject(h_NMinus1_TruthMatched_[criterion]);
    outputFile_->WriteTObject(h_NMinus1_ETEfficiencies_TruthMatched_[criterion]);
    outputFile_->WriteTObject(h_global1D_TruthMatched_[criterion]);
    outputFile_->WriteTObject(h_global1D_ETEfficiencies_TruthMatched_[criterion]);
  }
  for (const auto& criterion: diphotonHLTCriteria_) {
    outputFile_->WriteTObject(h_NMinus1_HLT_TruthMatched_[criterion]);
    outputFile_->WriteTObject(h_global1D_HLT_TruthMatched_[criterion]);
  }
  outputFile_->WriteTObject(h_HLTEfficiency_TruthMatched_);
  for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria_.size()); ++criterion1Index) {
    std::string& criterion1 = photonIDCriteria_.at(criterion1Index);
    for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria_.size(); ++criterion2Index) {
      std::string& criterion2 = photonIDCriteria_.at(criterion2Index);
      outputFile_->WriteTObject(h_NMinus2_TruthMatched_[criterion1][criterion2]);
      outputFile_->WriteTObject(h_global2D_TruthMatched_[criterion1][criterion2]);
    }
  }

  for (unsigned int sequenceIndex = 0; sequenceIndex < stepByStepSequences_.size(); ++sequenceIndex) {
    for (unsigned int stepIndex = 0; stepIndex < photonIDCriteria_.size(); ++stepIndex) {
      outputFile_->WriteTObject(h_stepByStep_TruthMatched_[sequenceIndex][stepIndex]);
    }
  }

  outputFile_->WriteTObject(h_nMediumPhotons_);
  outputFile_->WriteTObject(h_nLoosePhotons_);
  outputFile_->WriteTObject(h_nMediumPhotons_TruthMatched_);
  outputFile_->WriteTObject(h_nLoosePhotons_TruthMatched_);
  outputFile_->WriteTObject(h_selectionRegion_);
  outputFile_->WriteTObject(h_photonType_);
  outputFile_->WriteTObject(h_nPU_);
  outputFile_->WriteTObject(h_rho_);
  outputFile_->WriteTObject(h_chIso_raw_TruthMatched_);
  outputFile_->WriteTObject(h_neutIso_raw_TruthMatched_);
  outputFile_->WriteTObject(h_phoIso_raw_TruthMatched_);
  outputFile_->WriteTObject(h_phoET_TruthMatched_);
  outputFile_->WriteTObject(h_phoET_passingID_TruthMatched_);
  outputFile_->WriteTObject(h_overallETEfficiency_TruthMatched_);
  outputFile_->WriteTObject(h_mediumFakeCriteria_);
  outputFile_->WriteTObject(h_mediumFakeCriteria_TruthMatched_);
  outputFile_->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MinimalMiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<PileupSummaryInfo> > pileupInfo;
  iEvent.getByToken(pileupSummaryToken_, pileupInfo);
  float nPU = 1.0*(pileupInfo->begin()->getTrueNumInteractions());
  h_nPU_->Fill(nPU);

  std::vector<std::pair<float, float> > truthPhotonsEtaPhi;
  bool passesMC = passesTruthBasedSelection(iEvent, truthPhotonsEtaPhi);
  if (!(passesMC)) {
    return;
  }
  assert(static_cast<int>(truthPhotonsEtaPhi.size()) == 2);
  if (verbosity_ >= 4) {
    for (auto& etaPhiElement: truthPhotonsEtaPhi) {
      std::cout << "Found true photon with eta = " << etaPhiElement.first << ", phi = " << etaPhiElement.second << std::endl;
    }
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rho = *(rhoHandle.product());
  h_rho_->Fill(rho);

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonCollection_, photonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("MinimalMiniAODAnalyzer") << "no pat::Photons in event";
    return;
  }

  int nMediumPhotons = 0;
  int nLoosePhotons = 0;
  int nMediumPhotons_TruthMatched = 0;
  int nLoosePhotons_TruthMatched = 0;
  std::vector<std::pair<float, float> > fillQueue_mediumFakeCriteria;
  std::vector<std::pair<float, float> > fillQueue_mediumFakeCriteria_TruthMatched;
  std::map<std::string, bool> eventHLTBits;
  std::map<std::string, float> eventHLTProperties;
  float leadingET = -1.;
  float subLeadingET = -1.;
  int nTruthMatchedPhotonsInKinematicRegion = 0;
  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {
    float ET = iPho->et();
    float absEta = std::fabs(iPho->eta());

    if ((ET >= 25.) &&
        (absEta < 1.442)) {
      float min_truth_deltaR = getMinDeltaR(iPho->eta(), iPho->phi(), truthPhotonsEtaPhi);
      assert(min_truth_deltaR > 0.);
      bool isTruthMatched = (min_truth_deltaR < 0.05);
      if (isTruthMatched) ++nTruthMatchedPhotonsInKinematicRegion;
      if (isTruthMatched && verbosity_ >= 4) std::cout << "Found truth-matched photon at eta = " << iPho->eta() << ", phi = " << iPho->phi() << std::endl;
      float hOverE = iPho->hadTowOverEm();
      float R9 = iPho->full5x5_r9();
      bool passes_hOverE = (hOverE < 0.02197);
      float sigmaIEtaIEta = iPho->full5x5_sigmaIetaIeta();
      bool passes_sigmaIEtaIEta = (sigmaIEtaIEta < 0.01015);
      bool passes_sigmaIEtaIEta_loose = (sigmaIEtaIEta < 0.02);
      float chargedHadronIsolation = iPho->userFloat("phoChargedIsolation");
      if (isTruthMatched) h_chIso_raw_TruthMatched_->Fill(chargedHadronIsolation);
      float rhoCorrectedChargedHadronIsolation = getRhoCorrectedIsolation(absEta, PFTypeForEA::chargedHadron, chargedHadronIsolation, rho);
      bool passes_chIso = (rhoCorrectedChargedHadronIsolation < 1.141);
      bool passes_chIso_loose = (rhoCorrectedChargedHadronIsolation < 6.0);
      float neutralHadronIsolation = iPho->userFloat("phoNeutralHadronIsolation");
      if (isTruthMatched) h_neutIso_raw_TruthMatched_->Fill(neutralHadronIsolation/(1.189 + 0.01512*ET + 0.00002259*ET*ET));
      float rhoCorrectedNeutralHadronIsolation = getRhoCorrectedIsolation(absEta, PFTypeForEA::neutralHadron, neutralHadronIsolation, rho);
      float rhoCorrectedNeutralHadronIsolation_scaled = rhoCorrectedNeutralHadronIsolation/(1.189 + 0.01512*ET + 0.00002259*ET*ET);
      bool passes_neutIso = (rhoCorrectedNeutralHadronIsolation_scaled < 1.0);
      float photonIsolation = iPho->userFloat("phoPhotonIsolation");
      if (isTruthMatched) h_phoIso_raw_TruthMatched_->Fill(photonIsolation/(2.08 + 0.004017*ET));
      float rhoCorrectedPhotonIsolation = getRhoCorrectedIsolation(absEta, PFTypeForEA::photon, photonIsolation, rho);
      float rhoCorrectedPhotonIsolation_scaled = rhoCorrectedPhotonIsolation/(2.08 + 0.004017*ET);
      bool passes_phoIso = (rhoCorrectedPhotonIsolation_scaled < 1.0);
      float clusIso_scaled_for_leading = (iPho->ecalPFClusterIso())/(6.0 + 0.012*ET + rho*0.29);
      float clusIso_scaled_for_subLeading = (iPho->ecalPFClusterIso())/(6.0 + 0.012*ET + rho*0.16544);
      float trkIso_scaled = (iPho->trkSumPtHollowConeDR03())/(6.0 + 0.002*ET);
      
      std::map<std::string, bool> photonIDBits;
      std::map<std::string, float> photonProperties;
      photonIDBits["hOverE"] = passes_hOverE;
      photonIDBits["sigmaIEtaIEta"] = passes_sigmaIEtaIEta;
      photonIDBits["chIso"] = passes_chIso;
      photonIDBits["neutIso"] = passes_neutIso;
      photonIDBits["phoIso"] = passes_phoIso;
      photonProperties["hOverE"] = hOverE;
      photonProperties["sigmaIEtaIEta"] = sigmaIEtaIEta;
      photonProperties["chIso"] = rhoCorrectedChargedHadronIsolation;
      photonProperties["neutIso"] = rhoCorrectedNeutralHadronIsolation_scaled;
      photonProperties["phoIso"] = rhoCorrectedPhotonIsolation_scaled;
      checkMapKeysAgainstVector(photonIDBits, photonIDCriteria_);
      checkMapKeysAgainstVector(photonProperties, photonIDCriteria_);
      fillGlobal1DAndNMinus1Histograms(photonIDBits, photonProperties, isTruthMatched, ET);
      if (isTruthMatched) {
        h_photonType_->Fill(5.0);
        h_phoET_TruthMatched_->Fill(ET);
        fillGlobal2DAndNMinus2Histograms(photonIDBits, photonProperties);
        fillStepByStepHistograms(photonIDBits, photonProperties);
        h_overallETEfficiency_TruthMatched_->Fill(passes_hOverE && passes_neutIso && passes_phoIso && passes_chIso && passes_sigmaIEtaIEta, ET);
      }
      if (passes_hOverE && passes_neutIso && passes_phoIso) {
        fillQueue_mediumFakeCriteria.push_back(std::make_pair(sigmaIEtaIEta, rhoCorrectedChargedHadronIsolation));
        if (isTruthMatched) fillQueue_mediumFakeCriteria_TruthMatched.push_back(std::make_pair(sigmaIEtaIEta, rhoCorrectedChargedHadronIsolation));

        if (passes_chIso && passes_sigmaIEtaIEta) {
          ++nMediumPhotons;
          h_photonType_->Fill(1.0);
          if (isTruthMatched) {
            ++nMediumPhotons_TruthMatched;
            h_photonType_->Fill(2.0);
            h_phoET_passingID_TruthMatched_->Fill(ET);
          }
        }
        else if (passes_chIso_loose && passes_sigmaIEtaIEta_loose) {
          ++nLoosePhotons;
          h_photonType_->Fill(3.0);
          if (isTruthMatched) {
            ++nLoosePhotons_TruthMatched;
            h_photonType_->Fill(4.0);
          }
        }
      } // ends special 2D plots
      if (isTruthMatched) {
        if ((leadingET < 0.) || (ET > leadingET)) { // leading photon is not set; assume this one is the leading photon, OR leading photon is set but ET of this photon is greater than ET of currently leading photon
          if ((leadingET > 0.) && (ET > leadingET)) { // in the latter case, first copy currently leading photon info to subleading photon
            eventHLTProperties["hOverE_subLeading"] = eventHLTProperties["hOverE_leading"];
            eventHLTProperties["R9_subLeading"] = eventHLTProperties["R9_leading"];
            eventHLTProperties["sigmaIEtaIEta_subLeading"] = eventHLTProperties["sigmaIEtaIEta_leading"];
            eventHLTProperties["clusIso_subLeading"] = eventHLTProperties["clusIso_leading"]*((6.0 + 0.012*leadingET + 0.29*rho)/(6.0 + 0.012*ET + 0.16544*rho)); // needs to be recalculated because criteria are different between leading and subleading legs
            eventHLTProperties["trkIso_subLeading"] = eventHLTProperties["trkIso_leading"];
            eventHLTBits["hOverE_subLeading"] = eventHLTBits["hOverE_leading"];
            eventHLTBits["R9_subLeading"] = eventHLTBits["R9_leading"];
            eventHLTBits["sigmaIEtaIEta_subLeading"] = eventHLTBits["sigmaIEtaIEta_leading"];
            eventHLTBits["clusIso_subLeading"] = (eventHLTProperties["clusIso_subLeading"] < 1.0); // needs to be recalculated because clus iso is redefined
            eventHLTBits["trkIso_subLeading"] = (eventHLTProperties["trkIso_subLeading"] < 1.0); // the trk iso is only applicable for the subleading photon, this bit is always set to true for the leading photon so it has to be recalculated
            subLeadingET = leadingET;
          }
          // set leading
          eventHLTProperties["hOverE_leading"] = hOverE;
          eventHLTProperties["R9_leading"] = R9;
          eventHLTProperties["sigmaIEtaIEta_leading"] = sigmaIEtaIEta;
          eventHLTProperties["clusIso_leading"] = clusIso_scaled_for_leading;
          eventHLTProperties["trkIso_leading"] = trkIso_scaled;
          eventHLTBits["hOverE_leading"] = (hOverE <= 0.12);
          eventHLTBits["trkIso_leading"] = true;
          eventHLTBits["R9_leading"] = true;
          if ((R9 >= 0.5) && (R9 < 0.85)) {
            eventHLTBits["sigmaIEtaIEta_leading"] = true;
            eventHLTBits["clusIso_leading"] = true;
          }
          else if (R9 >= 0.85) {
            eventHLTBits["sigmaIEtaIEta_leading"] = (sigmaIEtaIEta < 0.015);
            eventHLTBits["clusIso_leading"] = (clusIso_scaled_for_leading <= 1.0);
          }
          else {
            eventHLTBits["R9_leading"] = false;
            eventHLTBits["sigmaIEtaIEta_leading"] = true;
            eventHLTBits["clusIso_leading"] = true;
          }
          leadingET = ET;
        }
        else { // leadingET > 0 and ET < leadingET; leading photon is set already
          if ((subLeadingET < 0.) || (ET > subLeadingET)) { // this photon is the subleading photon
            // set subleading
            eventHLTProperties["hOverE_subLeading"] = hOverE;
            eventHLTProperties["R9_subLeading"] = R9;
            eventHLTProperties["sigmaIEtaIEta_subLeading"] = sigmaIEtaIEta;
            eventHLTProperties["clusIso_subLeading"] = clusIso_scaled_for_subLeading;
            eventHLTProperties["trkIso_subLeading"] = trkIso_scaled;
            eventHLTBits["hOverE_subLeading"] = (hOverE <= 0.12);
            eventHLTBits["R9_subLeading"] = true;
            if ((R9 >= 0.5) && (R9 < 0.85)) {
              eventHLTBits["sigmaIEtaIEta_subLeading"] = true;
              eventHLTBits["clusIso_subLeading"] = true;
              eventHLTBits["trkIso_subLeading"] = true;
            }
            else if (R9 >= 0.85) {
              eventHLTBits["sigmaIEtaIEta_subLeading"] = (sigmaIEtaIEta < 0.015);
              eventHLTBits["clusIso_subLeading"] = (clusIso_scaled_for_subLeading <= 1.0);
              eventHLTBits["trkIso_subLeading"] = (trkIso_scaled <= 1.0);
            }
            else {
              eventHLTBits["R9_subLeading"] = false;
              eventHLTBits["sigmaIEtaIEta_subLeading"] = true;
              eventHLTBits["clusIso_subLeading"] = true;
              eventHLTBits["trkIso_subLeading"] = true;
            }
            subLeadingET = ET;
          }
        }
      }
    } // ends kinematic check
  } // ends loop over photon objects
  
  if (nTruthMatchedPhotonsInKinematicRegion == 2) {
    checkMapKeysAgainstVector(eventHLTBits, diphotonHLTCriteria_);
    checkMapKeysAgainstVector(eventHLTProperties, diphotonHLTCriteria_);
    fillGlobal1DAndNMinus1HLTHistograms(eventHLTBits, eventHLTProperties, leadingET, subLeadingET);
  }

  h_nMediumPhotons_->Fill(1.0*nMediumPhotons);
  h_nLoosePhotons_->Fill(1.0*nLoosePhotons);
  h_nMediumPhotons_TruthMatched_->Fill(1.0*nMediumPhotons_TruthMatched);
  h_nLoosePhotons_TruthMatched_->Fill(1.0*nLoosePhotons_TruthMatched);
  float selectionRegion = 0.;
  if (nMediumPhotons == 2) {
    selectionRegion = 1.;
  }
  else if ((nMediumPhotons == 1) && (nLoosePhotons >= 1)) {
    selectionRegion = 2.;
  }
  else if ((nMediumPhotons == 0) && (nLoosePhotons >= 2)) {
    selectionRegion = 3.;
  }
  h_selectionRegion_->Fill(selectionRegion);
  
  if ((nMediumPhotons == 2) ||
      ((nMediumPhotons == 1) && (nLoosePhotons >= 1)) ||
      ((nMediumPhotons == 0) && (nLoosePhotons >= 2))) {
    for (const auto& valuePair: fillQueue_mediumFakeCriteria) {
      h_mediumFakeCriteria_->Fill(valuePair.first, valuePair.second);
    }
    for (const auto& valuePair: fillQueue_mediumFakeCriteria_TruthMatched) {
      h_mediumFakeCriteria_TruthMatched_->Fill(valuePair.first, valuePair.second);
    }
  }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

float
MinimalMiniAODAnalyzer::getRhoCorrectedIsolation(const float& absEta, const PFTypeForEA& type, const float& isolation, const double& rho) {
  float EA = 0.;
  if (absEta < 1.) EA = EAValues_barrelCenter_.at(type);
  else if (absEta < 1.479) EA = EAValues_barrelEdge_.at(type);
  else {
    edm::LogError("MinimalMiniAODAnalyzer") << "Called for unsupported absEta = " << absEta;
    std::exit(EXIT_FAILURE);
  }
  float corrected = std::max(isolation - EA*(static_cast<float>(rho)), 0.0f);
  if (verbosity_ >= 5) std::cout << "Called for absEta = " << absEta << ", type = " << PFTypeNames_.at(type) << ", rho = " << rho << ". Isolation values: raw = " << isolation << ", corrected: " << corrected << std::endl;
  return corrected;
}

bool
MinimalMiniAODAnalyzer::passesTruthBasedSelection(const edm::Event& iEvent, std::vector<std::pair<float, float> >& truthPhotonsEtaPhi) {
  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  if (!genParticlesHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no reco::GenParticles in event";
    return false;
  }

  int nPhotonsWithDesiredMom = 0;
  for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    if ((ip->fromHardProcessFinalState()) &&
        (ip->isPromptFinalState()) &&
        (ip->isHardProcess())) {
      int pid = ip->pdgId();
      if (pid == 22) {
        const reco::Candidate * mom = ip->mother();
        int mom_pid = mom->pdgId();
        if (verbosity_ >= 4) std::cout << "Found final state photon with mom ID = " << mom_pid << std::endl;
        if (mom_pid == desired_parent_pid_) {
          truthPhotonsEtaPhi.push_back(std::make_pair(static_cast<float>(ip->eta()), static_cast<float>(ip->phi())));
          ++nPhotonsWithDesiredMom;
        }
      }
    }
  }
  return (nPhotonsWithDesiredMom == 2);
}

float
MinimalMiniAODAnalyzer::get_deltaR(const float& eta, const float& phi, const float& target_eta, const float& target_phi) {
  float deltaEta = target_eta - eta;
  float phi1 = target_phi;
  float phi2 = phi;
  if (phi2 > phi1) { // make sure phi1 > phi2
    phi1 = phi;
    phi2 = target_phi;
  }
  float deltaPhi = std::min((phi1 - phi2), (static_cast<float>(2.0*TMath::Pi()) - (phi1 - phi2)));
  return std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

float
MinimalMiniAODAnalyzer::getMinDeltaR(const float& eta, const float& phi, std::vector<std::pair<float, float> >& targetEtaPhiList) {
  float min_dR = -0.005;
  for (auto& angularVariablesToCompare: targetEtaPhiList) {
    float dR = get_deltaR(eta, phi, angularVariablesToCompare.first, angularVariablesToCompare.second);
    if ((min_dR < 0) || (dR < min_dR)) min_dR = dR;
  }
  return min_dR;
}

void
MinimalMiniAODAnalyzer::fillGlobal1DAndNMinus1Histograms(const std::map<std::string, bool>& photonIDBits, const std::map<std::string, float>& photonProperties, const bool& isTruthMatched, const float& ET) {
  std::map<std::string, bool> otherCriteriaAreMet;
  for (unsigned int criterionIndex = 0; criterionIndex < photonIDCriteria_.size(); ++criterionIndex) {
    std::string& criterion = photonIDCriteria_.at(criterionIndex);
    bool passesOtherRequirements = true;
    for (unsigned int otherCriterionIndex = 0; otherCriterionIndex < photonIDCriteria_.size(); ++otherCriterionIndex) {
      if (criterionIndex == otherCriterionIndex) continue;
      std::string& otherCriterion = photonIDCriteria_.at(otherCriterionIndex);
      if (!(photonIDBits.at(otherCriterion))) {
        passesOtherRequirements = false;
        break;
      }
    }
    otherCriteriaAreMet[criterion] = passesOtherRequirements;
  }
  for (const std::string& criterion: photonIDCriteria_) {
    if (otherCriteriaAreMet.at(criterion)) (h_NMinus1_.at(criterion))->Fill(photonProperties.at(criterion));
    if (isTruthMatched) {
      (h_global1D_TruthMatched_.at(criterion))->Fill(photonProperties.at(criterion));
      (h_global1D_ETEfficiencies_TruthMatched_.at(criterion))->Fill(photonIDBits.at(criterion), ET);
      if (otherCriteriaAreMet.at(criterion)) {
        (h_NMinus1_TruthMatched_.at(criterion))->Fill(photonProperties.at(criterion));
        (h_NMinus1_ETEfficiencies_TruthMatched_.at(criterion))->Fill(photonIDBits.at(criterion), ET);
      }
    }
  }
}

void
MinimalMiniAODAnalyzer::fillGlobal1DAndNMinus1HLTHistograms(const std::map<std::string, bool>& eventHLTBits, const std::map<std::string, float>& eventHLTProperties, const float& leadingET, const float& subLeadingET) {
  std::map<std::string, bool> otherCriteriaAreMet;
  for (unsigned int criterionIndex = 0; criterionIndex < diphotonHLTCriteria_.size(); ++criterionIndex) {
    std::string& criterion = diphotonHLTCriteria_.at(criterionIndex);
    bool passesOtherRequirements = true;
    for (unsigned int otherCriterionIndex = 0; otherCriterionIndex < diphotonHLTCriteria_.size(); ++otherCriterionIndex) {
      if (criterionIndex == otherCriterionIndex) continue;
      std::string& otherCriterion = diphotonHLTCriteria_.at(otherCriterionIndex);
      if (!(eventHLTBits.at(otherCriterion))) {
        passesOtherRequirements = false;
        break;
      }
    }
    otherCriteriaAreMet[criterion] = passesOtherRequirements;
  }

  for (const std::string& criterion: diphotonHLTCriteria_) {
    (h_global1D_HLT_TruthMatched_[criterion])->Fill(eventHLTProperties.at(criterion));
    if (otherCriteriaAreMet.at(criterion)) (h_NMinus1_HLT_TruthMatched_[criterion])->Fill(eventHLTProperties.at(criterion));
  }

  h_HLTEfficiency_TruthMatched_->Fill((otherCriteriaAreMet.at(diphotonHLTCriteria_.at(0)) && eventHLTBits.at(diphotonHLTCriteria_.at(0))), leadingET, subLeadingET);
}

void
MinimalMiniAODAnalyzer::fillGlobal2DAndNMinus2Histograms(const std::map<std::string, bool>& photonIDBits, const std::map<std::string, float>& photonProperties) {
  for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria_.size()); ++criterion1Index) {
    std::string& criterion1 = photonIDCriteria_.at(criterion1Index);
    for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria_.size(); ++criterion2Index) {
      std::string& criterion2 = photonIDCriteria_.at(criterion2Index);
      h_global2D_TruthMatched_[criterion1][criterion2]->Fill(photonProperties.at(criterion1), photonProperties.at(criterion2));
      bool passesOtherRequirements = true;
      for (unsigned int otherCriterionIndex = 0; otherCriterionIndex < photonIDCriteria_.size(); ++otherCriterionIndex) {
        if ((otherCriterionIndex == criterion1Index) || (otherCriterionIndex == criterion2Index)) continue;
        std::string& otherCriterion = photonIDCriteria_.at(otherCriterionIndex);
        if (verbosity_ >= 5) std::cout << "For criterion1 = " << criterion1 << ", criterion2 = " << criterion2 << ", checking bit value at: " << otherCriterion << std::endl;
        if (!(photonIDBits.at(otherCriterion))) {
          passesOtherRequirements = false;
          break;
        }
      }
      if (passesOtherRequirements) {
        h_NMinus2_TruthMatched_[criterion1][criterion2]->Fill(photonProperties.at(criterion1), photonProperties.at(criterion2));
      }
    }
  }
}

void
MinimalMiniAODAnalyzer::fillStepByStepHistograms(const std::map<std::string, bool>& photonIDBits, const std::map<std::string, float>& photonProperties) {
  for (unsigned int sequenceIndex = 0; sequenceIndex < stepByStepSequences_.size(); ++sequenceIndex) {
    const std::vector<std::string>& sequence = stepByStepSequences_.at(sequenceIndex);
    for (unsigned int stepIndex = 0; stepIndex < photonIDCriteria_.size(); ++stepIndex) {
      const std::string& criterion = sequence.at(stepIndex);
      h_stepByStep_TruthMatched_[sequenceIndex][stepIndex]->Fill(photonProperties.at(criterion));
      if (!(photonIDBits.at(criterion))) break;
    }
  }
}

template<typename valueType>
void MinimalMiniAODAnalyzer::checkMapKeysAgainstVector(const std::map<std::string, valueType>& mapToCheck, const std::vector<std::string>& allowedValues) {
  assert(mapToCheck.size() == allowedValues.size());
  for (auto const& mapElement : mapToCheck) {
    assert(std::find(allowedValues.begin(), allowedValues.end(), mapElement.first) != allowedValues.end());
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
MinimalMiniAODAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MinimalMiniAODAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MinimalMiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MinimalMiniAODAnalyzer);
