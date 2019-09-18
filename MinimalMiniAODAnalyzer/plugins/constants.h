#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

namespace constants{
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
  std::vector<std::string> diphotonHLTCriteria_ = {"R9_leading", "hOverE_leading", "sigmaIEtaIEta_leading", "clusIso_leading", "trkIso_leading", "R9_subLeading", "hOverE_subLeading", "sigmaIEtaIEta_subLeading", "clusIso_subLeading", "trkIso_subLeading"};

  float deltaR_photonTruthMatching = 0.05;
  float cut_hOverE = 0.02197;
  float cut_sigmaIEtaIEta = 0.01015;
  float cutLoose_sigmaIEtaIEta = 0.02;
  float cut_chargedIsolation = 1.141;
  float cutLoose_chargedIsolation = 6.0;
  float neutralHadronScale_constant = 1.189;
  float neutralHadronScale_linear = 0.01512;
  float neutralHadronScale_quadratic = 0.00002259;
  float photonIsolationScale_constant = 2.08;
  float photonIsolationScale_linear = 0.004017;

  float clusIsoScale_HLT_leading_constant = 6.0;
  float clusIsoScale_HLT_leading_linear = 0.012;
  float clusIsoScale_HLT_leading_rhoTerm = 0.29;
  float clusIsoScale_HLT_subLeading_constant = 6.0;
  float clusIsoScale_HLT_subLeading_linear = 0.012;
  float clusIsoScale_HLT_subLeading_rhoTerm = 0.16544;
  float trkIsoScale_HLT_constant = 6.0;
  float trkIsoScale_HLT_linear = 0.002;
  float cut_HLT_hOverE = 0.12;
  float cut_HLT_sigmaIEtaIEta = 0.015;
}
