#ifndef H_DATASTRUCTURES
#define H_DATASTRUCTURES

#include <cmath>

#include "constants.h"

struct angularVariablesStruct{
  float eta, phi;

  angularVariablesStruct () : eta(-999.), phi(-1.) {}

  angularVariablesStruct (float eta_, float phi_) : eta(eta_), phi(phi_) {}

  float get_deltaR(const angularVariablesStruct& angularVariables) {
    float deltaEta = angularVariables.eta - this->eta;
    float phi1 = angularVariables.phi;
    float phi2 = this->phi;
    if (phi2 > phi1) { // make sure phi1 > phi2
      phi1 = this->phi;
      phi2 = angularVariables.phi;
    }
    float deltaPhi_direction1 = phi1 - phi2;
    float deltaPhi_direction2 = constants::TWOPI - deltaPhi_direction1;
    float deltaPhi = (deltaPhi_direction1 < deltaPhi_direction2) ? deltaPhi_direction1 : deltaPhi_direction2;
    return std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
  }
};

#endif
