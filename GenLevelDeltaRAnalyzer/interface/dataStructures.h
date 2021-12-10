#ifndef H_DATASTRUCTURES
#define H_DATASTRUCTURES

#include <cmath>

#include "constants.h"

struct indexAndDeltaRStruct{
  int index;
  float deltaR;

  indexAndDeltaRStruct () : index(-1), deltaR(-0.1) {}
  indexAndDeltaRStruct (int index_, float deltaR_) : index(index_), deltaR(deltaR_) {}
};

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

  std::pair<indexAndDeltaRStruct, indexAndDeltaRStruct> get_indices_and_deltaRs_closest_and_second_closest_from_list_of_angles(const std::vector<angularVariablesStruct>& list_of_angles) {
    int index_closest = -1;
    float deltaR_closest = -0.1;
    int index_second_closest = -1;
    float deltaR_second_closest = -0.1;
    bool closest_is_set = false;
    bool second_closest_is_set = false;

    if (list_of_angles.size() > 0) {
      for (size_t angle_index = 0; angle_index < list_of_angles.size(); ++angle_index) {
	float deltaR = this->get_deltaR(list_of_angles.at(angle_index));
	if (!(closest_is_set)) { // closest is unset: set it to whatever values we get from the first entry
	  index_closest = angle_index;
	  deltaR_closest = deltaR;
	  closest_is_set = true;
	  assert(!(second_closest_is_set));
	}
	else if (!(second_closest_is_set)) { // closest is set, but second closest is unset
	  if (deltaR < deltaR_closest) { // first check with deltaR_closest and swap if needed
	    index_second_closest = index_closest;
	    deltaR_second_closest = deltaR_closest;
	    second_closest_is_set = true;
	    index_closest = angle_index;
	    deltaR_closest = deltaR;
	  }
	  else { // else set the second closest to current entry
	    index_second_closest = angle_index;
	    deltaR_second_closest = deltaR;
	    second_closest_is_set = true;
	  }
	}
	else { // both closest and second-closest are set
	  if (deltaR < deltaR_closest) { // first check with deltaR_closest and swap if needed
	    index_second_closest = index_closest;
	    deltaR_second_closest = deltaR_closest;
	    index_closest = angle_index;
	    deltaR_closest = deltaR;
	  }
	  else if (deltaR < deltaR_second_closest) { // deltaR is larger than closest but smaller than second closest
	    index_second_closest = angle_index;
	    deltaR_second_closest = deltaR;
	  }
	}
      }
    }

    indexAndDeltaRStruct closest = indexAndDeltaRStruct(index_closest, deltaR_closest);
    indexAndDeltaRStruct second_closest = indexAndDeltaRStruct(index_second_closest, deltaR_second_closest);
    return std::make_pair(closest, second_closest);
  }
};

struct PTEtaPhiStruct{
  float pt, eta, phi;
  PTEtaPhiStruct () : pt(0.), eta(0.), phi(0.) {}
  PTEtaPhiStruct (float pt_, float eta_, float phi_) : pt(pt_), eta(eta_), phi(phi_) {}

  friend std::ostream& operator<< (std::ostream& out, const PTEtaPhiStruct& pt_eta_phi_) {
    out << "(pt: " << pt_eta_phi_.pt << ", eta: " << pt_eta_phi_.eta << ", phi: " << pt_eta_phi_.phi << ")";
    return out;
  }
};

#endif
