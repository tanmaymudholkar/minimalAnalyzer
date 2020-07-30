#ifndef H_CONSTANTS
#define H_CONSTANTS

#include "TMath.h"

namespace constants{ // for readability
  const float PI = static_cast<float>(TMath::Pi());
  const float TWOPI = static_cast<float>(2.0*TMath::Pi());
}

namespace EAValues{
  const float neutIso_central = 0.0668f;
  const float neutIso_peripheral = 0.1054f;
  const float phoIso_central = 0.1113f;
  const float phoIso_peripheral = 0.0953f;
  const float chIso_central = 0.0112f;
  const float chIso_peripheral = 0.0108f;
}

namespace photonCuts{
  const float pTLeading = 35.0f;
  const float pTSubleading = 25.0f;
  const float eta = 1.442f;
  const float sigma_ieta_ieta = 0.01015f;
  const float hOverE = 0.02197f;
  const float neutIso_const = 1.189f;
  const float neutIso_lin = 0.01512f;
  const float neutIso_quad = 0.00002259f;
  const float phoIso_const = 2.08f;
  const float phoIso_lin = 0.004017f;
  const float chIso = 1.141f;
}

#endif
