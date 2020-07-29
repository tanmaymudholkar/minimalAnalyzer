#ifndef H_CONSTANTS
#define H_CONSTANTS

#include "TMath.h"

namespace constants{ // for readability
  const float PI = static_cast<float>(TMath::Pi());
  const float TWOPI = static_cast<float>(2.0*TMath::Pi());
}

namespace photonCuts{
  const float pTLeading = 25.0f;
  const float pTSubleading = 35.0f;
}

#endif
