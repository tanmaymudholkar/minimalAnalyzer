#ifndef H_TRIGGERS
#define H_TRIGGERS

#include <iostream>
#include <string>
#include <array>

namespace triggerPatterns{
  const std::array<std::string, 6> patternsToSave { {
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v*",
        "HLT_Diphoton30PV_*_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v*",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_*_Mass55_v*",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v*",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
        "HLT_Diphoton30_*_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"
    } };
}

#endif
