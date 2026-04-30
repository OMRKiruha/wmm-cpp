//	WMM Subroutine library was tested in the following environments
// 1. Red Hat Linux with GCC Compiler
// 2. MS Windows XP with CodeGear C++ compiler
// 3. Sun Solaris with GCC Compiler
// Revision Number : $Revision : 1437 $
// Last changed by : $Author : Li-Yin Young $
// Last changed on : $Date : 2024 - 11 - 08 10 : 49 : 40 - 0700 $
// Edited by Kiryuhin Viacheslav on 04/21/26

#pragma once

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "CoordGeodetic.h"
#include "CoordSpherical.h"
#include "Date.h"
#include "Ellipsoid.h"
#include "GeoMagneticElements.h"
#include "Geoid.h"
#include "Gradient.h"
#include "LegendreFunction.h"
#include "MagneticConstants.h"
#include "MagneticModel.h"
#include "MagneticResults.h"
#include "SphericalHarmonicVariables.h"
#include "UTMParameters.h"

namespace wmm {

    struct CoordGeodeticStr {
        std::string Longitude;
        std::string Latitude;
    };

    enum PARAMS {
        SHDF,
        MODELNAME,
        PUBLISHER,
        RELEASEDATE,
        DATACUTOFF,
        MODELSTARTYEAR,
        MODELENDYEAR,
        EPOCH,
        INTSTATICDEG,
        INTSECVARDEG,
        EXTSTATICDEG,
        EXTSECVARDEG,
        GEOMAGREFRAD,
        NORMALIZATION,
        SPATBASFUNC
    };

    enum YYYYMMDD { YEAR, MONTH, DAY };

    // Conversions, Transformations, and other Calculations
    void BaseErrors(double declCoef, double declBaseline, double inclOffset, double fOffset, double multiplier, double H,
                    double *declErr, double *inclErr, double *fErr);

    void DegreeToDMSstring(double degreesOfArc, int unitDepth, std::string &out);

    void DMSstringToDegree(std::string_view DMSstring, double *degreesOfArc);

}  // namespace wmm