//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "SphericalHarmonicVariables.h"

#include "CoordSpherical.h"
#include "Ellipsoid.h"
#include "MagneticUtils.h"

namespace wmm {

    SphericalHarmonicVariables::SphericalHarmonicVariables(const Ellipsoid &ellip, const CoordSpherical &coordSpherical,
                                                           const int nMax) {
        relativeRadiusPower.resize(nMax + 1);
        cos_mlambda.resize(nMax + 1);
        sin_mlambda.resize(nMax + 1);

        // Computes Spherical variables. Variables computed are (a/r)^(n+2), cos_m(lambdag) and sin_m(lambdag) for spherical
        // harmonic summations. (Equations 10-12 in the WMM Technical Report)
        const double cos_lambda = cos(deg2Rad(coordSpherical.lambdag));
        const double sin_lambda = sin(deg2Rad(coordSpherical.lambdag));
        // for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2) for n  1..nMax-1 (this is much
        // faster than calling pow MAX_N+1 times).
        relativeRadiusPower.at(0) = (ellip.re / coordSpherical.r) * (ellip.re / coordSpherical.r);
        for(int n = 1; n <= nMax; n++) {
            relativeRadiusPower.at(n) = relativeRadiusPower.at(n - 1) * (ellip.re / coordSpherical.r);
        }

        // Compute cos(m*lambdag), sin(m*lambdag) for m = 0 ... nMax
        // cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
        // sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
        cos_mlambda.at(0) = 1.0;  // The size if cos_mlambda and sin_mlambda is nMax+1
        sin_mlambda.at(0) = 0.0;
        if(nMax + 1 >= 2) {
            cos_mlambda.at(1) = cos_lambda;
            sin_mlambda.at(1) = sin_lambda;
        }

        if(nMax + 1 >= 3) {
            for(int m = 2; m <= nMax; m++) {
                cos_mlambda.at(m) = (cos_mlambda.at(m - 1) * cos_lambda) - (sin_mlambda.at(m - 1) * sin_lambda);
                sin_mlambda.at(m) = (cos_mlambda.at(m - 1) * sin_lambda) + (sin_mlambda.at(m - 1) * cos_lambda);
            }
        }
    }
}  // namespace wmm