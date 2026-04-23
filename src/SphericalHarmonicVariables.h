//
// Created by Professional on 22.04.2026.
//

#pragma once

#include <vector>

namespace wmm {
    struct Ellipsoid;
    struct CoordSpherical;

    struct SphericalHarmonicVariables {
        SphericalHarmonicVariables() = delete;
        SphericalHarmonicVariables(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, int nMax);
        ~SphericalHarmonicVariables() = default;

        std::vector<double> RelativeRadiusPower;  // [earth_reference_radius_km / sph. radius ]^n
        std::vector<double> cos_mlambda;          // cp(m)  - cosine of (m*spherical coord. longitude)
        std::vector<double> sin_mlambda;          // sp(m)  - sine of (m*spherical coord. longitude)
    };

}  // namespace wmm