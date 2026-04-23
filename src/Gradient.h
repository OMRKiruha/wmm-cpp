//
// Created by Professional on 22.04.2026.
//

#pragma once

#include "GeoMagneticElements.h"

namespace wmm{
    struct Ellipsoid;
    struct CoordGeodetic;
    struct MagneticModel;

    struct Gradient {
        void Calculate(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel);

        int UseGradient{};
        GeoMagneticElements GradPhi;     // phi
        GeoMagneticElements GradLambda;  // lambdag
        GeoMagneticElements GradZ;
    };
}