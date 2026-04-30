//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

#include "GeoMagneticElements.h"

namespace wmm {
    struct Ellipsoid;
    struct CoordGeodetic;
    struct MagneticModel;

    struct Gradient {
        void calculate(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel);

        int useGradient{};
        GeoMagneticElements gradPhi;     // phi
        GeoMagneticElements gradLambda;  // lambdag
        GeoMagneticElements gradZ;
    };
}  // namespace wmm