//
// Created by Professional on 22.04.2026.
//

#pragma once

#include "Geoid.h"

namespace wmm {
    struct Ellipsoid;
    struct CoordSpherical;

    struct CoordGeodetic {
        void CheckGeographicPole();
        void fromSpherical(const Ellipsoid &ellip, const CoordSpherical &coordSpherical);
        void fromCartesian(const Ellipsoid &ellip, double x, double y, double z);
        static void equivalentLatLon(double lat, double lon, double *repairedLat, double *repairedLon);
        int convertGeoidToEllipsoidHeight(const Geoid &geoid);

        double lambda{};                // geodetic longitude
        double phi{};                   // geodetic latitude
        double HeightAboveEllipsoid{};  // height above the ellipsoid (HaE)
        double HeightAboveGeoid{};      // (height above the EGM96 geoid model)
        int UseGeoid{1};
    };
}  // namespace wmm
