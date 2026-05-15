//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

namespace wmm {
    struct Geoid;
    struct Ellipsoid;
    struct CoordSpherical;

    struct CoordGeodetic {
        [[maybe_unused]] void fromSpherical(const Ellipsoid &ellip, const CoordSpherical &coordSpherical);
        void fromCartesian(const Ellipsoid &ellip, double x, double y, double z);
        [[maybe_unused]] void checkGeographicPole();
        static void equivalentLatLon(double lat, double lon, double &repairedLat, double &repairedLon);
        bool convertGeoidToEllipsoidHeight(const Geoid &geoid);

        double lambda{0.0};                // geodetic longitude
        double phi{0.0};                   // geodetic latitude
        double heightAboveEllipsoid{0.0};  // height above the ellipsoid (HaE)
        double heightAboveGeoid{0.0};      // height above the EGM96 geoid model
        bool isUseGeoid{true};             // by default use geoid in calculations
    };
}  // namespace wmm
