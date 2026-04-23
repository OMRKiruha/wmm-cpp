//
// Created by Professional on 22.04.2026.
//

#pragma once

#include <array>
#include <memory>

namespace wmm {

    static constexpr int USE_GEOID    = 1;     // 1 Geoid - Ellipsoid difference should be corrected, 0 otherwise
    static constexpr int GeoidCols    = 1441;  // 360 degrees of longitude at 15 minute spacing
    static constexpr int GeoidRows    = 721;   // 180 degrees of latitude  at 15 minute spacing
    static constexpr int GeoidHeights = GeoidCols * GeoidRows;

    struct CoordGeodetic;

    struct Geoid {
        Geoid();
        ~Geoid() = default;

        int GetGeoidHeight(double Latitude, double Longitude, double *DeltaHeight) const;

        using GeoidHeightArray_t = std::array<float, GeoidHeights>;
        // Defaults is EGM-96 model file parameters
        int NumbGeoidCols{GeoidCols};  // 360 degrees of longitude at 15 minute spacing
        int NumbGeoidRows{GeoidRows};  // 180 degrees of latitude  at 15 minute spacing
        int NumbHeaderItems{6};        // min, max lat, min, max long, lat, long spacing
        int ScaleFactor{4};            // 4 grid cells per degree at 15 minute spacing
        std::unique_ptr<GeoidHeightArray_t> GeoidHeightBuffer{nullptr};
        int NumbGeoidElevs{GeoidHeights};
        int Geoid_Initialized{0};      // indicates successful initialization
        int UseGeoid{USE_GEOID};       // Is the Geoid being used?
    };
}  // namespace wmm