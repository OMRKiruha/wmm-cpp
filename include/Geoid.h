//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

#include <array>
#include <memory>

namespace wmm {

    static constexpr bool useGeoid    = true;  // 1 Geoid - Ellipsoid difference should be corrected, 0 otherwise
    static constexpr int geoidCols    = 1441;  // 360 degrees of longitude at 15 minute spacing
    static constexpr int geoidRows    = 721;   // 180 degrees of latitude  at 15 minute spacing
    static constexpr int geoidHeights = geoidCols * geoidRows;

    struct CoordGeodetic;

    struct Geoid {
        Geoid();

        bool getGeoidHeight(double latitude, double longitude, double &deltaHeight) const;

        using GeoidHeightArray_t = std::array<float, geoidHeights>;
        // Defaults is EGM-96 model file parameters
        int numbGeoidCols{geoidCols};    // 360 degrees of longitude at 15 minute spacing
        int numbGeoidRows{geoidRows};    // 180 degrees of latitude  at 15 minute spacing
        int numbHeaderItems{6};          // min, max lat, min, max long, lat, long spacing
        int scaleFactor{4};              // 4 grid cells per degree at 15 minute spacing
        std::unique_ptr<GeoidHeightArray_t> geoidHeightBuffer{nullptr};
        int numbGeoidElevs{geoidHeights};
        bool isGeoidInitialized{false};  // indicates successful initialization
        bool isUseGeoid{useGeoid};       // Is the Geoid being used?
    };
}  // namespace wmm