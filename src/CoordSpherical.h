//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

namespace wmm {
    struct Ellipsoid;
    struct CoordGeodetic;

    struct CoordSpherical {
        static CoordSpherical sphericalFromGeodetic(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic);
        void fromGeodetic(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic);
        void toCartesian(double &x, double &y, double &z) const;

        double lambdag{};  // geocentric longitude
        double phig{};     // geocentric latitude
        double r{};        // distance from the center of the ellipsoid
    };
}  // namespace wmm