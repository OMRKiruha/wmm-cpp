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

        double lambdag{0.0};  // geocentric longitude
        double phig{0.0};     // geocentric latitude
        double r{0.0};        // distance from the center of the ellipsoid
    };
}  // namespace wmm