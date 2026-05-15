//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "CoordSpherical.h"

#include "CoordGeodetic.h"
#include "Ellipsoid.h"
#include "MagneticUtils.h"

#include <cmath>

namespace wmm {

    void CoordSpherical::toCartesian(double &x, double &y, double &z) const {
        const double radphi    = phig * (std::numbers::pi / 180);
        const double radlambda = lambdag * (std::numbers::pi / 180);

        x = r * cos(radphi) * cos(radlambda);
        y = r * cos(radphi) * sin(radlambda);
        z = r * sin(radphi);
    }

    /** @brief Convert geodetic coordinates, (defined by the WGS-84 reference ellipsoid), to Earth Centered Earth Fixed
     * Cartesian coordinates, and then to spherical coordinates.
     */
    void CoordSpherical::fromGeodetic(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic) {
        const double cosLat = std::cos(deg2Rad(coordGeodetic.phi));
        const double sinLat = std::sin(deg2Rad(coordGeodetic.phi));

        // compute the local radius of curvature on the WGS-84 reference ellipsoid
        const double rc = ellip.a / std::sqrt(1.0 - (ellip.epssq * sinLat * sinLat));

        // compute ECEF Cartesian coordinates of specified point (for longitude=0)
        const double xp = (rc + coordGeodetic.heightAboveEllipsoid) * cosLat;
        const double zp = ((rc * (1.0 - ellip.epssq)) + coordGeodetic.heightAboveEllipsoid) * sinLat;

        // compute spherical radius and angle lambdag and phi of specified point
        r = std::sqrt((xp * xp) + (zp * zp));
        // geocentric latitude
        phig = rad2Deg(std::asin(zp / r));
        // longitude
        lambdag = coordGeodetic.lambda;
    }

    CoordSpherical CoordSpherical::sphericalFromGeodetic(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic) {
        CoordSpherical sph;
        sph.fromGeodetic(ellip, coordGeodetic);
        return sph;
    }

}  // namespace wmm