//
// Created by Professional on 22.04.2026.
//

#include "CoordSpherical.h"

#include "CoordGeodetic.h"
#include "Ellipsoid.h"
#include "MagneticUtils.h"

#include <cmath>

namespace wmm {

    void CoordSpherical::toCartesian(double &x, double &y, double &z) const {
        double radphi;
        double radlambda;

        radphi    = phig * (M_PI / 180);
        radlambda = lambdag * (M_PI / 180);

        x = r * cos(radphi) * cos(radlambda);
        y = r * cos(radphi) * sin(radlambda);
        z = r * sin(radphi);
    }

    /** Convert geodetic coordinates, (defined by the WGS-84 reference ellipsoid), to Earth Centered Earth Fixed
     * Cartesian coordinates, and then to spherical coordinates.
     * INPUT   Ellip  data  structure with the following elements
     *         CoordGeodetic  Pointer to the  data  structure with the following elements updates
     * OUTPUT  CoordSpherical 	Pointer to the data structure with the following elements
     */
    int CoordSpherical::fromGeodetic(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic) {
        const double CosLat = std::cos(Deg2Rad(coordGeodetic.phi));
        const double SinLat = std::sin(Deg2Rad(coordGeodetic.phi));

        /* compute the local radius of curvature on the WGS-84 reference ellipsoid */
        const double rc = ellip.a / std::sqrt(1.0 - ellip.epssq * SinLat * SinLat);

        /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */
        const double xp = (rc + coordGeodetic.HeightAboveEllipsoid) * CosLat;
        const double zp = (rc * (1.0 - ellip.epssq) + coordGeodetic.HeightAboveEllipsoid) * SinLat;

        /* compute spherical radius and angle lambdag and phi of specified point */
        r = std::sqrt(xp * xp + zp * zp);
        /* geocentric latitude */
        phig = Rad2Deg(std::asin(zp / r));
        /* longitude */
        lambdag = coordGeodetic.lambda;
        return true;
    }

}  // namespace wmm