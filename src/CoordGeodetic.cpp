//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "CoordGeodetic.h"

#include "CoordSpherical.h"
#include "Ellipsoid.h"
#include "Geoid.h"
#include "MagneticConstants.h"
#include "MagneticUtils.h"

#include <algorithm>
#include <cmath>

namespace wmm {

    /** @brief This converts spherical coordinates back to geodetic coordinates.  It is not used in
     * the WMM but may be necessary for some applications, such as geomagnetic coordinates
     */
    [[maybe_unused]] void CoordGeodetic::fromSpherical(const Ellipsoid &ellip, const CoordSpherical &coordSpherical) {
        double x{0.0};
        double y{0.0};
        double z{0.0};

        coordSpherical.toCartesian(x, y, z);
        this->fromCartesian(ellip, x, y, z);
    }

    /** @brief This converts the Cartesian x, y, and z coordinates to Geodetic Coordinates
     * x is defined as the direction pointing out of the core toward the point defined
     * by 0 degrees latitude and longitude.
     * y is defined as the direction from the core toward 90 degrees east longitude along
     * the equator
     * z is defined as the direction from the core out the geographic north pole
     */
    void CoordGeodetic::fromCartesian(const Ellipsoid &ellip, const double x, const double y, const double z) {
        // 1.0 compute semi-minor axis and set sign to that of z in order to get sign of Phi correct
        double modified_b{0.0};
        if(z < 0.0) {
            modified_b = -ellip.b;
        } else {
            modified_b = ellip.b;
        }

        // 2.0 compute intermediate values for latitude
        const double r = std::sqrt((x * x) + (y * y));
        const double e = ((modified_b * z) - ((ellip.a * ellip.a) - (modified_b * modified_b))) / (ellip.a * r);
        const double f = ((modified_b * z) + ((ellip.a * ellip.a) - (modified_b * modified_b))) / (ellip.a * r);

        // 3.0 find solution to: t^4 + 2*E*t^3 + 2*F*t - 1 = 0
        const double p = (4.0 / 3.0) * ((e * f) + 1.0);
        const double q = 2.0 * ((e * e) - (f * f));
        const double d = (p * p * p) + (q * q);

        double v{};
        if(d >= 0.0) {
            v = pow((sqrt(d) - q), (1.0 / 3.0)) - pow((sqrt(d) + q), (1.0 / 3.0));
        } else {
            v = 2.0 * sqrt(-p) * cos(acos(q / (p * sqrt(-p))) / 3.0);
        }

        // 4.0 improve v NOTE: not really necessary unless point is near pole
        if(v * v < fabs(p)) {
            v = -((v * v * v) + (2.0 * q)) / (3.0 * p);
        }
        const double g    = (sqrt((e * e) + v) + e) / 2.0;
        const double t    = sqrt((g * g) + ((f - (v * g)) / ((2.0 * g) - e))) - g;
        const double rlat = atan((ellip.a * (1.0 - (t * t))) / (2.0 * modified_b * t));

        phi = rad2Deg(rlat);

        // 5.0 compute height above ellipsoid
        heightAboveEllipsoid = ((r - (ellip.a * t)) * cos(rlat)) + ((z - modified_b) * sin(rlat));

        // 6.0 compute longitude east of Greenwich
        double zlong = atan2(y, x);
        if(zlong < 0.0) {
            zlong = zlong + (2 * M_PI);
        }

        lambda = rad2Deg(zlong);
        while(lambda > 180) {
            lambda -= 360;
        }
    }

    /** @brief Check if the latitude is equal to -90 or 90. If it is, offset it by 1e-5 to avoid
     * division by zero. This is not currently used in the Geomagnetic main function. This may be
     * used to avoid calling SummationSpecial.
     */
    [[maybe_unused]] void CoordGeodetic::checkGeographicPole() {
        phi = std::clamp(phi, -90.0 + GEO_POLE_TOLERANCE, 90.0 - GEO_POLE_TOLERANCE);
    }

    /** @brief The function converts the specified WGS84 geoid height at the specified geodetic
     * coordinates to the equivalent ellipsoid height, using the EGM96 gravity model.
     */
    bool CoordGeodetic::convertGeoidToEllipsoidHeight(const Geoid &geoid) {
        bool isError{};

        if(isUseGeoid) {
            double lon{};
            double lat{};
            double deltaHeight{};
            // Geoid correction required
            // To ensure that latitude is less than 90 call EquivalentLatLon()
            equivalentLatLon(phi, lambda, lat, lon);
            isError              = geoid.getGeoidHeight(lat, lon, deltaHeight);
            heightAboveEllipsoid = heightAboveGeoid + (deltaHeight / 1000);  //  Input and output
            // should be kilometers, However getGeoidHeight returns Geoid height in meters - Hence division by 1000
        } else  // Geoid correction not required, copy the MSL height to Ellipsoid height
        {
            heightAboveEllipsoid = heightAboveGeoid;
            isError              = true;
        }
        return isError;
    }

    /** @brief This function takes a latitude and longitude that are ordinarily out of range and gives
     * in range values that are equivalent on the Earth's surface. This is required to get correct values
     * for the geoid function.
     **/
    void CoordGeodetic::equivalentLatLon(const double lat, const double lon, double &repairedLat, double &repairedLon) {
        double colat = 90 - lat;
        repairedLon  = lon;

        if(colat < 0) {
            colat = -colat;
        }

        while(colat > 360) {
            colat -= 360;
        }

        if(colat > 180) {
            colat -= 180;
            repairedLon += 180;
        }

        repairedLat = 90 - colat;

        if(repairedLon > 360) {
            repairedLon -= 360;
        }

        if(repairedLon < -180) {
            repairedLon += 360;
        }
    }
}  // namespace wmm