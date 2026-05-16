//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "Gradient.h"

#include "CoordGeodetic.h"
#include "CoordSpherical.h"

#include <cmath>

namespace wmm {

    /** @brief It should be noted that the x[2], y[2], and z[2] variables are NOT the same
     * coordinate system as the directions in which the gradients are taken.  These
     * variables represent a Cartesian coordinate system where the Earth's center is
     * the origin, 'z' points up toward the North (rotational) pole and 'x' points toward
     * the prime meridian.  'y' points toward longitude = 90 degrees East.
     * The gradient is preformed along a local Cartesian coordinate system with the
     * origin at CoordGeodetic.  'z' points down toward the Earth's core, x points
     * North, tangent to the local longitude line, and 'y' points East, tangent to
     * the local latitude line.
     **/
    void Gradient::calculate(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel) {
        constexpr double phiDelta = 0.01;
        constexpr double hDelta   = -1;
        double x1{};
        double y1{};
        double z1{};

        double x2{};
        double y2{};
        double z2{};

        // Initialization
        CoordSpherical adjCoordSpherical;
        adjCoordSpherical.fromGeodetic(ellip, coordGeodetic);

        GeoMagneticElements geomagneticElements;
        geomagneticElements.calculate(ellip, adjCoordSpherical, coordGeodetic, timedMagneticModel);
        CoordGeodetic adjCoordGeodetic{coordGeodetic};

        // Gradient along x
        adjCoordGeodetic.phi = coordGeodetic.phi + phiDelta;
        adjCoordSpherical.fromGeodetic(ellip, adjCoordGeodetic);

        GeoMagneticElements adjGeoMagneticElements1;
        adjGeoMagneticElements1.calculate(ellip, adjCoordSpherical, adjCoordGeodetic, timedMagneticModel);
        adjCoordSpherical.toCartesian(x1, y1, z1);

        adjCoordGeodetic.phi = coordGeodetic.phi - phiDelta;
        adjCoordSpherical.fromGeodetic(ellip, adjCoordGeodetic);

        GeoMagneticElements adjGeoMagneticElements2;
        adjGeoMagneticElements2.calculate(ellip, adjCoordSpherical, adjCoordGeodetic, timedMagneticModel);
        adjCoordSpherical.toCartesian(x2, y2, z2);

        double distance = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2)));
        gradPhi         = adjGeoMagneticElements1 - adjGeoMagneticElements2;
        gradPhi.scale(1 / distance);
        adjCoordGeodetic = coordGeodetic;

        // Gradient along y

        /*It is perhaps noticeable that the method here for calculation is substantially
         different than that for the gradient along x.  As we near the North pole
         the longitude lines approach each other, and the calculation that works well
         for latitude lines becomes unstable when 0.01 degrees represents sufficiently
         small numbers, and fails to function correctly at all at the North Pole */

        adjCoordSpherical.fromGeodetic(ellip, coordGeodetic);
        gradLambda.gradY(ellip, adjCoordSpherical, coordGeodetic, timedMagneticModel, geomagneticElements);

        // Gradient along z
        adjCoordGeodetic.heightAboveEllipsoid = coordGeodetic.heightAboveEllipsoid + hDelta;
        adjCoordGeodetic.heightAboveGeoid     = coordGeodetic.heightAboveGeoid + hDelta;
        adjCoordSpherical.fromGeodetic(ellip, adjCoordGeodetic);

        adjGeoMagneticElements1.calculate(ellip, adjCoordSpherical, adjCoordGeodetic, timedMagneticModel);
        adjCoordSpherical.toCartesian(x1, y1, z1);

        adjCoordGeodetic.heightAboveEllipsoid = coordGeodetic.heightAboveEllipsoid - hDelta;
        adjCoordGeodetic.heightAboveGeoid     = coordGeodetic.heightAboveGeoid - hDelta;
        adjCoordSpherical.fromGeodetic(ellip, adjCoordGeodetic);

        adjGeoMagneticElements2.calculate(ellip, adjCoordSpherical, adjCoordGeodetic, timedMagneticModel);
        adjCoordSpherical.toCartesian(x2, y2, z2);

        distance = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2)));
        gradZ    = adjGeoMagneticElements1 - adjGeoMagneticElements2;
        gradZ.scale(1 / distance);
    }
}  // namespace wmm