//
// Created by Professional on 22.04.2026.
//

#include "Gradient.h"

#include "CoordGeodetic.h"
#include "CoordSpherical.h"

#include <cmath>

namespace wmm {

    /** It should be noted that the x[2], y[2], and z[2] variables are NOT the same
     * coordinate system as the directions in which the gradients are taken.  These
     * variables represent a Cartesian coordinate system where the Earth's center is
     * the origin, 'z' points up toward the North (rotational) pole and 'x' points toward
     * the prime meridian.  'y' points toward longitude = 90 degrees East.
     * The gradient is preformed along a local Cartesian coordinate system with the
     * origin at CoordGeodetic.  'z' points down toward the Earth's core, x points
     * North, tangent to the local longitude line, and 'y' points East, tangent to
     * the local latitude line.
     **/
    void Gradient::Calculate(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel) {
        double phiDelta = 0.01, /* DeltaY = 0.01, */ hDelta = -1, x[2], y[2], z[2], distance;

        CoordSpherical AdjCoordSpherical;
        CoordGeodetic AdjCoordGeodetic;
        GeoMagneticElements GeomagneticElements, AdjGeoMagneticElements[2];

        /* Initialization */
        AdjCoordSpherical.fromGeodetic(ellip, coordGeodetic);
        GeomagneticElements.geomag(ellip, AdjCoordSpherical, coordGeodetic, timedMagneticModel);
        AdjCoordGeodetic = coordGeodetic;


        /* Gradient along x */
        AdjCoordGeodetic.phi = coordGeodetic.phi + phiDelta;
        AdjCoordSpherical.fromGeodetic(ellip, AdjCoordGeodetic);
        AdjGeoMagneticElements[0].geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel);
        AdjCoordSpherical.toCartesian(x[0], y[0], z[0]);
        AdjCoordGeodetic.phi = coordGeodetic.phi - phiDelta;
        AdjCoordSpherical.fromGeodetic(ellip, AdjCoordGeodetic);
        AdjGeoMagneticElements[1].geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel);
        AdjCoordSpherical.toCartesian(x[1], y[1], z[1]);


        distance = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) + (z[0] - z[1]) * (z[0] - z[1]));
        GradPhi  = AdjGeoMagneticElements[0] - AdjGeoMagneticElements[1];
        GradPhi.scale(1 / distance);
        AdjCoordGeodetic = coordGeodetic;

        /*Gradient along y*/

        /*It is perhaps noticeable that the method here for calculation is substantially
         different than that for the gradient along x.  As we near the North pole
         the longitude lines approach each other, and the calculation that works well
         for latitude lines becomes unstable when 0.01 degrees represents sufficiently
         small numbers, and fails to function correctly at all at the North Pole */

        AdjCoordSpherical.fromGeodetic(ellip, coordGeodetic);
        GradLambda.gradY(ellip, AdjCoordSpherical, coordGeodetic, timedMagneticModel, GeomagneticElements);

        /*Gradient along z*/
        AdjCoordGeodetic.HeightAboveEllipsoid = coordGeodetic.HeightAboveEllipsoid + hDelta;
        AdjCoordGeodetic.HeightAboveGeoid     = coordGeodetic.HeightAboveGeoid + hDelta;
        AdjCoordSpherical.fromGeodetic(ellip, AdjCoordGeodetic);
        AdjGeoMagneticElements[0].geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel);
        AdjCoordSpherical.toCartesian(x[0], y[0], z[0]);
        AdjCoordGeodetic.HeightAboveEllipsoid = coordGeodetic.HeightAboveEllipsoid - hDelta;
        AdjCoordGeodetic.HeightAboveGeoid     = coordGeodetic.HeightAboveGeoid - hDelta;
        AdjCoordSpherical.fromGeodetic(ellip, AdjCoordGeodetic);
        AdjGeoMagneticElements[1].geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel);
        AdjCoordSpherical.toCartesian(x[1], y[1], z[1]);

        distance = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) + (z[0] - z[1]) * (z[0] - z[1]));
        GradZ    = AdjGeoMagneticElements[0] - AdjGeoMagneticElements[1];
        GradZ.scale(1 / distance);
        AdjCoordGeodetic = coordGeodetic;
    }
}  // namespace wmm