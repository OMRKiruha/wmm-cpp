//
// Created by Professional on 22.04.2026.
//

#include "Geoid.h"

#include "MagneticUtils.h"
#include "EGM9615.h"

#include <cmath>
#include <memory>

namespace wmm {

    Geoid::Geoid() {
        GeoidHeightBuffer = std::make_unique<Geoid::GeoidHeightArray_t>(GeoidHeightsArray);
        Geoid_Initialized = 1;
    }

    /**
     * The  function GetGeoidHeight returns the height of the
     * EGM96 geiod above or below the WGS84 ellipsoid,
     * at the specified geodetic coordinates,
     * using a grid of height adjustments from the EGM96 gravity model.
     *
     *    Latitude            : Geodetic latitude in radians           (input)
     *    Longitude           : Geodetic longitude in radians          (input)
     *    DeltaHeight         : Height Adjustment, in meters.          (output)
     *    Geoid				  : Geoid with Geoid grid		   (input)
     */
    int Geoid::GetGeoidHeight(double Latitude, double Longitude, double *DeltaHeight) const {
        long Index;
        double ElevationSE, ElevationSW, ElevationNE, ElevationNW;
        double OffsetX, OffsetY;
        double PostX, PostY;
        double UpperY, LowerY;
        int Error_Code = 0;

        if(!Geoid_Initialized) {
            PrintError(5);
            return (false);
        }
        if((Latitude < -90) || (Latitude > 90)) {     /* Latitude out of range */
            Error_Code |= 1;
        }
        if((Longitude < -180) || (Longitude > 360)) { /* Longitude out of range */
            Error_Code |= 1;
        }

        if(!Error_Code) { /* no errors */
            /*  Compute X and Y Offsets into Geoid Height Array:                          */

            if(Longitude < 0.0) {
                OffsetX = (Longitude + 360.0) * ScaleFactor;
            } else {
                OffsetX = Longitude * ScaleFactor;
            }
            OffsetY = (90.0 - Latitude) * ScaleFactor;

            /*  Find Four Nearest Geoid Height Cells for specified Latitude, Longitude;   */
            /*  Assumes that (0,0) of Geoid Height Array is at Northwest corner:          */

            PostX = floor(OffsetX);
            if((PostX + 1) == NumbGeoidCols) {
                PostX--;
            }
            PostY = floor(OffsetY);
            if((PostY + 1) == NumbGeoidRows) {
                PostY--;
            }

            Index       = (long)(PostY * NumbGeoidCols + PostX);
            ElevationNW = (double)GeoidHeightBuffer->at(Index);
            ElevationNE = (double)GeoidHeightBuffer->at(Index + 1);

            Index       = (long)((PostY + 1) * NumbGeoidCols + PostX);
            ElevationSW = (double)GeoidHeightBuffer->at(Index);
            ElevationSE = (double)GeoidHeightBuffer->at(Index + 1);

            /*  Perform Bi-Linear Interpolation to compute Height above Ellipsoid:        */

            const double DeltaX = OffsetX - PostX;
            const double DeltaY = OffsetY - PostY;

            UpperY = ElevationNW + DeltaX * (ElevationNE - ElevationNW);
            LowerY = ElevationSW + DeltaX * (ElevationSE - ElevationSW);

            *DeltaHeight = UpperY + DeltaY * (LowerY - UpperY);
        } else {
            PrintError(17);
            return (false);
        }
        return true;
    }
}  // namespace wmm