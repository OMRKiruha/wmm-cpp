//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "Geoid.h"

#include "EGM9615.h"
#include "MagneticUtils.h"

#include <cmath>

namespace wmm {

    Geoid::Geoid() {
        geoidHeightBuffer = std::make_unique<Geoid::GeoidHeightArray_t>(GeoidHeightsArray);
        isGeoidInitialized = true;
    }

    /** @brief The function returns the height of the EGM96 geiod above or below
     * the WGS84 ellipsoid, at the specified geodetic coordinates, using a grid of height
     * adjustments from the EGM96 gravity model.
     */
    bool Geoid::getGeoidHeight(double latitude, double longitude, double &deltaHeight) const {
        bool isError{false};

        if(!isGeoidInitialized) {
            printError(5);
            return false;
        }
        if((latitude < -90) || (latitude > 90)) {      // latitude out of range
            isError = true;
        }
        if((longitude < -180) || (longitude > 360)) {  // longitude out of range
            isError = true;
        }

        if(!isError) {
            double offsetX{};
            double offsetY{};

            //  Compute X and Y Offsets into Geoid Height Array:
            if(longitude < 0.0) {
                offsetX = (longitude + 360.0) * scaleFactor;
            } else {
                offsetX = longitude * scaleFactor;
            }
            offsetY = (90.0 - latitude) * scaleFactor;

            //  Find Four Nearest Geoid Height Cells for specified latitude, longitude;
            //  Assumes that (0,0) of Geoid Height Array is at Northwest corner:
            double postX = floor(offsetX);
            if((postX + 1) == numbGeoidCols) {
                postX--;
            }
            double postY = floor(offsetY);
            if((postY + 1) == numbGeoidRows) {
                postY--;
            }

            long index               = static_cast<long>((postY * numbGeoidCols) + postX);
            const double elevationNW = geoidHeightBuffer->at(index);
            const double elevationNE = geoidHeightBuffer->at(index + 1);

            index                    = static_cast<long>(((postY + 1) * numbGeoidCols) + postX);
            const double elevationSW = geoidHeightBuffer->at(index);
            const double elevationSE = geoidHeightBuffer->at(index + 1);

            //  Perform Bi-Linear Interpolation to compute Height above Ellipsoid:
            const double deltaX = offsetX - postX;
            const double deltaY = offsetY - postY;

            const double upperY = elevationNW + (deltaX * (elevationNE - elevationNW));
            const double lowerY = elevationSW + (deltaX * (elevationSE - elevationSW));

            deltaHeight = upperY + deltaY * (lowerY - upperY);
        } else {
            printError(17);
            return false;
        }
        return true;
    }
}  // namespace wmm