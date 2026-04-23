//
// Created by Li-Yin Young on 8/19/22.
// Edited by Kiryuhin Viacheslav on 04/21/26
//
#include "magcalc.h"

#include "GeoMagneticElements.h"
#include "MagneticModel.h"

namespace wmm {

    /* The function is for point calculation for ewmm_point. The input of hight is already determined whether to
     * covert to Ellipsoid height
     */
    void point_calc(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, const CoordSpherical &coordSpherical,
                    const Date &userDate, const MagneticModel &magneticModel, GeoMagneticElements *geoMagneticElements,
                    GeoMagneticElements *errors) {
        // Time adjust the coefficients, Equation 19, WMM Technical report
        MagneticModel timedMagneticModel{magneticModel.applyDate(userDate)};

        // Computes the geoMagnetic field elements and their time change
        geoMagneticElements->geomag(ellip, coordSpherical, coordGeodetic, timedMagneticModel);

        geoMagneticElements->CalculateGridVariation(coordGeodetic);
#ifdef WMMHR
        errors->WMMHRErrorCalc(geoMagneticElements->H);
#else
        errors->WMMErrorCalc(geoMagneticElements->H);
#endif
    }

}  // namespace wmm