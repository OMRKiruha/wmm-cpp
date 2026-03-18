//
// Created by Li-Yin Young on 8/19/22.
//
#include "magcalc.h"

namespace wmm {

    /*
     * The function is for point calculation for ewmm_point. The input of hight is already determined whether to
     * covert to Ellipsoid height
     */
    void point_calc(Ellipsoid ellip, CoordGeodetic coordGeodetic, CoordSpherical *coordSpherical, Date userDate,
                    MagneticModel *magneticModel, MagneticModel *timedMagneticModel,
                    GeoMagneticElements *geoMagneticElements, GeoMagneticElements *errors) {
        TimelyModifyMagneticModel(userDate, *magneticModel,
                                  timedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */
        Geomag(ellip, *coordSpherical, coordGeodetic, timedMagneticModel,
               geoMagneticElements);                   /* Computes the geoMagnetic field elements and their time change*/
        CalculateGridVariation(coordGeodetic, geoMagneticElements);
#ifdef WMMHR
        WMMHRErrorCalc(geoMagneticElements->H, errors);
#else
        WMMErrorCalc(geoMagneticElements->H, errors);
#endif
    }

}  // namespace wmm