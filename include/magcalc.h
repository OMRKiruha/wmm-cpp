//
// Created by Li-Yin Young on 8/19/22.
//

#pragma once

#include "GeomagnetismHeader.h"

namespace wmm {

    void point_calc(Ellipsoid ellip, CoordGeodetic coordGeodetic, CoordSpherical *coordSpherical, Date userDate,
                    MagneticModel *magneticModel, MagneticModel *timedMagneticModel,
                    GeoMagneticElements *geoMagneticElements, GeoMagneticElements *errors);

}