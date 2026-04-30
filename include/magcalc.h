//
// Created by Li-Yin Young on 8/19/22.
// Edited by Kiryuhin Viacheslav on 04/21/26
//

#pragma once

namespace wmm {
    struct Ellipsoid;
    struct CoordGeodetic;
    struct CoordSpherical;
    struct Date;
    struct MagneticModel;
    struct GeoMagneticElements;

    void point_calc(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, const CoordSpherical &coordSpherical,
                    const Date &userDate, const MagneticModel &magneticModel, GeoMagneticElements *geoMagneticElements,
                    GeoMagneticElements *errors);

    void point_calc(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, const Date &userDate,
                    const MagneticModel &magneticModel, GeoMagneticElements *geoMagneticElements,
                    GeoMagneticElements *errors);

}  // namespace wmm