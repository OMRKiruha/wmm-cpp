//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

#include <array>
#include <cmath>

namespace wmm {
    struct Ellipsoid {
        // Defaults is WGS-84 parameters
        double a{6378.137};                         // semi-major axis of the ellipsoid
        double b{6356.7523142};                     // semi-minor axis of the ellipsoid
        double fla{1 / 298.257223563};              // flattening
        double eps{sqrt(1 - ((b * b) / (a * a)))};  // first eccentricity
        double epssq{eps * eps};                    // first eccentricity squared
        double re{6371.2};                          // mean radius of  ellipsoid
        double K0     = 0.9996;
        double K0R4   = 6367449.1458234153093 * K0;
        double K0R4oa = K0R4 / 6378137;
        std::array<double, 8> Acoeff{8.37731820624469723600E-04, 7.60852777357248641400E-07, 1.19764550324249124400E-09,
                                     2.42917068039708917100E-12, 5.71181837042801392800E-15, 1.47999793137966169400E-17,
                                     4.10762410937071532000E-20, 1.21078503892257704200E-22};
    };
}  // namespace wmm
