//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

namespace wmm {
    struct CoordGeodetic;
    struct Ellipsoid;

    struct UTMParameters {
        void getTransverseMercator(const CoordGeodetic &coordGeodetic);

        double easting{0.0};   // (X) in meters
        double northing{0.0};  // (Y) in meters
        int zone{0};           // UTM zone
        char hemiSphere{};
        double centralMeridian{0.0};
        double convergenceOfMeridians{0.0};
        double pointScale{0.0};

    private:
        static bool getUtmParameters(double latitude, double longitude, int &zone, char &hemisphere,
                                     double &centralMeridian);

        static void TMfwd4(const Ellipsoid &ellip, double Lam0, double falseE, double falseN, int XYonly, double Lambda,
                           double Phi, double &X, double &Y, double &pscale, double &CoM);
    };
}  // namespace wmm