//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

namespace wmm {
    struct CoordGeodetic;
    struct Ellipsoid;

    struct UTMParameters {
        void getTransverseMercator(const CoordGeodetic &coordGeodetic);

        double easting{};   // (X) in meters
        double northing{};  // (Y) in meters
        int zone{};         // UTM zone
        char hemiSphere{};
        double centralMeridian{};
        double convergenceOfMeridians{};
        double pointScale{};

    private:
        static bool getUtmParameters(double latitude, double longitude, int *zone, char *hemisphere, double *centralMeridian);

        static void TMfwd4(const Ellipsoid &ellip, double Lam0, double falseE, double falseN, int XYonly, double Lambda,
                           double Phi, double *X, double *Y, double *pscale, double *CoM);
    };
}  // namespace wmm