//
// Created by Professional on 22.04.2026.
//

#pragma once

namespace wmm {
    struct CoordGeodetic;
    struct Ellipsoid;

    struct UTMParameters {
        void GetTransverseMercator(const CoordGeodetic &coordGeodetic);

        double Easting{};   // (X) in meters
        double Northing{};  // (Y) in meters
        int Zone{};         // UTM Zone
        char HemiSphere{};
        double CentralMeridian{};
        double ConvergenceOfMeridians{};
        double PointScale{};

    private:
        static int GetUTMParameters(double latitude, double longitude, int *zone, char *hemisphere, double *centralMeridian);

        static void TMfwd4(const Ellipsoid &ellip, double Lam0, double falseE, double falseN, int XYonly, double Lambda,
                           double Phi, double *X, double *Y, double *pscale, double *CoM);
    };
}  // namespace wmm