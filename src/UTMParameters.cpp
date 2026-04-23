//
// Created by Professional on 22.04.2026.
//

#include "UTMParameters.h"

#include "CoordGeodetic.h"
#include "Ellipsoid.h"
#include "MagneticConstants.h"
#include "MagneticUtils.h"

namespace wmm {
    /* Gets the UTM Parameters for a given Latitude and Longitude.
     *INPUT: CoordGeodetic
     *OUTPUT : UTMParameters
     */
    void UTMParameters::GetTransverseMercator(const CoordGeodetic &coordGeodetic) {
        double Lam0{};
        double falseN{};
        double X{};
        double Y{};
        double pscale{};
        double CoM{};
        int zone{};
        char hemisphere{};

        /** Get the map projection  parameters */
        const double lambda = Deg2Rad(coordGeodetic.lambda);
        const double phi    = Deg2Rad(coordGeodetic.phi);

        GetUTMParameters(phi, lambda, &zone, &hemisphere, &Lam0);

        if(hemisphere == 'n' || hemisphere == 'N') {
            falseN = 0;
        }
        if(hemisphere == 's' || hemisphere == 'S') {
            falseN = 10000000;
        }
        constexpr double falseE = 500000;

        /** Execution of the forward T.M. algorithm */
        constexpr int XYonly = 0;
        Ellipsoid wgs86{};
        TMfwd4(wgs86, Lam0, falseE, falseN, XYonly, lambda, phi, &X, &Y, &pscale, &CoM);

        /** Report results */
        Easting                = X;             /* UTM Easting (X) in meters */
        Northing               = Y;             /* UTM Northing (Y) in meters */
        Zone                   = zone;          /* UTM Zone */
        HemiSphere             = hemisphere;
        CentralMeridian        = Rad2Deg(Lam0); /* Central Meridian of the UTM Zone */
        ConvergenceOfMeridians = Rad2Deg(CoM);  /* Convergence of meridians of the UTM Zone and location */
        PointScale             = pscale;
    }

    /** The function GetUTMParameters converts geodetic (latitude and
     * longitude) coordinates to UTM projection parameters (zone, hemisphere and central meridian)
     * If any errors occur, the error code(s) are returned
     * by the function, otherwise true is returned.
     *
     *    latitude          : latitude in radians                 (input)
     *    longitude         : longitude in radians                (input)
     *    zone              : UTM zone                            (output)
     *    hemisphere        : North or South hemisphere           (output)
     *    centralMeridian	: Central Meridian of the UTM zone in radians	   (output)
     */
    int UTMParameters::GetUTMParameters(double latitude, double longitude, int *zone, char *hemisphere,
                                        double *centralMeridian) {
        long Lat_Degrees;
        long Long_Degrees;
        long temp_zone;
        int Error_Code = 0;

        if((latitude < Deg2Rad(UTM_MIN_LAT_DEGREE)) ||
           (latitude > Deg2Rad(UTM_MAX_LAT_DEGREE))) { /* latitude out of range */
            PrintError(23);
            Error_Code = 1;
        }
        if((longitude < -M_PI) || (longitude > (2 * M_PI))) { /* longitude out of range */
            PrintError(24);
            Error_Code = 1;
        }
        if(!Error_Code) { /* no errors */
            if(longitude < 0) {
                longitude += (2 * M_PI) + 1.0e-10;
            }
            Lat_Degrees  = (long)(latitude * 180.0 / M_PI);
            Long_Degrees = (long)(longitude * 180.0 / M_PI);

            if(longitude < M_PI) {
                temp_zone = (long)(31 + ((longitude * 180.0 / M_PI) / 6.0));
            } else {
                temp_zone = (long)(((longitude * 180.0 / M_PI) / 6.0) - 29);
            }
            if(temp_zone > 60) {
                temp_zone = 1;
            }
            /* UTM special cases */
            if((Lat_Degrees > 55) && (Lat_Degrees < 64) && (Long_Degrees > -1) && (Long_Degrees < 3)) {
                temp_zone = 31;
            }
            if((Lat_Degrees > 55) && (Lat_Degrees < 64) && (Long_Degrees > 2) && (Long_Degrees < 12)) {
                temp_zone = 32;
            }
            if((Lat_Degrees > 71) && (Long_Degrees > -1) && (Long_Degrees < 9)) {
                temp_zone = 31;
            }
            if((Lat_Degrees > 71) && (Long_Degrees > 8) && (Long_Degrees < 21)) {
                temp_zone = 33;
            }
            if((Lat_Degrees > 71) && (Long_Degrees > 20) && (Long_Degrees < 33)) {
                temp_zone = 35;
            }
            if((Lat_Degrees > 71) && (Long_Degrees > 32) && (Long_Degrees < 42)) {
                temp_zone = 37;
            }

            if(!Error_Code) {
                if(temp_zone >= 31) {
                    *centralMeridian = (6 * temp_zone - 183) * M_PI / 180.0;
                } else {
                    *centralMeridian = (6 * temp_zone + 177) * M_PI / 180.0;
                }
                *zone = temp_zone;
                if(latitude < 0) {
                    *hemisphere = 'S';
                } else {
                    *hemisphere = 'N';
                }
            }
        } /* END OF if (!Error_Code) */
        return (Error_Code);
    } /* GetUTMParameters */

    /**  Transverse Mercator forward equations including point-scale and CoM
     *   Algorithm developed by: C. Rollins   August 7, 2006
     *   C software written by:  K. Robins
     *
     *        Constants fixed by choice of ellipsoid and choice of projection parameters
     *          Eps          Eccentricity (epsilon) of the ellipsoid
     *          Epssq        Eccentricity squared
     *        ( R4           Meridional isoperimetric radius   )
     *        ( K0           Central scale factor              )
     *          K0R4         K0 times R4
     *          K0R4oa       K0 times Ratio of R4 over semi-major axis
     *          Acoeff       Trig series coefficients, omega as a function of chi
     *          Lam0         Longitude of the central meridian in radians
     *          K0           Central scale factor, for example, 0.9996 for UTM
     *          falseE       False easting, for example, 500000 for UTM
     *          falseN       False northing
     *
     *   Processing option
     *          XYonly       If one (1), then only X and Y will be properly computed.  Values returned for point-scale
     *                       and CoM will merely be the trivial values for points on the central meridian
     *
     *   Input items that identify the point to be converted
     *          Lambda       Longitude (from Greenwich) in radians
     *          Phi          Latitude in radians
     *
     *   Output items
     *          X            X coordinate (Easting) in meters
     *          Y            Y coordinate (Northing) in meters
     *          pscale       point-scale (dimensionless)
     *          CoM          Convergence-of-meridians in radians
     */
    void UTMParameters::TMfwd4(const Ellipsoid &ellip, const double Lam0, const double falseE, const double falseN,
                               const int XYonly, const double Lambda, const double Phi, double *X, double *Y, double *pscale,
                               double *CoM) {
        /** Ellipsoid to sphere
         *  Convert longitude (Greenwhich) to longitude from the central meridian
         *  It is unnecessary to find the (-Pi, Pi] equivalent of the result.
         *  Compute its cosine and sine.         */
        const double Lam  = Lambda - Lam0;
        const double CLam = cos(Lam);
        const double SLam = sin(Lam);

        /** Latitude  */
        const double CPhi = cos(Phi);
        const double SPhi = sin(Phi);

        /** Convert geodetic latitude, Phi, to conformal latitude, Chi
         *  Only the cosine and sine of Chi are actually needed.        */
        const double P     = exp(ellip.eps * std::atanh(ellip.eps * SPhi));
        const double part1 = (1 + SPhi) / P;
        const double part2 = (1 - SPhi) * P;
        const double denom = 1 / (part1 + part2);
        const double CChi  = 2 * CPhi * denom;
        const double SChi  = (part1 - part2) * denom;

        /** Sphere to first plane
         *  Apply spherical theory of transverse Mercator to get (u,v) coordinates
         *  Note the order of the arguments in Fortran's version of ArcTan, i.e.
         *            atan2(y, x) = ATan(y/x)
         *  The two argument form of ArcTan is needed here. */
        const double T = CChi * SLam;
        const double U = std::atanh(T);
        const double V = std::atan2(SChi, CChi * CLam);

        /** Trigonometric multiple angles
         *  Compute Cosh of even multiples of U
         *  Compute Sinh of even multiples of U
         *  Compute Cos  of even multiples of V
         *  Compute Sin  of even multiples of V */
        const double Tsq    = T * T;
        const double denom2 = 1 / (1 - Tsq);
        const double c2u    = (1 + Tsq) * denom2;
        const double s2u    = 2 * T * denom2;
        const double c2v    = (-1 + CChi * CChi * (1 + CLam * CLam)) * denom2;
        const double s2v    = 2 * CLam * CChi * SChi * denom2;

        const double c4u = 1 + 2 * s2u * s2u;
        const double s4u = 2 * c2u * s2u;
        const double c4v = 1 - 2 * s2v * s2v;
        const double s4v = 2 * c2v * s2v;

        const double c6u = c4u * c2u + s4u * s2u;
        const double s6u = s4u * c2u + c4u * s2u;
        const double c6v = c4v * c2v - s4v * s2v;
        const double s6v = s4v * c2v + c4v * s2v;

        const double c8u = 1 + 2 * s4u * s4u;
        const double s8u = 2 * c4u * s4u;
        const double c8v = 1 - 2 * s4v * s4v;
        const double s8v = 2 * c4v * s4v;

        /** First plane to second plane
         *  Accumulate terms for X and Y */
        double Xstar = ellip.Acoeff[3] * s8u * c8v;
        Xstar        = Xstar + ellip.Acoeff[2] * s6u * c6v;
        Xstar        = Xstar + ellip.Acoeff[1] * s4u * c4v;
        Xstar        = Xstar + ellip.Acoeff.at(0) * s2u * c2v;
        Xstar        = Xstar + U;

        double Ystar = ellip.Acoeff[3] * c8u * s8v;
        Ystar        = Ystar + ellip.Acoeff[2] * c6u * s6v;
        Ystar        = Ystar + ellip.Acoeff[1] * c4u * s4v;
        Ystar        = Ystar + ellip.Acoeff.at(0) * c2u * s2v;
        Ystar        = Ystar + V;

        /** Apply isoperimetric radius, scale adjustment, and offsets  */
        *X = ellip.K0R4 * Xstar + falseE;
        *Y = ellip.K0R4 * Ystar + falseN;

        /** Point-scale and CoM */
        if(XYonly == 1) {
            *pscale = ellip.K0;
            *CoM    = 0;
        } else {
            double sig1 = 8 * ellip.Acoeff[3] * c8u * c8v;
            sig1        = sig1 + 6 * ellip.Acoeff[2] * c6u * c6v;
            sig1        = sig1 + 4 * ellip.Acoeff[1] * c4u * c4v;
            sig1        = sig1 + 2 * ellip.Acoeff.at(0) * c2u * c2v;
            sig1        = sig1 + 1;

            double sig2 = 8 * ellip.Acoeff[3] * s8u * s8v;
            sig2        = sig2 + 6 * ellip.Acoeff[2] * s6u * s6v;
            sig2        = sig2 + 4 * ellip.Acoeff[1] * s4u * s4v;
            sig2        = sig2 + 2 * ellip.Acoeff.at(0) * s2u * s2v;

            /*    Combined square roots  */
            const double comroo = sqrt((1 - ellip.epssq * SPhi * SPhi) * denom2 * (sig1 * sig1 + sig2 * sig2));

            *pscale = ellip.K0R4oa * 2 * denom * comroo;
            *CoM    = atan2(SChi * SLam, CLam) + atan2(sig2, sig1);
        }
    }

}  // namespace wmm