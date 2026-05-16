//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "UTMParameters.h"

#include "CoordGeodetic.h"
#include "Ellipsoid.h"
#include "MagneticConstants.h"
#include "MagneticUtils.h"

#include <iostream>

namespace wmm {

    /** @brief Gets the UTM Parameters for a given Latitude and Longitude.
     */
    void UTMParameters::getTransverseMercator(const CoordGeodetic &coordGeodetic) {
        // Get the map projection  parameters
        const double lambda = deg2Rad(coordGeodetic.lambda);
        const double phi    = deg2Rad(coordGeodetic.phi);

        char hemisphere{};
        double Lam0{};
        getUtmParameters(phi, lambda, zone, hemisphere, Lam0);

        double falseN{};
        if(hemisphere == 'n' || hemisphere == 'N') {
            falseN = 0;
        }
        if(hemisphere == 's' || hemisphere == 'S') {
            falseN = 10000000;
        }
        constexpr double falseE = 500000;

        // Execution of the forward T.M. algorithm
        const Ellipsoid wgs86{};
        constexpr int XYonly = 0;
        double CoM{};
        TMfwd4(wgs86, Lam0, falseE, falseN, XYonly, lambda, phi, easting, northing, pointScale, CoM);

        // Report results
        hemiSphere             = hemisphere;
        centralMeridian        = rad2Deg(Lam0);  // Central Meridian of the UTM zone
        convergenceOfMeridians = rad2Deg(CoM);   // Convergence of meridians of the UTM zone and location
    }

    /** @brief The function getUtmParameters converts geodetic (latitude and longitude) coordinates
     * to UTM projection parameters (zone, hemisphere and central meridian)
     * If any errors occur, the error code(s) are returned by the function, otherwise true is returned.
     */
    bool UTMParameters::getUtmParameters(const double latitude, double longitude, int &zone, char &hemisphere,
                                         double &centralMeridian) {
        // latitude out of range
        if((latitude < deg2Rad(UTM_MIN_LAT_DEGREE)) || (latitude > deg2Rad(UTM_MAX_LAT_DEGREE))) {
            std::cerr << "\nError: Latitude out of range in getUtmParameters\n";
            return true;
        }

        // longitude out of range
        if((longitude < -std::numbers::pi) || (longitude > (2 * std::numbers::pi))) {
            std::cerr << "\nError: Longitude out of range in getUtmParameters\n";
            return true;
        }

        if(longitude < 0) {
            longitude += (2 * std::numbers::pi) + 1.0e-10;
        }
        const long Lat_Degrees  = static_cast<long>(rad2Deg(latitude));
        const long Long_Degrees = static_cast<long>(rad2Deg(longitude));

        long temp_zone{};
        if(longitude < std::numbers::pi) {
            temp_zone = static_cast<long>(31 + (rad2Deg(longitude) / 6.0));
        } else {
            temp_zone = static_cast<long>((rad2Deg(longitude) / 6.0) - 29);
        }

        if(temp_zone > 60) {
            temp_zone = 1;
        }

        // UTM special cases
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

        if(temp_zone >= 31) {
            centralMeridian = deg2Rad((6 * temp_zone) - 183);
        } else {
            centralMeridian = deg2Rad((6 * temp_zone) + 177);
        }

        zone = temp_zone;

        if(latitude < 0) {
            hemisphere = 'S';
        } else {
            hemisphere = 'N';
        }

        return false;
    }

    /** @brief Transverse Mercator forward equations including point-scale and CoM
     *   Algorithm developed by: C. Rollins   August 7, 2006
     *   C software written by:  K. Robins
     */
    void UTMParameters::TMfwd4(const Ellipsoid &ellip, const double Lam0, const double falseE, const double falseN,
                               const int XYonly, const double Lambda, const double Phi, double &X, double &Y, double &pscale,
                               double &CoM) {
        // Ellipsoid to sphere
        // Convert longitude (Greenwhich) to longitude from the central meridian
        // It is unnecessary to find the (-Pi, Pi] equivalent of the result.
        // Compute its cosine and sine.
        const double Lam    = Lambda - Lam0;
        const double cosLam = cos(Lam);
        const double sinLam = sin(Lam);

        // Latitude
        const double cosPhi = cos(Phi);
        const double sinPhi = sin(Phi);

        // Convert geodetic latitude, Phi, to conformal latitude, Chi
        // Only the cosine and sine of Chi are actually needed.
        const double P     = exp(ellip.eps * std::atanh(ellip.eps * sinPhi));
        const double part1 = (1 + sinPhi) / P;
        const double part2 = (1 - sinPhi) * P;
        const double denom = 1 / (part1 + part2);
        const double CChi  = 2 * cosPhi * denom;
        const double SChi  = (part1 - part2) * denom;

        // Sphere to first plane
        //  Apply spherical theory of transverse Mercator to get (u,v) coordinates
        //  Note the order of the arguments in Fortran's version of ArcTan, i.e.
        //            atan2(y, x) = ATan(y/x)
        //  The two argument form of ArcTan is needed here.
        const double T = CChi * sinLam;
        const double U = std::atanh(T);
        const double V = std::atan2(SChi, CChi * cosLam);

        // Trigonometric multiple angles
        // Compute Cosh of even multiples of U
        // Compute Sinh of even multiples of U
        // Compute Cos  of even multiples of V
        // Compute Sin  of even multiples of V
        const double Tsq    = T * T;
        const double denom2 = 1 / (1 - Tsq);
        const double c2u    = (1 + Tsq) * denom2;
        const double s2u    = 2 * T * denom2;
        const double c2v    = (-1 + (CChi * CChi * (1 + (cosLam * cosLam)))) * denom2;
        const double s2v    = 2 * cosLam * CChi * SChi * denom2;

        const double c4u = 1 + (2 * s2u * s2u);
        const double s4u = 2 * c2u * s2u;
        const double c4v = 1 - (2 * s2v * s2v);
        const double s4v = 2 * c2v * s2v;

        const double c6u = (c4u * c2u) + (s4u * s2u);
        const double s6u = (s4u * c2u) + (c4u * s2u);
        const double c6v = (c4v * c2v) - (s4v * s2v);
        const double s6v = (s4v * c2v) + (c4v * s2v);

        const double c8u = 1 + (2 * s4u * s4u);
        const double s8u = 2 * c4u * s4u;
        const double c8v = 1 - (2 * s4v * s4v);
        const double s8v = 2 * c4v * s4v;

        // First plane to second plane
        // Accumulate terms for X and Y
        double Xstar = ellip.Acoeff.at(3) * s8u * c8v;
        Xstar        = Xstar + (ellip.Acoeff.at(2) * s6u * c6v);
        Xstar        = Xstar + (ellip.Acoeff.at(1) * s4u * c4v);
        Xstar        = Xstar + (ellip.Acoeff.at(0) * s2u * c2v);
        Xstar        = Xstar + U;

        double Ystar = ellip.Acoeff.at(3) * c8u * s8v;
        Ystar        = Ystar + (ellip.Acoeff.at(2) * c6u * s6v);
        Ystar        = Ystar + (ellip.Acoeff.at(1) * c4u * s4v);
        Ystar        = Ystar + (ellip.Acoeff.at(0) * c2u * s2v);
        Ystar        = Ystar + V;

        // Apply isoperimetric radius, scale adjustment, and offsets
        X = (ellip.K0R4 * Xstar) + falseE;
        Y = (ellip.K0R4 * Ystar) + falseN;

        // Point-scale and CoM
        if(XYonly == 1) {
            pscale = ellip.K0;
            CoM    = 0;
        } else {
            double sig1 = 8 * ellip.Acoeff.at(3) * c8u * c8v;
            sig1        = sig1 + (6 * ellip.Acoeff.at(2) * c6u * c6v);
            sig1        = sig1 + (4 * ellip.Acoeff.at(1) * c4u * c4v);
            sig1        = sig1 + (2 * ellip.Acoeff.at(0) * c2u * c2v);
            sig1        = sig1 + 1;

            double sig2 = 8 * ellip.Acoeff.at(3) * s8u * s8v;
            sig2        = sig2 + (6 * ellip.Acoeff.at(2) * s6u * s6v);
            sig2        = sig2 + (4 * ellip.Acoeff.at(1) * s4u * s4v);
            sig2        = sig2 + (2 * ellip.Acoeff.at(0) * s2u * s2v);

            // Combined square roots
            const double comroo = sqrt((1 - (ellip.epssq * sinPhi * sinPhi)) * denom2 * ((sig1 * sig1) + (sig2 * sig2)));

            pscale = ellip.K0R4oa * 2 * denom * comroo;
            CoM    = atan2(SChi * sinLam, cosLam) + atan2(sig2, sig1);
        }
    }

}  // namespace wmm