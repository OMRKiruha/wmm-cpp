//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "LegendreFunction.h"

#include "CoordSpherical.h"
#include "MagneticUtils.h"

#include <cmath>
#include <iostream>

namespace wmm {
    LegendreFunction::LegendreFunction(const CoordSpherical &coordSpherical, const int nMax) : nMax(nMax) {
        const int numTerms = ((nMax + 1) * (nMax + 2) / 2);
        Pcup.resize(numTerms + 1);
        dPcup.resize(numTerms + 1);

        // Computes all the Schmidt-semi normalized associated Legendre functions up to degree nMax.
        // If nMax <= 16, function PcupLow is used. Otherwise, PcupHigh is called.
        const double sin_phi = sin(deg2Rad(coordSpherical.phig));  // sin  (geocentric latitude)

        if(nMax <= 16 || (1 - fabs(sin_phi)) < 1.0e-10) {          // If nMax is less tha 16 or at the poles
            PcupLow(sin_phi);
        } else {
            try {
                PcupHigh(sin_phi);
            } catch(const std::exception &e) {
                std::cerr << e.what();
            }
        }
    }

    /** @brief This function evaluates all of the Schmidt-semi normalized associated Legendre functions up to degree nMax.
     * The functions are initially scaled by 10^280 sin^m in order to minimize the effects of underflow at large m near the
     * poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). Note that this function performs the same operation
     * as PcupLow. However this function also can be used for high degree (large nMax) models.
     *
     * Notes:
     *  Adopted from the FORTRAN code written by Mark Wieczorek September 25, 2005.
     *  Manoj Nair, Nov, 2009 Manoj.C.Nair\@Noaa.Gov
     *
     *  Change from the previous version
     *  The prevous version computes the derivatives as dP(n,m)(x)/dx, where x = sin(latitude) (or cos(colatitude) ).
     *  However, the WMM Geomagnetic routines requires dP(n,m)(x)/dlatitude. Hence the derivatives are multiplied by
     * sin(latitude). Removed the options for CS phase and normalizations.
     *
     *  Note: In geomagnetism, the derivatives of ALF are usually found with respect to the colatitudes. Here the derivatives
     * are found with respect to the latitude. The difference is a sign reversal for the derivative of the Associated
     * Legendre Functions.
     *
     *  The derivatives can't be computed for latitude = |90| degrees.
     */
    void LegendreFunction::PcupHigh(const double x) {
        constexpr double scalef = 1.0e-280;

        if(fabs(x) == 1.0) {
            throw std::logic_error("Error in PcupHigh: derivative cannot be calculated at poles\n");
        }

        const double z = sqrt((1.0 - x) * (1.0 + x));
        if(z == 0) {
            throw std::logic_error("Error in PcupHigh: derivative cannot be calculated. z == 0 \n");
        }

        const int numTerms = ((nMax + 1) * (nMax + 2) / 2);
        std::vector<double> f1(numTerms + 1);
        std::vector<double> preSqr(numTerms + 1);
        std::vector<double> f2(numTerms + 1);

        for(int n = 0; n <= 2 * nMax + 1; ++n) {
            preSqr.at(n) = std::sqrt(n);
        }

        int k = 2;
        for(int n = 2; n <= nMax; n++) {
            ++k;
            f1[k] = (2. * n - 1) / static_cast<double>(n);
            f2[k] = static_cast<double>(n - 1) / static_cast<double>(n);
            for(int m = 1; m <= n - 2; m++) {
                ++k;
                f1[k] = (2. * n - 1) / preSqr[n + m] / preSqr[n - m];
                f2[k] = preSqr[n - m - 1] * preSqr[n + m - 1] / preSqr[n + m] / preSqr[n - m];
            }
            k += 2;
        }

        // z = sin (geocentric latitude)
        double pm2  = 1.0;
        Pcup.at(0)  = 1.0;
        dPcup.at(0) = 0.0;
        if(nMax == 0) {
            throw std::logic_error("Error in PcupHigh: derivative cannot be calculated. nMax == 0 \n");
        }
        double pm1  = x;
        Pcup.at(1)  = pm1;
        dPcup.at(1) = z;
        k           = 1;

        for(int n = 2; n <= nMax; n++) {
            k                = k + n;
            const double plm = (f1.at(k) * x * pm1) - (f2.at(k) * pm2);
            Pcup.at(k)       = plm;
            dPcup.at(k)      = static_cast<double>(n) * (pm1 - x * plm) / z;
            pm2              = pm1;
            pm1              = plm;
        }

        double pmm      = preSqr[2] * scalef;
        double rescalem = 1.0 / scalef;
        int kstart      = 0;

        int m{};
        for(m = 1; m <= nMax - 1; ++m) {
            rescalem = rescalem * z;

            // Calculate Pcup(m, m)
            kstart           = kstart + m + 1;
            pmm              = pmm * preSqr.at((2 * m) + 1) / preSqr.at(2 * m);
            Pcup.at(kstart)  = pmm * rescalem / preSqr.at((2 * m) + 1);
            dPcup.at(kstart) = -(static_cast<double>(m) * x * Pcup.at(kstart) / z);
            pm2              = pmm / preSqr.at((2 * m) + 1);

            // Calculate Pcup(m + 1, m)
            k           = kstart + m + 1;
            pm1         = x * preSqr.at(2 * m + 1) * pm2;
            Pcup.at(k)  = pm1 * rescalem;
            dPcup.at(k) = ((pm2 * rescalem) * preSqr.at(2 * m + 1) - x * static_cast<double>(m + 1) * Pcup.at(k)) / z;

            // Calculate Pcup(n, m)
            for(int n = m + 2; n <= nMax; ++n) {
                k                = k + n;
                const double plm = (x * f1.at(k) * pm1) - (f2.at(k) * pm2);
                Pcup.at(k)       = plm * rescalem;
                dPcup.at(k) =
                    (preSqr.at(n + m) * preSqr.at(n - m) * (pm1 * rescalem) - static_cast<double>(n) * x * Pcup.at(k)) / z;
                pm2 = pm1;
                pm1 = plm;
            }
        }

        // Calculate Pcup(nMax,nMax)*/
        rescalem         = rescalem * z;
        kstart           = kstart + m + 1;
        pmm              = pmm / preSqr.at(2 * nMax);
        Pcup.at(kstart)  = pmm * rescalem;
        dPcup.at(kstart) = -static_cast<double>(nMax) * x * Pcup.at(kstart) / z;
    }

    /**  This function evaluates all of the Schmidt-semi normalized associated Legendre
     *   functions up to degree nMax.
     *
     *   Notes: Overflow may occur if nMax > 20 , especially for high-latitudes. Use PcupHigh for large nMax.
     *   Written by Manoj Nair, June, 2009 . Manoj.C.Nair@Noaa.Gov.
     *
     *   Note: In geomagnetism, the derivatives of ALF are usually found with respect to the colatitudes. Here the
     *   derivatives are found with respect to the latitude. The difference is a sign reversal for the derivative of
     *   the Associated Legendre Functions.
     */
    void LegendreFunction::PcupLow(const double x) {
        Pcup.at(0)  = 1.0;
        dPcup.at(0) = 0.0;
        // sin (geocentric latitude) - sin_phi
        const double z = sqrt((1.0 - x) * (1.0 + x));

        const int NumTerms = ((nMax + 1) * (nMax + 2) / 2);
        std::vector<double> schmidtQuasiNorm(NumTerms + 1);

        //	 First,	Compute the Gauss-normalized associated Legendre functions
        for(int n = 1; n <= nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = ((n * (n + 1) / 2) + m);
                int index1{0};
                if(n == m) {
                    index1          = (n - 1) * n / 2 + m - 1;
                    Pcup.at(index)  = z * Pcup.at(index1);
                    dPcup.at(index) = z * dPcup.at(index1) + x * Pcup.at(index1);
                } else if(n == 1 && m == 0) {
                    index1          = (n - 1) * n / 2 + m;
                    Pcup.at(index)  = x * Pcup.at(index1);
                    dPcup.at(index) = x * dPcup.at(index1) - z * Pcup.at(index1);
                } else if(n > 1 && n != m) {
                    index1           = (n - 2) * (n - 1) / 2 + m;
                    const int index2 = ((n - 1) * n / 2) + m;
                    if(m > n - 2) {
                        Pcup.at(index)  = x * Pcup.at(index2);
                        dPcup.at(index) = x * dPcup.at(index2) - z * Pcup.at(index2);
                    } else {
                        const double k  = static_cast<double>(((n - 1) * (n - 1)) - (m * m)) / ((2. * n - 1) * (2. * n - 3));
                        Pcup.at(index)  = x * Pcup[index2] - k * Pcup[index1];
                        dPcup.at(index) = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];
                    }
                }
            }
        }

        // Compute the ration between the the Schmidt quasi-normalized associated Legendre
        // functions and the Gauss-normalized version.
        schmidtQuasiNorm.at(0) = 1.0;
        for(int n = 1; n <= nMax; n++) {
            int index  = (n * (n + 1) / 2);
            int index1 = (n - 1) * n / 2;

            // for m = 0
            schmidtQuasiNorm.at(index) = schmidtQuasiNorm.at(index1) * (2. * n - 1) / static_cast<double>(n);

            for(int m = 1; m <= n; m++) {
                index  = (n * (n + 1) / 2 + m);
                index1 = (n * (n + 1) / 2 + m - 1);
                schmidtQuasiNorm.at(index) =
                    schmidtQuasiNorm.at(index1) * sqrt(((n - m + 1) * (m == 1 ? 2. : 1.)) / static_cast<double>(n + m));
            }
        }

        // Converts the  Gauss-normalized associated Legendre functions to the Schmidt quasi-normalized version using
        // pre-computed relation stored in the variable schmidtQuasiNorm
        for(int n = 1; n <= nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = ((n * (n + 1) / 2) + m);
                Pcup.at(index)  = Pcup.at(index) * schmidtQuasiNorm.at(index);
                dPcup.at(index) = -dPcup.at(index) * schmidtQuasiNorm.at(index);
                // The sign is changed since the new WMM routines use derivative with respect to latitude instead of
                // co-latitude
            }
        }
    }

}  // namespace wmm