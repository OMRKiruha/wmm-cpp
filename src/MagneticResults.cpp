//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "MagneticResults.h"

#include "CoordGeodetic.h"
#include "CoordSpherical.h"
#include "LegendreFunction.h"
#include "MagneticModel.h"
#include "MagneticUtils.h"
#include "SphericalHarmonicVariables.h"

#include <cmath>

namespace wmm {
    /** Computes Geomagnetic Field Elements X, Y and Z in Spherical coordinate system using spherical harmonic summation.
     *
     * The vector Magnetic field is given by -grad V, where V is Geomagnetic scalar potential
     * The gradient in spherical coordinates is given by:
     *
     *          dV ^     1 dV ^        1     dV ^
     * grad V = -- r  +  - -- t  +  -------- -- p
     *          dr       r dt       r sin(t) dp
     *
     * Manoj Nair, June, 2009 Manoj.C.Nair\@Noaa.Gov
     */
    void MagneticResults::summation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                                    const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical) {
        Bz = 0.0;
        By = 0.0;
        Bx = 0.0;

        for(int n = 1; n <= magneticModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = ((n * (n + 1) / 2) + m);

                //		    nMax  	(n+2) 	  n     m            m           m
                //  Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                //            n=1      	      m=0   n            n           n
                // Equation 12 in the WMM Technical report.  Derivative with respect to radius.
                Bz -= sphVariables.relativeRadiusPower.at(n) *
                      (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.main_Field_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      (n + 1) * legendreFunction.Pcup.at(index);

                //		  1 nMax  (n+2)    n     m            m           m
                // By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                //         n=1             m=0   n            n           n
                // Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.
                By += sphVariables.relativeRadiusPower.at(n) *
                      (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.sin_mlambda.at(m) -
                       magneticModel.main_Field_Coeff_H.at(index) * sphVariables.cos_mlambda.at(m)) *
                      m * legendreFunction.Pcup.at(index);

                //		   nMax  (n+2) n     m            m           m
                //  Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                //         n=1         m=0   n            n           n
                // Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius.
                Bx -= sphVariables.relativeRadiusPower.at(n) *
                      (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.main_Field_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      legendreFunction.dPcup.at(index);
            }
        }

        const double cos_phi = cos(deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            By = By / cos_phi;
        } else
        // Special calculation for component - By - at Geographic poles.
        // If the user wants to avoid using this function,  please make sure that the latitude
        // is not exactly +/-90. An option is to make use the function CheckGeographicPoles.
        {
            summationSpecial(magneticModel, sphVariables, coordSpherical);
        }
    }

    /** @brief Special calculation for the component By at Geographic poles.
     * Manoj Nair, June, 2009 manoj.c.nair\@noaa.gov
     * See Section 1.4, "SINGULARITIES AT THE GEOGRAPHIC POLES", WMM Technical report
     */
    void MagneticResults::summationSpecial(const MagneticModel &magneticModel,
                                           const SphericalHarmonicVariables &sphVariables,
                                           const CoordSpherical &coordSpherical) {
        std::vector<double> PcupS(magneticModel.nMax + 1);

        PcupS.at(0)              = 1;
        double schmidtQuasiNorm1 = 1.0;

        By                   = 0.0;
        const double sin_phi = sin(deg2Rad(coordSpherical.phig));

        for(int n = 1; n <= magneticModel.nMax; n++) {
            // Compute the ration between the Gauss-normalized associated Legendre functions
            // and the Schmidt quasi-normalized version. This is equivalent to
            // sqrt( (m == 0 ? 1 : 2) * ( n - m )!/( n + m!) ) * ( 2n - 1 )!!/(n - m)!

            const int index                = ((n * (n + 1) / 2) + 1);
            const double schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2. * n - 1) / static_cast<double>(n);
            const double schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((2. * n) / static_cast<double>(n + 1));
            schmidtQuasiNorm1              = schmidtQuasiNorm2;
            if(n == 1) {
                PcupS.at(n) = PcupS.at(n - 1);
            } else {
                const double k = static_cast<double>(((n - 1) * (n - 1)) - 1) / ((2. * n - 1) * (2. * n - 3));
                PcupS.at(n)    = sin_phi * PcupS.at(n - 1) - k * PcupS.at(n - 2);
            }

            //		  1 nMax  (n+2)    n     m            m           m
            // By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
            //         n=1             m=0   n            n           n
            // Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.
            By += sphVariables.relativeRadiusPower.at(n) *
                  (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.sin_mlambda.at(1) -
                   magneticModel.main_Field_Coeff_H.at(index) * sphVariables.cos_mlambda.at(1)) *
                  PcupS.at(n) * schmidtQuasiNorm3;
        }
    }

    /** @brief This Function sums the secular variation coefficients to get the secular variation of the Magnetic vector.
     */
    void MagneticResults::secVarSummation(const LegendreFunction &legendreFunction, MagneticModel &magneticModel,
                                          const SphericalHarmonicVariables &sphVariables,
                                          const CoordSpherical &coordSpherical) {
        magneticModel.secularVariationUsed = true;
        Bz                                 = 0.0;
        By                                 = 0.0;
        Bx                                 = 0.0;
        for(int n = 1; n <= magneticModel.nMaxSecVar; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = ((n * (n + 1) / 2) + m);

                //		    nMax  	(n+2) 	  n     m            m           m
                // Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                //          n=1      	      m=0   n            n           n
                //  Derivative with respect to radius.
                Bz -= sphVariables.relativeRadiusPower.at(n) *
                      (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      (n + 1) * legendreFunction.Pcup.at(index);

                //		  1 nMax  (n+2)    n     m            m           m
                //  By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                //         n=1             m=0   n            n           n
                // Derivative with respect to longitude, divided by radius.
                By += sphVariables.relativeRadiusPower.at(n) *
                      (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.sin_mlambda.at(m) -
                       magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.cos_mlambda.at(m)) *
                      m * legendreFunction.Pcup.at(index);

                //		   nMax  (n+2) n     m            m           m
                // Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                //         n=1         m=0   n            n           n
                // Derivative with respect to latitude, divided by radius.
                Bx -= sphVariables.relativeRadiusPower.at(n) *
                      (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      legendreFunction.dPcup.at(index);
            }
        }
        const double cos_phi = cos(deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            By = By / cos_phi;
        } else  // Special calculation for component By at Geographic poles
        {
            secVarSummationSpecial(magneticModel, sphVariables, coordSpherical);
        }
    }

    /** @brief Special calculation for the secular variation summation at the poles.
     */
    void MagneticResults::secVarSummationSpecial(const MagneticModel &magneticModel,
                                                 const SphericalHarmonicVariables &sphVariables,
                                                 const CoordSpherical &coordSpherical) {
        std::vector<double> PcupS(magneticModel.nMaxSecVar + 1);

        PcupS[0]                 = 1;
        double schmidtQuasiNorm1 = 1.0;

        By                   = 0.0;
        const double sin_phi = sin(deg2Rad(coordSpherical.phig));

        for(int n = 1; n <= magneticModel.nMaxSecVar; n++) {
            const int index                = ((n * (n + 1) / 2) + 1);
            const double schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2. * n - 1) / static_cast<double>(n);
            const double schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((2. * n) / static_cast<double>(n + 1));
            schmidtQuasiNorm1              = schmidtQuasiNorm2;
            if(n == 1) {
                PcupS.at(n) = PcupS.at(n - 1);
            } else {
                const double k = static_cast<double>(((n - 1) * (n - 1)) - 1) / ((2. * n - 1) * (2. * n - 3));
                PcupS.at(n)    = sin_phi * PcupS.at(n - 1) - k * PcupS.at(n - 2);
            }

            //		  1 nMax  (n+2)    n     m            m           m
            // By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
            //         n=1             m=0   n            n           n
            // Derivative with respect to longitude, divided by radius.
            By += sphVariables.relativeRadiusPower.at(n) *
                  (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.sin_mlambda.at(1) -
                   magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.cos_mlambda.at(1)) *
                  PcupS.at(n) * schmidtQuasiNorm3;
        }
    }

    /** @brief Rotate the Magnetic Vectors to Geodetic Coordinates
     * Manoj Nair, June, 2009 Manoj.C.Nair\@Noaa.Gov
     * Equation 16, WMM Technical report
     **/
    void MagneticResults::rotateMagneticVector(const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
                                               const MagneticResults &magneticResultsSph) {
        // Difference between the spherical and Geodetic latitudes
        const double psi = (std::numbers::pi / 180) * (coordSpherical.phig - coordGeodetic.phi);

        // Rotate spherical field components to the Geodetic system
        Bz = magneticResultsSph.Bx * sin(psi) + magneticResultsSph.Bz * cos(psi);
        Bx = magneticResultsSph.Bx * cos(psi) - magneticResultsSph.Bz * sin(psi);
        By = magneticResultsSph.By;
    }

    void MagneticResults::GradYSummation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                                         const SphericalHarmonicVariables &sphVariables,
                                         const CoordSpherical &coordSpherical) {
        Bz = 0.0;
        By = 0.0;
        Bx = 0.0;
        for(int n = 1; n <= magneticModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = ((n * (n + 1) / 2) + m);

                Bz -= sphVariables.relativeRadiusPower[n] *
                      (-1 * magneticModel.main_Field_Coeff_G[index] * sphVariables.sin_mlambda[m] +
                       magneticModel.main_Field_Coeff_H[index] * sphVariables.cos_mlambda[m]) *
                      (double)(n + 1) * (double)(m)*legendreFunction.Pcup[index] * (1 / coordSpherical.r);
                By += sphVariables.relativeRadiusPower[n] *
                      (magneticModel.main_Field_Coeff_G[index] * sphVariables.cos_mlambda[m] +
                       magneticModel.main_Field_Coeff_H[index] * sphVariables.sin_mlambda[m]) *
                      (double)(m * m) * legendreFunction.Pcup[index] * (1 / coordSpherical.r);
                Bx -= sphVariables.relativeRadiusPower[n] *
                      (-1 * magneticModel.main_Field_Coeff_G[index] * sphVariables.sin_mlambda[m] +
                       magneticModel.main_Field_Coeff_H[index] * sphVariables.cos_mlambda[m]) *
                      (double)(m)*legendreFunction.dPcup[index] * (1 / coordSpherical.r);
            }
        }

        const double cos_phi = std::cos(deg2Rad(coordSpherical.phig));
        if(std::fabs(cos_phi) > 1.0e-10) {
            By = By / (cos_phi * cos_phi);
            Bx = Bx / (cos_phi);
            Bz = Bz / (cos_phi);
        } else
        // Special calculation for component - By - at Geographic poles.
        // If the user wants to avoid using this function,  please make sure that
        // the latitude is not exactly +/-90. An option is to make use the function
        // CheckGeographicPoles.
        {
            // summationSpecial(magneticModel, sphVariables, coordSpherical, gradY);
        }
    }
}  // namespace wmm