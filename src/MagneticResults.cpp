//
// Created by Professional on 22.04.2026.
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
     * INPUT :  LegendreFunction
     *          MagneticModel
     *          SphVariables
     *          CoordSpherical
     * OUTPUT : MagneticResults
     *Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
     */
    void MagneticResults::Summation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                                    const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical) {
        Bz = 0.0;
        By = 0.0;
        Bx = 0.0;

        for(int n = 1; n <= magneticModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = (n * (n + 1) / 2 + m);

                /*		    nMax  	(n+2) 	  n     m            m           m
                        Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                        n=1      	      m=0   n            n           n  */
                /* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
                Bz -= sphVariables.RelativeRadiusPower.at(n) *
                      (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.main_Field_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      static_cast<double>(n + 1) * legendreFunction.Pcup.at(index);

                /*		  1 nMax  (n+2)    n     m            m           m
                        By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1             m=0   n            n           n  */
                /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
                By += sphVariables.RelativeRadiusPower.at(n) *
                      (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.sin_mlambda.at(m) -
                       magneticModel.main_Field_Coeff_H.at(index) * sphVariables.cos_mlambda.at(m)) *
                      static_cast<double>(m) * legendreFunction.Pcup.at(index);
                /*		   nMax  (n+2) n     m            m           m
                        Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1         m=0   n            n           n  */
                /* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */

                Bx -= sphVariables.RelativeRadiusPower.at(n) *
                      (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.main_Field_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      legendreFunction.dPcup.at(index);
            }
        }

        const double cos_phi = cos(Deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            By = By / cos_phi;
        } else
        /* Special calculation for component - By - at Geographic poles.
         * If the user wants to avoid using this function,  please make sure that
         * the latitude is not exactly +/-90. An option is to make use the function
         * CheckGeographicPoles.
         */
        {
            SummationSpecial(magneticModel, sphVariables, coordSpherical);
        }
    }

    /** Special calculation for the component By at Geographic poles.
     * Manoj Nair, June, 2009 manoj.c.nair@noaa.gov
     * INPUT:  MagneticModel
     *         SphVariables
     *         CoordSpherical
     * OUTPUT: MagneticResults
     * See Section 1.4, "SINGULARITIES AT THE GEOGRAPHIC POLES", WMM Technical report
     */
    void MagneticResults::SummationSpecial(const MagneticModel &magneticModel,
                                           const SphericalHarmonicVariables &sphVariables,
                                           const CoordSpherical &coordSpherical) {
        std::vector<double> PcupS(magneticModel.nMax + 1);

        PcupS.at(0)              = 1;
        double schmidtQuasiNorm1 = 1.0;

        By                   = 0.0;
        const double sin_phi = sin(Deg2Rad(coordSpherical.phig));

        for(int n = 1; n <= magneticModel.nMax; n++) {
            /* Compute the ration between the Gauss-normalized associated Legendre functions and the Schmidt quasi-normalized
      version. This is equivalent to sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)! */

            const int index                = (n * (n + 1) / 2 + 1);
            const double schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2. * n - 1) / static_cast<double>(n);
            const double schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((2. * n) / static_cast<double>(n + 1));
            schmidtQuasiNorm1              = schmidtQuasiNorm2;
            if(n == 1) {
                PcupS.at(n) = PcupS.at(n - 1);
            } else {
                const double k = static_cast<double>(((n - 1) * (n - 1)) - 1) / ((2. * n - 1) * (2. * n - 3));
                PcupS.at(n)    = sin_phi * PcupS.at(n - 1) - k * PcupS.at(n - 2);
            }

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
            /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */

            By += sphVariables.RelativeRadiusPower.at(n) *
                  (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.sin_mlambda.at(1) -
                   magneticModel.main_Field_Coeff_H.at(index) * sphVariables.cos_mlambda.at(1)) *
                  PcupS.at(n) * schmidtQuasiNorm3;
        }
    }

    /** This Function sums the secular variation coefficients to get the secular variation of the Magnetic vector.
     * INPUT :  LegendreFunction
     *          MagneticModel
     *          SphVariables
     *          CoordSpherical
     * OUTPUT : MagneticResults
     */
    void MagneticResults::SecVarSummation(const LegendreFunction &legendreFunction, MagneticModel &magneticModel,
                                          const SphericalHarmonicVariables &sphVariables,
                                          const CoordSpherical &coordSpherical) {
        magneticModel.secularVariationUsed = true;
        Bz                                 = 0.0;
        By                                 = 0.0;
        Bx                                 = 0.0;
        for(int n = 1; n <= magneticModel.nMaxSecVar; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = (n * (n + 1) / 2 + m);

                /*		    nMax  	(n+2) 	  n     m            m           m
                        Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                        n=1      	      m=0   n            n           n  */
                /*  Derivative with respect to radius.*/
                Bz -= sphVariables.RelativeRadiusPower.at(n) *
                      (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      static_cast<double>(n + 1) * legendreFunction.Pcup.at(index);

                /*		  1 nMax  (n+2)    n     m            m           m
                        By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1             m=0   n            n           n  */
                /* Derivative with respect to longitude, divided by radius. */
                By += sphVariables.RelativeRadiusPower.at(n) *
                      (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.sin_mlambda.at(m) -
                       magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.cos_mlambda.at(m)) *
                      static_cast<double>(m) * legendreFunction.Pcup.at(index);

                /*		   nMax  (n+2) n     m            m           m
                        Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1         m=0   n            n           n  */
                /* Derivative with respect to latitude, divided by radius. */
                Bx -= sphVariables.RelativeRadiusPower.at(n) *
                      (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                       magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                      legendreFunction.dPcup.at(index);
            }
        }
        double cos_phi = cos(Deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            By = By / cos_phi;
        } else
        /* Special calculation for component By at Geographic poles */
        {
            SecVarSummationSpecial(magneticModel, sphVariables, coordSpherical);
        }
    }

    /** Special calculation for the secular variation summation at the poles.
     * INPUT:  MagneticModel
     *         SphVariables
     *         CoordSpherical
     * OUTPUT: MagneticResults
     */
    void MagneticResults::SecVarSummationSpecial(const MagneticModel &magneticModel,
                                                 const SphericalHarmonicVariables &sphVariables,
                                                 const CoordSpherical &coordSpherical) {
        std::vector<double> PcupS(magneticModel.nMaxSecVar + 1);

        PcupS[0]                 = 1;
        double schmidtQuasiNorm1 = 1.0;

        By                   = 0.0;
        const double sin_phi = sin(Deg2Rad(coordSpherical.phig));

        for(int n = 1; n <= magneticModel.nMaxSecVar; n++) {
            const int index                = (n * (n + 1) / 2 + 1);
            const double schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2. * n - 1) / static_cast<double>(n);
            const double schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((2. * n) / static_cast<double>(n + 1));
            schmidtQuasiNorm1              = schmidtQuasiNorm2;
            if(n == 1) {
                PcupS.at(n) = PcupS.at(n - 1);
            } else {
                const double k = static_cast<double>(((n - 1) * (n - 1)) - 1) / ((2. * n - 1) * (2. * n - 3));
                PcupS.at(n)    = sin_phi * PcupS.at(n - 1) - k * PcupS.at(n - 2);
            }

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
            /* Derivative with respect to longitude, divided by radius. */
            By += sphVariables.RelativeRadiusPower.at(n) *
                  (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.sin_mlambda.at(1) -
                   magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.cos_mlambda.at(1)) *
                  PcupS.at(n) * schmidtQuasiNorm3;
        }
    }

    /** Rotate the Magnetic Vectors to Geodetic Coordinates
     * Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
     * Equation 16, WMM Technical report
     * INPUT : CoordSpherical
     *         CoordGeodetic
     *         MagneticResultsSph
     **/
    void MagneticResults::RotateMagneticVector(const CoordSpherical &CoordSpherical, const CoordGeodetic &coordGeodetic,
                                               const MagneticResults &magneticResultsSph) {
        // Difference between the spherical and Geodetic latitudes
        const double Psi = (M_PI / 180) * (CoordSpherical.phig - coordGeodetic.phi);

        // Rotate spherical field components to the Geodetic system
        Bz = magneticResultsSph.Bx * sin(Psi) + magneticResultsSph.Bz * cos(Psi);
        Bx = magneticResultsSph.Bx * cos(Psi) - magneticResultsSph.Bz * sin(Psi);
        By = magneticResultsSph.By;
    }

    void MagneticResults::GradYSummation(const LegendreFunction &LegendreFunction, const MagneticModel &magneticModel,
                                         const SphericalHarmonicVariables &SphVariables,
                                         const CoordSpherical &coordSpherical) {
        int m, n, index;
        double cos_phi;
        Bz = 0.0;
        By = 0.0;
        Bx = 0.0;
        for(n = 1; n <= magneticModel.nMax; n++) {
            for(m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);

                Bz -= SphVariables.RelativeRadiusPower[n] *
                      (-1 * magneticModel.main_Field_Coeff_G[index] * SphVariables.sin_mlambda[m] +
                       magneticModel.main_Field_Coeff_H[index] * SphVariables.cos_mlambda[m]) *
                      (double)(n + 1) * (double)(m)*LegendreFunction.Pcup[index] * (1 / coordSpherical.r);
                By += SphVariables.RelativeRadiusPower[n] *
                      (magneticModel.main_Field_Coeff_G[index] * SphVariables.cos_mlambda[m] +
                       magneticModel.main_Field_Coeff_H[index] * SphVariables.sin_mlambda[m]) *
                      (double)(m * m) * LegendreFunction.Pcup[index] * (1 / coordSpherical.r);
                Bx -= SphVariables.RelativeRadiusPower[n] *
                      (-1 * magneticModel.main_Field_Coeff_G[index] * SphVariables.sin_mlambda[m] +
                       magneticModel.main_Field_Coeff_H[index] * SphVariables.cos_mlambda[m]) *
                      (double)(m)*LegendreFunction.dPcup[index] * (1 / coordSpherical.r);
            }
        }

        cos_phi = cos(Deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            By = By / (cos_phi * cos_phi);
            Bx = Bx / (cos_phi);
            Bz = Bz / (cos_phi);
        } else
        /* Special calculation for component - By - at Geographic poles.
         * If the user wants to avoid using this function,  please make sure that
         * the latitude is not exactly +/-90. An option is to make use the function
         * CheckGeographicPoles.
         */
        {
            /* SummationSpecial(MagneticModel, SphVariables, coordSpherical, gradY); */
        }
    }

}  // namespace wmm