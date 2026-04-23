//
// Created by Professional on 22.04.2026.
//

#include "GeoMagneticElements.h"

#include "CoordGeodetic.h"
#include "CoordSpherical.h"
#include "Ellipsoid.h"
#include "LegendreFunction.h"
#include "MagneticConstants.h"
#include "MagneticModel.h"
#include "MagneticResults.h"
#include "MagneticUtils.h"
#include "SphericalHarmonicVariables.h"
#include "UTMParameters.h"

#include <cmath>

namespace wmm {

    /** Calculate all the Geomagnetic elements from X,Y and Z components
     * INPUT     MagneticResultsGeo
     */
    void GeoMagneticElements::calculate(const MagneticResults &magneticResultsGeo) {
        X = magneticResultsGeo.Bx;
        Y = magneticResultsGeo.By;
        Z = magneticResultsGeo.Bz;

        H    = sqrt(magneticResultsGeo.Bx * magneticResultsGeo.Bx + magneticResultsGeo.By * magneticResultsGeo.By);
        F    = sqrt(H * H + magneticResultsGeo.Bz * magneticResultsGeo.Bz);
        Decl = Rad2Deg(atan2(Y, X));
        Incl = Rad2Deg(atan2(Z, H));
    }

    /**This takes the Magnetic Variation in x, y, and z and uses it to calculate the secular variation of each of the
     * Geomagnetic elements.
     * INPUT : MagneticVariation
     */
    void GeoMagneticElements::calculateSecularVariation(const MagneticResults &magneticVariation) {
        Xdot    = magneticVariation.Bx;
        Ydot    = magneticVariation.By;
        Zdot    = magneticVariation.Bz;
        Hdot    = (X * Xdot + Y * Ydot) / H; /* See equation 19 in the WMM technical report */
        Fdot    = (X * Xdot + Y * Ydot + Z * Zdot) / F;
        Decldot = 180.0 / M_PI * (X * Ydot - Y * Xdot) / (H * H);
        Incldot = 180.0 / M_PI * (H * Zdot - Z * Hdot) / (F * F);
        GVdot   = Decldot;
    }

    void GeoMagneticElements::calculateGradientElements(const MagneticResults &gradResults,
                                                        const GeoMagneticElements &magneticElements) {
        X = gradResults.Bx;
        Y = gradResults.By;
        Z = gradResults.Bz;

        H    = (X * magneticElements.X + Y * magneticElements.Y) / magneticElements.H;
        F    = (X * magneticElements.X + Y * magneticElements.Y + Z * magneticElements.Z) / magneticElements.F;
        Decl = 180.0 / M_PI * (magneticElements.X * Y - magneticElements.Y * X) / (magneticElements.H * magneticElements.H);
        Incl = 180.0 / M_PI * (magneticElements.H * Z - magneticElements.Z * H) / (magneticElements.F * magneticElements.F);
        GV   = Decl;
    }

    /**
     * The main subroutine that calls a sequence of WMM sub-functions to calculate the magnetic field elements for a single
     * point. The function expects the model coefficients and point coordinates as input and returns the magnetic field
     * elements and their rate of change. Though, this subroutine can be called successively to calculate a time series,
     * profile or grid of magnetic field, these are better achieved by the subroutine Grid.
     *
     * INPUT: ellip
     *        CoordSpherical
     *        CoordGeodetic
     *        TimedMagneticModel
     *
     * OUTPUT : GeoMagneticElements
     */
    void GeoMagneticElements::geomag(const Ellipsoid &ellip, const CoordSpherical &coordSpherical,
                                     const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel) {
        MagneticResults MagneticResultsSph{};
        MagneticResults MagneticResultsGeo{};
        MagneticResults MagneticResultsSphVar{};
        MagneticResults MagneticResultsGeoVar{};

        // Create and compute ALF functions
        LegendreFunction legendreFunction(coordSpherical, timedMagneticModel.nMax);

        // Create and compute Spherical Harmonic variables
        SphericalHarmonicVariables sphVariables(ellip, coordSpherical, timedMagneticModel.nMax);

        // Accumulate the spherical harmonic coefficients
        MagneticResultsSph.Summation(legendreFunction, timedMagneticModel, sphVariables, coordSpherical);

        // Sum the Secular Variation Coefficients
        MagneticResultsSphVar.SecVarSummation(legendreFunction, timedMagneticModel, sphVariables, coordSpherical);

        // Map the computed Magnetic fields to Geodeitic coordinates
        MagneticResultsGeo.RotateMagneticVector(coordSpherical, coordGeodetic, MagneticResultsSph);

        // Map the secular variation field components to Geodetic coordinates
        MagneticResultsGeoVar.RotateMagneticVector(coordSpherical, coordGeodetic, MagneticResultsSphVar);

        // Calculate the Geomagnetic elements, Equation 19, WMM Technical report
        this->calculate(MagneticResultsGeo);

        // Calculate the secular variation of each of the Geomagnetic elements
        this->calculateSecularVariation(MagneticResultsGeoVar);
    }

    void GeoMagneticElements::gradY(const Ellipsoid &ellip, const CoordSpherical &coordSpherical,
                                    const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel,
                                    const GeoMagneticElements &geoMagneticElements) {
        // Create and compute ALF functions
        LegendreFunction legendreFunction(coordSpherical, timedMagneticModel.nMax);

        // Create and compute Spherical Harmonic variables
        SphericalHarmonicVariables sphVariables(ellip, coordSpherical, timedMagneticModel.nMax);

        // Accumulate the spherical harmonic coefficients
        MagneticResults GradYResultsSph;
        GradYResultsSph.GradYSummation(legendreFunction, timedMagneticModel, sphVariables, coordSpherical);

        // Map the computed Magnetic fields to Geodetic coordinates
        MagneticResults GradYResultsGeo;
        GradYResultsGeo.RotateMagneticVector(coordSpherical, coordGeodetic, GradYResultsSph);

        // Calculate the Geomagnetic elements, Equation 18 , WMM Technical report
        this->calculateGradientElements(GradYResultsGeo, geoMagneticElements);
    }

    /** Computes the grid variation for |latitudes| > MAX_LAT_DEGREE
     * Grivation (or grid variation) is the angle between grid north and magnetic north. This routine calculates Grivation
     * for the Polar Stereographic projection for polar locations (Latitude => |55| deg). Otherwise, it computes the grid
     * variation in UTM projection system. However, the UTM projection codes may be used to compute the grid variation at any
     * latitudes.
     *
     * INPUT :  CoordGeodetic location
     * OUTPUT:  GeoMagneticElements elements
     **/
    int GeoMagneticElements::CalculateGridVariation(const CoordGeodetic &location) {
        UTMParameters UTMParameters;

        if(location.phi >= PS_MAX_LAT_DEGREE) {
            GV = Decl - location.lambda;
            return 1;
        }

        if(location.phi <= PS_MIN_LAT_DEGREE) {
            GV = Decl + location.lambda;
            return 1;
        }

        UTMParameters.GetTransverseMercator(location);
        GV = Decl - UTMParameters.ConvergenceOfMeridians;
        return 0;
    }

    void GeoMagneticElements::WMMErrorCalc(double H_) {
        double decl_variable, decl_constant;
        F             = WMM_UNCERTAINTY_F;
        H             = WMM_UNCERTAINTY_H;
        X             = WMM_UNCERTAINTY_X;
        Z             = WMM_UNCERTAINTY_Z;
        Incl          = WMM_UNCERTAINTY_I;
        Y             = WMM_UNCERTAINTY_Y;
        decl_variable = (WMM_UNCERTAINTY_D_COEF / H_);
        decl_constant = (WMM_UNCERTAINTY_D_OFFSET);
        Decl          = sqrt(decl_constant * decl_constant + decl_variable * decl_variable);
        if(Decl > 180) {
            Decl = 180;
        }
    }

    void GeoMagneticElements::WMMHRErrorCalc(double H_) {
        double decl_variable, decl_constant;
        F             = WMMHR_UNCERTAINTY_F;
        H             = WMMHR_UNCERTAINTY_H;
        X             = WMMHR_UNCERTAINTY_X;
        Z             = WMMHR_UNCERTAINTY_Z;
        Incl          = WMMHR_UNCERTAINTY_I;
        Y             = WMMHR_UNCERTAINTY_Y;
        decl_variable = (WMMHR_UNCERTAINTY_D_COEF / H_);
        decl_constant = (WMMHR_UNCERTAINTY_D_OFFSET);
        Decl          = sqrt(decl_constant * decl_constant + decl_variable * decl_variable);
        if(Decl > 180) {
            Decl = 180;
        }
    }

    // This function scales all the geomagnetic elements to scale a vector use MagneticResultsScale
    GeoMagneticElements GeoMagneticElements::scaled(const double factor) const {
        GeoMagneticElements product;
        product.X       = X * factor;
        product.Y       = Y * factor;
        product.Z       = Z * factor;
        product.H       = H * factor;
        product.F       = F * factor;
        product.Incl    = Incl * factor;
        product.Decl    = Decl * factor;
        product.GV      = GV * factor;
        product.Xdot    = Xdot * factor;
        product.Ydot    = Ydot * factor;
        product.Zdot    = Zdot * factor;
        product.Hdot    = Hdot * factor;
        product.Fdot    = Fdot * factor;
        product.Incldot = Incldot * factor;
        product.Decldot = Decldot * factor;
        product.GVdot   = GVdot * factor;
        return product;
    }

    // This function scales all the geomagnetic elements to scale a vector use MagneticResultsScale
    void GeoMagneticElements::scale(const double factor) {
        X *= factor;
        Y *= factor;
        Z *= factor;
        H *= factor;
        F *= factor;
        Incl *= factor;
        Decl *= factor;
        GV *= factor;
        Xdot *= factor;
        Ydot *= factor;
        Zdot *= factor;
        Hdot *= factor;
        Fdot *= factor;
        Incldot *= factor;
        Decldot *= factor;
        GVdot *= factor;
    }

    /** This algorithm does not result in the difference of F being derived from
     * the Pythagorean theorem.  This function should be used for computing residuals
     * or changes in elements.*/
    GeoMagneticElements GeoMagneticElements::operator-(const GeoMagneticElements &subtrahend) const {
        GeoMagneticElements difference;
        difference.X = X - subtrahend.X;
        difference.Y = Y - subtrahend.Y;
        difference.Z = Z - subtrahend.Z;

        difference.H    = H - subtrahend.H;
        difference.F    = F - subtrahend.F;
        difference.Decl = Decl - subtrahend.Decl;
        difference.Incl = Incl - subtrahend.Incl;

        difference.Xdot = Xdot - subtrahend.Xdot;
        difference.Ydot = Ydot - subtrahend.Ydot;
        difference.Zdot = Zdot - subtrahend.Zdot;

        difference.Hdot    = Hdot - subtrahend.Hdot;
        difference.Fdot    = Fdot - subtrahend.Fdot;
        difference.Decldot = Decldot - subtrahend.Decldot;
        difference.Incldot = Incldot - subtrahend.Incldot;

        difference.GV    = GV - subtrahend.GV;
        difference.GVdot = GVdot - subtrahend.GVdot;

        return difference;
    }

    // Errors.Decl, Errors.Incl, Errors.F are all assumed to exist
    void GeoMagneticElements::ErrorCalc(GeoMagneticElements B) {
        const double cos2D = cos(Deg2Rad(B.Decl)) * cos(Deg2Rad(B.Decl));
        const double cos2I = cos(Deg2Rad(B.Incl)) * cos(Deg2Rad(B.Incl));
        const double sin2D = sin(Deg2Rad(B.Decl)) * sin(Deg2Rad(B.Decl));
        const double sin2I = sin(Deg2Rad(B.Incl)) * sin(Deg2Rad(B.Incl));
        const double eD    = Deg2Rad(Decl);
        const double eI    = Deg2Rad(Incl);
        const double EDSq  = eD * eD;
        const double EISq  = eI * eI;
        X = sqrt(cos2D * cos2I * F * F + B.F * B.F * sin2D * cos2I * EDSq + B.F * B.F * cos2D * sin2I * EISq);
        Y = sqrt(sin2D * cos2I * F * F + B.F * B.F * cos2D * cos2I * EDSq + B.F * B.F * sin2D * sin2I * EISq);
        Z = sqrt(sin2I * F * F + B.F * B.F * cos2I * EISq);
        H = sqrt(cos2I * F * F + B.F * B.F * sin2I * EISq);
    }

}  // namespace wmm