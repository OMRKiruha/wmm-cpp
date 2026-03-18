
#include "GeomagnetismHeader.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>

/* $Id: GeomagnetismLibrary.c 1521 2017-01-24 17:52:41Z awoods $
 *
 * ABSTRACT
 *
 * The purpose of Geomagnetism Library is primarily to support the World Magnetic Model (WMM) 2015-2020.
 * It however is built to be used for spherical harmonic models of the Earth's magnetic field
 * generally and supports models even with a large (>>12) number of degrees.  It is also used in many
 * other geomagnetic models distributed by NCEI.
 *
 * REUSE NOTES
 *
 * Geomagnetism Library is intended for reuse by any application that requires
 * Computation of Geomagnetic field from a spherical harmonic model.
 *
 * REFERENCES
 *
 *    Further information on Geoid can be found in the WMM Technical Documents.
 *
 *
 * LICENSES
 *
 *  The WMM source code is in the public domain and not licensed or under copyright.
 *	The information and software may be used freely by the public. As required by 17 U.S.C. 403,
 *	third parties producing copyrighted works consisting predominantly of the material produced by
 *	U.S. government agencies must provide notice with such work(s) identifying the U.S. Government material
 *	incorporated and stating that such material is not subject to copyright protection.
 *
 * RESTRICTIONS
 *
 *    Geomagnetism library has no restrictions.
 *
 * ENVIRONMENT
 *
 *    Geomagnetism library was tested in the following environments
 *
 *    1. Red Hat Linux  with GCC Compiler
 *    2. MS Windows 7 with MinGW compiler
 *    3. Sun Solaris with GCC Compiler
 *
 *


 *  National Centers for Environmental Information
 *  NOAA E/NE42, 325 Broadway
 *  Boulder, CO 80305 USA
 *  Attn: Arnaud Chulliat
 *  Phone:  (303) 497-6522
 *  Email:  Arnaud.Chulliat@noaa.gov

 *  Software and Model Support
 *  National Centers for Environmental Information
 *  NOAA E/NE42
 *  325 Broadway
 *  Boulder, CO 80305 USA
 *  Attn: Adam Woods or Manoj Nair
 *  Phone:  (303) 497-6640 or -4642
 *  Email:  geomag.models@noaa.gov
 *  URL: http://www.ngdc.noaa.gov/Geomagnetic/WMM/DoDWMM.shtml


 *  For more details on the subroutines, please consult the WMM
 *  Technical Documentations at
 *  http://www.ngdc.noaa.gov/Geomagnetic/WMM/DoDWMM.shtml

 *  Nov 23, 2009
 *  Written by Manoj C Nair and Adam Woods
 *  Manoj.C.Nair@noaa.Gov
 *  Adam.Woods@noaa.gov
 */

namespace wmm {

    MagneticModel::MagneticModel(const int numTerms) {
        main_Field_Coeff_G.resize(numTerms + 1);
        main_Field_Coeff_H.resize(numTerms + 1);
        secular_Var_Coeff_G.resize(numTerms + 1);
        secular_Var_Coeff_H.resize(numTerms + 1);
    }

    Date::Date() {
        const std::time_t t = std::time(nullptr);
        const std::tm *now  = std::localtime(&t);
        Year                = now->tm_year + 1900;
        Month               = now->tm_mon + 1;
        Day                 = now->tm_mday;

        DecimalYear = date_to_decYear(Year, Month, Day);
    }

    LegendreFunction::LegendreFunction(const int numTerms) {
        Pcup.resize(numTerms + 1);
        dPcup.resize(numTerms + 1);
    }

    SphericalHarmonicVariables::SphericalHarmonicVariables(const int numTerms) {
        RelativeRadiusPower.resize(numTerms + 1);
        cos_mlambda.resize(numTerms + 1);
        sin_mlambda.resize(numTerms + 1);
    }

    /******************************************************************************
     ************************************Wrapper***********************************
     * This grouping consists of functions call groups of other functions to do a
     * complete calculation of some sort.  For example, the Geomag function
     * does everything necessary to compute the geomagnetic elements from a given
     * geodetic point in space and magnetic model adjusted for the appropriate
     * date. These functions are the external functions necessary to create a
     * program that uses or calculates the magnetic field.
     ******************************************************************************
     ******************************************************************************/

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
     *
     * CALLS: ComputeSphericalHarmonicVariables(ellip, CoordSpherical, timedMagneticModel->nMax, &SphVariables); (Compute
     * Spherical Harmonic variables  )
     *       AssociatedLegendreFunction(CoordSpherical, timedMagneticModel->nMax, LegendreFunction); Compute ALF
     *       Summation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical, &MagneticResultsSph);
     * Accumulate the spherical harmonic coefficients
     *       SecVarSummation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical, &MagneticResultsSphVar);
     * Sum the Secular Variation Coefficients
     *       RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSph, &MagneticResultsGeo); Map the
     * computed Magnetic fields to Geodetic coordinates
     *       CalculateGeoMagneticElements(&MagneticResultsGeo, GeoMagneticElements); Calculate the Geomagnetic elements
     *       CalculateSecularVariationElements(MagneticResultsGeoVar, GeoMagneticElements); Calculate the secular
     * variation of each of the Geomagnetic elements
     */
    int Geomag(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
               MagneticModel *timedMagneticModel, GeoMagneticElements *geoMagneticElements) {
        MagneticResults MagneticResultsSph{};
        MagneticResults MagneticResultsGeo{};
        MagneticResults MagneticResultsSphVar{};
        MagneticResults MagneticResultsGeoVar{};

        const int NumTerms = ((timedMagneticModel->nMax + 1) * (timedMagneticModel->nMax + 2) / 2);

        /* For storing the ALF functions */
        LegendreFunction legendreFunction(NumTerms);
        SphericalHarmonicVariables sphVariables(timedMagneticModel->nMax);

        /* Compute Spherical Harmonic variables  */
        ComputeSphericalHarmonicVariables(ellip, coordSpherical, timedMagneticModel->nMax, &sphVariables);

        /* Compute ALF  */
        AssociatedLegendreFunction(coordSpherical, timedMagneticModel->nMax, &legendreFunction);

        /* Accumulate the spherical harmonic coefficients*/
        Summation(legendreFunction, *timedMagneticModel, sphVariables, coordSpherical, &MagneticResultsSph);

        /*Sum the Secular Variation Coefficients  */
        SecVarSummation(legendreFunction, timedMagneticModel, sphVariables, coordSpherical, &MagneticResultsSphVar);

        /* Map the computed Magnetic fields to Geodeitic coordinates  */
        RotateMagneticVector(coordSpherical, coordGeodetic, MagneticResultsSph, &MagneticResultsGeo);

        /* Map the secular variation field components to Geodetic coordinates*/
        RotateMagneticVector(coordSpherical, coordGeodetic, MagneticResultsSphVar, &MagneticResultsGeoVar);

        /* Calculate the Geomagnetic elements, Equation 19 , WMM Technical report */
        CalculateGeoMagneticElements(MagneticResultsGeo, geoMagneticElements);

        /*Calculate the secular variation of each of the Geomagnetic elements*/
        CalculateSecularVariationElements(MagneticResultsGeoVar, geoMagneticElements);

        return true;
    }

    /*It should be noted that the x[2], y[2], and z[2] variables are NOT the same
     coordinate system as the directions in which the gradients are taken.  These
     variables represent a Cartesian coordinate system where the Earth's center is
     the origin, 'z' points up toward the North (rotational) pole and 'x' points toward
     the prime meridian.  'y' points toward longitude = 90 degrees East.
     The gradient is preformed along a local Cartesian coordinate system with the
     origin at CoordGeodetic.  'z' points down toward the Earth's core, x points
     North, tangent to the local longitude line, and 'y' points East, tangent to
     the local latitude line.*/
    void CalcGradient(Ellipsoid ellip, CoordGeodetic coordGeodetic, MagneticModel *timedMagneticModel, Gradient *gradient) {
        double phiDelta = 0.01, /* DeltaY = 0.01, */ hDelta = -1, x[2], y[2], z[2], distance;

        CoordSpherical AdjCoordSpherical;
        CoordGeodetic AdjCoordGeodetic;
        GeoMagneticElements GeomagneticElements, AdjGeoMagneticElements[2];


        /* Initialization */
        GeodeticToSpherical(ellip, coordGeodetic, &AdjCoordSpherical);
        Geomag(ellip, AdjCoordSpherical, coordGeodetic, timedMagneticModel, &GeomagneticElements);
        AdjCoordGeodetic = CoordGeodeticAssign(coordGeodetic);


        /* Gradient along x */

        AdjCoordGeodetic.phi = coordGeodetic.phi + phiDelta;
        GeodeticToSpherical(ellip, AdjCoordGeodetic, &AdjCoordSpherical);
        Geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel, &AdjGeoMagneticElements[0]);
        SphericalToCartesian(AdjCoordSpherical, &x[0], &y[0], &z[0]);
        AdjCoordGeodetic.phi = coordGeodetic.phi - phiDelta;
        GeodeticToSpherical(ellip, AdjCoordGeodetic, &AdjCoordSpherical);
        Geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel, &AdjGeoMagneticElements[1]);
        SphericalToCartesian(AdjCoordSpherical, &x[1], &y[1], &z[1]);


        distance = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) + (z[0] - z[1]) * (z[0] - z[1]));
        gradient->GradPhi = GeoMagneticElementsSubtract(AdjGeoMagneticElements[0], AdjGeoMagneticElements[1]);
        gradient->GradPhi = GeoMagneticElementsScale(gradient->GradPhi, 1 / distance);
        AdjCoordGeodetic  = CoordGeodeticAssign(coordGeodetic);

        /*Gradient along y*/

        /*It is perhaps noticeable that the method here for calculation is substantially
         different than that for the gradient along x.  As we near the North pole
         the longitude lines approach each other, and the calculation that works well
         for latitude lines becomes unstable when 0.01 degrees represents sufficiently
         small numbers, and fails to function correctly at all at the North Pole */

        GeodeticToSpherical(ellip, coordGeodetic, &AdjCoordSpherical);
        GradY(ellip, AdjCoordSpherical, coordGeodetic, timedMagneticModel, GeomagneticElements, &(gradient->GradLambda));

        /*Gradient along z*/
        AdjCoordGeodetic.HeightAboveEllipsoid = coordGeodetic.HeightAboveEllipsoid + hDelta;
        AdjCoordGeodetic.HeightAboveGeoid     = coordGeodetic.HeightAboveGeoid + hDelta;
        GeodeticToSpherical(ellip, AdjCoordGeodetic, &AdjCoordSpherical);
        Geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel, &AdjGeoMagneticElements[0]);
        SphericalToCartesian(AdjCoordSpherical, &x[0], &y[0], &z[0]);
        AdjCoordGeodetic.HeightAboveEllipsoid = coordGeodetic.HeightAboveEllipsoid - hDelta;
        AdjCoordGeodetic.HeightAboveGeoid     = coordGeodetic.HeightAboveGeoid - hDelta;
        GeodeticToSpherical(ellip, AdjCoordGeodetic, &AdjCoordSpherical);
        Geomag(ellip, AdjCoordSpherical, AdjCoordGeodetic, timedMagneticModel, &AdjGeoMagneticElements[1]);
        SphericalToCartesian(AdjCoordSpherical, &x[1], &y[1], &z[1]);

        distance = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) + (z[0] - z[1]) * (z[0] - z[1]));
        gradient->GradZ  = GeoMagneticElementsSubtract(AdjGeoMagneticElements[0], AdjGeoMagneticElements[1]);
        gradient->GradZ  = GeoMagneticElementsScale(gradient->GradZ, 1 / distance);
        AdjCoordGeodetic = CoordGeodeticAssign(coordGeodetic);
    }

    //    int robustReadMagneticModel_Large(const std::string_view filename, char *filenameSV, std::vector<MagneticModel>&
    //    magneticModels) {
    //        char line[MAXLINELENGTH], ModelName[] = "Enhanced Magnetic Model"; /*Model Name must be no longer than 31
    //        characters*/ int n, nMax = 0, nMaxSV = 0, num_terms, a, epochlength = 5, i;
    //
    //        std::ifstream file{filename.data()};
    //        if(!file.is_open()) {
    //            return 0;
    //        }
    //
    //        if (NULL == fgets(line, MAXLINELENGTH, MODELFILE)) {
    //            return 0;
    //        }
    //        do {
    //            if (NULL == fgets(line, MAXLINELENGTH, MODELFILE))
    //                break;
    //            a = sscanf(line, "%d", &n);
    //            if (n > nMax && (n < 99999 && a == 1 && n > 0))
    //                nMax = n;
    //        } while (n < 99999 && a == 1);
    //        fclose(MODELFILE);
    //        MODELFILE = fopen(filenameSV, "r");
    //        if (MODELFILE == 0) {
    //            return 0;
    //        }
    //        n = 0;
    //        if (NULL == fgets(line, MAXLINELENGTH, MODELFILE))
    //            return 0;
    //        do {
    //            if (NULL == fgets(line, MAXLINELENGTH, MODELFILE))
    //                break;
    //            a = sscanf(line, "%d", &n);
    //            if (n > nMaxSV && (n < 99999 && a == 1 && n > 0))
    //                nMaxSV = n;
    //        } while (n < 99999 && a == 1);
    //        fclose(MODELFILE);
    //        num_terms = CALCULATE_NUMTERMS(nMax);
    //        *magneticModel = AllocateModelMemory(num_terms);
    //        (*magneticModel)->nMax = nMax;
    //        (*magneticModel)->nMaxSecVar = nMaxSV;
    //        if (nMaxSV > 0) (*magneticModel)->SecularVariationUsed = true;
    //        for (i = 0; i < num_terms; i++) {
    //            (*magneticModel)->Main_Field_Coeff_G[i] = 0;
    //            (*magneticModel)->Main_Field_Coeff_H[i] = 0;
    //            (*magneticModel)->Secular_Var_Coeff_G[i] = 0;
    //            (*magneticModel)->Secular_Var_Coeff_H[i] = 0;
    //        }
    //        readMagneticModel_Large(filename, filenameSV, *magneticModel);
    //        (*magneticModel)->CoefficientFileEndDate = (*magneticModel)->epoch + epochlength;
    //        strlcpy_equivalent((*magneticModel)->ModelName, ModelName, sizeof((*magneticModel)->ModelName));
    //        magneticModel.EditionDate = (*magneticModel)->epoch;
    //        return 1;
    //    } /*robustReadMagneticModel_Large*/

    int robustReadMagModels(const std::string_view filename, MagneticModel &magneticModel, const int array_size) {
        std::string str;
        int n{};

        std::ifstream file{filename.data()};
        if(!file.is_open()) {
            return 0;
        }

        // First string in file is not empty
        if(std::getline(file, str); str.empty()) {
            return 0;
        }

        if(array_size == 1) {
            int nMax = 0;
            do {
                if(std::getline(file, str); str.empty()) {
                    break;
                }
                n = std::stoi(str.substr(0, 3));  // Extract first number in line
                if(n > nMax && (n < 999 && n > 0)) {
                    nMax = n;
                }
            } while(n < 999);

            const int num_terms      = (nMax * (nMax + 1) / 2 + nMax);
            magneticModel            = MagneticModel(num_terms);
            magneticModel.nMax       = nMax;
            magneticModel.nMaxSecVar = nMax;
            ReadMagneticModel(filename, &magneticModel);
            magneticModel.coefficientFileEndDate = magneticModel.epoch + 5;

        } else {
            return 0;
        }
        return 1;
    }

    /*End of Wrapper Functions*/

    /******************************************************************************
     ********************************User Interface********************************
     * This grouping consists of functions which interact with the directly with
     * the user and are generally specific to the XXX_point.c, XXX_grid.c, and
     * XXX_file.c programs. They deal with input from and output to the user.
     ******************************************************************************/

    /** This prints WMM errors.
     * INPUT : control Error look up number
     * OUTPUT: none
     * CALLS : none */
    void PrintError(const int control) {
        switch(control) {
            case 1:
                printf("\nError allocating in LegendreFunctionMemory.\n");
                break;
            case 2:
                printf("\nError allocating in AllocateModelMemory.\n");
                break;
            case 3:
                printf("\nError allocating in InitializeGeoid\n");
                break;
            case 4:
                printf("\nError in setting default values.\n");
                break;
            case 5:
                printf("\nError initializing Geoid.\n");
                break;
            case 6:
                printf("\nError opening wmmhr.cof\n.");
                break;
            case 7:
                printf("\nError opening WMMSV.COF\n.");
                break;
            case 8:
                printf("\nError reading Magnetic Model.\n");
                break;
            case 9:
                printf("\nError printing Command Prompt introduction.\n");
                break;
            case 10:
                printf("\nError converting from geodetic co-ordinates to spherical co-ordinates.\n");
                break;
            case 11:
                printf("\nError in time modifying the Magnetic model\n");
                break;
            case 12:
                printf("\nError in Geomagnetic\n");
                break;
            case 13:
                printf("\nError printing user data\n");
                break;
            case 14:
                printf("\nError allocating in SummationSpecial\n");
                break;
            case 15:
                printf("\nError allocating in SecVarSummationSpecial\n");
                break;
            case 16:
                printf("\nError in opening EGM9615.BIN file\n");
                break;
            case 17:
                printf("\nError: Latitude OR Longitude out of range in GetGeoidHeight\n");
                break;
            case 18:
                printf("\nError allocating in PcupHigh\n");
                break;
            case 19:
                printf("\nError allocating in PcupLow\n");
                break;
            case 20:
                printf("\nError opening coefficient file\n");
                break;
            case 21:
                printf("\nError: UnitDepth too large\n");
                break;
            case 22:
                printf("\nYour system needs Big endian version of EGM9615.BIN.  \n");
                printf("Please download this file from http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml.  \n");
                printf("Replace the existing EGM9615.BIN file with the downloaded one\n");
                break;
            default:
                printf("\nError: Unknown error\n");
        }
    }

    void PrintGradient(const Gradient &gradient) {
        printf("\nGradient\n");
        printf("\n                 Northward       Eastward        Downward\n");
        printf("X:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.X, gradient.GradLambda.X,
               gradient.GradZ.X);
        printf("Y:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.Y, gradient.GradLambda.Y,
               gradient.GradZ.Y);
        printf("Z:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.Z, gradient.GradLambda.Z,
               gradient.GradZ.Z);
        printf("H:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.H, gradient.GradLambda.H,
               gradient.GradZ.H);
        printf("F:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.F, gradient.GradLambda.F,
               gradient.GradZ.F);
        printf("Declination: %7.2f min/km %8.2f min/km %8.2f min/km \n", gradient.GradPhi.Decl * 60,
               gradient.GradLambda.Decl * 60, gradient.GradZ.Decl * 60);
        printf("Inclination: %7.2f min/km %8.2f min/km %8.2f min/km \n", gradient.GradPhi.Incl * 60,
               gradient.GradLambda.Incl * 60, gradient.GradZ.Incl * 60);
    }

    /** This function prints the results in  Geomagnetic Elements for a point calculation. It takes the calculated
     *  Geomagnetic elements "GeomagElements" as input.
     *  As well as the coordinates, date, and Magnetic Model.
     * INPUT :  GeomagElements
     *          CoordGeodetic
     *          Date
     *          MagneticModel
     * OUTPUT : none
     */
    void PrintUserData(GeoMagneticElements geomagElements, CoordGeodetic spaceInput, Date timeInput,
                       MagneticModel *magneticModel, Geoid *geoid) {
        std::string DeclString;
        std::string InclString;
        DegreeToDMSstring(geomagElements.Incl, 2, InclString);
        if(geomagElements.H < 6000 && geomagElements.H > 2000) {
            Warnings(1, geomagElements.H, *magneticModel);
        }
        if(geomagElements.H < 2000) {
            Warnings(2, geomagElements.H, *magneticModel);
        }
        if(magneticModel->secularVariationUsed == true) {
            DegreeToDMSstring(geomagElements.Decl, 2, DeclString);
            printf("\n Results For \n\n");
            if(spaceInput.phi < 0) {
                printf("Latitude	%.2fS\n", -spaceInput.phi);
            } else {
                printf("Latitude	%.2fN\n", spaceInput.phi);
            }
            if(spaceInput.lambda < 0) {
                printf("Longitude	%.2fW\n", -spaceInput.lambda);
            } else {
                printf("Longitude	%.2fE\n", spaceInput.lambda);
            }
            if(geoid->UseGeoid == 1) {
                printf("Altitude:	%.2f Kilometers above mean sea level\n", spaceInput.HeightAboveGeoid);
            } else {
                printf("Altitude:	%.2f Kilometers above the WGS-84 ellipsoid\n", spaceInput.HeightAboveEllipsoid);
            }
            printf("Date:		%.1f\n", timeInput.DecimalYear);
            printf("\n		Main Field\t\t\tSecular Change\n");
            printf("F	=	%-9.1f nT\t\t  Fdot = %.1f\tnT/yr\n", geomagElements.F, geomagElements.Fdot);
            printf("H	=	%-9.1f nT\t\t  Hdot = %.1f\tnT/yr\n", geomagElements.H, geomagElements.Hdot);
            printf("X	=	%-9.1f nT\t\t  Xdot = %.1f\tnT/yr\n", geomagElements.X, geomagElements.Xdot);
            printf("Y	=	%-9.1f nT\t\t  Ydot = %.1f\tnT/yr\n", geomagElements.Y, geomagElements.Ydot);
            printf("Z	=	%-9.1f nT\t\t  Zdot = %.1f\tnT/yr\n", geomagElements.Z, geomagElements.Zdot);
            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST)\t  Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * geomagElements.Decldot);
            } else {
                printf("Decl	=%20s  (EAST)\t  Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * geomagElements.Decldot);
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)\t  Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * geomagElements.Incldot);
            } else {
                printf("Incl	=%20s  (DOWN)\t  Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * geomagElements.Incldot);
            }
        } else {
            DegreeToDMSstring(geomagElements.Decl, 2, DeclString);
            printf("\n Results For \n\n");
            if(spaceInput.phi < 0) {
                printf("Latitude	%.2fS\n", -spaceInput.phi);
            } else {
                printf("Latitude	%.2fN\n", spaceInput.phi);
            }
            if(spaceInput.lambda < 0) {
                printf("Longitude	%.2fW\n", -spaceInput.lambda);
            } else {
                printf("Longitude	%.2fE\n", spaceInput.lambda);
            }
            if(geoid->UseGeoid == 1) {
                printf("Altitude:	%.2f Kilometers above MSL\n", spaceInput.HeightAboveGeoid);
            } else {
                printf("Altitude:	%.2f Kilometers above WGS-84 Ellipsoid\n", spaceInput.HeightAboveEllipsoid);
            }
            printf("Date:		%.1f\n", timeInput.DecimalYear);
            printf("\n	Main Field\n");
            printf("F	=	%-9.1f nT\n", geomagElements.F);
            printf("H	=	%-9.1f nT\n", geomagElements.H);
            printf("X	=	%-9.1f nT\n", geomagElements.X);
            printf("Y	=	%-9.1f nT\n", geomagElements.Y);
            printf("Z	=	%-9.1f nT\n", geomagElements.Z);
            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST)\n", DeclString.c_str());
            } else {
                printf("Decl	=%20s  (EAST)\n", DeclString.c_str());
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)\n", InclString.c_str());
            } else {
                printf("Incl	=%20s  (DOWN)\n", InclString.c_str());
            }
        }

        if(spaceInput.phi <= -55 || spaceInput.phi >= 55)
        /* Print Grid Variation */
        {
            DegreeToDMSstring(geomagElements.GV, 2, InclString);
            printf("\n\n Grid variation =%20s\n", InclString.c_str());
        }

    } /*PrintUserData*/

    /** Return value 0 means end program, Return value 1 means get new data, Return value 2 means continue.
     * This prints a warning to the screen determined by the control integer. It also takes the value of the parameter
     * causing the warning as a double.  This is unnecessary for some warnings. It requires the MagneticModel to determine
     * the current epoch.
     *
     * INPUT control :int : (Warning number)
     *       value   : double: Magnetic field strength
     *       MagneticModel
     *       OUTPUT : none
     *       CALLS : none
     */
    int Warnings(int control, double value, const MagneticModel &magneticModel) {
        std::string ans{" "};

        switch(control) {
            case 1: /* Horizontal Field strength low */
                do {
                    printf("\nCaution: location is approaching the blackout zone around the magnetic pole as\n");
                    printf("      defined by the WMM military specification \n");
                    printf("      (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass\n");
                    printf("      accuracy may be degraded in this region.\n");
                    printf("Press enter to continue...\n");
                } while(nullptr == fgets(ans.data(), sizeof(ans), stdin));
                break;
            case 2: /* Horizontal Field strength very low */
                do {
                    printf("\nWarning: location is in the blackout zone around the magnetic pole as defined\n");
                    printf("      by the WMM military specification \n");
                    printf("      (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass\n");
                    printf("      accuracy is highly degraded in this region.\n");
                } while(nullptr == fgets(ans.data(), sizeof(ans), stdin));
                break;
            case 3: /* Elevation outside the recommended range */
                printf("\nWarning: The value you have entered of %.1f km for the elevation is outside of the recommended "
                       "range.\n Elevations above -10.0 km are recommended for accurate results. \n",
                       value);
                while(true) {
                    printf("\nPlease press 'C' to continue, 'G' to get new data or 'X' to exit...\n");
                    while(nullptr == fgets(ans.data(), sizeof(ans), stdin)) {
                        printf("\nInvalid input\n");
                    }
                    switch(ans[0]) {
                        case 'X':
                        case 'x':
                            return 0;
                        case 'G':
                        case 'g':
                            return 1;
                        case 'C':
                        case 'c':
                            return 2;
                        default:
                            printf("\nInvalid input %c\n", ans[0]);
                            break;
                    }
                }
                break;

            case 4: /*Date outside the recommended range*/
                printf("\nWARNING - TIME EXTENDS BEYOND INTENDED USAGE RANGE\n CONTACT NCEI FOR PRODUCT UPDATES:\n");
                printf("	National Centers for Environmental Information\n");
                printf("	NOAA E/NE42\n");
                printf("	325 Broadway\n");
                printf("\n	Boulder, CO 80305 USA");
                printf("	Attn: Manoj Nair or Arnaud Chulliat\n");
                printf("	Phone:	(303) 497-4642 or -6522\n");
                printf("	Email:	geomag.models@noaa.gov\n");
                printf("	Web: http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml\n");
                printf("\n VALID RANGE  = %d - %d\n", (int)magneticModel.min_year,
                       (int)magneticModel.coefficientFileEndDate);
                printf(" TIME   = %f\n", value);
                while(true) {
                    printf("\nPlease press 'C' to continue, 'N' to enter new data or 'X' to exit...\n");
                    while(nullptr == fgets(ans.data(), sizeof(ans), stdin)) {
                        printf("\nInvalid input\n");
                    }
                    switch(ans[0]) {
                        case 'X':
                        case 'x':
                            return 0;
                        case 'N':
                        case 'n':
                            return 1;
                        case 'C':
                        case 'c':
                            return 2;
                        default:
                            printf("\nInvalid input %c\n", ans[0]);
                            break;
                    }
                }
                break;
            case 5: /*Elevation outside the allowable range*/
                printf("\nError: The value you have entered of %f km for the elevation is outside of the recommended "
                       "range.\n Elevations above -10.0 km are recommended for accurate results. \n",
                       value);
                while(true) {
                    printf("\nPlease press 'C' to continue, 'G' to get new data or 'X' to exit...\n");
                    while(nullptr == fgets(ans.data(), sizeof(ans), stdin)) {
                        printf("\nInvalid input\n");
                    }
                    switch(ans[0]) {
                        case 'X':
                        case 'x':
                            return 0;
                        case 'G':
                        case 'g':
                            return 1;
                        case 'C':
                        case 'c':
                            return 2;
                        default:
                            printf("\nInvalid input %c\n", ans[0]);
                            break;
                    }
                }
                break;
        }
        return 2;
    } /*Warnings*/

    /*End of User Interface functions*/


    /******************************************************************************
     ********************************Memory and File Processing********************
     * This grouping consists of functions that read coefficient files into the
     * memory, allocate memory, free memory or print models into coefficient files.
     ******************************************************************************/

    /* This function assigns the first nMax degrees of the Source model to the Assignee model, leaving the other coefficients
 untouched*/
    void AssignMagneticModelCoeffs(MagneticModel *Assignee, MagneticModel *Source, int nMax, int nMaxSecVar) {
        int n, m, index;
        assert(nMax <= Source->nMax);
        assert(nMax <= Assignee->nMax);
        assert(nMaxSecVar <= Source->nMaxSecVar);
        assert(nMaxSecVar <= Assignee->nMaxSecVar);
        for(n = 1; n <= nMaxSecVar; n++) {
            for(m = 0; m <= n; m++) {
                index                                = (n * (n + 1) / 2 + m);
                Assignee->main_Field_Coeff_G[index]  = Source->main_Field_Coeff_G[index];
                Assignee->main_Field_Coeff_H[index]  = Source->main_Field_Coeff_H[index];
                Assignee->secular_Var_Coeff_G[index] = Source->secular_Var_Coeff_G[index];
                Assignee->secular_Var_Coeff_H[index] = Source->secular_Var_Coeff_H[index];
            }
        }
        for(n = nMaxSecVar + 1; n <= nMax; n++) {
            for(m = 0; m <= n; m++) {
                index                               = (n * (n + 1) / 2 + m);
                Assignee->main_Field_Coeff_G[index] = Source->main_Field_Coeff_G[index];
                Assignee->main_Field_Coeff_H[index] = Source->main_Field_Coeff_H[index];
            }
        }
    } /*AssignMagneticModelCoeffs*/

    void PrintWMMFormat(char *filename, const MagneticModel &magneticModel) {
        int index;
        FILE *OUT;
        Date Date;
        char Datestring[11];

        Date.DecimalYear = magneticModel.editionDate;
        YearToDate(&Date);
        sprintf(Datestring, "%d/%d/%d", Date.Month, Date.Day, Date.Year);
        OUT = fopen(filename, "w");
        fprintf(OUT, "    %.1f               %s              %s\n", magneticModel.epoch, magneticModel.modelName.c_str(),
                Datestring);
        for(int n = 1; n <= magneticModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);
                if(m != 0) {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f  %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index],
                            magneticModel.main_Field_Coeff_H[index], magneticModel.secular_Var_Coeff_G[index],
                            magneticModel.secular_Var_Coeff_H[index]);
                } else {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f  %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index], 0.0,
                            magneticModel.secular_Var_Coeff_G[index], 0.0);
                }
            }
        }
        fclose(OUT);
    } /*PrintWMMFormat*/

    void PrintEMMFormat(char *filename, char *filenameSV, const MagneticModel &magneticModel) {
        int index, n, m;
        FILE *OUT;
        Date Date;
        char Datestring[11];

        Date.DecimalYear = magneticModel.editionDate;
        YearToDate(&Date);
        sprintf(Datestring, "%d/%d/%d", Date.Month, Date.Day, Date.Year);
        OUT = fopen(filename, "w");
        fprintf(OUT, "    %.1f               %s              %s\n", magneticModel.epoch, magneticModel.modelName.c_str(),
                Datestring);
        for(n = 1; n <= magneticModel.nMax; n++) {
            for(m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);
                if(m != 0) {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index],
                            magneticModel.main_Field_Coeff_H[index]);
                } else {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index], 0.0);
                }
            }
        }
        fclose(OUT);
        OUT = fopen(filenameSV, "w");
        for(n = 1; n <= magneticModel.nMaxSecVar; n++) {
            for(m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);
                if(m != 0) {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.secular_Var_Coeff_G[index],
                            magneticModel.secular_Var_Coeff_H[index]);
                } else {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.secular_Var_Coeff_G[index], 0.0);
                }
            }
        }
        fclose(OUT);
    } /*PrintEMMFormat*/

    void PrintSHDFFormat(char *filename, MagneticModel *(*magneticModel)[], int epochs) {
        int i, n, m, index, epochRange;
        FILE *SHDF_file;
        SHDF_file = fopen(filename, "w");
        /*lines = (int)(UFM_DEGREE / 2.0 * (UFM_DEGREE + 3));*/
        for(i = 0; i < epochs; i++) {
            if(i < epochs - 1) {
                epochRange = (*magneticModel)[i + 1]->epoch - (*magneticModel)[i]->epoch;
            } else {
                epochRange = (*magneticModel)[i]->epoch - (*magneticModel)[i - 1]->epoch;
            }
            fprintf(SHDF_file, "%%SHDF 16695 Definitive Geomagnetic Reference Field Model Coefficient File\n");
            fprintf(SHDF_file, "%%ModelName: %s\n", (*magneticModel)[i]->modelName.c_str());
            fprintf(SHDF_file,
                    "%%Publisher: International Association of Geomagnetism and Aeronomy (IAGA), Working Group V-Mod\n");
            fprintf(SHDF_file, "%%ReleaseDate: Some Number\n");
            fprintf(SHDF_file, "%%DataCutOFF: Some Other Number\n");
            fprintf(SHDF_file, "%%ModelStartYear: %d\n", (int)(*magneticModel)[i]->epoch);
            fprintf(SHDF_file, "%%ModelEndYear: %d\n", (int)(*magneticModel)[i]->epoch + epochRange);
            fprintf(SHDF_file, "%%Epoch: %.0f\n", (*magneticModel)[i]->epoch);
            fprintf(SHDF_file, "%%IntStaticDeg: %d\n", (*magneticModel)[i]->nMax);
            fprintf(SHDF_file, "%%IntSecVarDeg: %d\n", (*magneticModel)[i]->nMaxSecVar);
            fprintf(SHDF_file, "%%ExtStaticDeg: 0\n");
            fprintf(SHDF_file, "%%ExtSecVarDeg: 0\n");
            fprintf(SHDF_file, "%%Normalization: Schmidt semi-normailized\n");
            fprintf(SHDF_file, "%%SpatBasFunc: spherical harmonics\n");
            fprintf(SHDF_file, "# To synthesize the field for a given date:\n");
            fprintf(SHDF_file, "# Use the sub-model of the epoch corresponding to each date\n");
            fprintf(SHDF_file, "#\n#\n#\n#\n# I/E, n, m, Gnm, Hnm, SV-Gnm, SV-Hnm\n#\n");
            n = 1;
            m = 0;
            for(n = 1; n <= (*magneticModel)[i]->nMax; n++) {
                for(m = 0; m <= n; m++) {
                    index = (n * (n + 1)) / 2 + m;
                    if(i < epochs - 1) {
                        if(m != 0) {
                            fprintf(SHDF_file, "I,%d,%d,%f,%f,%f,%f\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->main_Field_Coeff_H[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_H[index]);
                        } else {
                            fprintf(SHDF_file, "I,%d,%d,%f,,%f,\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index]);
                        }
                    } else {
                        if(m != 0) {
                            fprintf(SHDF_file, "I,%d,%d,%f,%f,%f,%f\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->main_Field_Coeff_H[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_H[index]);
                        } else {
                            fprintf(SHDF_file, "I,%d,%d,%f,,%f,\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index]);
                        }
                    }
                }
            }
        }
    } /*PrintSHDFFormat*/

    /** READ WORLD Magnetic MODEL SPHERICAL HARMONIC COEFFICIENTS (WMM.cof)
     *  INPUT :   filename
     *  UPDATES : MagneticModel
     *  CALLS : none */
    int ReadMagneticModel(const std::string_view filename, MagneticModel *magneticModel) {
        std::vector<std::string> header;
        std::vector<std::string> values;
        std::string str;

        std::ifstream file{filename.data()};
        if(!file.is_open()) {
            PrintError(20);
            return false;
        }

        magneticModel->main_Field_Coeff_H.at(0)  = 0.0;
        magneticModel->main_Field_Coeff_G.at(0)  = 0.0;
        magneticModel->secular_Var_Coeff_H.at(0) = 0.0;
        magneticModel->secular_Var_Coeff_G.at(0) = 0.0;

        for(std::string line; std::getline(file, line, ' ');) {
            if(!line.empty()) {
                auto &val = header.emplace_back(line);
                if(val.back() == '\n') {
                    val.pop_back();
                    break;
                }
            }
        }
        double epoch             = std::stod(header.at(0));
        magneticModel->modelName = header.at(1);
        std::string edit_date{header.at(2)};

        magneticModel->min_year = dateStr_to_decYear(edit_date);
        if(magneticModel->min_year == -1) {
            magneticModel->min_year = epoch;
        }

        magneticModel->epoch = epoch;
        bool isEnd           = false;
        while(!isEnd) {
            values.clear();

            for(std::string line; std::getline(file, line, ' ');) {
                if(!line.empty()) {
                    if(line.find('\n') != std::string::npos) {
                        isEnd = line.find("999999999") != std::string::npos;
                        line  = line.substr(0, line.find('\n'));
                        values.emplace_back(line);
                        break;
                    }
                    values.emplace_back(line);
                }
            }

            /* END OF FILE NOT ENCOUNTERED, GET VALUES */
            int n       = std::stoi(values.at(N));
            int m       = std::stoi(values.at(M));
            double gnm  = std::stod(values.at(GNM));
            double hnm  = std::stod(values.at(HNM));
            double dgnm = std::stod(values.at(DGNM));
            double dhnm = std::stod(values.at(DHNM));
            if(m <= n) {
                int index                                    = (n * (n + 1) / 2 + m);
                magneticModel->main_Field_Coeff_G.at(index)  = gnm;
                magneticModel->secular_Var_Coeff_G.at(index) = dgnm;
                magneticModel->main_Field_Coeff_H.at(index)  = hnm;
                magneticModel->secular_Var_Coeff_H.at(index) = dhnm;
            }
        }

        return true;
    }

    /**  To read the high-degree model coefficients (for example, NGDC 720)
     * INPUT :  filename   file name for static coefficients
     *                    filenameSV file name for secular variation coefficients
     *
     *    MagneticModel : Pointer to the data structure with the following fields required as inputs
     *                        nMaxSecVar : Number of secular variation coefficients
     *                                         nMax : 	Number of static coefficients
     *                                                    UPDATES : MagneticModel : Pointer to the data structure with the
     * following fields populated double epoch;       Base time of Geomagnetic model epoch (yrs) double *Main_Field_Coeff_G;
     * C - Gauss coefficients of main geomagnetic model (nT) double *Main_Field_Coeff_H;          C - Gauss coefficients of
     * main geomagnetic model (nT) double *Secular_Var_Coeff_G;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
     *             double *Secular_Var_Coeff_H;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
     *                 CALLS : none
     */
    //    int readMagneticModel_Large(char *filename, char *filenameSV, MagneticModel *magneticModel){
    //        FILE *COF_File;
    //        FILE *COFSV_File;
    //        char c_str[81], c_str2[81]; /* these strings are used to read a line from coefficient file */
    //        int i, m, n, index, a, b;
    //        double epoch, gnm, hnm, dgnm, dhnm;
    //        COF_File   = fopen(filename, "r");
    //        COFSV_File = fopen(filenameSV, "r");
    //        if(COF_File == NULL || COFSV_File == NULL) {
    //            PrintError(20);
    //            return false;
    //        }
    //        magneticModel->Main_Field_Coeff_H[0]  = 0.0;
    //        magneticModel->Main_Field_Coeff_G[0]  = 0.0;
    //        magneticModel->Secular_Var_Coeff_H[0] = 0.0;
    //        magneticModel->Secular_Var_Coeff_G[0] = 0.0;
    //        if(NULL == fgets(c_str, sizeof(c_str), COF_File)) {
    //            fclose(COF_File);
    //            fclose(COFSV_File);
    //            return false;
    //        }
    //        snprintf(c_str, sizeof(c_str), "%lf%s", epoch, magneticModel->ModelName);
    //        magneticModel->epoch = epoch;
    //        a                    = CALCULATE_NUMTERMS(magneticModel->nMaxSecVar);
    //        b                    = CALCULATE_NUMTERMS(magneticModel->nMax);
    //        for(i = 0; i < a; i++) {
    //            if(NULL == fgets(c_str, sizeof(c_str), COF_File)) {
    //                fclose(COF_File);
    //                fclose(COFSV_File);
    //                return false;
    //            }
    //            sscanf(c_str, "%d%d%lf%lf", &n, &m, &gnm, &hnm);
    //            if(NULL == fgets(c_str2, sizeof(c_str2), COFSV_File)) {
    //                fclose(COF_File);
    //                fclose(COFSV_File);
    //                return false;
    //            }
    //            sscanf(c_str2, "%d%d%lf%lf", &n, &m, &dgnm, &dhnm);
    //            if(m <= n) {
    //                index                                     = (n * (n + 1) / 2 + m);
    //                magneticModel->Main_Field_Coeff_G[index]  = gnm;
    //                magneticModel->Secular_Var_Coeff_G[index] = dgnm;
    //                magneticModel->Main_Field_Coeff_H[index]  = hnm;
    //                magneticModel->Secular_Var_Coeff_H[index] = dhnm;
    //            }
    //        }
    //        for(i = a; i < b; i++) {
    //            if(NULL == fgets(c_str, sizeof(c_str), COF_File)) {
    //                fclose(COF_File);
    //                fclose(COFSV_File);
    //                return false;
    //            }
    //            sscanf(c_str, "%d%d%lf%lf", &n, &m, &gnm, &hnm);
    //            if(m <= n) {
    //                index                                    = (n * (n + 1) / 2 + m);
    //                magneticModel->Main_Field_Coeff_G[index] = gnm;
    //                magneticModel->Main_Field_Coeff_H[index] = hnm;
    //            }
    //        }
    //        if(COF_File != NULL && COFSV_File != NULL) {
    //            fclose(COF_File);
    //            fclose(COFSV_File);
    //        }
    //
    //        return true;
    //    } /*readMagneticModel_Large*/


    /** readMagneticModels - Read the Magnetic Models from an SHDF format file
     *
     * Input:
     *  filename - Path to the SHDF format model file to be read
     *  array_size - Max No of models to be read from the file
     *
     * Output:
     *  magneticmodels[] - Array of magnetic models read from the file
     *
     * Return value:
     *  Returns the number of models read from the file.
     *  -2 implies that internal or external static degree was not found in the file, hence memory cannot be allocated
     *  -1 implies some error during file processing (I/O)
     *  0 implies no models were read from the file
     *  if ReturnValue > array_size then there were too many models in model file but only <array_size> number were read .
     *  if ReturnValue <= array_size then the function execution was successful.
     */
    //    int readMagneticModel_SHDF(char *filename, MagneticModel *(*magneticmodels)[], int array_size)    {
    //        char paramkeys[NOOFPARAMS][MAXLINELENGTH] = {
    //            "SHDF ",          "ModelName: ",     "Publisher: ",    "ReleaseDate: ",  "DataCutOff: ",   "ModelStartYear:
    //            ", "ModelEndYear: ", "Epoch: ",         "IntStaticDeg: ", "IntSecVarDeg: ", "ExtStaticDeg: ",
    //            "ExtSecVarDeg: ", "GeoMagRefRad: ", "Normalization: ", "SpatBasFunc: "};
    //
    //        char paramvalues[NOOFPARAMS][MAXLINELENGTH];
    //        char *line = (char *)malloc(MAXLINELENGTH);
    //        char *ptrreset;
    //        char paramvalue[MAXLINELENGTH];
    //        int paramvaluelength = 0;
    //        int paramkeylength   = 0;
    //        int i = 0, j = 0;
    //        int newrecord    = 1;
    //        int header_index = -1;
    //        int numterms;
    //        int tempint;
    //        int allocationflag = 0;
    //        char coefftype; /* Internal or External (I/E) */
    //
    //        /* For reading coefficients */
    //        int n, m;
    //        double gnm, hnm, dgnm, dhnm, cutoff;
    //        int index;
    //
    //        FILE *stream;
    //        ptrreset = line;
    //        stream   = fopen(filename, READONLYMODE);
    //        if(stream == NULL) {
    //            perror("File open error");
    //            return header_index;
    //        }
    //
    //        /* Read records from the model file and store header information. */
    //        while(fgets(line, MAXLINELENGTH, stream) != NULL) {
    //            j++;
    //            if(strlen(Trim(line)) == 0) {
    //                continue;
    //            }
    //            if(*line == '%') {
    //                line++;
    //                if(newrecord) {
    //                    if(header_index > -1) {
    //                        AssignHeaderValues((*magneticmodels)[header_index], paramvalues);
    //                    }
    //                    header_index++;
    //                    if(header_index >= array_size) {
    //                        fprintf(stderr, "Header limit exceeded - too many models in model file. (%d)\n", header_index);
    //                        return array_size + 1;
    //                    }
    //                    newrecord      = 0;
    //                    allocationflag = 0;
    //                }
    //                for(i = 0; i < NOOFPARAMS; i++) {
    //                    paramkeylength = strlen(paramkeys[i]);
    //                    if(!strncmp(line, paramkeys[i], paramkeylength)) {
    //                        paramvaluelength = strlen(line) - paramkeylength;
    //                        memset(paramvalues, '\0', paramvaluelength);
    //                        strlcpy_equivalent(paramvalue, line + paramkeylength, paramvaluelength);
    //                        paramvalue[paramvaluelength] = '\0';
    //                        strlcpy_equivalent(paramvalues[i], paramvalue, 1);
    //                        if(!strcmp(paramkeys[i], paramkeys[INTSTATICDEG]) ||
    //                           !strcmp(paramkeys[i], paramkeys[EXTSTATICDEG])) {
    //                            tempint = atoi(paramvalues[i]);
    //                            if(tempint > 0 && allocationflag == 0) {
    //                                numterms                        = CALCULATE_NUMTERMS(tempint);
    //                                (*magneticmodels)[header_index] = AllocateModelMemory(numterms);
    //                                /* model = (*magneticmodels)[header_index]; */
    //                                allocationflag = 1;
    //                            }
    //                        }
    //                        break;
    //                    }
    //                }
    //                line--;
    //            } else if(*line == '#') {
    //                /* process comments */
    //
    //            } else if(sscanf(line, "%c,%d,%d", &coefftype, &n, &m) == 3) {
    //                if(m == 0) {
    //                    sscanf(line, "%c,%d,%d,%lf,,%lf,", &coefftype, &n, &m, &gnm, &dgnm);
    //                    hnm  = 0;
    //                    dhnm = 0;
    //                } else {
    //                    sscanf(line, "%c,%d,%d,%lf,%lf,%lf,%lf", &coefftype, &n, &m, &gnm, &hnm, &dgnm, &dhnm);
    //                }
    //                newrecord = 1;
    //                if(!allocationflag) {
    //                    fprintf(stderr, "Degree not found in model. Memory cannot be allocated.\n");
    //                    return _DEGREE_NOT_FOUND;
    //                }
    //                if(m <= n) {
    //                    index                                                       = (n * (n + 1) / 2 + m);
    //                    (*magneticmodels)[header_index]->Main_Field_Coeff_G[index]  = gnm;
    //                    (*magneticmodels)[header_index]->Secular_Var_Coeff_G[index] = dgnm;
    //                    (*magneticmodels)[header_index]->Main_Field_Coeff_H[index]  = hnm;
    //                    (*magneticmodels)[header_index]->Secular_Var_Coeff_H[index] = dhnm;
    //                }
    //            }
    //        }
    //        if(header_index > -1) {
    //            AssignHeaderValues((*magneticmodels)[header_index], paramvalues);
    //        }
    //        fclose(stream);
    //
    //        cutoff = (*magneticmodels)[array_size - 1]->CoefficientFileEndDate;
    //
    //        for(i = 0; i < array_size; i++) {
    //            (*magneticmodels)[i]->CoefficientFileEndDate = cutoff;
    //        }
    //
    //        free(ptrreset);
    //        line     = NULL;
    //        ptrreset = NULL;
    //        return header_index + 1;
    //    } /*readMagneticModel_SHDF*/

    //    char *Trim(char *str) {
    //        char *end;
    //
    //        while(isspace(*str)) {
    //            str++;
    //        }
    //
    //        if(*str == 0) {
    //            return str;
    //        }
    //
    //        end = str + strlen(str) - 1;
    //        while(end > str && isspace(*end)) {
    //            end--;
    //        }
    //
    //        *(end + 1) = 0;
    //
    //        return str;
    //    }

    /*End of Memory and File Processing functions*/


    /******************************************************************************
     *************Conversions, Transformations, and other Calculations**************
     * This grouping consists of functions that perform unit conversions, coordinate
     * transformations and other simple or straightforward calculations that are
     * usually easily replicable with a typical scientific calculator.
     ******************************************************************************/


    void BaseErrors(double DeclCoef, double DeclBaseline, double InclOffset, double FOffset, double Multiplier, double H,
                    double *DeclErr, double *InclErr, double *FErr) {
        double declHorizontalAdjustmentSq;
        declHorizontalAdjustmentSq = (DeclCoef / H) * (DeclCoef / H);
        *DeclErr                   = sqrt(declHorizontalAdjustmentSq + DeclBaseline * DeclBaseline) * Multiplier;
        *InclErr                   = InclOffset * Multiplier;
        *FErr                      = FOffset * Multiplier;
    }

    /** Calculate all the Geomagnetic elements from X,Y and Z components
     * INPUT     MagneticResultsGeo
     * OUTPUT    GeoMagneticElements
     * CALLS : none
     */
    int CalculateGeoMagneticElements(const MagneticResults &magneticResultsGeo, GeoMagneticElements *geoMagneticElements) {
        geoMagneticElements->X = magneticResultsGeo.Bx;
        geoMagneticElements->Y = magneticResultsGeo.By;
        geoMagneticElements->Z = magneticResultsGeo.Bz;

        geoMagneticElements->H =
            sqrt(magneticResultsGeo.Bx * magneticResultsGeo.Bx + magneticResultsGeo.By * magneticResultsGeo.By);
        geoMagneticElements->F =
            sqrt(geoMagneticElements->H * geoMagneticElements->H + magneticResultsGeo.Bz * magneticResultsGeo.Bz);
        geoMagneticElements->Decl = Rad2Deg(atan2(geoMagneticElements->Y, geoMagneticElements->X));
        geoMagneticElements->Incl = Rad2Deg(atan2(geoMagneticElements->Z, geoMagneticElements->H));

        return true;
    }

    /** Computes the grid variation for |latitudes| > MAX_LAT_DEGREE
     * Grivation (or grid variation) is the angle between grid north and magnetic north. This routine calculates Grivation
     * for the Polar Stereographic projection for polar locations (Latitude => |55| deg). Otherwise, it computes the grid
     * variation in UTM projection system. However, the UTM projection codes may be used to compute the grid variation at any
     * latitudes.
     *
     * INPUT :  CoordGeodetic location
     * OUTPUT:  GeoMagneticElements elements
     * CALLS :  GetTransverseMercator     */
    int CalculateGridVariation(const CoordGeodetic &location, GeoMagneticElements *elements) {
        UTMParameters UTMParameters;

        if(location.phi >= PS_MAX_LAT_DEGREE) {
            elements->GV = elements->Decl - location.lambda;
            return 1;
        }

        if(location.phi <= PS_MIN_LAT_DEGREE) {
            elements->GV = elements->Decl + location.lambda;
            return 1;
        }

        GetTransverseMercator(location, &UTMParameters);
        elements->GV = elements->Decl - UTMParameters.ConvergenceOfMeridians;
        return 0;
    }

    void CalculateGradientElements(const MagneticResults &gradResults, const GeoMagneticElements &magneticElements,
                                   GeoMagneticElements *gradElements) {
        gradElements->X = gradResults.Bx;
        gradElements->Y = gradResults.By;
        gradElements->Z = gradResults.Bz;

        gradElements->H = (gradElements->X * magneticElements.X + gradElements->Y * magneticElements.Y) / magneticElements.H;
        gradElements->F = (gradElements->X * magneticElements.X + gradElements->Y * magneticElements.Y +
                           gradElements->Z * magneticElements.Z) /
                          magneticElements.F;
        gradElements->Decl = 180.0 / M_PI * (magneticElements.X * gradElements->Y - magneticElements.Y * gradElements->X) /
                             (magneticElements.H * magneticElements.H);
        gradElements->Incl = 180.0 / M_PI * (magneticElements.H * gradElements->Z - magneticElements.Z * gradElements->H) /
                             (magneticElements.F * magneticElements.F);
        gradElements->GV = gradElements->Decl;
    }

    /**This takes the Magnetic Variation in x, y, and z and uses it to calculate the secular variation of each of the
     * Geomagnetic elements.
     * INPUT : MagneticVariation   Data structure with the following elements double Bx;    ( North )
     * OUTPUT: MagneticElements    Pointer to the data  structure with the following elements updated
     * CALLS : none
     */
    int CalculateSecularVariationElements(const MagneticResults &magneticVariation, GeoMagneticElements *magneticElements) {
        magneticElements->Xdot = magneticVariation.Bx;
        magneticElements->Ydot = magneticVariation.By;
        magneticElements->Zdot = magneticVariation.Bz;
        magneticElements->Hdot =
            (magneticElements->X * magneticElements->Xdot + magneticElements->Y * magneticElements->Ydot) /
            magneticElements->H; /* See equation 19 in the WMM technical report */
        magneticElements->Fdot =
            (magneticElements->X * magneticElements->Xdot + magneticElements->Y * magneticElements->Ydot +
             magneticElements->Z * magneticElements->Zdot) /
            magneticElements->F;
        magneticElements->Decldot =
            180.0 / M_PI * (magneticElements->X * magneticElements->Ydot - magneticElements->Y * magneticElements->Xdot) /
            (magneticElements->H * magneticElements->H);
        magneticElements->Incldot =
            180.0 / M_PI * (magneticElements->H * magneticElements->Zdot - magneticElements->Z * magneticElements->Hdot) /
            (magneticElements->F * magneticElements->F);
        magneticElements->GVdot = magneticElements->Decldot;
        return true;
    }

    /*This converts the Cartesian x, y, and z coordinates to Geodetic Coordinates
     * x is defined as the direction pointing out of the core toward the point defined
     * by 0 degrees latitude and longitude.
     * y is defined as the direction from the core toward 90 degrees east longitude along
     * the equator
     * z is defined as the direction from the core out the geographic north pole
     */
    void CartesianToGeodetic(const Ellipsoid &ellip, double x, double y, double z, CoordGeodetic *coordGeodetic) {
        double modified_b, r, e, f, p, q, d, v, g, t, zlong, rlat;

        /*   1.0 compute semi-minor axis and set sign to that of z in order
         *   to get sign of Phi correct
         */

        if(z < 0.0) {
            modified_b = -ellip.b;
        } else {
            modified_b = ellip.b;
        }

        /*   2.0 compute intermediate values for latitude
         */
        r = sqrt(x * x + y * y);
        e = (modified_b * z - (ellip.a * ellip.a - modified_b * modified_b)) / (ellip.a * r);
        f = (modified_b * z + (ellip.a * ellip.a - modified_b * modified_b)) / (ellip.a * r);

        /*   3.0 find solution to:
         *       t^4 + 2*E*t^3 + 2*F*t - 1 = 0
         */
        p = (4.0 / 3.0) * (e * f + 1.0);
        q = 2.0 * (e * e - f * f);
        d = p * p * p + q * q;

        if(d >= 0.0) {
            v = pow((sqrt(d) - q), (1.0 / 3.0)) - pow((sqrt(d) + q), (1.0 / 3.0));
        } else {
            v = 2.0 * sqrt(-p) * cos(acos(q / (p * sqrt(-p))) / 3.0);
        }

        /*   4.0 improve v
         *   NOTE: not really necessary unless point is near pole
         */
        if(v * v < fabs(p)) {
            v = -(v * v * v + 2.0 * q) / (3.0 * p);
        }
        g = (sqrt(e * e + v) + e) / 2.0;
        t = sqrt(g * g + (f - v * g) / (2.0 * g - e)) - g;

        rlat               = atan((ellip.a * (1.0 - t * t)) / (2.0 * modified_b * t));
        coordGeodetic->phi = Rad2Deg(rlat);

        /*   5.0 compute height above ellipsoid
         */
        coordGeodetic->HeightAboveEllipsoid = (r - ellip.a * t) * cos(rlat) + (z - modified_b) * sin(rlat);

        /*   6.0 compute longitude east of Greenwich
         */
        zlong = atan2(y, x);
        if(zlong < 0.0) {
            zlong = zlong + 2 * M_PI;
        }

        coordGeodetic->lambda = Rad2Deg(zlong);
        while(coordGeodetic->lambda > 180) {
            coordGeodetic->lambda -= 360;
        }
    }

    CoordGeodetic CoordGeodeticAssign(const CoordGeodetic &coordGeodetic) {
        CoordGeodetic Assignee;
        Assignee.phi                  = coordGeodetic.phi;
        Assignee.lambda               = coordGeodetic.lambda;
        Assignee.HeightAboveEllipsoid = coordGeodetic.HeightAboveEllipsoid;
        Assignee.HeightAboveGeoid     = coordGeodetic.HeightAboveGeoid;
        Assignee.UseGeoid             = coordGeodetic.UseGeoid;
        return Assignee;
    }

    /** Converts a given calendar date into a decimal year,
     * it also outputs an error string if there is a problem
     * INPUT  CalendarDate  Pointer to the  data  structure with the following elements
     * OUTPUT CalendarDate  Pointer to the  data  structure with the following elements updated
     *        Error	        pointer to an error string
     * CALLS : none
     */
    int DateToYear(Date *calendarDate, std::string &error) {
        int totalDays = 0; /*Total number of days */
        static std::array<int, 13> days{0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        int extraDay = 0;

        if(calendarDate->Month == 0) {
            calendarDate->DecimalYear = calendarDate->Year;
            return true;
        }

        if((calendarDate->Year % 4 == 0 && calendarDate->Year % 100 != 0) || calendarDate->Year % 400 == 0) {
            extraDay = 1;
        }
        days[2] = 28 + extraDay;

        /******************Validation********************************/
        if(calendarDate->Month <= 0 || calendarDate->Month > 12) {
            // The Error is passed as pointer and its size if defined outside of the function. The functions used to pass
            // Error to DateToYear() defines the size of DMSstring are all 255.
            error = "\nError: The Month entered is invalid, valid months are '1 to 12'\n";
            return 0;
        }
        if(calendarDate->Day <= 0 || calendarDate->Day > days[calendarDate->Month]) {
            printf("\nThe number of days in month %d is %d\n", calendarDate->Month, days[calendarDate->Month]);
            error = "\nError: The day entered is invalid\n";
            return 0;
        }
        /****************Calculation of t***************************/
        for(int i = 1; i <= calendarDate->Month; i++) {
            totalDays += days[i - 1];
        }

        totalDays += calendarDate->Day;
        calendarDate->DecimalYear = calendarDate->Year + (totalDays - 1) / (365.0 + extraDay);
        return true;
    }

    /** This converts a given decimal degree into a DMS string.
     * INPUT  DegreesOfArc   decimal degree
     *           UnitDepth	How many iterations should be printed,
     *                        1 = Degrees
     *                        2 = Degrees, Minutes
     *                        3 = Degrees, Minutes, Seconds
     * OUPUT  DMSstring 	 pointer to DMSString.  Must be at least 30 characters.
     * CALLS : none
     */
    void DegreeToDMSstring(double DegreesOfArc, int UnitDepth, std::string &out) {
        int DMS;
        double temp = DegreesOfArc;

        if(UnitDepth > 3) {
            PrintError(21);
        }

        for(int i = 0; i < UnitDepth; i++) {
            DMS  = static_cast<int>(temp);
            temp = (temp - DMS) * 60;

            if(i == UnitDepth - 1 && temp >= 30) {
                DMS++;
            } else if(i == UnitDepth - 1 && temp <= -30) {
                DMS--;
            }

            out.append(std::to_string(DMS));

            switch(i) {
                case 0:
                    out.append(" Deg ");
                    break;
                case 1:
                    out.append(" Min ");
                    break;
                case 2:
                    out.append(" Sec ");
                    break;
                default:;
            }
        }
    }

    /** This converts a given DMS string into decimal degrees.
     * INPUT  DMSstring 	 pointer to DMSString
     * OUTPUT  DegreesOfArc   decimal degree
     * CALLS : none
     */
    void DMSstringToDegree(std::string_view DMSstring, double *DegreesOfArc) {
        int second, minute, degree, sign = 1, j = 0;
        j = sscanf(DMSstring.data(), "%d, %d, %d", &degree, &minute, &second);
        if(j != 3) {
            sscanf(DMSstring.data(), "%d %d %d", &degree, &minute, &second);
        }
        if(degree < 0) {
            sign = -1;
        }
        degree        = degree * sign;
        *DegreesOfArc = sign * (degree + minute / 60.0 + second / 3600.0);
    } /*DMSstringToDegree*/

    void ErrorCalc(GeoMagneticElements B, GeoMagneticElements *Errors) {
        /*Errors.Decl, Errors.Incl, Errors.F are all assumed to exist*/
        double cos2D, cos2I, sin2D, sin2I, EDSq, EISq, eD, eI;
        cos2D     = cos(Deg2Rad(B.Decl)) * cos(Deg2Rad(B.Decl));
        cos2I     = cos(Deg2Rad(B.Incl)) * cos(Deg2Rad(B.Incl));
        sin2D     = sin(Deg2Rad(B.Decl)) * sin(Deg2Rad(B.Decl));
        sin2I     = sin(Deg2Rad(B.Incl)) * sin(Deg2Rad(B.Incl));
        eD        = Deg2Rad(Errors->Decl);
        eI        = Deg2Rad(Errors->Incl);
        EDSq      = eD * eD;
        EISq      = eI * eI;
        Errors->X = sqrt(cos2D * cos2I * Errors->F * Errors->F + B.F * B.F * sin2D * cos2I * EDSq +
                         B.F * B.F * cos2D * sin2I * EISq);
        Errors->Y = sqrt(sin2D * cos2I * Errors->F * Errors->F + B.F * B.F * cos2D * cos2I * EDSq +
                         B.F * B.F * sin2D * sin2I * EISq);
        Errors->Z = sqrt(sin2I * Errors->F * Errors->F + B.F * B.F * cos2I * EISq);
        Errors->H = sqrt(cos2I * Errors->F * Errors->F + B.F * B.F * sin2I * EISq);
    }

    /** Convert geodetic coordinates, (defined by the WGS-84 reference ellipsoid), to Earth Centered Earth Fixed
     * Cartesian coordinates, and then to spherical coordinates.
     * INPUT   Ellip  data  structure with the following elements
     *         CoordGeodetic  Pointer to the  data  structure with the following elements updates
     * OUTPUT  CoordSpherical 	Pointer to the data structure with the following elements
     * CALLS : none
     */
    int GeodeticToSpherical(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, CoordSpherical *coordSpherical) {
        const double CosLat = std::cos(Deg2Rad(coordGeodetic.phi));
        const double SinLat = std::sin(Deg2Rad(coordGeodetic.phi));

        /* compute the local radius of curvature on the WGS-84 reference ellipsoid */
        const double rc = ellip.a / std::sqrt(1.0 - ellip.epssq * SinLat * SinLat);

        /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */
        const double xp = (rc + coordGeodetic.HeightAboveEllipsoid) * CosLat;
        const double zp = (rc * (1.0 - ellip.epssq) + coordGeodetic.HeightAboveEllipsoid) * SinLat;

        /* compute spherical radius and angle lambda and phi of specified point */
        coordSpherical->r = std::sqrt(xp * xp + zp * zp);
        /* geocentric latitude */
        coordSpherical->phig = Rad2Deg(std::asin(zp / coordSpherical->r));
        /* longitude */
        coordSpherical->lambda = coordGeodetic.lambda;
        return true;
    }

    GeoMagneticElements GeoMagneticElementsAssign(GeoMagneticElements Elements) {
        GeoMagneticElements Assignee;
        Assignee.X       = Elements.X;
        Assignee.Y       = Elements.Y;
        Assignee.Z       = Elements.Z;
        Assignee.H       = Elements.H;
        Assignee.F       = Elements.F;
        Assignee.Decl    = Elements.Decl;
        Assignee.Incl    = Elements.Incl;
        Assignee.GV      = Elements.GV;
        Assignee.Xdot    = Elements.Xdot;
        Assignee.Ydot    = Elements.Ydot;
        Assignee.Zdot    = Elements.Zdot;
        Assignee.Hdot    = Elements.Hdot;
        Assignee.Fdot    = Elements.Fdot;
        Assignee.Decldot = Elements.Decldot;
        Assignee.Incldot = Elements.Incldot;
        Assignee.GVdot   = Elements.GVdot;
        return Assignee;
    }

    GeoMagneticElements GeoMagneticElementsScale(GeoMagneticElements Elements, double factor) {
        /*This function scales all the geomagnetic elements to scale a vector use
     MagneticResultsScale*/
        GeoMagneticElements product;
        product.X       = Elements.X * factor;
        product.Y       = Elements.Y * factor;
        product.Z       = Elements.Z * factor;
        product.H       = Elements.H * factor;
        product.F       = Elements.F * factor;
        product.Incl    = Elements.Incl * factor;
        product.Decl    = Elements.Decl * factor;
        product.GV      = Elements.GV * factor;
        product.Xdot    = Elements.Xdot * factor;
        product.Ydot    = Elements.Ydot * factor;
        product.Zdot    = Elements.Zdot * factor;
        product.Hdot    = Elements.Hdot * factor;
        product.Fdot    = Elements.Fdot * factor;
        product.Incldot = Elements.Incldot * factor;
        product.Decldot = Elements.Decldot * factor;
        product.GVdot   = Elements.GVdot * factor;
        return product;
    }

    GeoMagneticElements GeoMagneticElementsSubtract(GeoMagneticElements minuend, GeoMagneticElements subtrahend) {
        /*This algorithm does not result in the difference of F being derived from
     the Pythagorean theorem.  This function should be used for computing residuals
     or changes in elements.*/
        GeoMagneticElements difference;
        difference.X = minuend.X - subtrahend.X;
        difference.Y = minuend.Y - subtrahend.Y;
        difference.Z = minuend.Z - subtrahend.Z;

        difference.H    = minuend.H - subtrahend.H;
        difference.F    = minuend.F - subtrahend.F;
        difference.Decl = minuend.Decl - subtrahend.Decl;
        difference.Incl = minuend.Incl - subtrahend.Incl;

        difference.Xdot = minuend.Xdot - subtrahend.Xdot;
        difference.Ydot = minuend.Ydot - subtrahend.Ydot;
        difference.Zdot = minuend.Zdot - subtrahend.Zdot;

        difference.Hdot    = minuend.Hdot - subtrahend.Hdot;
        difference.Fdot    = minuend.Fdot - subtrahend.Fdot;
        difference.Decldot = minuend.Decldot - subtrahend.Decldot;
        difference.Incldot = minuend.Incldot - subtrahend.Incldot;

        difference.GV    = minuend.GV - subtrahend.GV;
        difference.GVdot = minuend.GVdot - subtrahend.GVdot;

        return difference;
    }

    /* Gets the UTM Parameters for a given Latitude and Longitude.
     *INPUT: CoordGeodetic
     *OUTPUT : UTMParameters
     */
    int GetTransverseMercator(const CoordGeodetic &coordGeodetic, UTMParameters *utmParameters) {
        /* WGS84 ellipsoid */
        constexpr double K0     = 0.9996;
        constexpr double Eps    = 0.081819190842621494335;
        constexpr double Epssq  = 0.0066943799901413169961;
        constexpr double K0R4   = 6367449.1458234153093 * K0;
        constexpr double K0R4oa = K0R4 / 6378137;
        constexpr std::array Acoeff{8.37731820624469723600E-04, 7.60852777357248641400E-07, 1.19764550324249124400E-09,
                                    2.42917068039708917100E-12, 5.71181837042801392800E-15, 1.47999793137966169400E-17,
                                    4.10762410937071532000E-20, 1.21078503892257704200E-22};

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
        TMfwd4(Eps, Epssq, K0R4, K0R4oa, Acoeff, Lam0, K0, falseE, falseN, XYonly, lambda, phi, &X, &Y, &pscale, &CoM);

        /** Report results */
        utmParameters->Easting                = X;             /* UTM Easting (X) in meters */
        utmParameters->Northing               = Y;             /* UTM Northing (Y) in meters */
        utmParameters->Zone                   = zone;          /* UTM Zone */
        utmParameters->HemiSphere             = hemisphere;
        utmParameters->CentralMeridian        = Rad2Deg(Lam0); /* Central Meridian of the UTM Zone */
        utmParameters->ConvergenceOfMeridians = Rad2Deg(CoM);  /* Convergence of meridians of the UTM Zone and location */
        utmParameters->PointScale             = pscale;

        return 0;
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
    int GetUTMParameters(double latitude, double longitude, int *zone, char *hemisphere, double *centralMeridian) {
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

    /** Rotate the Magnetic Vectors to Geodetic Coordinates
     * Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
     * Equation 16, WMM Technical report
     * INPUT : CoordSpherical
     *         CoordGeodetic
     *         MagneticResultsSph
     * OUTPUT: MagneticResultsGeo
     * CALLS : none     */
    int RotateMagneticVector(const CoordSpherical &CoordSpherical, const CoordGeodetic &coordGeodetic,
                             const MagneticResults &magneticResultsSph, MagneticResults *magneticResultsGeo) {
        /* Difference between the spherical and Geodetic latitudes */
        const double Psi = (M_PI / 180) * (CoordSpherical.phig - coordGeodetic.phi);

        /* Rotate spherical field components to the Geodetic system */
        magneticResultsGeo->Bz = magneticResultsSph.Bx * sin(Psi) + magneticResultsSph.Bz * cos(Psi);
        magneticResultsGeo->Bx = magneticResultsSph.Bx * cos(Psi) - magneticResultsSph.Bz * sin(Psi);
        magneticResultsGeo->By = magneticResultsSph.By;
        return true;
    }

    void SphericalToCartesian(CoordSpherical coordSpherical, double *x, double *y, double *z) {
        double radphi;
        double radlambda;

        radphi    = coordSpherical.phig * (M_PI / 180);
        radlambda = coordSpherical.lambda * (M_PI / 180);

        *x = coordSpherical.r * cos(radphi) * cos(radlambda);
        *y = coordSpherical.r * cos(radphi) * sin(radlambda);
        *z = coordSpherical.r * sin(radphi);
    }

    void SphericalToGeodetic(Ellipsoid ellip, CoordSpherical coordSpherical, CoordGeodetic *coordGeodetic) {
        /*This converts spherical coordinates back to geodetic coordinates.  It is not used in the WMM but
     may be necessary for some applications, such as geomagnetic coordinates*/
        double x, y, z;

        SphericalToCartesian(coordSpherical, &x, &y, &z);
        CartesianToGeodetic(ellip, x, y, z, coordGeodetic);
    }

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
    void TMfwd4(const double Eps, const double Epssq, const double K0R4, const double K0R4oa,
                const std::array<double, 8> &Acoeff, const double Lam0, const double K0, const double falseE,
                const double falseN, const int XYonly, const double Lambda, const double Phi, double *X, double *Y,
                double *pscale, double *CoM) {
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
        const double P     = exp(Eps * std::atanh(Eps * SPhi));
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
        double Xstar = Acoeff[3] * s8u * c8v;
        Xstar        = Xstar + Acoeff[2] * s6u * c6v;
        Xstar        = Xstar + Acoeff[1] * s4u * c4v;
        Xstar        = Xstar + Acoeff.at(0) * s2u * c2v;
        Xstar        = Xstar + U;

        double Ystar = Acoeff[3] * c8u * s8v;
        Ystar        = Ystar + Acoeff[2] * c6u * s6v;
        Ystar        = Ystar + Acoeff[1] * c4u * s4v;
        Ystar        = Ystar + Acoeff.at(0) * c2u * s2v;
        Ystar        = Ystar + V;

        /** Apply isoperimetric radius, scale adjustment, and offsets  */
        *X = K0R4 * Xstar + falseE;
        *Y = K0R4 * Ystar + falseN;

        /** Point-scale and CoM */
        if(XYonly == 1) {
            *pscale = K0;
            *CoM    = 0;
        } else {
            double sig1 = 8 * Acoeff[3] * c8u * c8v;
            sig1        = sig1 + 6 * Acoeff[2] * c6u * c6v;
            sig1        = sig1 + 4 * Acoeff[1] * c4u * c4v;
            sig1        = sig1 + 2 * Acoeff.at(0) * c2u * c2v;
            sig1        = sig1 + 1;

            double sig2 = 8 * Acoeff[3] * s8u * s8v;
            sig2        = sig2 + 6 * Acoeff[2] * s6u * s6v;
            sig2        = sig2 + 4 * Acoeff[1] * s4u * s4v;
            sig2        = sig2 + 2 * Acoeff.at(0) * s2u * s2v;

            /*    Combined square roots  */
            const double comroo = sqrt((1 - Epssq * SPhi * SPhi) * denom2 * (sig1 * sig1 + sig2 * sig2));

            *pscale = K0R4oa * 2 * denom * comroo;
            *CoM    = atan2(SChi * SLam, CLam) + atan2(sig2, sig1);
        }
    }

    /* Converts a given Decimal year into a Year, Month and Date
it also outputs an error string if there is a problem
INPUT  CalendarDate  Pointer to the  data  structure with the following elements
                double DecimalYear;      decimal years
OUTPUT  CalendarDate  Pointer to the  data  structure with the following elements updated
* int Year
* int Month
* int Day
           Error    pointer to an error string
CALLS : none

    */
    int YearToDate(Date *CalendarDate) {
        int MonthDays[13], CumulativeDays = 0;
        int ExtraDay = 0;
        int i, DayOfTheYear;


        if(CalendarDate->DecimalYear == 0) {
            CalendarDate->Year  = 0;
            CalendarDate->Month = 0;
            CalendarDate->Day   = 0;
            return false;
        }

        CalendarDate->Year = (int)floor(CalendarDate->DecimalYear);


        if((CalendarDate->Year % 4 == 0 && CalendarDate->Year % 100 != 0) || CalendarDate->Year % 400 == 0) {
            ExtraDay = 1;
        }

        DayOfTheYear =
            floor((CalendarDate->DecimalYear - (double)CalendarDate->Year) * (365.0 + (double)ExtraDay) + 0.5) + 1;
        /*The above floor is used for rounding, this only works for positive integers*/


        MonthDays[0]  = 0;
        MonthDays[1]  = 31;
        MonthDays[2]  = 28 + ExtraDay;
        MonthDays[3]  = 31;
        MonthDays[4]  = 30;
        MonthDays[5]  = 31;
        MonthDays[6]  = 30;
        MonthDays[7]  = 31;
        MonthDays[8]  = 31;
        MonthDays[9]  = 30;
        MonthDays[10] = 31;
        MonthDays[11] = 30;
        MonthDays[12] = 31;


        for(i = 1; i <= 12; i++) {
            CumulativeDays = CumulativeDays + MonthDays[i];

            if(DayOfTheYear <= CumulativeDays) {
                CalendarDate->Month = i;
                CalendarDate->Day   = MonthDays[i] - (CumulativeDays - DayOfTheYear);
                break;
            }
        }


        return true;

    } /*YearToDate*/

    /******************************************************************************
     ********************************Spherical Harmonics***************************
     * This grouping consists of functions that together take gauss coefficients
     * and return a magnetic vector for an input location in spherical coordinates
     ******************************************************************************/

    /** Computes  all of the Schmidt-semi normalized associated Legendre
     * functions up to degree nMax. If nMax <= 16, function PcupLow is used.
     * Otherwise PcupHigh is called.
     * INPUT   CoordSpherical
     *         nMax        	integer 	 ( Maxumum degree of spherical harmonic secular model)
     *         LegendreFunction pointer
     * OUTPUT  LegendreFunction  Calculated Legendre variables in the data structure
     */
    int AssociatedLegendreFunction(const CoordSpherical &coordSpherical, const int nMax,
                                   LegendreFunction *legendreFunction) {
        bool FLAG;

        const double sin_phi = sin(Deg2Rad(coordSpherical.phig)); /* sin  (geocentric latitude) */

        if(nMax <= 16 || (1 - fabs(sin_phi)) < 1.0e-10) {         /* If nMax is less tha 16 or at the poles */
            FLAG = PcupLow(*legendreFunction, sin_phi, nMax);
        } else {
            FLAG = PcupHigh(*legendreFunction, sin_phi, nMax);
        }

        return FLAG;
    }

    /** Check if the latitude is equal to -90 or 90. If it is, offset it by 1e-5 to avoid division by zero. This is not
     * currently used in the Geomagnetic main function. This may be used to avoid calling SummationSpecial. The function
     * updates the input data structure.
     *
     * INPUT   CoordGeodetic Pointer to the  data  structure with the following elements
     *                double lambda; (longitude)
     *                double phi; ( geodetic latitude)
     *                double HeightAboveEllipsoid; (height above the ellipsoid (HaE) )
     *                double HeightAboveGeoid;(height above the Geoid )
     * OUTPUT  CoordGeodetic  Pointer to the  data  structure with the following elements updates
     *                double phi; ( geodetic latitude)
     * CALLS : none
     */
    int CheckGeographicPole(CoordGeodetic *coordGeodetic) {
        coordGeodetic->phi = std::clamp(coordGeodetic->phi, -90.0 + GEO_POLE_TOLERANCE, 90.0 - GEO_POLE_TOLERANCE);
        return true;
    }

    /** Computes Spherical variables. Variables computed are (a/r)^(n+2), cos_m(lamda) and sin_m(lambda) for spherical
     * harmonic summations. (Equations 10-12 in the WMM Technical Report)
     * INPUT   Ellip
     *         CoordSpherical
     *         nMax   integer 	 ( Maxumum degree of spherical harmonic secular model)
     * OUTPUT  SphVariables  Pointer to the   data structure with the following elements
     *         double RelativeRadiusPower[MAX_MODEL_DEGREES+1];   [earth_reference_radius_km  sph. radius ]^n
     *         double cos_mlambda[MAX_MODEL_DEGREES+1]; cp(m)  - cosine of (mspherical coord. longitude)
     *         double sin_mlambda[MAX_MODEL_DEGREES+1];  sp(m)  - sine of (mspherical coord. longitude)
     * CALLS : none
     */
    int ComputeSphericalHarmonicVariables(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const int nMax,
                                          SphericalHarmonicVariables *sphVariables) {
        const double cos_lambda = cos(Deg2Rad(coordSpherical.lambda));
        const double sin_lambda = sin(Deg2Rad(coordSpherical.lambda));
        /* for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2) for n  1..nMax-1 (this is much
         * faster than calling pow MAX_N+1 times).      */
        sphVariables->RelativeRadiusPower.at(0) = (ellip.re / coordSpherical.r) * (ellip.re / coordSpherical.r);
        for(int n = 1; n <= nMax; n++) {
            sphVariables->RelativeRadiusPower.at(n) =
                sphVariables->RelativeRadiusPower.at(n - 1) * (ellip.re / coordSpherical.r);
        }

        /*
         Compute cos(m*lambda), sin(m*lambda) for m = 0 ... nMax
               cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
               sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
         */

        sphVariables->cos_mlambda.at(0) = 1.0;  // The size if cos_mlambda and sin_mlambda is nMax+1
        sphVariables->sin_mlambda.at(0) = 0.0;
        if(nMax + 1 >= 2) {
            sphVariables->cos_mlambda.at(1) = cos_lambda;
            sphVariables->sin_mlambda.at(1) = sin_lambda;
        }

        if(nMax + 1 >= 3) {
            for(int m = 2; m <= nMax; m++) {
                sphVariables->cos_mlambda.at(m) =
                    sphVariables->cos_mlambda.at(m - 1) * cos_lambda - sphVariables->sin_mlambda.at(m - 1) * sin_lambda;
                sphVariables->sin_mlambda.at(m) =
                    sphVariables->cos_mlambda.at(m - 1) * sin_lambda + sphVariables->sin_mlambda.at(m - 1) * cos_lambda;
            }
        }
        return true;
    }

    void GradY(Ellipsoid ellip, CoordSpherical coordSpherical, CoordGeodetic coordGeodetic,
               MagneticModel *timedMagneticModel, GeoMagneticElements geoMagneticElements,
               GeoMagneticElements *GradYElements) {
        int NumTerms;
        MagneticResults GradYResultsSph, GradYResultsGeo;

        NumTerms = ((timedMagneticModel->nMax + 1) * (timedMagneticModel->nMax + 2) / 2);
        LegendreFunction legendreFunction(NumTerms);      /* For storing the ALF functions */
        SphericalHarmonicVariables SphVariables(timedMagneticModel->nMax);
        ComputeSphericalHarmonicVariables(ellip, coordSpherical, timedMagneticModel->nMax,
                                          &SphVariables); /* Compute Spherical Harmonic variables  */
        AssociatedLegendreFunction(coordSpherical, timedMagneticModel->nMax, &legendreFunction); /* Compute ALF  */
        GradYSummation(&legendreFunction, timedMagneticModel, SphVariables, coordSpherical,
                       &GradYResultsSph);       /* Accumulate the spherical harmonic coefficients*/
        RotateMagneticVector(coordSpherical, coordGeodetic, GradYResultsSph,
                             &GradYResultsGeo); /* Map the computed Magnetic fields to Geodetic coordinates  */
        CalculateGradientElements(
            GradYResultsGeo, geoMagneticElements,
            GradYElements);                     /* Calculate the Geomagnetic elements, Equation 18 , WMM Technical report */
    }

    void GradYSummation(LegendreFunction *LegendreFunction, MagneticModel *magneticModel,
                        SphericalHarmonicVariables SphVariables, CoordSpherical coordSpherical, MagneticResults *GradY) {
        int m, n, index;
        double cos_phi;
        GradY->Bz = 0.0;
        GradY->By = 0.0;
        GradY->Bx = 0.0;
        for(n = 1; n <= magneticModel->nMax; n++) {
            for(m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);

                GradY->Bz -= SphVariables.RelativeRadiusPower[n] *
                             (-1 * magneticModel->main_Field_Coeff_G[index] * SphVariables.sin_mlambda[m] +
                              magneticModel->main_Field_Coeff_H[index] * SphVariables.cos_mlambda[m]) *
                             (double)(n + 1) * (double)(m)*LegendreFunction->Pcup[index] * (1 / coordSpherical.r);
                GradY->By += SphVariables.RelativeRadiusPower[n] *
                             (magneticModel->main_Field_Coeff_G[index] * SphVariables.cos_mlambda[m] +
                              magneticModel->main_Field_Coeff_H[index] * SphVariables.sin_mlambda[m]) *
                             (double)(m * m) * LegendreFunction->Pcup[index] * (1 / coordSpherical.r);
                GradY->Bx -= SphVariables.RelativeRadiusPower[n] *
                             (-1 * magneticModel->main_Field_Coeff_G[index] * SphVariables.sin_mlambda[m] +
                              magneticModel->main_Field_Coeff_H[index] * SphVariables.cos_mlambda[m]) *
                             (double)(m)*LegendreFunction->dPcup[index] * (1 / coordSpherical.r);
            }
        }

        cos_phi = cos(Deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            GradY->By = GradY->By / (cos_phi * cos_phi);
            GradY->Bx = GradY->Bx / (cos_phi);
            GradY->Bz = GradY->Bz / (cos_phi);
        } else
        /* Special calculation for component - By - at Geographic poles.
         * If the user wants to avoid using this function,  please make sure that
         * the latitude is not exactly +/-90. An option is to make use the function
         * CheckGeographicPoles.
         */
        {
            /* SummationSpecial(MagneticModel, SphVariables, coordSpherical, GradY); */
        }
    }

    /** This function evaluates all of the Schmidt-semi normalized associated Legendre functions up to degree nMax. The
     * functions are initially scaled by 10^280 sin^m in order to minimize the effects of underflow at large m near the poles
     * (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). Note that this function performs the same operation as
     * PcupLow. However this function also can be used for high degree (large nMax) models.
     *
     * Calling Parameters:
     * INPUT nMax:	 Maximum spherical harmonic degree to compute.
     *       x:		 cos(colatitude) or sin(latitude).
     *
     * OUTPUT Pcup:	A vector of all associated Legendgre polynomials evaluated at x up to nMax. The lenght must by
     * greater or equal to (nMax+1)*(nMax+2)/2.
     *        dPcup:   Derivative of Pcup(x) with respect to latitude
     *
     * CALLS : none
     *
     * Notes:
     *  Adopted from the FORTRAN code written by Mark Wieczorek September 25, 2005.
     *
     *  Manoj Nair, Nov, 2009 Manoj.C.Nair@Noaa.Gov
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
    int PcupHigh(LegendreFunction &legendreFunction, const double x, const int nMax) {
        auto &Pcup  = legendreFunction.Pcup;
        auto &dPcup = legendreFunction.dPcup;

        double plm;
        int m, n;

        const int NumTerms = ((nMax + 1) * (nMax + 2) / 2);

        const double z = sqrt((1.0 - x) * (1.0 + x));

        if(z == 0) {
            return false;
        }

        if(fabs(x) == 1.0) {
            printf("Error in PcupHigh: derivative cannot be calculated at poles\n");
            return false;
        }

        std::vector<double> f1(NumTerms + 1);
        std::vector<double> PreSqr(NumTerms + 1);
        std::vector<double> f2(NumTerms + 1);

        constexpr double scalef = 1.0e-280;

        for(n = 0; n <= 2 * nMax + 1; ++n) {
            PreSqr.at(n) = std::sqrt(n);
        }

        int k = 2;

        for(n = 2; n <= nMax; n++) {
            k     = k + 1;
            f1[k] = (2. * n - 1) / static_cast<double>(n);
            f2[k] = static_cast<double>(n - 1) / static_cast<double>(n);
            for(m = 1; m <= n - 2; m++) {
                k     = k + 1;
                f1[k] = (2. * n - 1) / PreSqr[n + m] / PreSqr[n - m];
                f2[k] = PreSqr[n - m - 1] * PreSqr[n + m - 1] / PreSqr[n + m] / PreSqr[n - m];
            }
            k = k + 2;
        }

        /* z = sin (geocentric latitude) */
        double pm2  = 1.0;
        Pcup.at(0)  = 1.0;
        dPcup.at(0) = 0.0;
        if(nMax == 0) {
            return false;
        }
        double pm1  = x;
        Pcup.at(1)  = pm1;
        dPcup.at(1) = z;
        k           = 1;

        for(n = 2; n <= nMax; n++) {
            k           = k + n;
            plm         = f1.at(k) * x * pm1 - f2.at(k) * pm2;
            Pcup.at(k)  = plm;
            dPcup.at(k) = static_cast<double>(n) * (pm1 - x * plm) / z;
            pm2         = pm1;
            pm1         = plm;
        }

        double pmm      = PreSqr[2] * scalef;
        double rescalem = 1.0 / scalef;
        int kstart      = 0;

        for(m = 1; m <= nMax - 1; ++m) {
            rescalem = rescalem * z;

            /* Calculate Pcup(m,m)*/
            kstart           = kstart + m + 1;
            pmm              = pmm * PreSqr.at(2 * m + 1) / PreSqr.at(2 * m);
            Pcup.at(kstart)  = pmm * rescalem / PreSqr.at(2 * m + 1);
            dPcup.at(kstart) = -(static_cast<double>(m) * x * Pcup.at(kstart) / z);
            pm2              = pmm / PreSqr.at(2 * m + 1);

            /* Calculate Pcup(m+1,m)*/
            k           = kstart + m + 1;
            pm1         = x * PreSqr.at(2 * m + 1) * pm2;
            Pcup.at(k)  = pm1 * rescalem;
            dPcup.at(k) = ((pm2 * rescalem) * PreSqr.at(2 * m + 1) - x * static_cast<double>(m + 1) * Pcup.at(k)) / z;

            /* Calculate Pcup(n,m)*/
            for(n = m + 2; n <= nMax; ++n) {
                k          = k + n;
                plm        = x * f1.at(k) * pm1 - f2.at(k) * pm2;
                Pcup.at(k) = plm * rescalem;
                dPcup.at(k) =
                    (PreSqr.at(n + m) * PreSqr.at(n - m) * (pm1 * rescalem) - static_cast<double>(n) * x * Pcup.at(k)) / z;
                pm2 = pm1;
                pm1 = plm;
            }
        }

        /* Calculate Pcup(nMax,nMax)*/
        rescalem         = rescalem * z;
        kstart           = kstart + m + 1;
        pmm              = pmm / PreSqr.at(2 * nMax);
        Pcup.at(kstart)  = pmm * rescalem;
        dPcup.at(kstart) = -static_cast<double>(nMax) * x * Pcup.at(kstart) / z;

        return true;
    }

    /**  This function evaluates all of the Schmidt-semi normalized associated Legendre
     *   functions up to degree nMax.
     *
     *   Calling Parameters:
     *   INPUT nMax:	 Maximum spherical harmonic degree to compute.
     *         x:		 cos(colatitude) or sin(latitude).
     *
     *   OUTPUT Pcup:	 A vector of all associated Legendgre polynomials evaluated at x up to nMax.
     *          dPcup: Derivative of Pcup(x) with respect to latitude
     *
     *   Notes: Overflow may occur if nMax > 20 , especially for high-latitudes. Use PcupHigh for large nMax.
     *   Written by Manoj Nair, June, 2009 . Manoj.C.Nair@Noaa.Gov.
     *
     *   Note: In geomagnetism, the derivatives of ALF are usually found with respect to the colatitudes. Here the
     *   derivatives are found with respect to the latitude. The difference is a sign reversal for the derivative of
     *   the Associated Legendre Functions.
     */
    int PcupLow(LegendreFunction &legendreFunction, const double x, const int nMax) {
        auto &Pcup  = legendreFunction.Pcup;
        auto &dPcup = legendreFunction.dPcup;

        Pcup.at(0)  = 1.0;
        dPcup.at(0) = 0.0;
        /* sin (geocentric latitude) - sin_phi */
        const double z = sqrt((1.0 - x) * (1.0 + x));

        const int NumTerms = ((nMax + 1) * (nMax + 2) / 2);
        std::vector<double> schmidtQuasiNorm(NumTerms + 1);

        /*	 First,	Compute the Gauss-normalized associated Legendre functions */
        for(int n = 1; n <= nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = (n * (n + 1) / 2 + m);
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
                    const int index2 = (n - 1) * n / 2 + m;
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

        /* Compute the ration between the the Schmidt quasi-normalized associated Legendre
         * functions and the Gauss-normalized version. */
        schmidtQuasiNorm.at(0) = 1.0;
        for(int n = 1; n <= nMax; n++) {
            int index  = (n * (n + 1) / 2);
            int index1 = (n - 1) * n / 2;

            /* for m = 0 */
            schmidtQuasiNorm.at(index) = schmidtQuasiNorm.at(index1) * (2. * n - 1) / static_cast<double>(n);

            for(int m = 1; m <= n; m++) {
                index  = (n * (n + 1) / 2 + m);
                index1 = (n * (n + 1) / 2 + m - 1);
                schmidtQuasiNorm.at(index) =
                    schmidtQuasiNorm.at(index1) * sqrt(((n - m + 1) * (m == 1 ? 2. : 1.)) / static_cast<double>(n + m));
            }
        }

        /* Converts the  Gauss-normalized associated Legendre functions to the Schmidt quasi-normalized version using
         * pre-computed relation stored in the variable schmidtQuasiNorm */
        for(int n = 1; n <= nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = (n * (n + 1) / 2 + m);
                Pcup.at(index)  = Pcup.at(index) * schmidtQuasiNorm.at(index);
                dPcup.at(index) = -dPcup.at(index) * schmidtQuasiNorm.at(index);
                /* The sign is changed since the new WMM routines use derivative with respect to latitude insted of
                 * co-latitude */
            }
        }

        return true;
    }

    /** This Function sums the secular variation coefficients to get the secular variation of the Magnetic vector.
     * INPUT :  LegendreFunction
     *          MagneticModel
     *          SphVariables
     *          CoordSpherical
     * OUTPUT : MagneticResults
     *
     * CALLS : SecVarSummationSpecial
     */
    int SecVarSummation(const LegendreFunction &legendreFunction, MagneticModel *magneticModel,
                        const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical,
                        MagneticResults *magneticResults) {
        magneticModel->secularVariationUsed = true;
        magneticResults->Bz                 = 0.0;
        magneticResults->By                 = 0.0;
        magneticResults->Bx                 = 0.0;
        for(int n = 1; n <= magneticModel->nMaxSecVar; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = (n * (n + 1) / 2 + m);

                /*		    nMax  	(n+2) 	  n     m            m           m
                        Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                        n=1      	      m=0   n            n           n  */
                /*  Derivative with respect to radius.*/
                magneticResults->Bz -= sphVariables.RelativeRadiusPower.at(n) *
                                       (magneticModel->secular_Var_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                                        magneticModel->secular_Var_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                                       static_cast<double>(n + 1) * legendreFunction.Pcup.at(index);

                /*		  1 nMax  (n+2)    n     m            m           m
                        By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1             m=0   n            n           n  */
                /* Derivative with respect to longitude, divided by radius. */
                magneticResults->By += sphVariables.RelativeRadiusPower.at(n) *
                                       (magneticModel->secular_Var_Coeff_G.at(index) * sphVariables.sin_mlambda.at(m) -
                                        magneticModel->secular_Var_Coeff_H.at(index) * sphVariables.cos_mlambda.at(m)) *
                                       static_cast<double>(m) * legendreFunction.Pcup.at(index);

                /*		   nMax  (n+2) n     m            m           m
                        Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1         m=0   n            n           n  */
                /* Derivative with respect to latitude, divided by radius. */
                magneticResults->Bx -= sphVariables.RelativeRadiusPower.at(n) *
                                       (magneticModel->secular_Var_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                                        magneticModel->secular_Var_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                                       legendreFunction.dPcup.at(index);
            }
        }
        double cos_phi = cos(Deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            magneticResults->By = magneticResults->By / cos_phi;
        } else
        /* Special calculation for component By at Geographic poles */
        {
            SecVarSummationSpecial(*magneticModel, sphVariables, coordSpherical, magneticResults);
        }
        return true;
    }

    /** Special calculation for the secular variation summation at the poles.
     * INPUT:  MagneticModel
     *         SphVariables
     *         CoordSpherical
     * OUTPUT: MagneticResults
     * CALLS : none
     */
    int SecVarSummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                               const CoordSpherical &coordSpherical, MagneticResults *magneticResults) {
        std::vector<double> PcupS(magneticModel.nMaxSecVar + 1);

        PcupS[0]                 = 1;
        double schmidtQuasiNorm1 = 1.0;

        magneticResults->By  = 0.0;
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
            magneticResults->By += sphVariables.RelativeRadiusPower.at(n) *
                                   (magneticModel.secular_Var_Coeff_G.at(index) * sphVariables.sin_mlambda.at(1) -
                                    magneticModel.secular_Var_Coeff_H.at(index) * sphVariables.cos_mlambda.at(1)) *
                                   PcupS.at(n) * schmidtQuasiNorm3;
        }
        return true;
    }

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
     *
     * CALLS : SummationSpecial
     *
     *Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
     */
    int Summation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                  const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical,
                  MagneticResults *magneticResults) {
        magneticResults->Bz = 0.0;
        magneticResults->By = 0.0;
        magneticResults->Bx = 0.0;

        for(int n = 1; n <= magneticModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                const int index = (n * (n + 1) / 2 + m);

                /*		    nMax  	(n+2) 	  n     m            m           m
                        Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                        n=1      	      m=0   n            n           n  */
                /* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
                magneticResults->Bz -= sphVariables.RelativeRadiusPower.at(n) *
                                       (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                                        magneticModel.main_Field_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                                       static_cast<double>(n + 1) * legendreFunction.Pcup.at(index);

                /*		  1 nMax  (n+2)    n     m            m           m
                        By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1             m=0   n            n           n  */
                /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
                magneticResults->By += sphVariables.RelativeRadiusPower.at(n) *
                                       (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.sin_mlambda.at(m) -
                                        magneticModel.main_Field_Coeff_H.at(index) * sphVariables.cos_mlambda.at(m)) *
                                       static_cast<double>(m) * legendreFunction.Pcup.at(index);
                /*		   nMax  (n+2) n     m            m           m
                        Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                                   n=1         m=0   n            n           n  */
                /* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */

                magneticResults->Bx -= sphVariables.RelativeRadiusPower.at(n) *
                                       (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.cos_mlambda.at(m) +
                                        magneticModel.main_Field_Coeff_H.at(index) * sphVariables.sin_mlambda.at(m)) *
                                       legendreFunction.dPcup.at(index);
            }
        }

        const double cos_phi = cos(Deg2Rad(coordSpherical.phig));
        if(fabs(cos_phi) > 1.0e-10) {
            magneticResults->By = magneticResults->By / cos_phi;
        } else
        /* Special calculation for component - By - at Geographic poles.
         * If the user wants to avoid using this function,  please make sure that
         * the latitude is not exactly +/-90. An option is to make use the function
         * CheckGeographicPoles.
         */
        {
            SummationSpecial(magneticModel, sphVariables, coordSpherical, magneticResults);
        }
        return true;
    }

    /** Special calculation for the component By at Geographic poles.
     * Manoj Nair, June, 2009 manoj.c.nair@noaa.gov
     * INPUT:  MagneticModel
     *         SphVariables
     *         CoordSpherical
     * OUTPUT: MagneticResults
     * CALLS : none
     * See Section 1.4, "SINGULARITIES AT THE GEOGRAPHIC POLES", WMM Technical report
     */
    int SummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                         const CoordSpherical &coordSpherical, MagneticResults *magneticResults) {
        std::vector<double> PcupS(magneticModel.nMax + 1);

        PcupS.at(0)              = 1;
        double schmidtQuasiNorm1 = 1.0;

        magneticResults->By  = 0.0;
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

            magneticResults->By += sphVariables.RelativeRadiusPower.at(n) *
                                   (magneticModel.main_Field_Coeff_G.at(index) * sphVariables.sin_mlambda.at(1) -
                                    magneticModel.main_Field_Coeff_H.at(index) * sphVariables.cos_mlambda.at(1)) *
                                   PcupS.at(n) * schmidtQuasiNorm3;
        }

        return true;
    }

    /**
     * Time change the Model coefficients from the base year of the model using secular variation coefficients.
     * Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
     * Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.
     * If the address of "TimedMagneticModel" is the same as the address of "MagneticModel", then this procedure overwrites
     * the given item "MagneticModel".
     * INPUT:  UserDate
     *         MagneticModel
     * OUTPUT: TimedMagneticModel
     * CALLS : none
     */
    int TimelyModifyMagneticModel(const Date &userDate, const MagneticModel &magModel, MagneticModel *timedMagModel) {
        timedMagModel->editionDate = magModel.editionDate;
        timedMagModel->epoch       = magModel.epoch;
        timedMagModel->nMax        = magModel.nMax;
        timedMagModel->nMaxSecVar  = magModel.nMaxSecVar;
        const int a                = timedMagModel->nMaxSecVar;
        const int b                = (a * (a + 1) / 2 + a);
        timedMagModel->modelName   = magModel.modelName;

        timedMagModel->main_Field_Coeff_G.resize(b + 1);
        timedMagModel->main_Field_Coeff_H.resize(b + 1);
        timedMagModel->secular_Var_Coeff_G.resize(b + 1);
        timedMagModel->secular_Var_Coeff_H.resize(b + 1);

        for(int n = 1; n <= magModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                if(const int index = (n * (n + 1) / 2 + m); index <= b) {
                    timedMagModel->main_Field_Coeff_H.at(index) =
                        magModel.main_Field_Coeff_H.at(index) +
                        (userDate.DecimalYear - magModel.epoch) * magModel.secular_Var_Coeff_H.at(index);
                    timedMagModel->main_Field_Coeff_G.at(index) =
                        magModel.main_Field_Coeff_G.at(index) +
                        (userDate.DecimalYear - magModel.epoch) * magModel.secular_Var_Coeff_G.at(index);
                    /* We need a copy of the secular var coef to calculate secular change */
                    timedMagModel->secular_Var_Coeff_H.at(index) = magModel.secular_Var_Coeff_H.at(index);
                    timedMagModel->secular_Var_Coeff_G.at(index) = magModel.secular_Var_Coeff_G.at(index);
                } else {
                    timedMagModel->main_Field_Coeff_H.at(index) = magModel.main_Field_Coeff_H.at(index);
                    timedMagModel->main_Field_Coeff_G.at(index) = magModel.main_Field_Coeff_G.at(index);
                }
            }
        }
        return true;
    }

    /*End of Spherical Harmonic Functions*/


    /******************************************************************************
     *************************************Geoid************************************
     * This grouping consists of functions that make calculations to adjust
     * ellipsoid height to height above the geoid (Height above MSL).
     ******************************************************************************
     ******************************************************************************/


    int ConvertGeoidToEllipsoidHeight(CoordGeodetic *coordGeodetic, Geoid *Geoid)

    /*
 * The function Convert_Geoid_To_Ellipsoid_Height converts the specified WGS84
 * Geoid height at the specified geodetic coordinates to the equivalent
 * ellipsoid height, using the EGM96 gravity model.
 *
 *   coordGeodetic->phi        : Geodetic latitude in degress           (input)
 *    coordGeodetic->lambda     : Geodetic longitude in degrees          (input)
 *    coordGeodetic->HeightAboveEllipsoid	     : Ellipsoid height, in kilometers         (output)
 *    coordGeodetic->HeightAboveGeoid: Geoid height, in kilometers           (input)
 *
        CALLS : GetGeoidHeight (

 */
    {
        double DeltaHeight;
        int Error_Code;
        double lat, lon;

        if(Geoid->UseGeoid == 1) { /* Geoid correction required */
                                   /* To ensure that latitude is less than 90 call EquivalentLatLon() */
            EquivalentLatLon(coordGeodetic->phi, coordGeodetic->lambda, &lat, &lon);
            Error_Code                          = GetGeoidHeight(lat, lon, &DeltaHeight, Geoid);
            coordGeodetic->HeightAboveEllipsoid = coordGeodetic->HeightAboveGeoid + DeltaHeight / 1000; /*  Input and output
            should be kilometers, However GetGeoidHeight returns Geoid height in meters - Hence division by 1000 */
        } else /* Geoid correction not required, copy the MSL height to Ellipsoid height */
        {
            coordGeodetic->HeightAboveEllipsoid = coordGeodetic->HeightAboveGeoid;
            Error_Code                          = true;
        }
        return (Error_Code);
    } /* ConvertGeoidToEllipsoidHeight*/

    int GetGeoidHeight(double Latitude, double Longitude, double *DeltaHeight, Geoid *Geoid)
    /*
 * The  function GetGeoidHeight returns the height of the
 * EGM96 geiod above or below the WGS84 ellipsoid,
 * at the specified geodetic coordinates,
 * using a grid of height adjustments from the EGM96 gravity model.
 *
 *    Latitude            : Geodetic latitude in radians           (input)
 *    Longitude           : Geodetic longitude in radians          (input)
 *    DeltaHeight         : Height Adjustment, in meters.          (output)
 *    Geoid				  : Geoid with Geoid grid		   (input)
        CALLS : none
 */
    {
        long Index;
        double DeltaX, DeltaY;
        double ElevationSE, ElevationSW, ElevationNE, ElevationNW;
        double OffsetX, OffsetY;
        double PostX, PostY;
        double UpperY, LowerY;
        int Error_Code = 0;

        if(!Geoid->Geoid_Initialized) {
            PrintError(5);
            return (false);
        }
        if((Latitude < -90) || (Latitude > 90)) {     /* Latitude out of range */
            Error_Code |= 1;
        }
        if((Longitude < -180) || (Longitude > 360)) { /* Longitude out of range */
            Error_Code |= 1;
        }

        if(!Error_Code) { /* no errors */
            /*  Compute X and Y Offsets into Geoid Height Array:                          */

            if(Longitude < 0.0) {
                OffsetX = (Longitude + 360.0) * Geoid->ScaleFactor;
            } else {
                OffsetX = Longitude * Geoid->ScaleFactor;
            }
            OffsetY = (90.0 - Latitude) * Geoid->ScaleFactor;

            /*  Find Four Nearest Geoid Height Cells for specified Latitude, Longitude;   */
            /*  Assumes that (0,0) of Geoid Height Array is at Northwest corner:          */

            PostX = floor(OffsetX);
            if((PostX + 1) == Geoid->NumbGeoidCols) {
                PostX--;
            }
            PostY = floor(OffsetY);
            if((PostY + 1) == Geoid->NumbGeoidRows) {
                PostY--;
            }

            Index       = (long)(PostY * Geoid->NumbGeoidCols + PostX);
            ElevationNW = (double)Geoid->GeoidHeightBuffer[Index];
            ElevationNE = (double)Geoid->GeoidHeightBuffer[Index + 1];

            Index       = (long)((PostY + 1) * Geoid->NumbGeoidCols + PostX);
            ElevationSW = (double)Geoid->GeoidHeightBuffer[Index];
            ElevationSE = (double)Geoid->GeoidHeightBuffer[Index + 1];

            /*  Perform Bi-Linear Interpolation to compute Height above Ellipsoid:        */

            DeltaX = OffsetX - PostX;
            DeltaY = OffsetY - PostY;

            UpperY = ElevationNW + DeltaX * (ElevationNE - ElevationNW);
            LowerY = ElevationSW + DeltaX * (ElevationSE - ElevationSW);

            *DeltaHeight = UpperY + DeltaY * (LowerY - UpperY);
        } else {
            PrintError(17);
            return (false);
        }
        return true;
    } /*GetGeoidHeight*/

    void EquivalentLatLon(double lat, double lon, double *repairedLat, double *repairedLon)
    /*This function takes a latitude and longitude that are ordinarily out of range
 and gives in range values that are equivalent on the Earth's surface.  This is
 required to get correct values for the geoid function.*/
    {
        double colat;
        colat        = 90 - lat;
        *repairedLon = lon;
        if(colat < 0) {
            colat = -colat;
        }
        while(colat > 360) {
            colat -= 360;
        }
        if(colat > 180) {
            colat -= 180;
            *repairedLon = *repairedLon + 180;
        }
        *repairedLat = 90 - colat;
        if(*repairedLon > 360) {
            *repairedLon -= 360;
        }
        if(*repairedLon < -180) {
            *repairedLon += 360;
        }
    }

    /*End of Geoid Functions*/


    /*New Error Functions*/

    void WMMErrorCalc(double H, GeoMagneticElements *Uncertainty) {
        double decl_variable, decl_constant;
        Uncertainty->F    = WMM_UNCERTAINTY_F;
        Uncertainty->H    = WMM_UNCERTAINTY_H;
        Uncertainty->X    = WMM_UNCERTAINTY_X;
        Uncertainty->Z    = WMM_UNCERTAINTY_Z;
        Uncertainty->Incl = WMM_UNCERTAINTY_I;
        Uncertainty->Y    = WMM_UNCERTAINTY_Y;
        decl_variable     = (WMM_UNCERTAINTY_D_COEF / H);
        decl_constant     = (WMM_UNCERTAINTY_D_OFFSET);
        Uncertainty->Decl = sqrt(decl_constant * decl_constant + decl_variable * decl_variable);
        if(Uncertainty->Decl > 180) {
            Uncertainty->Decl = 180;
        }
    }

    void WMMHRErrorCalc(double H, GeoMagneticElements *Uncertainty) {
        double decl_variable, decl_constant;
        Uncertainty->F    = WMMHR_UNCERTAINTY_F;
        Uncertainty->H    = WMMHR_UNCERTAINTY_H;
        Uncertainty->X    = WMMHR_UNCERTAINTY_X;
        Uncertainty->Z    = WMMHR_UNCERTAINTY_Z;
        Uncertainty->Incl = WMMHR_UNCERTAINTY_I;
        Uncertainty->Y    = WMMHR_UNCERTAINTY_Y;
        decl_variable     = (WMMHR_UNCERTAINTY_D_COEF / H);
        decl_constant     = (WMMHR_UNCERTAINTY_D_OFFSET);
        Uncertainty->Decl = sqrt(decl_constant * decl_constant + decl_variable * decl_variable);
        if(Uncertainty->Decl > 180) {
            Uncertainty->Decl = 180;
        }
    }

    void PrintUserDataWithUncertainty(GeoMagneticElements geomagElements, GeoMagneticElements Errors,
                                      CoordGeodetic spaceInput, Date TimeInput, MagneticModel *magneticModel, Geoid *Geoid) {
        int dms_size = 150;
        std::string DeclString;
        std::string InclString;
        std::string GVString;

        DegreeToDMSstring(geomagElements.Incl, 2, InclString);
        if(geomagElements.H < 6000 && geomagElements.H > 2000) {
            Warnings(1, geomagElements.H, *magneticModel);
        }
        if(geomagElements.H < 2000) {
            Warnings(2, geomagElements.H, *magneticModel);
        }
        if(magneticModel->secularVariationUsed == true) {
            DegreeToDMSstring(geomagElements.Decl, 2, DeclString);
            printf("\n Results For \n\n");
            if(spaceInput.phi < 0) {
                printf("Latitude	%.2fS\n", -spaceInput.phi);
            } else {
                printf("Latitude	%.2fN\n", spaceInput.phi);
            }
            if(spaceInput.lambda < 0) {
                printf("Longitude	%.2fW\n", -spaceInput.lambda);
            } else {
                printf("Longitude	%.2fE\n", spaceInput.lambda);
            }
            if(Geoid->UseGeoid == 1) {
                printf("Altitude:	%.2f Kilometers above mean sea level\n", spaceInput.HeightAboveGeoid);
            } else {
                printf("Altitude:	%.2f Kilometers above the WGS-84 ellipsoid\n", spaceInput.HeightAboveEllipsoid);
            }
            printf("Date:		%.1f\n", TimeInput.DecimalYear);
            printf("\n		Main Field\t\t\tSecular Change\n");
            printf("F	=	%9.1f +/- %5.1f nT\t\t Fdot = %5.1f\tnT/yr\n", geomagElements.F, Errors.F, geomagElements.Fdot);
            printf("H	=	%9.1f +/- %5.1f nT\t\t Hdot = %5.1f\tnT/yr\n", geomagElements.H, Errors.H, geomagElements.Hdot);
            printf("X	=	%9.1f +/- %5.1f nT\t\t Xdot = %5.1f\tnT/yr\n", geomagElements.X, Errors.X, geomagElements.Xdot);
            printf("Y	=	%9.1f +/- %5.1f nT\t\t Ydot = %5.1f\tnT/yr\n", geomagElements.Y, Errors.Y, geomagElements.Ydot);
            printf("Z	=	%9.1f +/- %5.1f nT\t\t Zdot = %5.1f\tnT/yr\n", geomagElements.Z, Errors.Z, geomagElements.Zdot);
            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST) +/-%3.0f Min Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * Errors.Decl,
                       60 * geomagElements.Decldot);
            } else {
                printf("Decl	=%20s  (EAST) +/-%3.0f Min Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * Errors.Decl,
                       60 * geomagElements.Decldot);
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)   +/-%3.0f Min Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * Errors.Incl,
                       60 * geomagElements.Incldot);
            } else {
                printf("Incl	=%20s  (DOWN) +/-%3.0f Min Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * Errors.Incl,
                       60 * geomagElements.Incldot);
            }
        } else {
            DegreeToDMSstring(geomagElements.Decl, 2, DeclString);
            printf("\n Results For \n\n");
            if(spaceInput.phi < 0) {
                printf("Latitude	%.2fS\n", -spaceInput.phi);
            } else {
                printf("Latitude	%.2fN\n", spaceInput.phi);
            }
            if(spaceInput.lambda < 0) {
                printf("Longitude	%.2fW\n", -spaceInput.lambda);
            } else {
                printf("Longitude	%.2fE\n", spaceInput.lambda);
            }
            if(Geoid->UseGeoid == 1) {
                printf("Altitude:	%.2f Kilometers above MSL\n", spaceInput.HeightAboveGeoid);
            } else {
                printf("Altitude:	%.2f Kilometers above WGS-84 Ellipsoid\n", spaceInput.HeightAboveEllipsoid);
            }
            printf("Date:		%.1f\n", TimeInput.DecimalYear);
            printf("\n	Main Field\n");
            printf("F	=	%-9.1f +/-%5.1f nT\n", geomagElements.F, Errors.F);
            printf("H	=	%-9.1f +/-%5.1f nT\n", geomagElements.H, Errors.H);
            printf("X	=	%-9.1f +/-%5.1f nT\n", geomagElements.X, Errors.X);
            printf("Y	=	%-9.1f +/-%5.1f nT\n", geomagElements.Y, Errors.Y);
            printf("Z	=	%-9.1f +/-%5.1f nT\n", geomagElements.Z, Errors.Z);
            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST)+/-%4f\n", DeclString.c_str(), 60 * Errors.Decl);
            } else {
                printf("Decl	=%20s  (EAST)+/-%4f\n", DeclString.c_str(), 60 * Errors.Decl);
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)+/-%4f\n", InclString.c_str(), 60 * Errors.Incl);
            } else {
                printf("Incl	=%20s  (DOWN)+/-%4f\n", InclString.c_str(), 60 * Errors.Incl);
            }
        }

        DegreeToDMSstring(geomagElements.GV, 2, GVString);

        /* Print Grid Variation */

        if(spaceInput.phi < -55) {
            printf("\n\n Grid variation (SOUTH) =%20s\n", GVString.c_str());
        } else if(spaceInput.phi > 55) {
            printf("\n\n Grid variation (NORTH) =%20s \n", GVString.c_str());
        }
    } /*PrintUserDataWithUncertainty*/

    // Parse the date string format as mm/dd/yyyy
    double dateStr_to_decYear(const std::string_view edit_date) {
        int day{}, month{}, year{};
        double extra_day{0};
        static std::array days = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        int total_days         = 0;

        if(sscanf(edit_date.data(), "%d/%d/%d", &month, &day, &year) != 3) {
            printf("Failed to parse the date string. Please use the format mm/dd/yyyy\n");
            return -1;
        }

        if(month > 12 || month < 1) {
            printf("Month out of range\n");
            return -1;
        }

        if(day > days.at(month) || day < 1) {
            printf("Day out of range\n");
            return -1;
        }

        if((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
            extra_day = 1;
        }

        const double total_year_days = 365.0 + extra_day;
        days.at(2) += static_cast<int>(extra_day);

        for(int i = 0; i < month; i++) {
            total_days += days.at(i);
        }
        total_days += day;

        return year + (total_days - 1) / total_year_days;
    }

    // Parse the date string format as mm/dd/yyyy
    double date_to_decYear(const int year, const int month, const int day) {
        double extra_day{0};
        static std::array days = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        int total_days         = 0;

        if((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
            extra_day = 1;
        }

        const double total_year_days = 365.0 + extra_day;
        days.at(2) += static_cast<int>(extra_day);

        for(int i = 0; i < month; i++) {
            total_days += days.at(i);
        }
        total_days += day;

        return year + (total_days - 1) / total_year_days;
    }

}  // namespace wmm