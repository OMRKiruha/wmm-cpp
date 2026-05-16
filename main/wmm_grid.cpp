#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "GeomagInterativeLib.h"
#include "GeomagnetismHeader.h"

#include "wmm_warn.h"

/*
WMMHR grid program.

The Geomagnetism Library is used to make a command prompt program. The program prompts
the user to enter a location, performs the computations and prints the results to the
standard output. The program expects the files GeomagnetismLibrary.c, GeomagnetismHeader.h,
EWMM.COF and EGM9615.h to be in the same directory.

Manoj.C.Nair@Noaa.Gov
April 21, 2011

liyin.young@noaa.gov
Updated April, 2023
 */

int Grid(wmm::CoordGeodetic minimum, wmm::CoordGeodetic maximum, double cord_step_size, double altitude_step_size,
         double time_step, wmm::MagneticModel *MagneticModel, wmm::Geoid *Geoid, wmm::Ellipsoid Ellip, wmm::Date StartDate,
         wmm::Date EndDate, int ElementOption, int UncertaintyOption, int PrintOption, char *OutputFile);

int main() {
    wmm::MagneticModel *MagneticModels[1];
    wmm::Ellipsoid Ellip;
    wmm::CoordGeodetic minimum, maximum;
    wmm::Geoid Geoid;
    wmm::Date startdate, enddate;
    int ElementOption, PrintOption, i, epochs = 1, UncertaintyOption = 1;
    double cord_step_size, altitude_step_size, time_step_size;
#ifdef WMMHR
    char filename[] = "WMMHR.COF";
#else
    char filename[] = "WMM.COF";
#endif
    char OutputFilename[32];
    char VersionDate[12];
    char ans[20];

    if(!robustReadMagModels(filename, &MagneticModels, 1)) {
        printf("\n %s not found.  Press enter to exit... \n ", filename);
        fgets(ans, 20, stdin);
        return 1;
    }
    strncpy(VersionDate, VERSION_DATE_LARGE + 39, 11);
    VersionDate[11] = '\0';

    SetDefaults(&Ellip, &Geoid);
    /* Set EGM96 Geoid parameters */
    Geoid.geoidHeightBuffer  = std::make_unique<wmm::Geoid::GeoidHeightArray_t>(GeoidHeightBuffer);
    Geoid.isGeoidInitialized = 1;
/* Set EGM96 Geoid parameters END */
#ifdef WMMHR
    printf("\n\n Welcome to the World Magnetic Model High-Resolution(WMMHR) %d C-Program\n", (int)MagneticModels[0]->epoch);
#else
    printf("\n\n Welcome to the World Magnetic Model (WMM) %d C-Program\n", (int)MagneticModels[0]->epoch);
#endif
    printf("of the US National Centers for Environmental Information\n\t\t--- Grid Calculation Program ----\n\t");
    printf("       --- Model Release Date: %s ---\n\t", MODEL_RELEASE_DATE);
    printf("      --- Software Release Date: %s ---\n", VersionDate);
    printf("This program may be used to generate a grid/volume of magnetic field values\n");
    printf("over latitude, longitude, altitude and time axes. To skip an axis,\n");
    printf("keep the start and end values the same and enter zero for the step size.\n");
    printf("Output values are space separated latitude, longitude, altitude, date,\n");
    printf("Magnetic Element, and Uncertainty in that element.  For more information\n");
    printf("see accompanying software manual.");

    printf("\n\n                     Enter grid parameters \n\n");


    /* Get the Lat/Long, Altitude, Time limits from a user interface and print the grid to screen */

    GetUserGrid(&minimum, &maximum, &cord_step_size, &altitude_step_size, &time_step_size, &startdate, &enddate,
                &ElementOption, &PrintOption, OutputFilename, &Geoid, MagneticModels[0]);
    Grid(minimum, maximum, cord_step_size, altitude_step_size, time_step_size, MagneticModels[0], &Geoid, Ellip, startdate,
         enddate, ElementOption, UncertaintyOption, PrintOption, OutputFilename);

    for(i = 0; i < epochs; i++) {
        FreeMagneticModelMemory(MagneticModels[i]);
    }


    printf("\nPress any key to exit...\n");
    getchar();

    return 0;
}

int Grid(CoordGeodetic minimum, CoordGeodetic maximum, double cord_step_size, double altitude_step_size, double time_step,
         MagneticModel *MagneticModel, Geoid *Geoid, Ellipsoid Ellip, Date StartDate, Date EndDate, int ElementOption,
         int UncertaintyOption, int PrintOption, char *OutputFile)

/*This function calls WMM subroutines to generate a grid as defined by the user. The function may be used
to generate a grid of magnetic field elements, time series or a profile. The selected geomagnetic element
is either printed to the file GridResults.txt or to the screen depending on user option.

INPUT: minimum :Data structure with the following elements (minimum limits of the grid)
                                double lambdag; (longitude)
                                double phi; ( geodetic latitude)
                                double heightAboveEllipsoid; (height above the ellipsoid (HaE) )
                                double heightAboveGeoid;(height above the Geoid )
                maximum : same as the above (maximum limist of the grid)
                step_size  : double  : spatial step size, in decimal degrees
                a_step_size : double  :  double altitude step size (km)
                step_time : double  : time step size (decimal years)
                StartDate :  data structure with the following elements used
                                        double decimalYear;     ( decimal years )
                EndDate :	Same as the above;
                MagneticModel :	 data structure with the following elements
                        double EditionDate;
                        double epoch;       Base time of Geomagnetic model epoch (yrs)
                        char  ModelName[20];
                        double *Main_Field_Coeff_G;          C - Gauss coefficients of main geomagnetic model (nT)
                        double *Main_Field_Coeff_H;          C - Gauss coefficients of main geomagnetic model (nT)
                        double *Secular_Var_Coeff_G;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
                        double *Secular_Var_Coeff_H;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
                        int nMax;  Maximum degree of spherical harmonic model
                        int nMaxSecVar; Maxumum degree of spherical harmonic secular model
                        int SecularVariationUsed; Whether or not the magnetic secular variation vector will be needed by
program Geoid :  data structure with the following elements Pointer to data structure Geoid with the following elements int
numbGeoidCols ;   ( 360 degrees of longitude at 15 minute spacing ) int numbGeoidRows ;   ( 180 degrees of latitude  at 15
minute spacing ) int numbHeaderItems ;    ( min, max lat, min, max long, lat, long spacing ) int	scaleFactor;    ( 4 grid
cells per degree at 15 minute spacing  ) float *geoidHeightBuffer;   (Pointer to the memory to store the Geoid elevation data
) int numbGeoidElevs;    (number of points in the gridded file ) int  isGeoidInitialized ;  ( indicates successful
initialization ) Ellip  data  structure with the following elements double a; semi-major axis of the ellipsoid double b;
semi-minor axis of the ellipsoid double fla;  flattening double epssq; first eccentricity squared double eps;  first
eccentricity double re; mean radius of  ellipsoid ElementOption : int : Geomagnetic Element to print
 *        UncertaintyOption: int: 1-Append uncertainties.  Otherwise do not append uncertainties.
          PrintOption : int : 1 Print to File, Otherwise, print to screen

   OUTPUT: none (prints the output to a file )

   CALLS : AllocateModelMemory To allocate memory for model coefficients
      TimelyModifyMagneticModel This modifies the Magnetic coefficients to the correct date.
                  ConvertGeoidToEllipsoidHeight (&CoordGeodetic, &Geoid);   Convert height above msl to height above WGS-84
ellipsoid GeodeticToSpherical Convert from geodeitic to Spherical Equations: 7-8, WMM Technical report
                  ComputeSphericalHarmonicVariables Compute Spherical Harmonic variables
                  AssociatedLegendreFunction Compute ALF  Equations 5-6, WMM Technical report
                  summation Accumulate the spherical harmonic coefficients Equations 10:12 , WMM Technical report
                  rotateMagneticVector Map the computed Magnetic fields to Geodeitic coordinates Equation 16 , WMM Technical
report CalculateGeoMagneticElements Calculate the geoMagnetic elements, Equation 18 , WMM Technical report

 */
{
    int NumTerms;
    double a, b, c, d, PrintElement, ErrorElement = 0;

    MagneticModel *TimedMagneticModel;
    CoordSpherical CoordSpherical;
    MagneticResults MagneticResultsSph, MagneticResultsGeo, MagneticResultsSphVar, MagneticResultsGeoVar;
    SphericalHarmonicVariables *SphVariables;
    GeoMagneticElements GeoMagneticElements, Errors;
    LegendreFunction *LegendreFunction;
    Gradient Gradient;
    int print_boz_warning_weak   = FALSE;
    int print_boz_warning_strong = FALSE;
    int print_alt_warning        = FALSE;
    double min_wgsalt            = -1;
    double max_wgsalt            = 1900;

    FILE *fileout = NULL;

    if(PrintOption == 1) {
        fileout = fopen(OutputFile, "w");
        if(!fileout) {
            printf("Error opening %s to write", OutputFile);
            return FALSE;
        }
    }


    if(fabs(cord_step_size) < 1.0e-10) {
        cord_step_size = 99999.0; /*checks to make sure that the step_size is not too small*/
    }
    if(fabs(altitude_step_size) < 1.0e-10) {
        altitude_step_size = 99999.0;
    }
    if(fabs(time_step) < 1.0e-10) {
        time_step = 99999.0;
    }


    NumTerms           = ((MagneticModel->nMax + 1) * (MagneticModel->nMax + 2) / 2);
    TimedMagneticModel = AllocateModelMemory(NumTerms);
    LegendreFunction   = AllocateLegendreFunctionMemory(NumTerms); /* For storing the ALF functions */
    SphVariables       = AllocateSphVarMemory(MagneticModel->nMax);
    a                  = minimum.heightAboveGeoid;                 /*sets the loop initialization values*/
    b                  = minimum.phi;
    c                  = minimum.lambdag;
    d                  = StartDate.decimalYear;
    double alt         = minimum.heightAboveGeoid;


    for(minimum.heightAboveGeoid = a; minimum.heightAboveGeoid <= maximum.heightAboveGeoid;
        minimum.heightAboveGeoid += altitude_step_size)                                 /* Altitude loop*/
    {
        for(minimum.phi = b; minimum.phi <= maximum.phi; minimum.phi += cord_step_size) /*Latitude loop*/
        {
            for(minimum.lambdag = c; minimum.lambdag <= maximum.lambdag;
                minimum.lambdag += cord_step_size)                                      /*Longitude loop*/
            {
                alt = minimum.heightAboveGeoid;
                if(Geoid->isUseGeoid == 1) {
                    ConvertGeoidToEllipsoidHeight(
                        &minimum,
                        Geoid); /* This converts the height above mean sea level to height above the WGS-84 ellipsoid */
                } else {
                    minimum.heightAboveEllipsoid = minimum.heightAboveGeoid;
                }
#ifndef WMMHR
                if(minimum.heightAboveEllipsoid < min_wgsalt || minimum.heightAboveEllipsoid > max_wgsalt) {
                    printf("\n Unrecognized height: %.2f. \n %s \n", alt, WMM_MileSpec_WARN);
                    print_alt_warning = 1;
                }
#endif
                GeodeticToSpherical(Ellip, minimum, &CoordSpherical);
                ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, MagneticModel->nMax,
                                                  SphVariables); /* Compute Spherical Harmonic variables  */
                AssociatedLegendreFunction(CoordSpherical, MagneticModel->nMax,
                                           LegendreFunction);    /* Compute ALF  Equations 5-6, WMM Technical report*/

                for(StartDate.decimalYear = d; StartDate.decimalYear <= EndDate.decimalYear;
                    StartDate.decimalYear += time_step)          /*year loop*/
                {
                    TimelyModifyMagneticModel(
                        StartDate, MagneticModel,
                        TimedMagneticModel);        /*This modifies the Magnetic coefficients to the correct date. */
                    Summation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical,
                              &MagneticResultsSph); /* Accumulate the spherical harmonic coefficients Equations 10:12 , WMM
                                                       Technical report*/
                    SecVarSummation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical,
                                    &MagneticResultsSphVar); /*Sum the Secular Variation Coefficients, Equations 13:15 , WMM
                                                                Technical report  */
                    RotateMagneticVector(CoordSpherical, minimum, MagneticResultsSph,
                                         &MagneticResultsGeo);    /* Map the computed Magnetic fields to Geodetic coordinates
                                                                     Equation 16 , WMM Technical report */
                    RotateMagneticVector(CoordSpherical, minimum, MagneticResultsSphVar,
                                         &MagneticResultsGeoVar); /* Map the secular variation field components to Geodetic
                                                                     coordinates, Equation 17 , WMM Technical report*/
                    CalculateGeoMagneticElements(
                        &MagneticResultsGeo,
                        &GeoMagneticElements); /* Calculate the Geomagnetic elements, Equation 18 , WMM Technical report */
                    CalculateGridVariation(minimum, &GeoMagneticElements);
                    CalculateSecularVariationElements(
                        MagneticResultsGeoVar,
                        &GeoMagneticElements); /*Calculate the secular variation of each of the Geomagnetic elements,
                                                  Equation 19, WMM Technical report*/
#if WMMHR
                    WMMHRErrorCalc(GeoMagneticElements.H, &Errors);
#else
                    WMMErrorCalc(GeoMagneticElements.H, &Errors);
#endif
                    if(GeoMagneticElements.H <= 2000.0) {
                        print_boz_warning_strong = TRUE;
                    } else if(GeoMagneticElements.H <= 6000.0) {
                        print_boz_warning_weak = TRUE;
                    }

                    if(ElementOption >= 17) {
                        Gradient(Ellip, minimum, TimedMagneticModel, &Gradient);
                    }

                    switch(ElementOption) {
                        case 1:
                            PrintElement =
                                GeoMagneticElements
                                    .Decl; /*1. Angle between the magnetic field vector and true north, positive east*/
                            ErrorElement = Errors.Decl;
                            break;
                        case 2:
                            PrintElement = GeoMagneticElements.Incl; /*2. Angle between the magnetic field vector and the
                                                                        horizontal plane, positive downward*/
                            ErrorElement = Errors.Incl;
                            break;
                        case 3:
                            PrintElement = GeoMagneticElements.F; /*3. Magnetic Field Strength*/
                            ErrorElement = Errors.F;
                            break;
                        case 4:
                            PrintElement = GeoMagneticElements.H; /*4. Horizontal Magnetic Field Strength*/
                            ErrorElement = Errors.H;
                            break;
                        case 5:
                            PrintElement = GeoMagneticElements.X; /*5. Northern component of the magnetic field vector*/
                            ErrorElement = Errors.X;
                            break;
                        case 6:
                            PrintElement = GeoMagneticElements.Y; /*6. Eastern component of the magnetic field vector*/
                            ErrorElement = Errors.Y;
                            break;
                        case 7:
                            PrintElement = GeoMagneticElements.Z; /*7. Downward component of the magnetic field vector*/
                            ErrorElement = Errors.Z;
                            break;
                        case 8:
                            PrintElement = GeoMagneticElements.GV; /*8. The Grid Variation*/
                            ErrorElement = Errors.Decl;
                            break;
                        case 9:
                            PrintElement      = GeoMagneticElements.Decldot * 60; /*9. Yearly Rate of change in declination*/
                            UncertaintyOption = 0;
                            break;
                        case 10:
                            PrintElement = GeoMagneticElements.Incldot * 60; /*10. Yearly Rate of change in inclination*/
                            UncertaintyOption = 0;
                            break;
                        case 11:
                            PrintElement = GeoMagneticElements.Fdot; /*11. Yearly rate of change in Magnetic field strength*/
                            UncertaintyOption = 0;
                            break;
                        case 12:
                            PrintElement =
                                GeoMagneticElements.Hdot; /*12. Yearly rate of change in horizontal field strength*/
                            UncertaintyOption = 0;
                            break;
                        case 13:
                            PrintElement = GeoMagneticElements.Xdot; /*13. Yearly rate of change in the northern component*/
                            UncertaintyOption = 0;
                            break;
                        case 14:
                            PrintElement = GeoMagneticElements.Ydot; /*14. Yearly rate of change in the eastern component*/
                            UncertaintyOption = 0;
                            break;
                        case 15:
                            PrintElement = GeoMagneticElements.Zdot; /*15. Yearly rate of change in the downward component*/
                            UncertaintyOption = 0;
                            break;
                        case 16:
                            PrintElement      = GeoMagneticElements.GVdot;
                            UncertaintyOption = 0;
                            /*16. Yearly rate of change in grid variation*/;
                            break;
                        case 17:
                            PrintElement      = Gradient.gradPhi.X;
                            UncertaintyOption = 0;
                            break;
                        case 18:
                            PrintElement      = Gradient.gradPhi.Y;
                            UncertaintyOption = 0;
                            break;
                        case 19:
                            PrintElement      = Gradient.gradPhi.Z;
                            UncertaintyOption = 0;
                            break;
                        case 20:
                            PrintElement      = Gradient.gradLambda.X;
                            UncertaintyOption = 0;
                            break;
                        case 21:
                            PrintElement      = Gradient.gradLambda.Y;
                            UncertaintyOption = 0;
                            break;
                        case 22:
                            PrintElement      = Gradient.gradLambda.Z;
                            UncertaintyOption = 0;
                            break;
                        case 23:
                            PrintElement      = Gradient.gradZ.X;
                            UncertaintyOption = 0;
                            break;
                        case 24:
                            PrintElement      = Gradient.gradZ.Y;
                            UncertaintyOption = 0;
                            break;
                        case 25:
                            PrintElement      = Gradient.gradZ.Z;
                            UncertaintyOption = 0;
                            break;
                        default:
                            PrintElement =
                                GeoMagneticElements
                                    .Decl; /* 1. Angle between the magnetic field vector and true north, positive east*/
                            ErrorElement = Errors.Decl;
                    }

                    if(Geoid->isUseGeoid == 1) {
                        if(PrintOption == 1) {
                            fprintf(fileout, "%5.2f %6.2f %8.4f %7.2f %10.2f", minimum.phi, minimum.lambdag,
                                    minimum.heightAboveGeoid, StartDate.decimalYear, PrintElement);
                        } else {
                            printf("%5.2f %6.2f %8.4f %7.2f %10.2f", minimum.phi, minimum.lambdag, minimum.heightAboveGeoid,
                                   StartDate.decimalYear, PrintElement);
                        }
                    } else {
                        if(PrintOption == 1) {
                            fprintf(fileout, "%5.2f %6.2f %8.4f %7.2f %10.2f", minimum.phi, minimum.lambdag,
                                    minimum.heightAboveEllipsoid, StartDate.decimalYear, PrintElement);
                        } else {
                            printf("%5.2f %6.2f %8.4f %7.2f %10.2f", minimum.phi, minimum.lambdag,
                                   minimum.heightAboveEllipsoid, StartDate.decimalYear, PrintElement);
                        }
                    }
                    if(UncertaintyOption == 1) {
                        if(PrintOption == 1) {
                            fprintf(fileout, " %7.2f", ErrorElement);
                        } else {
                            printf(" %7.2f", ErrorElement);
                        }
                    }
                    if(PrintOption == 1) {
                        fprintf(fileout, "\n");
                    } else {
                        printf("\n"); /* Complete line */
                    }

                                      /**Below can be used for XYZ Printing format (longitude latitude output_data)
                                       *  fprintf(fileout, "%5.2f %6.2f %10.4f\n", minimum.lambdag, minimum.phi, PrintElement); **/

                } /* year loop */

            } /*Longitude Loop */

        } /* Latitude Loop */

    } /* Altitude Loop */

    if(PrintOption == 1) {
        if(print_boz_warning_strong) {
            fprintf(fileout, "%s", BOZ_WARN_TEXT_STRONG);
        } else if(print_boz_warning_weak) {
            fprintf(fileout, "%s", BOZ_WARN_TEXT_WEAK);
        }
#ifndef WMMHR
        if(!print_alt_warning) {
            fprintf(fileout, "%s\n", WMM_MileSpec_INFO);
        } else {
            fprintf(fileout, "%s\n", WMM_MileSpec_WARN);
        }
#endif
    } else {
        if(print_boz_warning_strong) {
            printf("%s\n", BOZ_WARN_TEXT_STRONG);
        } else if(print_boz_warning_weak) {
            printf("%s\n", BOZ_WARN_TEXT_WEAK);
        }
#ifndef WMMHR
        if(!print_alt_warning) {
            printf("%s\n", WMM_MileSpec_INFO);
        } else {
            printf("%s\n", WMM_MileSpec_WARN);
        }
#endif
    }

    return TRUE;
} /*Grid*/
