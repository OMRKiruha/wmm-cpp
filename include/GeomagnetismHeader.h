/*	WMM Subroutine library was tested in the following environments
 *
 *	1. Red Hat Linux  with GCC Compiler
 *	2. MS Windows XP with CodeGear C++ compiler
 *	3. Sun Solaris with GCC Compiler
 *
 *      Revision Number: $Revision: 1437 $
 *      Last changed by: $Author: Li-Yin Young $
 *      Last changed on: $Date: 2024-11-08 10:49:40 -0700 $
 */

#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#define WGS84ON 1
#define MSLON 2

namespace wmm {
    // Constants
    /*These error values come from the ISCWSA error model:
     *http://www.copsegrove.com/Pages/MWDGeomagneticModels.aspx
     */
    static constexpr double INCL_ERROR_BASE        = 0.20;
    static constexpr double DECL_ERROR_OFFSET_BASE = 0.36;
    static constexpr int F_ERROR_BASE              = 130;
    static constexpr int DECL_ERROR_SLOPE_BASE     = 5000;
    static constexpr double WMM_ERROR_MULTIPLIER   = 1.21;
    static constexpr double IGRF_ERROR_MULTIPLIER  = 1.21;

    static constexpr int WMMHR_UNCERTAINTY_F           = 134;
    static constexpr int WMMHR_UNCERTAINTY_H           = 130;
    static constexpr int WMMHR_UNCERTAINTY_X           = 135;
    static constexpr int WMMHR_UNCERTAINTY_Y           = 85;
    static constexpr int WMMHR_UNCERTAINTY_Z           = 134;
    static constexpr double WMMHR_UNCERTAINTY_I        = 0.19;
    static constexpr double WMMHR_UNCERTAINTY_D_OFFSET = 0.25;
    static constexpr int WMMHR_UNCERTAINTY_D_COEF      = 5205;

    /*These error values are the NCEI error model */
    static constexpr int WMM_UNCERTAINTY_F           = 138;
    static constexpr int WMM_UNCERTAINTY_H           = 133;
    static constexpr int WMM_UNCERTAINTY_X           = 137;
    static constexpr int WMM_UNCERTAINTY_Y           = 89;
    static constexpr int WMM_UNCERTAINTY_Z           = 141;
    static constexpr double WMM_UNCERTAINTY_I        = 0.20;
    static constexpr double WMM_UNCERTAINTY_D_OFFSET = 0.26;
    static constexpr int WMM_UNCERTAINTY_D_COEF      = 5417;

    static double Rad2Deg(const double rad) {
        return rad * (180.0 / M_PI);
    }

    static double Deg2Rad(const double deg) {
        return deg * (M_PI / 180.0);
    }

    static constexpr int PS_MIN_LAT_DEGREE     = -55;   /* Minimum Latitude for  Polar Stereographic projection in degrees */
    static constexpr int PS_MAX_LAT_DEGREE     = 55;    /* Maximum Latitude for Polar Stereographic projection in degrees  */
    static constexpr double UTM_MIN_LAT_DEGREE = -80.5; /* Minimum Latitude for UTM projection in degrees     */
    static constexpr double UTM_MAX_LAT_DEGREE = 84.5;  /* Maximum Latitude for UTM projection in degrees     */

    static constexpr double GEO_POLE_TOLERANCE = 1e-5;
    static constexpr int USE_GEOID             = 1;     /* 1 Geoid - Ellipsoid difference should be corrected, 0 otherwise */

    static constexpr int LAT_BOUND_MIN         = -90;
    static constexpr int LAT_BOUND_MAX         = 90;
    static constexpr int LON_BOUND_MIN         = -180;
    static constexpr int LON_BOUND_MAX         = 360;
    static constexpr int ALT_BOUND_MIN         = -10;
    static constexpr int NO_ALT_MAX            = -99999;
    static constexpr int USER_GAVE_UP          = -1;
    static constexpr double DEC_YEAR_BOUND_MIN = 2024.866;
    static constexpr double DEC_YEAR_BOUND_MAX = 2030;

    /*
Data types and prototype declaration for
World Magnetic Model (WMM) subroutines.

July 28, 2009

manoj.c.nair@noaa.gov*/

    struct MagneticModel {
        MagneticModel() = default;
        explicit MagneticModel(int numTerms);
        ~MagneticModel() = default;

        double editionDate{};
        double epoch{}; /*Base time of Geomagnetic model epoch (yrs)*/
        double min_year{};
        std::string modelName;
        std::vector<double>
            main_Field_Coeff_G; /* C - Gauss coefficients of main geomagnetic model (nT) Index is (n * (n + 1) / 2 + m) */
        std::vector<double> main_Field_Coeff_H;  /* C - Gauss coefficients of main geomagnetic model (nT) */
        std::vector<double> secular_Var_Coeff_G; /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */
        std::vector<double> secular_Var_Coeff_H; /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */
        int nMax{};                              /* Maximum degree of spherical harmonic model */
        int nMaxSecVar{};                        /* Maximum degree of spherical harmonic secular model */
        int secularVariationUsed{}; /* Whether or not the magnetic secular variation vector will be needed by program*/
        double coefficientFileEndDate{};
    };

    struct Ellipsoid {
        /* Sets WGS-84 parameters */
        double a{6378.137};                      /*semi-major axis of the ellipsoid*/
        double b{6356.7523142};                  /*semi-minor axis of the ellipsoid*/
        double fla{1 / 298.257223563};           /* flattening */
        double epssq{eps * eps};                 /*first eccentricity squared */
        double eps{sqrt(1 - (b * b) / (a * a))}; /* first eccentricity */
        double re{6371.2};                       /* mean radius of  ellipsoid*/
    };

    struct CoordGeodetic {
        double lambda{};               /* geodetic longitude */
        double phi{};                  /* geodetic latitude */
        double HeightAboveEllipsoid{}; /* height above the ellipsoid (HaE) */
        double HeightAboveGeoid{};     /* (height above the EGM96 geoid model ) */
        int UseGeoid{};
    };

    struct CoordSpherical {
        double lambda{}; /* geocentric longitude*/
        double phig{};   /* geocentric latitude*/
        double r{};      /* distance from the center of the ellipsoid*/
    };

    struct Date {
        Date();
        int Year{};
        int Month{};
        int Day{};
        double DecimalYear{}; /* decimal years */
    };

    struct LegendreFunction {
        explicit LegendreFunction(int numTerms);
        ~LegendreFunction() = default;

        std::vector<double> Pcup;  /* Legendre Function */
        std::vector<double> dPcup; /* Derivative of Legendre fcn */
    };

    struct MagneticResults {
        double Bx{}; /* North */
        double By{}; /* East */
        double Bz{}; /* Down */
    };

    struct SphericalHarmonicVariables {
        explicit SphericalHarmonicVariables(int numTerms);
        ~SphericalHarmonicVariables() = default;

        std::vector<double> RelativeRadiusPower; /* [earth_reference_radius_km / sph. radius ]^n      */
        std::vector<double> cos_mlambda;         /* cp(m)  - cosine of (m*spherical coord. longitude) */
        std::vector<double> sin_mlambda;         /* sp(m)  - sine of (m*spherical coord. longitude)   */
    };

    struct GeoMagneticElements {
        double Decl{};    /* 1. Angle between the magnetic field vector and true north, positive east*/
        double Incl{};    /*2. Angle between the magnetic field vector and the horizontal plane, positive down*/
        double F{};       /*3. Magnetic Field Strength*/
        double H{};       /*4. Horizontal Magnetic Field Strength*/
        double X{};       /*5. Northern component of the magnetic field vector*/
        double Y{};       /*6. Eastern component of the magnetic field vector*/
        double Z{};       /*7. Downward component of the magnetic field vector*/
        double GV{};      /*8. The Grid Variation*/
        double Decldot{}; /*9. Yearly Rate of change in declination*/
        double Incldot{}; /*10. Yearly Rate of change in inclination*/
        double Fdot{};    /*11. Yearly rate of change in Magnetic field strength*/
        double Hdot{};    /*12. Yearly rate of change in horizontal field strength*/
        double Xdot{};    /*13. Yearly rate of change in the northern component*/
        double Ydot{};    /*14. Yearly rate of change in the eastern component*/
        double Zdot{};    /*15. Yearly rate of change in the downward component*/
        double GVdot{};   /*16. Yearly rate of change in grid variation*/
    };

    struct Geoid {
        /* Sets EGM-96 model file parameters */
        int NumbGeoidCols{1441};  /* 360 degrees of longitude at 15 minute spacing */
        int NumbGeoidRows{721};   /* 180 degrees of latitude  at 15 minute spacing */
        int NumbHeaderItems{6};   /* min, max lat, min, max long, lat, long spacing*/
        int ScaleFactor{4};       /* 4 grid cells per degree at 15 minute spacing  */
        std::array<float, 1038961> GeoidHeightBuffer{};
        int NumbGeoidElevs{NumbGeoidCols * NumbGeoidRows};
        int Geoid_Initialized{0}; /* indicates successful initialization */
        int UseGeoid{USE_GEOID};  /*Is the Geoid being used?*/
    };

    struct Gradient {
        int UseGradient{};
        GeoMagneticElements GradPhi;    /* phi */
        GeoMagneticElements GradLambda; /* lambda */
        GeoMagneticElements GradZ;
    };

    struct CoordGeodeticStr {
        std::string Longitude;
        std::string Latitude;
    };

    struct UTMParameters {
        double Easting{};  /* (X) in meters*/
        double Northing{}; /* (Y) in meters */
        int Zone{};        /*UTM Zone*/
        char HemiSphere{};
        double CentralMeridian{};
        double ConvergenceOfMeridians{};
        double PointScale{};
    };

    enum PARAMS {
        SHDF,
        MODELNAME,
        PUBLISHER,
        RELEASEDATE,
        DATACUTOFF,
        MODELSTARTYEAR,
        MODELENDYEAR,
        EPOCH,
        INTSTATICDEG,
        INTSECVARDEG,
        EXTSTATICDEG,
        EXTSECVARDEG,
        GEOMAGREFRAD,
        NORMALIZATION,
        SPATBASFUNC
    };

    enum COEFFICIENTS { IE, N, M, GNM, HNM, DGNM, DHNM };

    enum YYYYMMDD { YEAR, MONTH, DAY };

    /*Prototypes */

    /*Functions that should be Magnetic Model member functions*/


    /*Wrapper Functions*/
    int Geomag(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
               MagneticModel *timedMagneticModel, GeoMagneticElements *geoMagneticElements);

    void CalcGradient(Ellipsoid ellip, CoordGeodetic coordGeodetic, MagneticModel *timedMagneticModel, Gradient *gradient);


    //    int robustReadMagneticModel_Large(std::string_view filename, char *filenameSV, std::vector<MagneticModel>&
    //    magneticModels);

    int robustReadMagModels(std::string_view filename, MagneticModel &magneticModel, int array_size);

    /*User Interface*/

    void PrintError(int control);

    void PrintGradient(const Gradient &gradient);

    void PrintUserData(GeoMagneticElements geomagElements, CoordGeodetic spaceInput, Date timeInput,
                       MagneticModel *magneticModel, Geoid *geoid);


    int Warnings(int control, double value, const MagneticModel &magneticModel);

    /*Memory and File Processing*/

    int ReadMagneticModel(std::string_view filename, MagneticModel *magneticModel);

    void AssignHeaderValues(MagneticModel *model, std::vector<std::string> values);

    void AssignMagneticModelCoeffs(MagneticModel *Assignee, MagneticModel *Source, int nMax, int nMaxSecVar);

    void PrintWMMFormat(char *filename, const MagneticModel &magneticModel);

    void PrintEMMFormat(char *filename, char *filenameSV, const MagneticModel &magneticModel);

    void PrintSHDFFormat(char *filename, MagneticModel *(*MagneticModel)[], int epochs);

    int readMagneticModel(char *filename, MagneticModel *MagneticModel);

    //    int readMagneticModel_Large(char *filename, char *filenameSV, MagneticModel *MagneticModel);

    //    int readMagneticModel_SHDF(char *filename, MagneticModel *(*magneticmodels)[], int array_size);

    //    char *Trim(char *str);

    /*Conversions, Transformations, and other Calculations*/
    void BaseErrors(double declCoef, double declBaseline, double inclOffset, double fOffset, double multiplier, double H,
                    double *declErr, double *inclErr, double *fErr);

    int CalculateGeoMagneticElements(const MagneticResults &magneticResultsGeo, GeoMagneticElements *geoMagneticElements);

    void CalculateGradientElements(const MagneticResults &gradResults, const GeoMagneticElements &magneticElements,
                                   GeoMagneticElements *gradElements);

    int CalculateSecularVariationElements(const MagneticResults &magneticVariation, GeoMagneticElements *magneticElements);

    int CalculateGridVariation(const CoordGeodetic &location, GeoMagneticElements *elements);

    void CartesianToGeodetic(const Ellipsoid &ellip, double x, double y, double z, CoordGeodetic *coordGeodetic);

    CoordGeodetic CoordGeodeticAssign(const CoordGeodetic &coordGeodetic);

    int DateToYear(Date *calendarDate, std::string &error);

    void DegreeToDMSstring(double degreesOfArc, int unitDepth, std::string &out);

    void DMSstringToDegree(std::string_view DMSstring, double *degreesOfArc);

    void ErrorCalc(GeoMagneticElements B, GeoMagneticElements *Errors);

    int GeodeticToSpherical(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, CoordSpherical *coordSpherical);

    GeoMagneticElements GeoMagneticElementsAssign(GeoMagneticElements elements);

    GeoMagneticElements GeoMagneticElementsScale(GeoMagneticElements elements, double factor);

    GeoMagneticElements GeoMagneticElementsSubtract(GeoMagneticElements minuend, GeoMagneticElements subtrahend);

    int GetTransverseMercator(const CoordGeodetic &coordGeodetic, UTMParameters *utmParameters);

    int GetUTMParameters(double latitude, double longitude, int *zone, char *hemisphere, double *centralMeridian);

    int RotateMagneticVector(const CoordSpherical &, const CoordGeodetic &coordGeodetic,
                             const MagneticResults &magneticResultsSph, MagneticResults *magneticResultsGeo);

    void SphericalToCartesian(CoordSpherical coordSpherical, double *x, double *y, double *z);

    void SphericalToGeodetic(Ellipsoid ellip, CoordSpherical coordSpherical, CoordGeodetic *CoordGeodetic);

    void TMfwd4(double Eps, double Epssq, double K0R4, double K0R4oa, const std::array<double, 8> &Acoeff, double Lam0,
                double K0, double falseE, double falseN, int XYonly, double Lambda, double Phi, double *X, double *Y,
                double *pscale, double *CoM);

    int YearToDate(Date *Date);


    /*Spherical Harmonics*/

    int AssociatedLegendreFunction(const CoordSpherical &coordSpherical, int nMax, LegendreFunction *legendreFunction);

    int CheckGeographicPole(CoordGeodetic *coordGeodetic);

    int ComputeSphericalHarmonicVariables(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, int nMax,
                                          SphericalHarmonicVariables *sphVariables);

    void GradY(Ellipsoid ellip, CoordSpherical coordSpherical, CoordGeodetic coordGeodetic,
               MagneticModel *timedMagneticModel, GeoMagneticElements geoMagneticElements,
               GeoMagneticElements *gradYElements);

    void GradYSummation(LegendreFunction *legendreFunction, MagneticModel *magneticModel,
                        SphericalHarmonicVariables sphVariables, CoordSpherical coordSpherical, MagneticResults *gradY);

    int PcupHigh(LegendreFunction &legendreFunction, double x, int nMax);

    int PcupLow(LegendreFunction &legendreFunction, double x, int nMax);

    int SecVarSummation(const LegendreFunction &legendreFunction, MagneticModel *magneticModel,
                        const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical,
                        MagneticResults *magneticResults);

    int SecVarSummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                               const CoordSpherical &coordSpherical, MagneticResults *magneticResults);

    int Summation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                  const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical,
                  MagneticResults *magneticResults);

    int SummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                         const CoordSpherical &coordSpherical, MagneticResults *magneticResults);

    int TimelyModifyMagneticModel(const Date &userDate, const MagneticModel &magModel, MagneticModel *timedMagModel);

    /*Geoid*/


    int ConvertGeoidToEllipsoidHeight(CoordGeodetic *coordGeodetic, Geoid *geoid);
    /*
     * The function Convert_Geoid_To_Ellipsoid_Height converts the specified WGS84
     * geoid height at the specified geodetic coordinates to the equivalent
     * ellipsoid height, using the EGM96 gravity model.
     *
     *    Latitude            : Geodetic latitude in radians           (input)
     *    Longitude           : Geodetic longitude in radians          (input)
     *    Geoid_Height        : Geoid height, in meters                (input)
     *    Ellipsoid_Height    : Ellipsoid height, in meters.           (output)
     *
     */

    int GetGeoidHeight(double latitude, double longitude, double *deltaHeight, Geoid *geoid);
    /*
     * The private function Get_Geoid_Height returns the height of the
     * WGS84 geiod above or below the WGS84 ellipsoid,
     * at the specified geodetic coordinates,
     * using a grid of height adjustments from the EGM96 gravity model.
     *
     *    Latitude            : Geodetic latitude in radians           (input)
     *    Longitude           : Geodetic longitude in radians          (input)
     *    DeltaHeight         : Height Adjustment, in meters.          (output)
     *
     */

    void EquivalentLatLon(double lat, double lon, double *repairedLat, double *repairedLon);

    void WMMErrorCalc(double H, GeoMagneticElements *Uncertainty);
    void WMMHRErrorCalc(double H, GeoMagneticElements *Uncertainty);
    void PrintUserDataWithUncertainty(GeoMagneticElements GeomagElements, GeoMagneticElements Errors,
                                      CoordGeodetic SpaceInput, Date TimeInput, MagneticModel *MagneticModel, Geoid *Geoid);

    double dateStr_to_decYear(std::string_view edit_date);

    double date_to_decYear(int year, int month, int day);

}  // namespace wmm