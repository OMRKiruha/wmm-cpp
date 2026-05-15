
#include "GeomagnetismHeader.h"
#include "MagneticUtils.h"

#include <cmath>
#include <cstdio>

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

    /******************************************************************************
     *************Conversions, Transformations, and other Calculations**************
     * This grouping consists of functions that perform unit conversions, coordinate
     * transformations and other simple or straightforward calculations that are
     * usually easily replicable with a typical scientific calculator.
     ******************************************************************************/


    void BaseErrors(const double DeclCoef, const double DeclBaseline, const double InclOffset, const double FOffset,
                    const double Multiplier, const double H,
                    double *DeclErr, double *InclErr, double *FErr) {
        double declHorizontalAdjustmentSq;
        declHorizontalAdjustmentSq = (DeclCoef / H) * (DeclCoef / H);
        *DeclErr                   = sqrt(declHorizontalAdjustmentSq + DeclBaseline * DeclBaseline) * Multiplier;
        *InclErr                   = InclOffset * Multiplier;
        *FErr                      = FOffset * Multiplier;
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
    void DegreeToDMSstring(const double DegreesOfArc, const int UnitDepth, std::string &out) {
        double temp = DegreesOfArc;

        if(UnitDepth > 3) {
            printError(21);
        }

        for(int i = 0; i < UnitDepth; i++) {
            int DMS  = static_cast<int>(temp);
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
    void DMSstringToDegree(const std::string_view DMSstring, double *DegreesOfArc) {
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

}  // namespace wmm