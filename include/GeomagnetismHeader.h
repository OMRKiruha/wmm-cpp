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

#pragma once

constexpr auto MODEL_RELEASE_DATE{"17 Dec 2024"};
constexpr auto VERSION_DATE_LARGE{"$Date: 2024-12-17 10:40:43 -0700 (Tue, 17 Dec 2024) $"};

// Easy use #include "GeomagnetismHeader.h" for all library
#include "CoordGeodetic.h"
#include "CoordSpherical.h"
#include "Date.h"
#include "Ellipsoid.h"
#include "GeoMagneticElements.h"
#include "Geoid.h"
#include "Gradient.h"
#include "LegendreFunction.h"
#include "MagneticConstants.h"
#include "MagneticModel.h"
#include "MagneticResults.h"
#include "MagneticUtils.h"
#include "PointCalc.h"
#include "SphericalHarmonicVariables.h"
#include "UTMParameters.h"