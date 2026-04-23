//
// Created by Professional on 23.04.2026.
//
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <string>

#include "GeomagnetismHeader.h"

#include "EGM9615.h"
 #include "GeomagInterativeLib.h"
#include "MagneticUtils.h"
#include "magcalc.h"
#include "version.h"

#include "wmm_warn.h"

/*---------------------------------------------------------------------------*/

/**
 * WMM Point Calculation Program.
 *
 * The Geomagnetism Library is used to make a command prompt program. The program prompts
 * the user to enter a location, performs the computations and prints the results to the
 * standard output. The program expects the files GeomagnetismLibrary.c, GeomagnetismHeader.h,
 * EWMM.COF and EGM9615.h to be in the same directory.
 *
 * Manoj.C.Nair@Noaa.Gov
 * April 21, 2011
 *
 * liyin.young@noaa.gov
 * Updated April, 2023
 *
 **/


char GeomagIntroduction_WMM(const wmm::MagneticModel &magneticModel, const std::string &modelDate);

void help_info(const wmm::MagneticModel &magneticModel, const std::string_view short_name);

int main() {
    char ans[20], b;

    const std::string filename{"WMM.COF"};  // Coefficients file must be placed in same directory as executable

    wmm::MagneticModel magneticModel;
    if(!magneticModel.readModel(filename)) {
        std::cerr << "\n " << filename << " not found.  Press enter to exit... \n ";
        std::cin.get();
        return 1;
    }

    wmm::Date userDate{};  // Default constructor Date using current system date
    wmm::MagneticModel timedMagneticModel{magneticModel.applyDate(userDate)};
    if(std::isnan(magneticModel.epoch) || std::isnan(timedMagneticModel.epoch)) {
        wmm::PrintError(2);
    }

    // Set EGM96 geoid parameters
    wmm::Geoid geoid{};

    // Set EGM96 geoid parameters END
    b = GeomagIntroduction_WMM(magneticModel, MODEL_RELEASE_DATE);

    wmm::Ellipsoid ellip{};
    wmm::CoordGeodetic coordGeodetic;
    wmm::CoordSpherical coordSpherical;
    wmm::GeoMagneticElements geoMagneticElements;
    wmm::GeoMagneticElements errors;

    while(b != 'x') {
        //        if(GetUserInput(magneticModel, &geoid, &coordGeodetic, &userDate) == 1) /*Get User Input */
        {
            coordGeodetic.lambda = 48.8;
            coordGeodetic.phi    = 54.4;
            // Convert from geodetic to Spherical Equations: 17-18, WMM Technical report
            coordSpherical.fromGeodetic(ellip, coordGeodetic);
            point_calc(ellip, coordGeodetic, coordSpherical, userDate, timedMagneticModel, &geoMagneticElements, &errors);

            if(geoMagneticElements.H <= 2000.0) {
                std::cout << std::endl << BOZ_WARN_TEXT_STRONG << std::endl;
            } else if(geoMagneticElements.H <= 6000.0) {
                std::cout << std::endl << BOZ_WARN_TEXT_WEAK << std::endl;
            }
#ifndef WMMHR
            if(coordGeodetic.HeightAboveEllipsoid < -1 || coordGeodetic.HeightAboveEllipsoid > 1900) {
                std::cout << std::endl << WMM_MileSpec_WARN << std::endl;
            } else {
                std::cout << std::endl << WMM_MileSpec_INFO << std::endl;
            }
#endif
            // Print the results
            PrintUserDataWithUncertainty(geoMagneticElements, errors, coordGeodetic, userDate, timedMagneticModel, geoid);
        }

        do {
            printf("\n\n Do you need more point data ? (y or n) \n ");
        } while(nullptr == fgets(ans, 20, stdin));
        switch(ans[0]) {
            case 'Y':
            case 'y':
                break;
            case 'N':
            case 'n':
                return 0;
            default:
                return 0;
        }
    }

    return 0;
}

// Prints the introduction to the Geomagnetic program.  It needs the Magnetic model for the epoch.
char GeomagIntroduction_WMM(const wmm::MagneticModel &magneticModel, const std::string &modelDate) {
    std::string versionDate{VERSIONDATE_LARGE};
    versionDate = versionDate.substr(39, 11);

    int ans_size = 10;
    char help[ans_size];

    std::string msg{};
#ifdef WMMHR
    std::string short_name{"WMMHR"};
    msg.append("\n\n Welcome to the World Magnetic Model High-Resolution (WMMHR) ");
#else
    std::string short_name{"WMM"};
    msg.append("\n\n Welcome to the World Magnetic Model ");
    msg.append(std::to_string(static_cast<int>(magneticModel.epoch)) + " C-Program\n\n");
#endif
    msg.append("              --- Model Release Date: " + modelDate + " ---\n");
    msg.append("            --- Software Release Date: " + versionDate + " ---\n\n");
    msg.append("\n This program estimates the strength and direction of ");
    msg.append("\n Earth's main Magnetic field for a given point/area.");
    msg.append("\n Enter h for help and contact information or c to continue.");
    msg.append("\n >");
    std::cout << msg;

    while(fgets(help, ans_size, stdin) == NULL || (help[0] != 'C' && help[0] != 'c')) {
        if(help[0] == 'h' || help[0] == 'H') {
            help_info(magneticModel, short_name);
        }
        std::cout << "\n Enter h for help and contact information or c to continue."
                  << "\n >";
    }


    return help[0];
} /*GeomagIntroduction_WMM*/

void help_info(const wmm::MagneticModel &magneticModel, const std::string_view short_name) {
    printf("\n Help information ");

#ifdef WMMHR
    printf("\n The World Magnetic Model High-Resolution (WMMHR) for %d", (int)magneticModel->epoch);
#else
    printf("\n The World Magnetic Model (WMM) for %d", (int)magneticModel.epoch);
#endif
    printf("\n is a model of Earth's main magnetic field. The %s", short_name.data());
    printf("\n is recomputed every five (5) years, in years divisible by ");
    printf("\n five (i.e. 2020, 2025). See the contact information below");
    printf("\n to obtain more information on the %s and associated software.", short_name.data());
    printf("\n ");
    printf("\n Input required is the location in geodetic latitude and");
    printf("\n longitude (positive for northern latitudes and eastern ");
    printf("\n longitudes), geodetic altitude in meters, and the date of ");
    printf("\n interest in years.");

    printf("\n\n\n The program computes the estimated Magnetic Declination");
    printf("\n (Decl) which is sometimes called MagneticVAR, Inclination (Incl), Total");
    printf("\n Intensity (F or TI), Horizontal Intensity (H or HI), Vertical");
    printf("\n Intensity (Z), and Grid Variation (GV). Declination and Grid");
    printf("\n Variation are measured in units of degrees and are considered");
    printf("\n positive when east or north. Inclination is measured in units");
    printf("\n of degrees and is considered positive when pointing down (into");
    printf("\n the Earth). The %s is referenced to the WGS-84 ellipsoid and", short_name.data());
    printf("\n is valid for 5 years after the base epoch. Uncertainties for the");
    printf("\n %s are one standard deviation uncertainties averaged over the globe.", short_name.data());
    printf("\n We represent the uncertainty as constant values in Incl, F, H, X,");
    printf("\n Y, and Z. Uncertainty in Declination varies depending on the strength");
    printf("\n of the horizontal field. For more information see the %s Technical", short_name.data());
    printf("\n Report.");

    printf("\n\n\n It is very important to note that a degree and order 133 model,");
    printf("\n such as %s, describes only the longest wavelength spatial magnetic ", short_name.data());
    printf("\n fluctuations due to Earth's core, mantle and crust. Not included in the %s series", short_name.data());
    printf("\n models are intermediate and short wavelength spatial fluctuations ");
    printf("\n that originate in Earth's mantle and crust. Consequently, isolated");
    printf("\n angular errors at various positions on the surface (primarily over");
    printf("\n land, along continental margins and over oceanic sea-mounts, ridges and");
    printf("\n trenches) of several degrees may be expected. Also not included in");
    printf("\n the model are temporal fluctuations of magnetospheric and ionospheric");
    printf("\n origin. On the days during and immediately following magnetic storms,");
    printf("\n temporal fluctuations can cause substantial deviations of the Geomagnetic");
    printf("\n field from model values. If the required declination accuracy is");
    printf("\n more stringent than the %s series of models provide, the user is", short_name.data());
    printf("\n advised to request special (regional or local) surveys be performed");
#ifdef WMMHR
    printf("\n and models prepared. The World Magnetic Model High-Resolution is a joint product of");
#else
    printf("\n and models prepared. The World Magnetic Model is a joint product of");
#endif
    printf("\n the United States' National Geospatial-Intelligence Agency (NGA) and");
    printf("\n the United Kingdom's Defence Geographic Centre (DGC). The %s was", short_name.data());
    printf("\n developed jointly by the National Centers for Environmental Information");
    printf("\n (NCEI, Boulder CO, USA) and the British Geological Survey (BGS, ");
    printf("\n Edinburgh, Scotland).");

    printf("\n\n\n Contact Information");

    printf("\n  Software and Model Support");
    printf("\n	National Centers for Environmental Information");
    printf("\n	NOAA E/NE42");
    printf("\n	325 Broadway");
    printf("\n	Boulder, CO 80305 USA");
    printf("\n	Attn: Manoj Nair or Arnaud Chulliat");
    printf("\n	Email:  geomag.Models@noaa.gov \n");
}
