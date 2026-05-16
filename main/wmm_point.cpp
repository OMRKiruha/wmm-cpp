/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <string>

#include "GeomagInterativeLib.h"
#include "GeomagnetismHeader.h"

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

void help_info(const wmm::MagneticModel &magneticModel, const std::string &shortName);

constexpr auto max_size = std::numeric_limits<std::streamsize>::max();

int main() {
#ifdef WMMHR
    const std::string filename{"WMMHR.COF"};
    const std::string program_name{"wmmhr_point"};
#else
    const std::string filename{"WMM.COF"};
    const std::string program_name{"wmm_point"};
#endif

    wmm::MagneticModel magneticModel;
    if(!magneticModel.readModel(filename)) {
        std::cerr << std::format("\n {} not found.  Press enter to exit... \n", filename);
        std::cin.get();
        return 1;
    }

    if(std::isnan(magneticModel.epoch)) {
        std::cerr << "\nError in MagneticModel.\n";
    }

    // Set EGM96 geoid parameters
    wmm::Date userDate{};
    const wmm::Geoid geoid{};
    const wmm::Ellipsoid ellip{};
    wmm::CoordGeodetic coordGeodetic;
    wmm::CoordSpherical coordSpherical;
    wmm::GeoMagneticElements geoMagneticElements;
    wmm::GeoMagneticElements errors;

    if(GeomagIntroduction_WMM(magneticModel, MODEL_RELEASE_DATE) != 'x') {
        //        if(GetUserInput(magneticModel, &geoid, &coordGeodetic, &userDate) == 1) /*Get User Input */
        {
            coordGeodetic.lambda = 48.8;
            coordGeodetic.phi    = 54.4;
            coordGeodetic.convertGeoidToEllipsoidHeight(geoid);
            // Convert from geodetic to Spherical Equations: 17-18, WMM Technical report
            coordSpherical.fromGeodetic(ellip, coordGeodetic);
            pointCalc(ellip, coordGeodetic, coordSpherical, userDate, magneticModel, geoMagneticElements, errors);

            if(geoMagneticElements.H <= 2000.0) {
                std::cout << std::endl << BOZ_WARN_TEXT_STRONG << std::endl;
            } else if(geoMagneticElements.H <= 6000.0) {
                std::cout << std::endl << BOZ_WARN_TEXT_WEAK << std::endl;
            }
#ifndef WMMHR
            if(coordGeodetic.heightAboveEllipsoid < -1 || coordGeodetic.heightAboveEllipsoid > 1900) {
                std::cout << std::endl << WMM_MileSpec_WARN << std::endl;
            } else {
                std::cout << std::endl << WMM_MileSpec_INFO << std::endl;
            }
#endif
            // Print the results
            PrintUserDataWithUncertainty(geoMagneticElements, errors, coordGeodetic, userDate, magneticModel, geoid);
        }

        std::cout << "\n\n Do you need more point data ? (y or n) \n ";
        std::cin.clear();
        std::cin.ignore(max_size, '\n');
        switch(std::cin.get()) {
            case 'Y':
            case 'y':
                break;
            case 'N':
            case 'n':
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

/** @brief Prints the introduction to the Geomagnetic program.  It needs the Magnetic model for the epoch.
 */
char GeomagIntroduction_WMM(const wmm::MagneticModel &magneticModel, const std::string &modelDate) {
    std::string versionDate{VERSION_DATE_LARGE};
    versionDate = versionDate.substr(39, 11);

    std::string msg{};
#ifdef WMMHR
    const std::string short_name{"WMMHR"};
    msg.append("\n\n Welcome to the World Magnetic Model High-Resolution (WMMHR) ");
#else
    const std::string short_name{"WMM"};
    msg.append("\n\n Welcome to the World Magnetic Model ");
#endif
    msg.append(std::to_string(static_cast<int>(magneticModel.epoch)) + " C++ Program\n\n");
    msg.append("              --- Model Release Date: " + modelDate + " ---\n");
    msg.append("            --- Software Release Date: " + versionDate + " ---\n\n");
    msg.append("\n This program estimates the strength and direction of ");
    msg.append("\n Earth's main Magnetic field for a given point/area.");
    msg.append("\n Enter h for help and contact information or c to continue.");
    msg.append("\n >");
    std::cout << msg;

    char c{};
    while(true) {
        c = std::cin.get();
        if(c == 'C' || c == 'c') {
            break;
        }
        if(c == 'H' || c == 'h') {
            help_info(magneticModel, short_name);
        }
        std::cin.clear();
        std::cin.ignore(max_size, '\n');
        std::cout << "\n Enter h for help and contact information or c to continue. \n >";
    }

    return c;
}

void help_info(const wmm::MagneticModel &magneticModel, const std::string &shortName) {
    std::string msg{};
    msg.append("\n Help information ");

#ifdef WMMHR
    msg.append("\n The World Magnetic Model High-Resolution (WMMHR) for ");
#else
    msg.append("\n The World Magnetic Model (WMM) for ");
#endif
    msg.append(std::to_string(static_cast<int>(magneticModel.epoch)));
    msg.append("\n is a model of Earth's main magnetic field. The " + shortName);
    msg.append("\n is recomputed every five (5) years, in years divisible by ");
    msg.append("\n five (i.e. 2020, 2025). See the contact information below");
    msg.append("\n to obtain more information on the " + shortName + " and associated software." + shortName);

    msg.append("\n\n Input required is the location in geodetic latitude and");
    msg.append("\n longitude (positive for northern latitudes and eastern ");
    msg.append("\n longitudes), geodetic altitude in meters, and the date of ");
    msg.append("\n interest in years.");

    msg.append("\n\n The program computes the estimated Magnetic Declination");
    msg.append("\n (Decl) which is sometimes called MagneticVAR, Inclination (Incl), Total");
    msg.append("\n Intensity (F or TI), Horizontal Intensity (H or HI), Vertical");
    msg.append("\n Intensity (Z), and Grid Variation (GV). Declination and Grid");
    msg.append("\n Variation are measured in units of degrees and are considered");
    msg.append("\n positive when east or north. Inclination is measured in units");
    msg.append("\n of degrees and is considered positive when pointing down (into");
    msg.append("\n the Earth). The " + shortName + " is referenced to the WGS-84 ellipsoid and");
    msg.append("\n is valid for 5 years after the base epoch. Uncertainties for the");
    msg.append("\n " + shortName + " are one standard deviation uncertainties averaged over the globe.");
    msg.append("\n We represent the uncertainty as constant values in Incl, F, H, X,");
    msg.append("\n Y, and Z. Uncertainty in Declination varies depending on the strength");
    msg.append("\n of the horizontal field. For more information see the " + shortName + " Technical");
    msg.append("\n Report.");

    msg.append("\n\n It is very important to note that a degree and order 133 model,");
    msg.append("\n such as " + shortName + ", describes only the longest wavelength spatial magnetic ");
    msg.append("\n fluctuations due to Earth's core, mantle and crust. Not included in the " + shortName + " series");
    msg.append("\n models are intermediate and short wavelength spatial fluctuations ");
    msg.append("\n that originate in Earth's mantle and crust. Consequently, isolated");
    msg.append("\n angular errors at various positions on the surface (primarily over");
    msg.append("\n land, along continental margins and over oceanic sea-mounts, ridges and");
    msg.append("\n trenches) of several degrees may be expected. Also not included in");
    msg.append("\n the model are temporal fluctuations of magnetospheric and ionospheric");
    msg.append("\n origin. On the days during and immediately following magnetic storms,");
    msg.append("\n temporal fluctuations can cause substantial deviations of the Geomagnetic");
    msg.append("\n field from model values. If the required declination accuracy is");
    msg.append("\n more stringent than the " + shortName + " series of models provide, the user is");
    msg.append("\n advised to request special (regional or local) surveys be performed");
#ifdef WMMHR
    msg.append("\n and models prepared. The World Magnetic Model High-Resolution is a joint product of");
#else
    msg.append("\n and models prepared. The World Magnetic Model is a joint product of");
#endif
    msg.append("\n the United States' National Geospatial-Intelligence Agency (NGA) and");
    msg.append("\n the United Kingdom's Defence Geographic Centre (DGC). The " + shortName + " was");
    msg.append("\n developed jointly by the National Centers for Environmental Information");
    msg.append("\n (NCEI, Boulder CO, USA) and the British Geological Survey (BGS, ");
    msg.append("\n Edinburgh, Scotland).");

    msg.append("\n\n Contact Information");
    msg.append("\n  Software and Model Support");
    msg.append("\n	National Centers for Environmental Information");
    msg.append("\n	NOAA E/NE42");
    msg.append("\n	325 Broadway");
    msg.append("\n	Boulder, CO 80305 USA");
    msg.append("\n	Attn: Manoj Nair or Arnaud Chulliat");
    msg.append("\n	Email:  geomag.Models@noaa.gov \n");

    std::cout << msg;
}
