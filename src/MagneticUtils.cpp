//
// Created by Professional on 21.04.2026.
//

#include "MagneticUtils.h"

#include "Date.h"
#include "MagneticModel.h"

#include <array>
#include <iostream>
#include <numeric>

namespace wmm {

    static std::array days = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    bool dateStr_to_ymd(const std::string_view str, int &year, int &month, int &day) {
        if(sscanf(str.data(), "%d/%d/%d", &month, &day, &year) != 3) {
            std::cerr << "Failed to parse the date string. Please use the format mm/dd/yyyy\n";
            return false;
        }

        if(month > 12 || month < 1) {
            std::cerr << "Month out of range\n";
            return false;
        }

        if((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
            days.at(2) = 29;
        } else {
            days.at(2) = 28;
        }

        if(day > days.at(month) || day < 1) {
            std::cerr << "Day out of range\n";
            return false;
        }
        return true;
    }

    // Parse the date string format as mm/dd/yyyy
    double dateStr_to_decYear(const std::string_view edit_date) {
        int day{}, month{}, year{};

        if(!dateStr_to_ymd(edit_date, year, month, day)) {
            return 0;
        }

        return date_to_decYear(year, month, day);
    }

    // Parse the date string format as mm/dd/yyyy
    double date_to_decYear(const int year, const int month, const int day) {
        if((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
            days.at(2) = 29;
        } else {
            days.at(2) = 28;
        }

        const double total_year_days = std::accumulate(days.begin(), days.end(), 0);

        int total_days = 0;
        for(int i = 0; i < month; i++) {
            total_days += days.at(i);
        }
        total_days += day;

        const auto out = year + (total_days - 1) / total_year_days;
        return out;
    }

    /** Converts a given Decimal year into a Year, Month and Date it also outputs an error string if there is a problem
     * INPUT  CalendarDate  Pointer to the  data  structure with the following elements
     *            double DecimalYear;      decimal years
     * OUTPUT  CalendarDate  Pointer to the  data  structure with the following elements updated
     * int Year
     * int Month
     * int Day
     *       Error    pointer to an error string
     * CALLS : none
     **/
    int decYear_to_date(const std::string_view edit_date, Date *CalendarDate) {  // TODO
        int MonthDays[13], CumulativeDays = 0;
        int ExtraDay = 0;
        int i, DayOfTheYear;


        if(CalendarDate->DecimalYear < 1900) {
            CalendarDate->Year  = 0;
            CalendarDate->Month = 0;
            CalendarDate->Day   = 0;
            return false;
        }

        CalendarDate->Year = floor(CalendarDate->DecimalYear);


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

    /** This prints WMM errors.
     * INPUT : errorNumer Error look up number
     * OUTPUT: none
     **/
    void PrintError(const int errorNumber) {
        switch(errorNumber) {
            case 1:
                std::cerr << "\nError allocating in LegendreFunctionMemory.\n";
                break;
            case 2:
                std::cerr << "\nError allocating in AllocateModelMemory.\n";
                break;
            case 3:
                std::cerr << "\nError allocating in InitializeGeoid\n";
                break;
            case 4:
                std::cerr << "\nError in setting default values.\n";
                break;
            case 5:
                std::cerr << "\nError initializing Geoid.\n";
                break;
            case 6:
                std::cerr << "\nError opening wmmhr.cof\n.";
                break;
            case 7:
                std::cerr << "\nError opening WMMSV.COF\n.";
                break;
            case 8:
                std::cerr << "\nError reading Magnetic Model.\n";
                break;
            case 9:
                std::cerr << "\nError printing Command Prompt introduction.\n";
                break;
            case 10:
                std::cerr << "\nError converting from geodetic co-ordinates to spherical co-ordinates.\n";
                break;
            case 11:
                std::cerr << "\nError in time modifying the Magnetic model\n";
                break;
            case 12:
                std::cerr << "\nError in Geomagnetic\n";
                break;
            case 13:
                std::cerr << "\nError printing user data\n";
                break;
            case 14:
                std::cerr << "\nError allocating in SummationSpecial\n";
                break;
            case 15:
                std::cerr << "\nError allocating in SecVarSummationSpecial\n";
                break;
            case 16:
                std::cerr << "\nError in opening EGM9615.BIN file\n";
                break;
            case 17:
                std::cerr << "\nError: Latitude OR Longitude out of range in GetGeoidHeight\n";
                break;
            case 18:
                std::cerr << "\nError allocating in PcupHigh\n";
                break;
            case 19:
                std::cerr << "\nError allocating in PcupLow\n";
                break;
            case 20:
                std::cerr << "\nError opening coefficient file\n";
                break;
            case 21:
                std::cerr << "\nError: UnitDepth too large\n";
                break;
            case 22:
                std::cerr << "\nYour system needs Big endian version of EGM9615.BIN.  \n"
                          << "Please download this file from https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml.  \n"
                          << "Replace the existing EGM9615.BIN file with the downloaded one\n";
                break;
            default:
                std::cerr << "\nError: Unknown error\n";
        }
    }

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
        std::string ans{"     "};

        switch(control) {
            case 1: /* Horizontal Field strength low */
                do {
                    std::cerr << "\nCaution: location is approaching the blackout zone around the magnetic pole as\n"
                              << "         defined by the WMM military specification \n"
                              << "         (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass\n"
                              << "         accuracy may be degraded in this region.\n"
                              << "Press enter to continue...\n";
                } while(nullptr == fgets(ans.data(), ans.size(), stdin));
                break;
            case 2: /* Horizontal Field strength very low */
                do {
                    std::cerr << "\nWarning: location is in the blackout zone around the magnetic pole as defined\n"
                              << "         by the WMM military specification \n"
                              << "         (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass\n"
                              << "         accuracy is highly degraded in this region.\n";
                } while(nullptr == fgets(ans.data(), ans.size(), stdin));
                break;
            case 3: /* Elevation outside the recommended range */
                printf("\nWarning: The value you have entered of %.1f km for the elevation is outside of the recommended "
                       "range.\n Elevations above -10.0 km are recommended for accurate results. \n",
                       value);
                while(true) {
                    std::cerr << "\nPlease press 'C' to continue, 'G' to get new data or 'X' to exit...\n";
                    while(nullptr == fgets(ans.data(), ans.size(), stdin)) {
                        std::cerr << "\nInvalid input\n";
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
                std::cerr << "\nWARNING - TIME EXTENDS BEYOND INTENDED USAGE RANGE\n CONTACT NCEI FOR PRODUCT UPDATES:\n"
                          << "	National Centers for Environmental Information\n"
                          << "	NOAA E/NE42\n"
                          << "	325 Broadway\n"
                          << "\n	Boulder, CO 80305 USA"
                          << "	Attn: Manoj Nair or Arnaud Chulliat\n"
                          << "	Phone:	(303) 497-4642 or -6522\n"
                          << "	Email:	geomag.models@noaa.gov\n"
                          << "	Web: https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml\n";
                printf("\n VALID RANGE  = %d - %d\n", (int)magneticModel.min_year,
                       (int)magneticModel.coefficientFileEndDate);
                printf(" TIME   = %f\n", value);
                while(true) {
                    std::cerr << "\nPlease press 'C' to continue, 'N' to enter new data or 'X' to exit...\n";
                    while(nullptr == fgets(ans.data(), ans.size(), stdin)) {
                        std::cerr << "\nInvalid input\n";
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
                    std::cerr << "\nPlease press 'C' to continue, 'G' to get new data or 'X' to exit...\n";
                    while(nullptr == fgets(ans.data(), sizeof(ans), stdin)) {
                        std::cerr << "\nInvalid input\n";
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
}  // namespace wmm