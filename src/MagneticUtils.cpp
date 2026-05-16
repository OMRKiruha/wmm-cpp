//
// Created by Kiryuhin Viacheslav on 21.04.2026.
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
            std::cerr << "month out of range\n";
            return false;
        }

        if((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
            days.at(2) = 29;
        } else {
            days.at(2) = 28;
        }

        if(day > days.at(month) || day < 1) {
            std::cerr << "day out of range\n";
            return false;
        }
        return true;
    }

    // Parse the date string format as mm/dd/yyyy
    double dateStr_to_decYear(const std::string_view edit_date) {
        int day{};
        int month{};
        int year{};

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

        const auto out = year + ((total_days - 1) / total_year_days);
        return out;
    }

    /** Converts a given Decimal year into a year, month and Date it also outputs an error string if there is a problem
     */
    int decYear_to_date(const std::string_view edit_date, Date *date) {  // TODO
        int MonthDays[13], CumulativeDays = 0;
        int ExtraDay = 0;
        int i, DayOfTheYear;


        if(date->decimalYear < 1900) {
            date->year  = 0;
            date->month = 0;
            date->day   = 0;
            return false;
        }

        date->year = floor(date->decimalYear);


        if((date->year % 4 == 0 && date->year % 100 != 0) || date->year % 400 == 0) {
            ExtraDay = 1;
        }

        DayOfTheYear = floor(((date->decimalYear - date->year) * (365.0 + (double)ExtraDay)) + 0.5) + 1;
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
                date->month = i;
                date->day   = MonthDays[i] - (CumulativeDays - DayOfTheYear);
                break;
            }
        }

        return true;
    }

    /** @brief Return value 0 means end program, Return value 1 means get new data, Return value 2 means continue.
     * This prints a warning to the screen determined by the control integer. It also takes the value of the parameter
     * causing the warning as a double.  This is unnecessary for some warnings. It requires the MagneticModel to determine
     * the current epoch.
     */
    int warnings(const int control, const double value, const MagneticModel &magneticModel) {
        switch(control) {
            case 1:  // Horizontal Field strength low
                std::cerr << "\nCaution: location is approaching the blackout zone around the magnetic pole as\n"
                          << "         defined by the WMM military specification \n"
                          << "         (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass\n"
                          << "         accuracy may be degraded in this region.\n"
                          << "Press enter to continue...\n";
                std::cin.get();
                break;
            case 2:  // Horizontal Field strength very low
                std::cerr << "\nWarning: location is in the blackout zone around the magnetic pole as defined\n"
                          << "         by the WMM military specification \n"
                          << "         (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass\n"
                          << "         accuracy is highly degraded in this region.\n";
                std::cin.get();
                break;
            case 3:  // Elevation outside the recommended range
                std::cerr << "\nWarning: The value you have entered of %.1f km for the elevation is outside of the "
                             "recommended range.\n Elevations above -10.0 km are recommended for accurate results. \n"
                          << value;
                while(true) {
                    std::cerr << "\nPlease press 'C' to continue, 'G' to get new data or 'X' to exit...\n";
                    switch(const auto c = std::cin.get()) {
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
                            std::cerr << std::format("\nInvalid input {}\n", static_cast<char>(c));
                            break;
                    }
                }
            case 4:  // Date outside the recommended range
                std::cerr << "\nWARNING - TIME EXTENDS BEYOND INTENDED USAGE RANGE\n CONTACT NCEI FOR PRODUCT UPDATES:\n"
                          << "	National Centers for Environmental Information\n"
                          << "	NOAA E/NE42\n"
                          << "	325 Broadway\n"
                          << "\n	Boulder, CO 80305 USA"
                          << "	Attn: Manoj Nair or Arnaud Chulliat\n"
                          << "	Phone:	(303) 497-4642 or -6522\n"
                          << "	Email:	geomag.models@noaa.gov\n"
                          << "	Web: https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml\n"
                          << "\n VALID RANGE  = " << static_cast<int>(magneticModel.min_year) << " - "
                          << static_cast<int>(magneticModel.coefficientFileEndDate) << "\n"
                          << " TIME   = " << value << "\n";
                while(true) {
                    std::cerr << "\nPlease press 'C' to continue, 'N' to enter new data or 'X' to exit...\n";
                    switch(const auto c = std::cin.get()) {
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
                            std::cerr << std::format("\nInvalid input {}\n", c);
                            break;
                    }
                }
            case 5:  // Elevation outside the allowable range
                std::cerr << "\nError: The value you have entered of " << value
                          << " km for the elevation is outside of the recommended range.\n"
                          << " Elevations above -10.0 km are recommended for accurate results. \n";
                while(true) {
                    std::cerr << "\nPlease press 'C' to continue, 'G' to get new data or 'X' to exit...\n";
                    switch(const auto c = std::cin.get()) {
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
                            std::cerr << std::format("\nInvalid input {}\n", c);
                            break;
                    }
                }
            default:
                return 2;
        }
        return 2;
    }
}  // namespace wmm