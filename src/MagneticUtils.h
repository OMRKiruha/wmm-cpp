//
// Created by Kiryuhin Viacheslav on 21.04.2026.
//

#pragma once

#include <cmath>
#include <string>

namespace wmm {

    struct MagneticModel;

    static double rad2Deg(const double rad) {
        return rad * (180.0 / M_PI);
    }

    static double deg2Rad(const double deg) {
        return deg * (M_PI / 180.0);
    }

    bool dateStr_to_ymd(std::string_view str, int &year, int &month, int &day);

    double dateStr_to_decYear(std::string_view edit_date);

    double date_to_decYear(int year, int month, int day);

    // User Interface
    void printError(int errorNumber);

    int warnings(int control, double value, const MagneticModel &magneticModel);

}  // namespace wmm