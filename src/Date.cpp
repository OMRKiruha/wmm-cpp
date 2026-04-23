//
// Created by Professional on 22.04.2026.
//

#include "Date.h"

#include "MagneticUtils.h"

#include <chrono>

namespace wmm{
    Date::Date() {
        const std::time_t t = std::time(nullptr);
        const std::tm *now  = std::localtime(&t);
        Year                = now->tm_year + 1900;
        Month               = now->tm_mon + 1;
        Day                 = now->tm_mday;

        calcDecYear();
    }

    Date::Date(std::string_view str) {
        if(dateStr_to_ymd(str, Year, Month, Day)) {
            calcDecYear();
        } else {
            Year        = 0;
            Month       = 0;
            Day         = 0;
            DecimalYear = 0;
        }
    }

    void Date::calcDecYear() {
        DecimalYear = date_to_decYear(Year, Month, Day);
    }
}