//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "Date.h"

#include "MagneticUtils.h"

#include <chrono>

namespace wmm{
    Date::Date() {
        const std::time_t t = std::time(nullptr);
        const std::tm *now  = std::localtime(&t);
        year                = now->tm_year + 1900;
        month               = now->tm_mon + 1;
        day                 = now->tm_mday;

        calcDecYear();
    }

    Date::Date(std::string_view str) {
        if(dateStr_to_ymd(str, year, month, day)) {
            calcDecYear();
        } else {
            year        = 0;
            month       = 0;
            day         = 0;
            decimalYear = 0;
        }
    }

    void Date::calcDecYear() {
        decimalYear = date_to_decYear(year, month, day);
    }
}