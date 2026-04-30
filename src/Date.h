//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

#include <string>

namespace wmm {
    struct Date {
        Date();
        Date(std::string_view str);

        void calcDecYear();

        int year{};
        int month{};
        int day{};
        double decimalYear{};  // decimal years
    };

}  // namespace wmm
