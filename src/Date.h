//
// Created by Professional on 22.04.2026.
//

#pragma once

#include <string>

namespace wmm {
    struct Date {
        Date();
        Date(std::string_view str);

        void calcDecYear();

        int Year{};
        int Month{};
        int Day{};
        double DecimalYear{};  // decimal years
    };

}  // namespace wmm
