//
// Created by Kiryuhin Viacheslav on 30.04.2026.
//

#pragma once

#include <string>

namespace wmm {
    struct DMS {
        DMS(std::string_view str);

        int degree{};
        int minute{};
        int second{};
        double angle{};

        bool isValid() const;
    };
}  // namespace wmm