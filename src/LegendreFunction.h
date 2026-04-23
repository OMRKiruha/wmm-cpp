//
// Created by Professional on 22.04.2026.
//

#pragma once

#include <vector>

namespace wmm {
    struct CoordSpherical;

    struct LegendreFunction {
        LegendreFunction() = delete;
        LegendreFunction(const CoordSpherical &coordSpherical, int nMax);
        ~LegendreFunction() = default;

        int nMax{};
        std::vector<double> Pcup;   // Legendre Function
        std::vector<double> dPcup;  // Derivative of Legendre fcn

    private:
        void PcupHigh(double x);
        void PcupLow(double x);
    };
}  // namespace wmm