//
// Created by Professional on 22.04.2026.
//

#pragma once

#include <string>
#include <vector>
#include <cmath>

namespace wmm {
    struct Date;

    struct MagneticModel {
        MagneticModel() = default;
        explicit MagneticModel(int numTerms);
        ~MagneticModel() = default;

        bool readModel(std::string_view filename, int array_size = 1);
        [[nodiscard]] MagneticModel applyDate(const Date &userDate) const;

        double editionDate{};
        double epoch{NAN};  // Base time of Geomagnetic model epoch (yrs)
        double min_year{};
        std::string modelName;
        std::vector<double>
            main_Field_Coeff_G;  // C - Gauss coefficients of main geomagnetic model (nT) Index is (n * (n + 1) / 2 + m)
        std::vector<double> main_Field_Coeff_H;   // C - Gauss coefficients of main geomagnetic model (nT)
        std::vector<double> secular_Var_Coeff_G;  // CD - Gauss coefficients of secular geomagnetic model (nT/yr)
        std::vector<double> secular_Var_Coeff_H;  // CD - Gauss coefficients of secular geomagnetic model (nT/yr)
        int nMax{};                               // Maximum degree of spherical harmonic model
        int nMaxSecVar{};                         // Maximum degree of spherical harmonic secular model
        int secularVariationUsed{};  // Whether or not the magnetic secular variation vector will be needed by program
        double coefficientFileEndDate{};

    private:
        bool readMagneticModelCoefficients(std::ifstream &file);
    };
}  // namespace wmm
