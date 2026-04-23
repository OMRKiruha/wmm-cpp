//
// Created by Professional on 22.04.2026.
//

#include "MagneticModel.h"

#include "Date.h"
#include "MagneticUtils.h"

#include <fstream>

namespace wmm {
    enum COEFFICIENTS : uint8_t { N = 0, M, GNM, HNM, DGNM, DHNM };

    MagneticModel::MagneticModel(const int numTerms) {
        main_Field_Coeff_G.resize(numTerms + 1);
        main_Field_Coeff_H.resize(numTerms + 1);
        secular_Var_Coeff_G.resize(numTerms + 1);
        secular_Var_Coeff_H.resize(numTerms + 1);
    }

    bool MagneticModel::readModel(std::string_view filename, int array_size) {
        std::string str;

        std::ifstream file{filename.data()};
        if(!file.is_open()) {
            return false;
        }

        // Read first line
        std::vector<std::string> header;
        for(std::string line; std::getline(file, line, ' ');) {
            if(!line.empty()) {
                auto &val = header.emplace_back(line);
                if(val.back() == '\n') {
                    val.pop_back();
                    break;
                }
            }
        }

        // First string in file is not empty
        if(header.empty()) {
            return false;
        }

        double fileEpoch = std::stod(header.at(0));
        modelName        = header.at(1);
        std::string edit_date{header.at(2)};

        min_year = dateStr_to_decYear(edit_date);
        if(min_year == -1) {
            min_year = fileEpoch;
        }

        epoch = fileEpoch;


        if(array_size == 1) {
            int _nMax{0};
            int n{0};
            do {
                if(std::getline(file, str); str.empty()) {
                    break;
                }
                n = std::stoi(str.substr(0, 3));  // Extract first number in line
                if(n > _nMax && (n < 999 && n > 0)) {
                    _nMax = n;
                }
            } while(n < 999);

            const int numTerms = (_nMax * (_nMax + 1) / 2 + _nMax);
            main_Field_Coeff_G.resize(numTerms + 1);
            main_Field_Coeff_H.resize(numTerms + 1);
            secular_Var_Coeff_G.resize(numTerms + 1);
            secular_Var_Coeff_H.resize(numTerms + 1);
            nMax       = _nMax;
            nMaxSecVar = _nMax;
            readMagneticModelCoefficients(file);
            coefficientFileEndDate = epoch + 5;

        } else {
            return false;
        }
        return true;
    }

    /** Read World Magnetic Model spherical harmonic coefficients (WMM.cof)
     *  INPUT :   filename
     *  UPDATES : MagneticModel
     *  CALLS : none */
    bool MagneticModel::readMagneticModelCoefficients(std::ifstream &file) {
        std::vector<std::string> values;
        std::string str;

        if(!file.is_open()) {
            PrintError(20);
            return false;
        }

        file.seekg(0, std::ios::beg);
        std::getline(file, str);

        main_Field_Coeff_H.at(0)  = 0.0;
        main_Field_Coeff_G.at(0)  = 0.0;
        secular_Var_Coeff_H.at(0) = 0.0;
        secular_Var_Coeff_G.at(0) = 0.0;


        bool isEnd = false;
        while(!isEnd) {
            values.clear();

            for(std::string line; std::getline(file, line, ' ');) {
                if(!line.empty()) {
                    if(line.find('\n') != std::string::npos) {
                        isEnd = line.find("999999999") != std::string::npos;
                        line  = line.substr(0, line.find('\n'));
                        values.emplace_back(line);
                        break;
                    }
                    values.emplace_back(line);
                }
            }

            /* END OF FILE NOT ENCOUNTERED, GET VALUES */
            int n       = std::stoi(values.at(N));
            int m       = std::stoi(values.at(M));
            double gnm  = std::stod(values.at(GNM));
            double hnm  = std::stod(values.at(HNM));
            double dgnm = std::stod(values.at(DGNM));
            double dhnm = std::stod(values.at(DHNM));
            if(m <= n) {
                int index                     = (n * (n + 1) / 2 + m);
                main_Field_Coeff_G.at(index)  = gnm;
                secular_Var_Coeff_G.at(index) = dgnm;
                main_Field_Coeff_H.at(index)  = hnm;
                secular_Var_Coeff_H.at(index) = dhnm;
            }
        }
        return true;
    }

    /**
     * Time change the Model coefficients from the base year of the model using secular variation coefficients.
     * Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
     * Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.
     * If the address of "TimedMagneticModel" is the same as the address of "MagneticModel", then this procedure overwrites
     * the given item "MagneticModel".
     * INPUT:  UserDate
     * OUTPUT: TimedMagneticModel
     */
    MagneticModel MagneticModel::applyDate(const Date &userDate) const {
        MagneticModel timedMagModel;
        timedMagModel.editionDate = editionDate;
        timedMagModel.epoch       = epoch;
        timedMagModel.nMax        = nMax;
        timedMagModel.nMaxSecVar  = nMaxSecVar;
        const int a               = timedMagModel.nMaxSecVar;
        const int b               = (a * (a + 1) / 2 + a);
        timedMagModel.modelName   = modelName;

        timedMagModel.main_Field_Coeff_G.resize(b + 1);
        timedMagModel.main_Field_Coeff_H.resize(b + 1);
        timedMagModel.secular_Var_Coeff_G.resize(b + 1);
        timedMagModel.secular_Var_Coeff_H.resize(b + 1);

        for(int n = 1; n <= nMax; n++) {
            for(int m = 0; m <= n; m++) {
                if(const int index = (n * (n + 1) / 2 + m); index <= b) {
                    const auto timeShift = userDate.DecimalYear - epoch;
                    timedMagModel.main_Field_Coeff_H.at(index) =
                        main_Field_Coeff_H.at(index) + timeShift * secular_Var_Coeff_H.at(index);
                    timedMagModel.main_Field_Coeff_G.at(index) =
                        main_Field_Coeff_G.at(index) + timeShift * secular_Var_Coeff_G.at(index);
                    // We need a copy of the secular var coef to calculate secular change
                    timedMagModel.secular_Var_Coeff_H.at(index) = secular_Var_Coeff_H.at(index);
                    timedMagModel.secular_Var_Coeff_G.at(index) = secular_Var_Coeff_G.at(index);
                } else {
                    timedMagModel.main_Field_Coeff_H.at(index) = main_Field_Coeff_H.at(index);
                    timedMagModel.main_Field_Coeff_G.at(index) = main_Field_Coeff_G.at(index);
                }
            }
        }
        return timedMagModel;
    }
}  // namespace wmm