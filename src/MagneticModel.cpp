//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#include "MagneticModel.h"

#include "Date.h"
#include "MagneticUtils.h"

#include <filesystem>
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
        std::ifstream file{std::filesystem::path{filename}};
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

        const double fileEpoch = std::stod(header.at(0));
        modelName              = header.at(1);
        const std::string edit_date{header.at(2)};

        min_year = dateStr_to_decYear(edit_date);
        if(min_year == -1) {
            min_year = fileEpoch;
        }

        epoch = fileEpoch;

        if(array_size == 1) {
            int _nMax{0};
            int n{0};
            std::string str;
            do {
                if(std::getline(file, str); str.empty()) {
                    break;
                }
                n = std::stoi(str.substr(0, 3));  // Extract first number in line
                if(n > _nMax && (n < 999 && n > 0)) {
                    _nMax = n;
                }
            } while(n < 999);

            const int numTerms = ((_nMax * (_nMax + 1) / 2) + _nMax);
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

    /** @brief Read World Magnetic Model spherical harmonic coefficients (WMM.cof)
     */
    bool MagneticModel::readMagneticModelCoefficients(std::ifstream &file) {
        std::vector<std::string> values;
        std::string str;

        if(!file.is_open()) {
            printError(20);
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
            const int n       = std::stoi(values.at(N));
            const int m       = std::stoi(values.at(M));
            const double gnm  = std::stod(values.at(GNM));
            const double hnm  = std::stod(values.at(HNM));
            const double dgnm = std::stod(values.at(DGNM));
            const double dhnm = std::stod(values.at(DHNM));
            if(m <= n) {
                const int index               = ((n * (n + 1) / 2) + m);
                main_Field_Coeff_G.at(index)  = gnm;
                secular_Var_Coeff_G.at(index) = dgnm;
                main_Field_Coeff_H.at(index)  = hnm;
                secular_Var_Coeff_H.at(index) = dhnm;
            }
        }
        return true;
    }

    /** @brief Time change the Model coefficients from the base year of the model using secular variation coefficients.
     * Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
     * Copy the SV coefficients.  If input "t1" is the same as "t0", then this is merely a copy operation.
     * If the address of "TimedMagneticModel" is the same as the address of "MagneticModel", then this procedure overwrites
     * the given item "MagneticModel".
     */
    MagneticModel MagneticModel::applyDate(const Date &userDate) const {
        MagneticModel timedMagModel;
        timedMagModel.editionDate = editionDate;
        timedMagModel.epoch       = epoch;
        timedMagModel.nMax        = nMax;
        timedMagModel.nMaxSecVar  = nMaxSecVar;
        const int a               = timedMagModel.nMaxSecVar;
        const int b               = ((a * (a + 1) / 2) + a);
        timedMagModel.modelName   = modelName;

        timedMagModel.main_Field_Coeff_G.resize(b + 1);
        timedMagModel.main_Field_Coeff_H.resize(b + 1);
        timedMagModel.secular_Var_Coeff_G.resize(b + 1);
        timedMagModel.secular_Var_Coeff_H.resize(b + 1);

        for(int n = 1; n <= nMax; n++) {
            for(int m = 0; m <= n; m++) {
                if(const int index = ((n * (n + 1) / 2) + m); index <= b) {
                    const auto timeShift = userDate.decimalYear - epoch;
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

    //    int robustReadMagneticModel_Large(const std::string_view filename, char *filenameSV, std::vector<MagneticModel>&
    //    magneticModels) {
    //        char line[MAXLINELENGTH], ModelName[] = "Enhanced Magnetic Model"; /*Model Name must be no longer than 31
    //        characters*/ int n, nMax = 0, nMaxSV = 0, num_terms, a, epochlength = 5, i;
    //
    //        std::ifstream file{filename.data()};
    //        if(!file.is_open()) {
    //            return 0;
    //        }
    //
    //        if (NULL == fgets(line, MAXLINELENGTH, MODELFILE)) {
    //            return 0;
    //        }
    //        do {
    //            if (NULL == fgets(line, MAXLINELENGTH, MODELFILE))
    //                break;
    //            a = sscanf(line, "%d", &n);
    //            if (n > nMax && (n < 99999 && a == 1 && n > 0))
    //                nMax = n;
    //        } while (n < 99999 && a == 1);
    //        fclose(MODELFILE);
    //        MODELFILE = fopen(filenameSV, "r");
    //        if (MODELFILE == 0) {
    //            return 0;
    //        }
    //        n = 0;
    //        if (NULL == fgets(line, MAXLINELENGTH, MODELFILE))
    //            return 0;
    //        do {
    //            if (NULL == fgets(line, MAXLINELENGTH, MODELFILE))
    //                break;
    //            a = sscanf(line, "%d", &n);
    //            if (n > nMaxSV && (n < 99999 && a == 1 && n > 0))
    //                nMaxSV = n;
    //        } while (n < 99999 && a == 1);
    //        fclose(MODELFILE);
    //        num_terms = CALCULATE_NUMTERMS(nMax);
    //        *magneticModel = AllocateModelMemory(num_terms);
    //        (*magneticModel)->nMax = nMax;
    //        (*magneticModel)->nMaxSecVar = nMaxSV;
    //        if (nMaxSV > 0) (*magneticModel)->SecularVariationUsed = true;
    //        for (i = 0; i < num_terms; i++) {
    //            (*magneticModel)->Main_Field_Coeff_G[i] = 0;
    //            (*magneticModel)->Main_Field_Coeff_H[i] = 0;
    //            (*magneticModel)->Secular_Var_Coeff_G[i] = 0;
    //            (*magneticModel)->Secular_Var_Coeff_H[i] = 0;
    //        }
    //        readMagneticModel_Large(filename, filenameSV, *magneticModel);
    //        (*magneticModel)->CoefficientFileEndDate = (*magneticModel)->epoch + epochlength;
    //        strlcpy_equivalent((*magneticModel)->ModelName, ModelName, sizeof((*magneticModel)->ModelName));
    //        magneticModel.EditionDate = (*magneticModel)->epoch;
    //        return 1;
    //    } /*robustReadMagneticModel_Large*/

    /*End of Wrapper Functions*/

    /******************************************************************************
     ********************************Memory and File Processing********************
     * This grouping consists of functions that read coefficient files into the
     * memory, allocate memory, free memory or print models into coefficient files.
     ******************************************************************************/

    /**  To read the high-degree model coefficients (for example, NGDC 720)
     * INPUT :  filename   file name for static coefficients
     *                    filenameSV file name for secular variation coefficients
     *
     *    MagneticModel : Pointer to the data structure with the following fields required as inputs
     *                        nMaxSecVar : Number of secular variation coefficients
     *                                         nMax : 	Number of static coefficients
     *                                                    UPDATES : MagneticModel : Pointer to the data structure with the
     * following fields populated double epoch;       Base time of Geomagnetic model epoch (yrs) double *Main_Field_Coeff_G;
     * C - Gauss coefficients of main geomagnetic model (nT) double *Main_Field_Coeff_H;          C - Gauss coefficients of
     * main geomagnetic model (nT) double *Secular_Var_Coeff_G;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
     *             double *Secular_Var_Coeff_H;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
     *                 CALLS : none
     */
    //    int readMagneticModel_Large(char *filename, char *filenameSV, MagneticModel *magneticModel){
    //        FILE *COF_File;
    //        FILE *COFSV_File;
    //        char c_str[81], c_str2[81]; /* these strings are used to read a line from coefficient file */
    //        int i, m, n, index, a, b;
    //        double epoch, gnm, hnm, dgnm, dhnm;
    //        COF_File   = fopen(filename, "r");
    //        COFSV_File = fopen(filenameSV, "r");
    //        if(COF_File == NULL || COFSV_File == NULL) {
    //            printError(20);
    //            return false;
    //        }
    //        magneticModel->Main_Field_Coeff_H[0]  = 0.0;
    //        magneticModel->Main_Field_Coeff_G[0]  = 0.0;
    //        magneticModel->Secular_Var_Coeff_H[0] = 0.0;
    //        magneticModel->Secular_Var_Coeff_G[0] = 0.0;
    //        if(NULL == fgets(c_str, sizeof(c_str), COF_File)) {
    //            fclose(COF_File);
    //            fclose(COFSV_File);
    //            return false;
    //        }
    //        snprintf(c_str, sizeof(c_str), "%lf%s", epoch, magneticModel->ModelName);
    //        magneticModel->epoch = epoch;
    //        a                    = CALCULATE_NUMTERMS(magneticModel->nMaxSecVar);
    //        b                    = CALCULATE_NUMTERMS(magneticModel->nMax);
    //        for(i = 0; i < a; i++) {
    //            if(NULL == fgets(c_str, sizeof(c_str), COF_File)) {
    //                fclose(COF_File);
    //                fclose(COFSV_File);
    //                return false;
    //            }
    //            sscanf(c_str, "%d%d%lf%lf", &n, &m, &gnm, &hnm);
    //            if(NULL == fgets(c_str2, sizeof(c_str2), COFSV_File)) {
    //                fclose(COF_File);
    //                fclose(COFSV_File);
    //                return false;
    //            }
    //            sscanf(c_str2, "%d%d%lf%lf", &n, &m, &dgnm, &dhnm);
    //            if(m <= n) {
    //                index                                     = (n * (n + 1) / 2 + m);
    //                magneticModel->Main_Field_Coeff_G[index]  = gnm;
    //                magneticModel->Secular_Var_Coeff_G[index] = dgnm;
    //                magneticModel->Main_Field_Coeff_H[index]  = hnm;
    //                magneticModel->Secular_Var_Coeff_H[index] = dhnm;
    //            }
    //        }
    //        for(i = a; i < b; i++) {
    //            if(NULL == fgets(c_str, sizeof(c_str), COF_File)) {
    //                fclose(COF_File);
    //                fclose(COFSV_File);
    //                return false;
    //            }
    //            sscanf(c_str, "%d%d%lf%lf", &n, &m, &gnm, &hnm);
    //            if(m <= n) {
    //                index                                    = (n * (n + 1) / 2 + m);
    //                magneticModel->Main_Field_Coeff_G[index] = gnm;
    //                magneticModel->Main_Field_Coeff_H[index] = hnm;
    //            }
    //        }
    //        if(COF_File != NULL && COFSV_File != NULL) {
    //            fclose(COF_File);
    //            fclose(COFSV_File);
    //        }
    //
    //        return true;
    //    } /*readMagneticModel_Large*/


    /** readMagneticModels - Read the Magnetic Models from an SHDF format file
     *
     * Input:
     *  filename - Path to the SHDF format model file to be read
     *  array_size - Max No of models to be read from the file
     *
     * Output:
     *  magneticmodels[] - Array of magnetic models read from the file
     *
     * Return value:
     *  Returns the number of models read from the file.
     *  -2 implies that internal or external static degree was not found in the file, hence memory cannot be allocated
     *  -1 implies some error during file processing (I/O)
     *  0 implies no models were read from the file
     *  if ReturnValue > array_size then there were too many models in model file but only <array_size> number were read .
     *  if ReturnValue <= array_size then the function execution was successful.
     */
    //    int readMagneticModel_SHDF(char *filename, MagneticModel *(*magneticmodels)[], int array_size)    {
    //        char paramkeys[NOOFPARAMS][MAXLINELENGTH] = {
    //            "SHDF ",          "ModelName: ",     "Publisher: ",    "ReleaseDate: ",  "DataCutOff: ",   "ModelStartYear:
    //            ", "ModelEndYear: ", "Epoch: ",         "IntStaticDeg: ", "IntSecVarDeg: ", "ExtStaticDeg: ",
    //            "ExtSecVarDeg: ", "GeoMagRefRad: ", "Normalization: ", "SpatBasFunc: "};
    //
    //        char paramvalues[NOOFPARAMS][MAXLINELENGTH];
    //        char *line = (char *)malloc(MAXLINELENGTH);
    //        char *ptrreset;
    //        char paramvalue[MAXLINELENGTH];
    //        int paramvaluelength = 0;
    //        int paramkeylength   = 0;
    //        int i = 0, j = 0;
    //        int newrecord    = 1;
    //        int header_index = -1;
    //        int numterms;
    //        int tempint;
    //        int allocationflag = 0;
    //        char coefftype; /* Internal or External (I/E) */
    //
    //        /* For reading coefficients */
    //        int n, m;
    //        double gnm, hnm, dgnm, dhnm, cutoff;
    //        int index;
    //
    //        FILE *stream;
    //        ptrreset = line;
    //        stream   = fopen(filename, READONLYMODE);
    //        if(stream == NULL) {
    //            perror("File open error");
    //            return header_index;
    //        }
    //
    //        /* Read records from the model file and store header information. */
    //        while(fgets(line, MAXLINELENGTH, stream) != NULL) {
    //            j++;
    //            if(strlen(Trim(line)) == 0) {
    //                continue;
    //            }
    //            if(*line == '%') {
    //                line++;
    //                if(newrecord) {
    //                    if(header_index > -1) {
    //                        AssignHeaderValues((*magneticmodels)[header_index], paramvalues);
    //                    }
    //                    header_index++;
    //                    if(header_index >= array_size) {
    //                        fprintf(stderr, "Header limit exceeded - too many models in model file. (%d)\n", header_index);
    //                        return array_size + 1;
    //                    }
    //                    newrecord      = 0;
    //                    allocationflag = 0;
    //                }
    //                for(i = 0; i < NOOFPARAMS; i++) {
    //                    paramkeylength = strlen(paramkeys[i]);
    //                    if(!strncmp(line, paramkeys[i], paramkeylength)) {
    //                        paramvaluelength = strlen(line) - paramkeylength;
    //                        memset(paramvalues, '\0', paramvaluelength);
    //                        strlcpy_equivalent(paramvalue, line + paramkeylength, paramvaluelength);
    //                        paramvalue[paramvaluelength] = '\0';
    //                        strlcpy_equivalent(paramvalues[i], paramvalue, 1);
    //                        if(!strcmp(paramkeys[i], paramkeys[INTSTATICDEG]) ||
    //                           !strcmp(paramkeys[i], paramkeys[EXTSTATICDEG])) {
    //                            tempint = atoi(paramvalues[i]);
    //                            if(tempint > 0 && allocationflag == 0) {
    //                                numterms                        = CALCULATE_NUMTERMS(tempint);
    //                                (*magneticmodels)[header_index] = AllocateModelMemory(numterms);
    //                                /* model = (*magneticmodels)[header_index]; */
    //                                allocationflag = 1;
    //                            }
    //                        }
    //                        break;
    //                    }
    //                }
    //                line--;
    //            } else if(*line == '#') {
    //                /* process comments */
    //
    //            } else if(sscanf(line, "%c,%d,%d", &coefftype, &n, &m) == 3) {
    //                if(m == 0) {
    //                    sscanf(line, "%c,%d,%d,%lf,,%lf,", &coefftype, &n, &m, &gnm, &dgnm);
    //                    hnm  = 0;
    //                    dhnm = 0;
    //                } else {
    //                    sscanf(line, "%c,%d,%d,%lf,%lf,%lf,%lf", &coefftype, &n, &m, &gnm, &hnm, &dgnm, &dhnm);
    //                }
    //                newrecord = 1;
    //                if(!allocationflag) {
    //                    fprintf(stderr, "Degree not found in model. Memory cannot be allocated.\n");
    //                    return _DEGREE_NOT_FOUND;
    //                }
    //                if(m <= n) {
    //                    index                                                       = (n * (n + 1) / 2 + m);
    //                    (*magneticmodels)[header_index]->Main_Field_Coeff_G[index]  = gnm;
    //                    (*magneticmodels)[header_index]->Secular_Var_Coeff_G[index] = dgnm;
    //                    (*magneticmodels)[header_index]->Main_Field_Coeff_H[index]  = hnm;
    //                    (*magneticmodels)[header_index]->Secular_Var_Coeff_H[index] = dhnm;
    //                }
    //            }
    //        }
    //        if(header_index > -1) {
    //            AssignHeaderValues((*magneticmodels)[header_index], paramvalues);
    //        }
    //        fclose(stream);
    //
    //        cutoff = (*magneticmodels)[array_size - 1]->CoefficientFileEndDate;
    //
    //        for(i = 0; i < array_size; i++) {
    //            (*magneticmodels)[i]->CoefficientFileEndDate = cutoff;
    //        }
    //
    //        free(ptrreset);
    //        line     = NULL;
    //        ptrreset = NULL;
    //        return header_index + 1;
    //    } /*readMagneticModel_SHDF*/

    //    char *Trim(char *str) {
    //        char *end;
    //
    //        while(isspace(*str)) {
    //            str++;
    //        }
    //
    //        if(*str == 0) {
    //            return str;
    //        }
    //
    //        end = str + strlen(str) - 1;
    //        while(end > str && isspace(*end)) {
    //            end--;
    //        }
    //
    //        *(end + 1) = 0;
    //
    //        return str;
    //    }
}  // namespace wmm