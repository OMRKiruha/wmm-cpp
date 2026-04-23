
#include "GeomagnetismHeader.h"
#include "MagneticUtils.h"

#include <cmath>
#include <cstdio>

/* $Id: GeomagnetismLibrary.c 1521 2017-01-24 17:52:41Z awoods $
 *
 * ABSTRACT
 *
 * The purpose of Geomagnetism Library is primarily to support the World Magnetic Model (WMM) 2015-2020.
 * It however is built to be used for spherical harmonic models of the Earth's magnetic field
 * generally and supports models even with a large (>>12) number of degrees.  It is also used in many
 * other geomagnetic models distributed by NCEI.
 *
 * REUSE NOTES
 *
 * Geomagnetism Library is intended for reuse by any application that requires
 * Computation of Geomagnetic field from a spherical harmonic model.
 *
 * REFERENCES
 *
 *    Further information on Geoid can be found in the WMM Technical Documents.
 *
 *
 * LICENSES
 *
 *  The WMM source code is in the public domain and not licensed or under copyright.
 *	The information and software may be used freely by the public. As required by 17 U.S.C. 403,
 *	third parties producing copyrighted works consisting predominantly of the material produced by
 *	U.S. government agencies must provide notice with such work(s) identifying the U.S. Government material
 *	incorporated and stating that such material is not subject to copyright protection.
 *
 * RESTRICTIONS
 *
 *    Geomagnetism library has no restrictions.
 *
 * ENVIRONMENT
 *
 *    Geomagnetism library was tested in the following environments
 *
 *    1. Red Hat Linux  with GCC Compiler
 *    2. MS Windows 7 with MinGW compiler
 *    3. Sun Solaris with GCC Compiler
 *
 *


 *  National Centers for Environmental Information
 *  NOAA E/NE42, 325 Broadway
 *  Boulder, CO 80305 USA
 *  Attn: Arnaud Chulliat
 *  Phone:  (303) 497-6522
 *  Email:  Arnaud.Chulliat@noaa.gov

 *  Software and Model Support
 *  National Centers for Environmental Information
 *  NOAA E/NE42
 *  325 Broadway
 *  Boulder, CO 80305 USA
 *  Attn: Adam Woods or Manoj Nair
 *  Phone:  (303) 497-6640 or -4642
 *  Email:  geomag.models@noaa.gov
 *  URL: http://www.ngdc.noaa.gov/Geomagnetic/WMM/DoDWMM.shtml


 *  For more details on the subroutines, please consult the WMM
 *  Technical Documentations at
 *  http://www.ngdc.noaa.gov/Geomagnetic/WMM/DoDWMM.shtml

 *  Nov 23, 2009
 *  Written by Manoj C Nair and Adam Woods
 *  Manoj.C.Nair@noaa.Gov
 *  Adam.Woods@noaa.gov
 */

namespace wmm {

    /******************************************************************************
     ************************************Wrapper***********************************
     * This grouping consists of functions call groups of other functions to do a
     * complete calculation of some sort.  For example, the geomag function
     * does everything necessary to compute the geomagnetic elements from a given
     * geodetic point in space and magnetic model adjusted for the appropriate
     * date. These functions are the external functions necessary to create a
     * program that uses or calculates the magnetic field.
     ******************************************************************************
     ******************************************************************************/

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
    //            PrintError(20);
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

    /*End of Memory and File Processing functions*/


    /******************************************************************************
     *************Conversions, Transformations, and other Calculations**************
     * This grouping consists of functions that perform unit conversions, coordinate
     * transformations and other simple or straightforward calculations that are
     * usually easily replicable with a typical scientific calculator.
     ******************************************************************************/


    void BaseErrors(double DeclCoef, double DeclBaseline, double InclOffset, double FOffset, double Multiplier, double H,
                    double *DeclErr, double *InclErr, double *FErr) {
        double declHorizontalAdjustmentSq;
        declHorizontalAdjustmentSq = (DeclCoef / H) * (DeclCoef / H);
        *DeclErr                   = sqrt(declHorizontalAdjustmentSq + DeclBaseline * DeclBaseline) * Multiplier;
        *InclErr                   = InclOffset * Multiplier;
        *FErr                      = FOffset * Multiplier;
    }

    /** This converts a given decimal degree into a DMS string.
     * INPUT  DegreesOfArc   decimal degree
     *           UnitDepth	How many iterations should be printed,
     *                        1 = Degrees
     *                        2 = Degrees, Minutes
     *                        3 = Degrees, Minutes, Seconds
     * OUPUT  DMSstring 	 pointer to DMSString.  Must be at least 30 characters.
     * CALLS : none
     */
    void DegreeToDMSstring(double DegreesOfArc, int UnitDepth, std::string &out) {
        int DMS;
        double temp = DegreesOfArc;

        if(UnitDepth > 3) {
            PrintError(21);
        }

        for(int i = 0; i < UnitDepth; i++) {
            DMS  = static_cast<int>(temp);
            temp = (temp - DMS) * 60;

            if(i == UnitDepth - 1 && temp >= 30) {
                DMS++;
            } else if(i == UnitDepth - 1 && temp <= -30) {
                DMS--;
            }

            out.append(std::to_string(DMS));

            switch(i) {
                case 0:
                    out.append(" Deg ");
                    break;
                case 1:
                    out.append(" Min ");
                    break;
                case 2:
                    out.append(" Sec ");
                    break;
                default:;
            }
        }
    }

    /** This converts a given DMS string into decimal degrees.
     * INPUT  DMSstring 	 pointer to DMSString
     * OUTPUT  DegreesOfArc   decimal degree
     * CALLS : none
     */
    void DMSstringToDegree(std::string_view DMSstring, double *DegreesOfArc) {
        int second, minute, degree, sign = 1, j = 0;
        j = sscanf(DMSstring.data(), "%d, %d, %d", &degree, &minute, &second);
        if(j != 3) {
            sscanf(DMSstring.data(), "%d %d %d", &degree, &minute, &second);
        }
        if(degree < 0) {
            sign = -1;
        }
        degree        = degree * sign;
        *DegreesOfArc = sign * (degree + minute / 60.0 + second / 3600.0);
    } /*DMSstringToDegree*/

}  // namespace wmm