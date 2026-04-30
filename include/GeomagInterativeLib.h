#pragma once

#include "GeomagnetismHeader.h"

#define WGS84ON 1
#define MSLON 2

namespace wmm {
    struct Geoid;
    struct CoordGeodetic;
    struct Date;
    struct MagneticModel;
    struct GeoMagneticElements;
    struct Gradient;

    bool ValidateDMSstring(std::string_view input, int min, int max, std::string &Error);
    void GetDeg(std::string_view msg, double *latitude, double min, double max);
    int GetAltitude(std::string_view Query_String, Geoid *geoid, CoordGeodetic *coords, int bounds[2], int AltitudeSetting);
    int GetUserGrid(CoordGeodetic *minimum, CoordGeodetic *maximum, double *step_size, double *a_step_size,
                    double *step_time, Date *StartDate, Date *EndDate, int *ElementOption, int *PrintOption,
                    char *OutputFile, Geoid *Geoid, MagneticModel *model);
    int GetUserInput(const MagneticModel &magneticModel, Geoid *geoid, CoordGeodetic *coordGeodetic, Date *magneticDate);

    void clear_input_buffer();
    void GetMinGridInput(double *coord, double *val_bound, std::string_view var_name);
    void GetMaxGridInput(double *coord, double *val_bound, std::string_view var_name);

    void GetMaxGridInputAlt(double *coord, double min_val);

    void GetMinGridInputDecYear(double *coord, double *val_bound);
    void GetMaxGridInputDecYear(double *coord, double *val_bound);

    void PrintUserDataWithUncertainty(const GeoMagneticElements &GeomagElements, const GeoMagneticElements &Errors,
                                      const CoordGeodetic &SpaceInput, const Date &TimeInput,
                                      const MagneticModel &MagneticModel, const Geoid &Geoid);

    void PrintWMMFormat(char *filename, const MagneticModel &magneticModel);

    void PrintEMMFormat(char *filename, char *filenameSV, const MagneticModel &magneticModel);

    void PrintSHDFFormat(char *filename, MagneticModel *(*MagneticModel)[], int epochs);

    void PrintGradient(const Gradient &gradient);

    void PrintUserData(const GeoMagneticElements &geomagElements, const CoordGeodetic &spaceInput, const Date &timeInput,
                       const MagneticModel &magneticModel, const Geoid &geoid);

}  // namespace wmm