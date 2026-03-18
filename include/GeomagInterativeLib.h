#pragma once

#include "GeomagnetismHeader.h"

namespace wmm {

    int ValidateDMSstring(const std::string_view input, int min, int max, std::string *Error);
    void GetDeg(std::string_view Query_String, double *latitude, double bounds[2]);
    int GetAltitude(std::string_view Query_String, Geoid *geoid, CoordGeodetic *coords, int bounds[2], int AltitudeSetting);
    int GetUserGrid(CoordGeodetic *minimum, CoordGeodetic *maximum, double *step_size, double *a_step_size,
                        double *step_time, Date *StartDate, Date *EndDate, int *ElementOption, int *PrintOption,
                        char *OutputFile, Geoid *Geoid, MagneticModel *model);
    int GetUserInput(MagneticModel *magneticModel, Geoid *geoid, CoordGeodetic *coordGeodetic, Date *magneticDate);

    void clear_input_buffer();
    void GetMinGridInput(double *coord, double *val_bound, char *var_name);
    void GetMaxGridInput(double *coord, double *val_bound, char *var_name);

    void GetMaxGridInputAlt(double *coord, double min_val);

    void GetMinGridInputDecYear(double *coord, double *val_bound);
    void GetMaxGridInputDecYear(double *coord, double *val_bound);

}  // namespace wmm