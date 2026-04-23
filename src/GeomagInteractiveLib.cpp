
#include <cmath>
#include <cstdio>
#include <string>

#include "CoordGeodetic.h"
#include "Date.h"
#include "GeoMagneticElements.h"
#include "Geoid.h"
#include "GeomagInterativeLib.h"
#include "Gradient.h"
#include "MagneticModel.h"
#include "MagneticUtils.h"

namespace wmm {

    void clear_input_buffer() {
        int c;
        while((c = getchar()) != '\n' && c != EOF) {
            // Simply read and discard each character
        }
    }

    /* Validates a latitude DMS string, and returns 1 for a success and returns 0 for a failure.
     * It copies an error message to the Error string in the event of a failure.
     *
     * INPUT : input (DMS string)
     * OUTPUT : Error : Error string
     * CALLS : none
     */
    int ValidateDMSstring(const std::string_view input, int min, int max, std::string *Error) {
        int degree, minute, second, j = 0, n, max_minute = 60, max_second = 60;
        int i;
        degree         = -1000;
        minute         = -1;
        second         = -1;
        int Error_size = 255;

        for(i = 0; i <= input.size() - 1; i++) /*tests for legal characters*/
        {
            if((input[i] < '0' || input[i] > '9') &&
               (input[i] != ',' && input[i] != ' ' && input[i] != '-' && input[i] != '\0' && input[i] != '\n')) {
                //  The Error is passed as char pointer from GetDeg(), so I can't use "sizeof()" to estimate its size.
                //  The size of Error is 255 and defined in GetDeg(). GetDeg() is the only function which will call
                //  ValidateDMSstring().
                *Error = "\nError: Input contains an illegal character, legal characters for Degree, Minute, Second format "
                         "are:\n '0-9' ',' '-' '[space]' '[Enter]'\n";
                return false;
            }
            if(input[i] == ',') {
                j++;
            }
        }
        if(j == 2) {
            j = sscanf(input.data(), "%d, %d, %d", &degree, &minute, &second); /*tests for legal formatting and range*/
        } else {
            j = sscanf(input.data(), "%d %d %d", &degree, &minute, &second);
        }
        if(j == 1) {
            minute = 0;
            second = 0;
            j      = 3;
        }
        if(j != 3) {
            *Error = "\nError: Not enough numbers used for Degrees, Minutes, Seconds format\n or they were "
                     "incorrectly formatted\n The legal format is DD,MM,SS or DD MM SS\n";
            return false;
        }
        if(degree > max || degree < min) {
            sprintf(Error->data(), "\nError: Degree input is outside legal range\n The legal range is from %d to %d\n", min,
                    max);
            return false;
        }
        if(degree == max || degree == min) {
            max_minute = 0;
        }
        if(minute > max_minute || minute < 0) {
            *Error = "\nError: Minute input is outside legal range\n The legal minute range is from 0 to 60\n";
            return false;
        }
        if(minute == max_minute) {
            max_second = 0;
        }
        if(second > max_second || second < 0) {
            *Error = "\nError: Second input is outside legal range\n The legal second range is from 0 to 60\n";
            return false;
        }
        return true;
    } /*ValidateDMSstring*/

    void GetDeg(std::string_view Query_String, double *latitude, double bounds[2]) {
        /*Gets a degree value from the user using the standard input*/
        std::string buffer;
        buffer.reserve(65);
        std::string Error_Message;
        int done, i, j;

        printf("%s", Query_String.data());
        while(nullptr == fgets(buffer.data(), 64, stdin)) {
            printf("%s", Query_String.data());
            if(buffer[64 - 1] != '\n') {
                clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read the
                                       // whole stdin
            }
        }

        for(i = 0, done = 0, j = 0; i < (int)64 && !done; i++) {
            if(buffer[i] == '.') {
                j = sscanf(buffer.c_str(), "%lf", latitude);
                if(j == 1) {
                    done = 1;
                } else {
                    done = -1;
                }
            }
            if(buffer[i] == ',') {
                if(ValidateDMSstring(buffer, bounds[0], bounds[1], &Error_Message)) {
                    DMSstringToDegree(buffer, latitude);
                    done = 1;
                } else {
                    done = -1;
                }
            }
            if(buffer[i] == ' ') /* This detects if there is a ' ' somewhere in the string,
             if there is the program tries to interpret the input as Degrees Minutes Seconds.*/
            {
                if(ValidateDMSstring(buffer, bounds[0], bounds[1], &Error_Message)) {
                    DMSstringToDegree(buffer, latitude);
                    done = 1;
                } else {
                    done = -1;
                }
            }
            if(buffer[i] == '\0' || done == -1) {
                if(ValidateDMSstring(buffer, bounds[0], bounds[1], &Error_Message) && done != -1) {
                    sscanf(buffer.c_str(), "%lf", latitude);
                    done = 1;
                } else {
                    printf("%s", &Error_Message);
                    buffer.clear();
                    printf("\nError encountered, please re-enter as '(-)DDD,MM,SS' or in Decimal Degrees DD.ddd:\n");
                    while(NULL == fgets(buffer.data(), 64, stdin)) {
                        printf("\nError encountered, please re-enter as '(-)DDD,MM,SS' or in Decimal Degrees DD.ddd:\n");
                        if(buffer[64 - 1] != '\n') {
                            clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to
                                                   // read the whole stdin
                        }
                    }

                    i    = -1;
                    done = 0;
                }
            }
        }
    }

    int GetAltitude(std::string_view Query_String, Geoid *geoid, CoordGeodetic *coords, int bounds[2], int AltitudeSetting) {
        int done, j, UpBoundOn;
        char tmp;
        std::string buffer;
        double value;
        done = 0;
        if(bounds[1] != NO_ALT_MAX) {
            UpBoundOn = true;
        } else {
            UpBoundOn = false;
        }
        printf("%s", Query_String.data());

        while(!done) {
            buffer.clear();
            while(NULL == fgets(buffer.data(), 64, stdin)) {
                printf("%s", Query_String.data());
            }
            j = 0;
            if((AltitudeSetting != MSLON) &&
               (buffer[0] == 'e' || buffer[0] == 'E' ||
                AltitudeSetting == WGS84ON)) /* User entered height above WGS-84 ellipsoid, copy it to
                                                CoordGeodetic->HeightAboveEllipsoid */
            {
                if(buffer[0] == 'e' || buffer[0] == 'E') {
                    j = sscanf(buffer.c_str(), "%c%lf", &tmp, &coords->HeightAboveEllipsoid);
                } else {
                    j = sscanf(buffer.c_str(), "%lf", &coords->HeightAboveEllipsoid);
                }
                if(j == 2) {
                    j = 1;
                }
                geoid->UseGeoid          = 0;
                coords->HeightAboveGeoid = coords->HeightAboveEllipsoid;
                value                    = coords->HeightAboveEllipsoid;
            } else /* User entered height above MSL, convert it to the height above WGS-84 ellipsoid */
            {
                geoid->UseGeoid = 1;
                j               = sscanf(buffer.c_str(), "%lf", &coords->HeightAboveGeoid);
                coords->convertGeoidToEllipsoidHeight(*geoid);
                value = coords->HeightAboveGeoid;
            }
            if(j == 1) {
                done = 1;
            } else {
                printf("\nIllegal Format, please re-enter as '(-)HHH.hhh:'\n");
            }
            if((value < bounds[0] || (value > bounds[1] && UpBoundOn)) && done == 1) {
                if(UpBoundOn) {
                    done = 0;
                    printf("\nWarning: The value you have entered of %f km for the elevation is outside of the required "
                           "range.\n",
                           value);
                    printf(" An elevation between %d km and %d km is needed. \n", bounds[0], bounds[1]);
                    if(AltitudeSetting == WGS84ON) {
                        printf("Please enter height above WGS-84 Ellipsoid (in kilometers):\n");
                    } else if(AltitudeSetting == MSLON) {
                        printf("Please enter height above mean sea level (in kilometers):\n");
                    } else {
                        printf("Please enter height in kilometers (prepend E for height above WGS-84 Ellipsoid):");
                    }
                } else {
                    switch(Warnings(3, value, {})) {
                        case 0:
                            return USER_GAVE_UP;
                        case 1:
                            done = 0;
                            printf("Please enter height above sea level (in kilometers):\n");
                            break;
                        case 2:
                            break;
                    }
                }
            }
        }
        return 0;
    }

    void GetMinGridInput(double *coord, double *val_bound, const std::string_view var_name) {
        std::string buffer;

        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
            *coord = 0;
            printf("Unrecognized input default %lf used\n", *coord);
        }


        while(*coord < val_bound[0] || *coord > val_bound[1]) {
            printf("Error: Degree input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the minimal %s:", var_name.data());
            if(NULL == fgets(buffer.data(), 20, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
                *coord = 0;
                printf("Unrecognized input default %lf used\n", *coord);
            }
        }
    }

    void GetMinGridInputDecYear(double *coord, double *val_bound) {
        std::string buffer;

        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
            *coord = 0;
        }

        while(*coord < val_bound[0] || *coord >= val_bound[1]) {
            printf("Error: Decimal year input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the minimal decimal year:");
            if(NULL == fgets(buffer.data(), 20, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
                *coord = 0;
                if(buffer[64 - 1] != '\n') {
                    clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read the
                                           // whole stdin
                }
            }
        }
    }

    void GetMaxGridInputDecYear(double *coord, double *val_bound) {
        std::string buffer;

        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
            *coord = 0;
        }


        while(*coord < val_bound[0] || *coord >= val_bound[1]) {
            printf("Error: Decimal year input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the maximal decimal year:");
            if(NULL == fgets(buffer.data(), 20, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
                *coord = 0;
                if(buffer[64 - 1] != '\n') {
                    clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read the
                                           // whole stdin
                }
            }
        }
    }

    void GetMaxGridInputAlt(double *coord, double min) {
        std::string buffer;

        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
            *coord = 0;
            printf("Unrecognized input default %lf used\n", *coord);
        }


        while(*coord < min) {
            printf("Error maximum altitude is less than minimum altitude\n");
            printf("Please re-enter the maximal altitude:");
            if(NULL == fgets(buffer.data(), 20, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
                *coord = 0;
                printf("Unrecognized input default %lf used\n", *coord);
            }
        }
    }

    void GetMaxGridInput(double *coord, double *val_bound, const std::string_view var_name) {
        std::string buffer;

        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
            *coord = 0;
            printf("Unrecognized input default %lf used\n", *coord);
        }


        while(*coord < val_bound[0] || *coord > val_bound[1]) {
            printf("Error: Degree input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the maximal %s:", var_name.data());
            if(NULL == fgets(buffer.data(), 20, stdin) || sscanf(buffer.data(), "%lf", coord) != 1) {
                *coord = 0;
                printf("Unrecognized input default %lf used\n", *coord);
            }
        }
    }

    int GetUserInput(const MagneticModel &magneticModel, Geoid *geoid, CoordGeodetic *coordGeodetic, Date *magneticDate)

    /*
    This prompts the user for coordinates, and accepts many entry formats.
    It takes the magneticModel and geoid as input and outputs the Geographic coordinates and Date as objects.
    Returns 0 when the user wants to exit and 1 if the user enters valid input data.
    INPUT :  magneticModel  : Data structure with the following elements used here
                            double epoch;       Base time of Geomagnetic model epoch (yrs)
                    : geoid Pointer to data structure geoid (used for converting HeightAboveGeoid to
    HeightABoveEllipsoid

    OUTPUT: coordGeodetic : Pointer to data structure. Following elements are updated
                            double lambdag; (longitude)
                            double phi; ( geodetic latitude)
                            double HeightAboveEllipsoid; (height above the ellipsoid (HaE) )
                            double HeightAboveGeoid;(height above the geoid )

                    magneticDate : Pointer to data structure Date with the following elements updated
                            int	Year; (If user directly enters decimal year this field is not populated)
                            int	Month;(If user directly enters decimal year this field is not populated)
                            int	Day; (If user directly enters decimal year this field is not populated)
                            double DecimalYear;      decimal years

    CALLS: 	DMSstringToDegree(buffer, &coordGeodetic->lambdag); (The program uses this to convert the string into a
    decimal longitude.) ValidateDMSstringlong(buffer, Error_Message) ValidateDMSstringlat(buffer, Error_Message)
    Warnings ConvertGeoidToEllipsoidHeight DateToYear

     */
    {
        std::string Error_Message;
        std::string buffer;
        buffer.reserve(64);
        int i, j, a, b, c, done = 0;
        double lat_bound[2] = {LAT_BOUND_MIN, LAT_BOUND_MAX};
        double lon_bound[2] = {LON_BOUND_MIN, LON_BOUND_MAX};
        int alt_bound[2]    = {ALT_BOUND_MIN, NO_ALT_MAX};

        std::string Qstring{"\nPlease enter latitude\nNorth latitude positive, For example:"
                            "\n30, 30, 30 (D,M,S) or 30.508 (Decimal Degrees) (both are north)\n"};
        GetDeg(Qstring, &coordGeodetic->phi, lat_bound);

        Qstring = "\nPlease enter longitude\nEast longitude positive, West negative.  For example:\n-100.5 or "
                  "-100, 30, 0 for 100.5 degrees west\n";
        GetDeg(Qstring, &coordGeodetic->lambda, lon_bound);

        Qstring = "\nPlease enter height above mean sea level (in kilometers):\n[For height above WGS-84 ellipsoid "
                  "prefix E, for example (E20.1)]\n";
        if(GetAltitude(Qstring, geoid, coordGeodetic, alt_bound, false) == USER_GAVE_UP) {
            return false;
        }

        printf("\nPlease enter the decimal year or calendar date\n (YYYY.yyy, MM DD YYYY or MM/DD/YYYY):\n");
        while(NULL == fgets(buffer.data(), 64, stdin)) {
            printf("\nPlease enter the decimal year or calendar date\n (YYYY.yyy, MM DD YYYY or MM/DD/YYYY):\n");
        }

        for(i = 0, done = 0; i < 64 && !done; i++) {
            if(buffer[i] == '.') {
                j = sscanf(buffer.data(), "%lf", &magneticDate->DecimalYear);
                if(j == 1) {
                    done = 1;
                } else {
                    buffer[i] = '\0';
                }
            }
            if(buffer[i] == '/') {
                if(!dateStr_to_ymd(buffer, magneticDate->Year, magneticDate->Month, magneticDate->Day)) {
                    printf("\nPlease re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a decimal year\n");
                    while(nullptr == fgets(buffer.data(), 64, stdin)) {
                        printf("\nPlease re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a decimal year\n");
                        if(buffer[64 - 1] != '\n') {
                            clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to
                                                   // read the whole stdin
                        }
                    }
                    i = 0;
                } else {
                    magneticDate->calcDecYear();
                    done = 1;
                }
            }
            if((i < 64 - 1 && buffer[i] == ' ' && buffer[i + 1] != '/') || buffer[i] == '\0') {
                if(3 == sscanf(buffer.c_str(), "%d %d %d", &a, &b, &c)) {
                    if(dateStr_to_ymd(buffer, magneticDate->Year, magneticDate->Month, magneticDate->Day)) {
                        magneticDate->calcDecYear();
                    }
                } else if(1 == sscanf(buffer.c_str(), "%d %d %d", &a, &b, &c)) {
                    magneticDate->DecimalYear = a;
                    done                      = 1;
                }
                if(!(magneticDate->DecimalYear == a)) {
                    if(false /*!DateToYear(magneticDate, Error_Message)*/) {
                        printf("%s", Error_Message.c_str());
                        buffer.clear();
                        printf(
                            "\nError encountered, please re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a decimal "
                            "year\n");
                        while(NULL == fgets(buffer.data(), 64, stdin)) {
                            printf("\nError encountered, please re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a "
                                   "decimal year\n");
                            if(buffer.back() != '\n') {
                                clear_input_buffer();
                            }
                        }

                        i = -1;
                    } else {
                        done = 1;
                    }
                }
            }
            if(buffer[i] == '\0' && i != -1 && done != 1) {
                buffer.clear();
                printf("\nError encountered, please re-enter as MM/DD/YYYY, MM DD YYYY, or as YYYY.yyy:\n");
                while(NULL == fgets(buffer.data(), 64, stdin)) {
                    printf("\nError encountered, please re-enter as MM/DD/YYYY, MM DD YYYY, or as YYYY.yyy:\n");
                    if(buffer[64 - 1] != '\n') {
                        clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read
                                               // the whole stdin
                    }
                }

                i = -1;
            }
            if(done) {
                if(magneticDate->DecimalYear > magneticModel.coefficientFileEndDate ||
                   magneticDate->DecimalYear < magneticModel.min_year) {
                    switch(Warnings(4, magneticDate->DecimalYear, magneticModel)) {
                        case 0:
                            return 0;
                        case 1:
                            done = 0;
                            i    = -1;
                            buffer.clear();
                            printf("\nPlease enter the decimal year or calendar date\n (YYYY.yyy, MM DD YYYY or "
                                   "MM/DD/YYYY):\n");
                            while(NULL == fgets(buffer.data(), 64, stdin)) {
                                printf("\nPlease enter the decimal year or calendar date\n (YYYY.yyy, MM DD YYYY or "
                                       "MM/DD/YYYY):\n");
                                if(buffer[64 - 1] != '\n') {
                                    clear_input_buffer();  // Remove the left characters from stdin if the buffer is not
                                                           // able to read the whole stdin
                                }
                            }

                            break;
                        case 2:
                            break;
                    }
                }
            }
        }
        return true;
    } /*GetUserInput*/

    /** Prompts user to enter parameters to compute a grid - for use with the grid function
Note: The user entries are not validated before here. The function populates the input variables & data structures.

    UPDATE : minimum Pointer to data structure with the following elements
             double lambdag; (longitude)
        double phi; ( geodetic latitude)
        double HeightAboveEllipsoid; (height above the ellipsoid (HaE) )
        double HeightAboveGeoid;(height above the Geoid )

            maximum   -same as the above -USE_GEOID
                                step_size  : double pointer : spatial step size, in decimal degrees
            a_step_size : double pointer :  double altitude step size (km)
                                               step_time : double pointer : time step size (decimal years)
                                                                                StartDate : pointer to data structure with
the following elements updates double DecimalYear;     ( decimal years ) EndDate :	Same as the above CALLS : none


    */
    int GetUserGrid(CoordGeodetic *minimum, CoordGeodetic *maximum, double *step_size, double *a_step_size,
                    double *step_time, Date *StartDate, Date *EndDate, int *ElementOption, int *PrintOption,
                    char *OutputFile, Geoid *Geoid, MagneticModel *model) {
        FILE *fileout;
        char filename[] = "GridProgramDirective.txt";
        std::string buffer;
        char default_outfile[16] = "GridResults.txt";

        int dummy;
        double lat_bound[2] = {LAT_BOUND_MIN, LAT_BOUND_MAX};
        double lon_bound[2] = {LON_BOUND_MIN, LON_BOUND_MAX};


        printf("Please Enter Minimum Latitude (in decimal degrees):\n");
        std::string var_name{"latitude"};
        GetMinGridInput(&minimum->phi, lat_bound, var_name);

        buffer.clear();
        lat_bound[0] = minimum->phi;
        printf("Please Enter Maximum Latitude (in decimal degrees):\n");
        GetMaxGridInput(&maximum->phi, lat_bound, var_name);

        buffer.clear();
        printf("Please Enter Minimum Longitude (in decimal degrees):\n");
        var_name = "lontitude";
        GetMinGridInput(&minimum->lambda, lon_bound, var_name);

        buffer.clear();
        lon_bound[0] = minimum->lambda;
        printf("Please Enter Maximum Longitude (in decimal degrees):\n");
        GetMaxGridInput(&maximum->lambda, lon_bound, var_name);

        buffer.clear();
        printf("Please Enter Step Size (in decimal degrees):\n");
        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", step_size) != 1) {
            *step_size = fmax(maximum->phi - minimum->phi, maximum->lambda - minimum->lambda);
            printf("Unrecognized input default %lf used\n", *step_size);

        } else {
            sscanf(buffer.data(), "%lf", step_size);
        }


        buffer.clear();

        printf("Select height (default : above MSL) \n1. Above Mean Sea Level\n2. Above WGS-84 Ellipsoid \n");
        if(NULL == fgets(buffer.data(), 64, stdin)) {
            Geoid->UseGeoid = 1;
            printf("Unrecognized option, height above MSL used.");

        } else {
            sscanf(buffer.data(), "%d", &dummy);
            if(dummy == 2) {
                Geoid->UseGeoid = 0;
            } else {
                Geoid->UseGeoid = 1;
            }
        }


        buffer.clear();
        if(Geoid->UseGeoid == 1) {
            printf("Please Enter Minimum Height above MSL (in km):\n");
            if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", &minimum->HeightAboveGeoid) != 1) {
                minimum->HeightAboveGeoid = 0;
                printf("Unrecognized input default %lf used\n", minimum->HeightAboveGeoid);
            } else {
                sscanf(buffer.data(), "%lf", &minimum->HeightAboveGeoid);
            }
            buffer.clear();
            printf("Please Enter Maximum Height above MSL (in km):\n");
            GetMaxGridInputAlt(&maximum->HeightAboveGeoid, minimum->HeightAboveGeoid);
            buffer.clear();

        } else {
            printf("Please Enter Minimum Height above the WGS-84 Ellipsoid (in km):\n");
            if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", &maximum->HeightAboveGeoid) != 1) {
                maximum->HeightAboveGeoid = 0;
                printf("Unrecognized input default %lf used\n", maximum->HeightAboveGeoid);
            } else {
                sscanf(buffer.data(), "%lf", &minimum->HeightAboveGeoid);
            }

            buffer.clear();
            printf("Please Enter Maximum Height above the WGS-84 Ellipsoid (in km):\n");
            GetMaxGridInputAlt(&maximum->HeightAboveGeoid, minimum->HeightAboveGeoid);
            buffer.clear();
        }
        printf("Please Enter height step size (in km):\n");
        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", a_step_size) != 1) {
            *a_step_size = maximum->HeightAboveGeoid - minimum->HeightAboveGeoid;
            printf("Unrecognized input default %lf used\n", *a_step_size);

        } else {
            sscanf(buffer.data(), "%lf", a_step_size);
        }

        buffer.clear();

        double dec_year_bound[2] = {model->min_year, model->coefficientFileEndDate};
        printf("\nPlease Enter the decimal year starting time:\n");
        GetMinGridInputDecYear(&StartDate->DecimalYear, dec_year_bound);

        buffer.clear();

        dec_year_bound[0] = StartDate->DecimalYear;
        printf("Please Enter the decimal year ending time:\n");
        GetMaxGridInputDecYear(&EndDate->DecimalYear, dec_year_bound);

        buffer.clear();
        printf("Please Enter the time step size:\n");
        if(NULL == fgets(buffer.data(), 64, stdin) || sscanf(buffer.data(), "%lf", step_time) != 1) {
            *step_time = EndDate->DecimalYear - StartDate->DecimalYear;
            printf("Unrecognized input, default of %lf used\n", *step_time);

        } else {
            sscanf(buffer.data(), "%lf", step_time);
        }


        buffer.clear();
        printf("Enter a geomagnetic element to print. Your options are:\n");
        printf(" 1. Declination	9.   Ddot\n 2. Inclination	10. Idot\n 3. F		11. Fdot\n 4. H		12. Hdot\n 5. X		13. "
               "Xdot\n 6. Y		14. Ydot\n 7. Z		15. Zdot\n 8. GV		16. GVdot\nFor gradients enter: 17\n");
        if(NULL == fgets(buffer.data(), 64, stdin)) {
            *ElementOption = 1;
            printf("Unrecognized input, default of %d used\n", *ElementOption);
        }
        sscanf(buffer.data(), "%d", ElementOption);


        buffer.clear();
        if(*ElementOption == 17) {
            printf("Enter a gradient element to print. Your options are:\n");
            printf(" 1. dX/dphi \t2. dY/dphi \t3. dZ/dphi\n");
            printf(" 4. dX/dlambda \t5. dY/dlambda \t6. dZ/dlambda\n");
            printf(" 7. dX/dz \t8. dY/dz \t9. dZ/dz\n");
            buffer.clear();
            if(NULL == fgets(buffer.data(), 20, stdin)) {
                *ElementOption = 1;
                printf("Unrecognized input, default of %d used\n", *ElementOption);
            } else {
                sscanf(buffer.data(), "%d", ElementOption);
            }
            buffer.clear();
            *ElementOption += 16;
        }
        printf("Select output :\n");
        printf(" 1. Print to a file \n 2. Print to Screen\n");
        if(NULL == fgets(buffer.data(), 64, stdin)) {
            *PrintOption = 2;
            printf("Unrecognized input, default of printing to screen\n");

        } else {
            sscanf(buffer.data(), "%d", PrintOption);
        }


        buffer.clear();
        fileout = fopen(filename, "a");
        if(*PrintOption == 1) {
            printf("Please enter output filename\nfor default ('%s') press enter:\n", default_outfile);
            if(NULL == fgets(buffer.data(), 64, stdin) || buffer.length() <= 1) {
                fprintf(fileout, "\nResults printed in: %s\n", default_outfile);
                // The size of OutputFile is 32 which is defined in #127 in GeomagnetismHeader.h
                OutputFile = "GridResults.txt";
            } else {
                sscanf(buffer.data(), "%s", OutputFile);
                fprintf(fileout, "\nResults printed in: %s\n", OutputFile);
            }
            /*strcpy(OutputFile, buffer);*/

            buffer.clear();
            /*sscanf(buffer, "%s", OutputFile);*/
        } else {
            fprintf(fileout, "\nResults printed in Console\n");
        }
        fprintf(fileout,
                "Minimum Latitude: %f\t\tMaximum Latitude: %f\t\tStep Size: %f\nMinimum Longitude: %f\t\tMaximum Longitude: "
                "%f\t\tStep Size: %f\n",
                minimum->phi, maximum->phi, *step_size, minimum->lambda, maximum->lambda, *step_size);
        if(Geoid->UseGeoid == 1) {
            fprintf(fileout, "Minimum Altitude above MSL: %f\tMaximum Altitude above MSL: %f\tStep Size: %f\n",
                    minimum->HeightAboveGeoid, maximum->HeightAboveGeoid, *a_step_size);
        } else {
            fprintf(
                fileout,
                "Minimum Altitude above WGS-84 Ellipsoid: %f\tMaximum Altitude above WGS-84 Ellipsoid: %f\tStep Size: %f\n",
                minimum->HeightAboveEllipsoid, maximum->HeightAboveEllipsoid, *a_step_size);
        }
        fprintf(fileout, "Starting Date: %f\t\tEnding Date: %f\t\tStep Time: %f\n\n\n", StartDate->DecimalYear,
                EndDate->DecimalYear, *step_time);
        fclose(fileout);
        return true;
    }

    void printCoord(const CoordGeodetic &coord, const int useGeoid) {
        printf("\n Results For \n\n");
        if(coord.phi < 0) {
            printf("Latitude	%.2fS\n", -coord.phi);
        } else {
            printf("Latitude	%.2fN\n", coord.phi);
        }
        if(coord.lambda < 0) {
            printf("Longitude	%.2fW\n", -coord.lambda);
        } else {
            printf("Longitude	%.2fE\n", coord.lambda);
        }
        if(useGeoid == 1) {
            printf("Altitude:	%.2f Kilometers above MSL\n", coord.HeightAboveGeoid);
        } else {
            printf("Altitude:	%.2f Kilometers above WGS-84 Ellipsoid\n", coord.HeightAboveEllipsoid);
        }
    }

    void PrintUserDataWithUncertainty(const GeoMagneticElements &geomagElements, const GeoMagneticElements &Errors,
                                      const CoordGeodetic &spaceInput, const Date &TimeInput,
                                      const MagneticModel &magneticModel, const Geoid &Geoid) {
        std::string DeclString;
        std::string InclString;
        std::string GVString;

        DegreeToDMSstring(geomagElements.Incl, 2, InclString);
        //        DegreeToDMSstring(geomagElements.Decl, 2, DeclString);
        DeclString = std::to_string(geomagElements.Decl);

        if(geomagElements.H < 6000 && geomagElements.H > 2000) {
            Warnings(1, geomagElements.H, magneticModel);
        }
        if(geomagElements.H < 2000) {
            Warnings(2, geomagElements.H, magneticModel);
        }

        printCoord(spaceInput, Geoid.UseGeoid);
        printf("Date:		%.1f \t\t %d/%d/%d\n", TimeInput.DecimalYear, TimeInput.Day, TimeInput.Month, TimeInput.Year);

        if(magneticModel.secularVariationUsed == true) {
            printf("\n		Main Field\t\t\tSecular Change\n");
            printf("F	=	%9.1f +/- %5.1f nT\t\t Fdot = %5.1f\tnT/yr\n", geomagElements.F, Errors.F, geomagElements.Fdot);
            printf("H	=	%9.1f +/- %5.1f nT\t\t Hdot = %5.1f\tnT/yr\n", geomagElements.H, Errors.H, geomagElements.Hdot);
            printf("X	=	%9.1f +/- %5.1f nT\t\t Xdot = %5.1f\tnT/yr\n", geomagElements.X, Errors.X, geomagElements.Xdot);
            printf("Y	=	%9.1f +/- %5.1f nT\t\t Ydot = %5.1f\tnT/yr\n", geomagElements.Y, Errors.Y, geomagElements.Ydot);
            printf("Z	=	%9.1f +/- %5.1f nT\t\t Zdot = %5.1f\tnT/yr\n", geomagElements.Z, Errors.Z, geomagElements.Zdot);

            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST) +/-%3.0f Min Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * Errors.Decl,
                       60 * geomagElements.Decldot);
            } else {
                printf("Decl	=%20s  (EAST) +/-%3.0f Min Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * Errors.Decl,
                       60 * geomagElements.Decldot);
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)   +/-%3.0f Min Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * Errors.Incl,
                       60 * geomagElements.Incldot);
            } else {
                printf("Incl	=%20s  (DOWN) +/-%3.0f Min Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * Errors.Incl,
                       60 * geomagElements.Incldot);
            }
        } else {
            printf("\n	Main Field\n");
            printf("F	=	%-9.1f +/-%5.1f nT\n", geomagElements.F, Errors.F);
            printf("H	=	%-9.1f +/-%5.1f nT\n", geomagElements.H, Errors.H);
            printf("X	=	%-9.1f +/-%5.1f nT\n", geomagElements.X, Errors.X);
            printf("Y	=	%-9.1f +/-%5.1f nT\n", geomagElements.Y, Errors.Y);
            printf("Z	=	%-9.1f +/-%5.1f nT\n", geomagElements.Z, Errors.Z);

            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST)+/-%4f\n", DeclString.c_str(), 60 * Errors.Decl);
            } else {
                printf("Decl	=%20s  (EAST)+/-%4f\n", DeclString.c_str(), 60 * Errors.Decl);
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)+/-%4f\n", InclString.c_str(), 60 * Errors.Incl);
            } else {
                printf("Incl	=%20s  (DOWN)+/-%4f\n", InclString.c_str(), 60 * Errors.Incl);
            }
        }

        DegreeToDMSstring(geomagElements.GV, 2, GVString);

        /* Print Grid Variation */

        if(spaceInput.phi < -55) {
            printf("\n\n Grid variation (SOUTH) =%20s\n", GVString.c_str());
        } else if(spaceInput.phi > 55) {
            printf("\n\n Grid variation (NORTH) =%20s \n", GVString.c_str());
        }
    } /*PrintUserDataWithUncertainty*/

    void PrintWMMFormat(char *filename, const MagneticModel &magneticModel) {
        int index;
        FILE *OUT;
        Date Date;
        char Datestring[11];

        Date.DecimalYear = magneticModel.editionDate;
        //        YearToDate(&Date);
        sprintf(Datestring, "%d/%d/%d", Date.Month, Date.Day, Date.Year);
        OUT = fopen(filename, "w");
        fprintf(OUT, "    %.1f               %s              %s\n", magneticModel.epoch, magneticModel.modelName.c_str(),
                Datestring);
        for(int n = 1; n <= magneticModel.nMax; n++) {
            for(int m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);
                if(m != 0) {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f  %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index],
                            magneticModel.main_Field_Coeff_H[index], magneticModel.secular_Var_Coeff_G[index],
                            magneticModel.secular_Var_Coeff_H[index]);
                } else {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f  %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index], 0.0,
                            magneticModel.secular_Var_Coeff_G[index], 0.0);
                }
            }
        }
        fclose(OUT);
    } /*PrintWMMFormat*/

    void PrintEMMFormat(char *filename, char *filenameSV, const MagneticModel &magneticModel) {
        int index, n, m;
        FILE *OUT;
        Date Date;
        char Datestring[11];

        Date.DecimalYear = magneticModel.editionDate;
        //        YearToDate(&Date);
        sprintf(Datestring, "%d/%d/%d", Date.Month, Date.Day, Date.Year);
        OUT = fopen(filename, "w");
        fprintf(OUT, "    %.1f               %s              %s\n", magneticModel.epoch, magneticModel.modelName.c_str(),
                Datestring);
        for(n = 1; n <= magneticModel.nMax; n++) {
            for(m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);
                if(m != 0) {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index],
                            magneticModel.main_Field_Coeff_H[index]);
                } else {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.main_Field_Coeff_G[index], 0.0);
                }
            }
        }
        fclose(OUT);
        OUT = fopen(filenameSV, "w");
        for(n = 1; n <= magneticModel.nMaxSecVar; n++) {
            for(m = 0; m <= n; m++) {
                index = (n * (n + 1) / 2 + m);
                if(m != 0) {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.secular_Var_Coeff_G[index],
                            magneticModel.secular_Var_Coeff_H[index]);
                } else {
                    fprintf(OUT, " %2d %2d %9.4f %9.4f\n", n, m, magneticModel.secular_Var_Coeff_G[index], 0.0);
                }
            }
        }
        fclose(OUT);
    } /*PrintEMMFormat*/

    void PrintSHDFFormat(char *filename, MagneticModel *(*magneticModel)[], int epochs) {
        int i, n, m, index, epochRange;
        FILE *SHDF_file;
        SHDF_file = fopen(filename, "w");
        /*lines = (int)(UFM_DEGREE / 2.0 * (UFM_DEGREE + 3));*/
        for(i = 0; i < epochs; i++) {
            if(i < epochs - 1) {
                epochRange = (*magneticModel)[i + 1]->epoch - (*magneticModel)[i]->epoch;
            } else {
                epochRange = (*magneticModel)[i]->epoch - (*magneticModel)[i - 1]->epoch;
            }
            fprintf(SHDF_file, "%%SHDF 16695 Definitive Geomagnetic Reference Field Model Coefficient File\n");
            fprintf(SHDF_file, "%%ModelName: %s\n", (*magneticModel)[i]->modelName.c_str());
            fprintf(SHDF_file,
                    "%%Publisher: International Association of Geomagnetism and Aeronomy (IAGA), Working Group V-Mod\n");
            fprintf(SHDF_file, "%%ReleaseDate: Some Number\n");
            fprintf(SHDF_file, "%%DataCutOFF: Some Other Number\n");
            fprintf(SHDF_file, "%%ModelStartYear: %d\n", (int)(*magneticModel)[i]->epoch);
            fprintf(SHDF_file, "%%ModelEndYear: %d\n", (int)(*magneticModel)[i]->epoch + epochRange);
            fprintf(SHDF_file, "%%Epoch: %.0f\n", (*magneticModel)[i]->epoch);
            fprintf(SHDF_file, "%%IntStaticDeg: %d\n", (*magneticModel)[i]->nMax);
            fprintf(SHDF_file, "%%IntSecVarDeg: %d\n", (*magneticModel)[i]->nMaxSecVar);
            fprintf(SHDF_file, "%%ExtStaticDeg: 0\n");
            fprintf(SHDF_file, "%%ExtSecVarDeg: 0\n");
            fprintf(SHDF_file, "%%Normalization: Schmidt semi-normailized\n");
            fprintf(SHDF_file, "%%SpatBasFunc: spherical harmonics\n");
            fprintf(SHDF_file, "# To synthesize the field for a given date:\n");
            fprintf(SHDF_file, "# Use the sub-model of the epoch corresponding to each date\n");
            fprintf(SHDF_file, "#\n#\n#\n#\n# I/E, n, m, Gnm, Hnm, SV-Gnm, SV-Hnm\n#\n");
            n = 1;
            m = 0;
            for(n = 1; n <= (*magneticModel)[i]->nMax; n++) {
                for(m = 0; m <= n; m++) {
                    index = (n * (n + 1)) / 2 + m;
                    if(i < epochs - 1) {
                        if(m != 0) {
                            fprintf(SHDF_file, "I,%d,%d,%f,%f,%f,%f\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->main_Field_Coeff_H[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_H[index]);
                        } else {
                            fprintf(SHDF_file, "I,%d,%d,%f,,%f,\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index]);
                        }
                    } else {
                        if(m != 0) {
                            fprintf(SHDF_file, "I,%d,%d,%f,%f,%f,%f\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->main_Field_Coeff_H[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_H[index]);
                        } else {
                            fprintf(SHDF_file, "I,%d,%d,%f,,%f,\n", n, m, (*magneticModel)[i]->main_Field_Coeff_G[index],
                                    (*magneticModel)[i]->secular_Var_Coeff_G[index]);
                        }
                    }
                }
            }
        }
    } /*PrintSHDFFormat*/

    void PrintGradient(const Gradient &gradient) {
        printf("\nGradient\n");
        printf("\n                 Northward       Eastward        Downward\n");
        printf("X:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.X, gradient.GradLambda.X,
               gradient.GradZ.X);
        printf("Y:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.Y, gradient.GradLambda.Y,
               gradient.GradZ.Y);
        printf("Z:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.Z, gradient.GradLambda.Z,
               gradient.GradZ.Z);
        printf("H:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.H, gradient.GradLambda.H,
               gradient.GradZ.H);
        printf("F:           %7.1f nT/km %9.1f nT/km %9.1f nT/km \n", gradient.GradPhi.F, gradient.GradLambda.F,
               gradient.GradZ.F);
        printf("Declination: %7.2f min/km %8.2f min/km %8.2f min/km \n", gradient.GradPhi.Decl * 60,
               gradient.GradLambda.Decl * 60, gradient.GradZ.Decl * 60);
        printf("Inclination: %7.2f min/km %8.2f min/km %8.2f min/km \n", gradient.GradPhi.Incl * 60,
               gradient.GradLambda.Incl * 60, gradient.GradZ.Incl * 60);
    }

    /** This function prints the results in  Geomagnetic Elements for a point calculation. It takes the calculated
     *  Geomagnetic elements "GeomagElements" as input.
     *  As well as the coordinates, date, and Magnetic Model.
     * INPUT :  GeomagElements
     *          CoordGeodetic
     *          Date
     *          MagneticModel
     * OUTPUT : none
     */
    void PrintUserData(const GeoMagneticElements &geomagElements, const CoordGeodetic &spaceInput, const Date &timeInput,
                       const MagneticModel &magneticModel, const Geoid &geoid) {
        std::string DeclString;
        std::string InclString;
        DegreeToDMSstring(geomagElements.Incl, 2, InclString);
        if(geomagElements.H < 6000 && geomagElements.H > 2000) {
            Warnings(1, geomagElements.H, magneticModel);
        }
        if(geomagElements.H < 2000) {
            Warnings(2, geomagElements.H, magneticModel);
        }

        printCoord(spaceInput, geoid.UseGeoid);

        DegreeToDMSstring(geomagElements.Decl, 2, DeclString);
        if(magneticModel.secularVariationUsed == true) {
            printf("Date:		%.1f\n", timeInput.DecimalYear);
            printf("\n		Main Field\t\t\tSecular Change\n");
            printf("F	=	%-9.1f nT\t\t  Fdot = %.1f\tnT/yr\n", geomagElements.F, geomagElements.Fdot);
            printf("H	=	%-9.1f nT\t\t  Hdot = %.1f\tnT/yr\n", geomagElements.H, geomagElements.Hdot);
            printf("X	=	%-9.1f nT\t\t  Xdot = %.1f\tnT/yr\n", geomagElements.X, geomagElements.Xdot);
            printf("Y	=	%-9.1f nT\t\t  Ydot = %.1f\tnT/yr\n", geomagElements.Y, geomagElements.Ydot);
            printf("Z	=	%-9.1f nT\t\t  Zdot = %.1f\tnT/yr\n", geomagElements.Z, geomagElements.Zdot);
            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST)\t  Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * geomagElements.Decldot);
            } else {
                printf("Decl	=%20s  (EAST)\t  Ddot = %.1f\tMin/yr\n", DeclString.c_str(), 60 * geomagElements.Decldot);
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)\t  Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * geomagElements.Incldot);
            } else {
                printf("Incl	=%20s  (DOWN)\t  Idot = %.1f\tMin/yr\n", InclString.c_str(), 60 * geomagElements.Incldot);
            }
        } else {
            printf("Date:		%.1f\n", timeInput.DecimalYear);
            printf("\n	Main Field\n");
            printf("F	=	%-9.1f nT\n", geomagElements.F);
            printf("H	=	%-9.1f nT\n", geomagElements.H);
            printf("X	=	%-9.1f nT\n", geomagElements.X);
            printf("Y	=	%-9.1f nT\n", geomagElements.Y);
            printf("Z	=	%-9.1f nT\n", geomagElements.Z);
            if(geomagElements.Decl < 0) {
                printf("Decl	=%20s  (WEST)\n", DeclString.c_str());
            } else {
                printf("Decl	=%20s  (EAST)\n", DeclString.c_str());
            }
            if(geomagElements.Incl < 0) {
                printf("Incl	=%20s  (UP)\n", InclString.c_str());
            } else {
                printf("Incl	=%20s  (DOWN)\n", InclString.c_str());
            }
        }

        if(spaceInput.phi <= -55 || spaceInput.phi >= 55)
        /* Print Grid Variation */
        {
            DegreeToDMSstring(geomagElements.GV, 2, InclString);
            printf("\n\n Grid variation =%20s\n", InclString.c_str());
        }

    } /*PrintUserData*/
}  // namespace wmm