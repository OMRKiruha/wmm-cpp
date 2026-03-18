#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "GeomagInterativeLib.h"

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
        std::string Error_Message;
        int done, i, j;

        printf("%s", Query_String.data());
        while(NULL == fgets(buffer.data(), 64, stdin)) {
            printf("%s", Query_String.data());
            if(buffer[sizeof(buffer) - 1] != '\n') {
                clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read the
                                       // whole stdin
            }
        }

        for(i = 0, done = 0, j = 0; i < (int)sizeof(buffer) && !done; i++) {
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
                        if(buffer[sizeof(buffer) - 1] != '\n') {
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
                ConvertGeoidToEllipsoidHeight(coords, geoid);
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

    void GetMinGridInput(double *coord, double *val_bound, char *var_name) {
        char buffer[20];

        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", coord) != 1) {
            *coord = 0;
            printf("Unrecognized input default %lf used\n", *coord);
        }


        while(*coord < val_bound[0] || *coord > val_bound[1]) {
            printf("Error: Degree input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the minimal %s:", var_name);
            if(NULL == fgets(buffer, 20, stdin) || sscanf(buffer, "%lf", coord) != 1) {
                *coord = 0;
                printf("Unrecognized input default %lf used\n", *coord);
            }
        }
    }

    void GetMinGridInputDecYear(double *coord, double *val_bound) {
        char buffer[20];

        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", coord) != 1) {
            *coord = 0;
        }

        while(*coord < val_bound[0] || *coord >= val_bound[1]) {
            printf("Error: Decimal year input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the minimal decimal year:");
            if(NULL == fgets(buffer, 20, stdin) || sscanf(buffer, "%lf", coord) != 1) {
                *coord = 0;
                if(buffer[sizeof(buffer) - 1] != '\n') {
                    clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read the
                                           // whole stdin
                }
            }
        }
    }

    void GetMaxGridInputDecYear(double *coord, double *val_bound) {
        char buffer[20];

        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", coord) != 1) {
            *coord = 0;
        }


        while(*coord < val_bound[0] || *coord >= val_bound[1]) {
            printf("Error: Decimal year input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the maximal decimal year:");
            if(NULL == fgets(buffer, 20, stdin) || sscanf(buffer, "%lf", coord) != 1) {
                *coord = 0;
                if(buffer[sizeof(buffer) - 1] != '\n') {
                    clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read the
                                           // whole stdin
                }
            }
        }
    }

    void GetMaxGridInputAlt(double *coord, double min) {
        char buffer[20];

        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", coord) != 1) {
            *coord = 0;
            printf("Unrecognized input default %lf used\n", *coord);
        }


        while(*coord < min) {
            printf("Error maximum altitude is less than minimum altitude\n");
            printf("Please re-enter the maximal altitude:");
            if(NULL == fgets(buffer, 20, stdin) || sscanf(buffer, "%lf", coord) != 1) {
                *coord = 0;
                printf("Unrecognized input default %lf used\n", *coord);
            }
        }
    }

    void GetMaxGridInput(double *coord, double *val_bound, char *var_name) {
        char buffer[20];

        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", coord) != 1) {
            *coord = 0;
            printf("Unrecognized input default %lf used\n", *coord);
        }


        while(*coord < val_bound[0] || *coord > val_bound[1]) {
            printf("Error: Degree input is outside range\n"
                   " The range is from %.2f to %.2f\n",
                   val_bound[0], val_bound[1]);
            printf("Please re-enter the maximal %s:", var_name);
            if(NULL == fgets(buffer, 20, stdin) || sscanf(buffer, "%lf", coord) != 1) {
                *coord = 0;
                printf("Unrecognized input default %lf used\n", *coord);
            }
        }
    }

    int GetUserInput(MagneticModel *magneticModel, Geoid *geoid, CoordGeodetic *coordGeodetic, Date *magneticDate)

    /*
    This prompts the user for coordinates, and accepts many entry formats.
    It takes the magneticModel and geoid as input and outputs the Geographic coordinates and Date as objects.
    Returns 0 when the user wants to exit and 1 if the user enters valid input data.
    INPUT :  magneticModel  : Data structure with the following elements used here
                            double epoch;       Base time of Geomagnetic model epoch (yrs)
                    : geoid Pointer to data structure geoid (used for converting HeightAboveGeoid to
    HeightABoveEllipsoid

    OUTPUT: coordGeodetic : Pointer to data structure. Following elements are updated
                            double lambda; (longitude)
                            double phi; ( geodetic latitude)
                            double HeightAboveEllipsoid; (height above the ellipsoid (HaE) )
                            double HeightAboveGeoid;(height above the geoid )

                    magneticDate : Pointer to data structure Date with the following elements updated
                            int	Year; (If user directly enters decimal year this field is not populated)
                            int	Month;(If user directly enters decimal year this field is not populated)
                            int	Day; (If user directly enters decimal year this field is not populated)
                            double DecimalYear;      decimal years

    CALLS: 	DMSstringToDegree(buffer, &coordGeodetic->lambda); (The program uses this to convert the string into a
    decimal longitude.) ValidateDMSstringlong(buffer, Error_Message) ValidateDMSstringlat(buffer, Error_Message)
    Warnings ConvertGeoidToEllipsoidHeight DateToYear

     */
    {
        std::string Error_Message;
        std::string buffer;
        int i, j, a, b, c, done = 0;
        double lat_bound[2] = {LAT_BOUND_MIN, LAT_BOUND_MAX};
        double lon_bound[2] = {LON_BOUND_MIN, LON_BOUND_MAX};
        int alt_bound[2]    = {ALT_BOUND_MIN, NO_ALT_MAX};
        int Qstring_size    = 1028;
        std::string Qstring;
        Qstring = "\nPlease enter latitude\nNorth latitude positive, For example:\n30, 30, 30 (D,M,S) or 30.508 "
                  "(Decimal Degrees) (both are north)\n";
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

        for(i = 0, done = 0; i < sizeof(buffer) && !done; i++) {
            if(buffer[i] == '.') {
                j = sscanf(buffer.data(), "%lf", &magneticDate->DecimalYear);
                if(j == 1) {
                    done = 1;
                } else {
                    buffer[i] = '\0';
                }
            }
            if(buffer[i] == '/') {
                sscanf(buffer.c_str(), "%d/%d/%d", &magneticDate->Month, &magneticDate->Day, &magneticDate->Year);
                if(!DateToYear(magneticDate, Error_Message)) {
                    printf("%s", Error_Message.c_str());
                    printf("\nPlease re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a decimal year\n");
                    while(NULL == fgets(buffer.data(), 64, stdin)) {
                        printf("\nPlease re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a decimal year\n");
                        if(buffer[sizeof(buffer) - 1] != '\n') {
                            clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to
                                                   // read the whole stdin
                        }
                    }

                    i = 0;
                } else {
                    done = 1;
                }
            }
            if((i < sizeof(buffer) - 1 && buffer[i] == ' ' && buffer[i + 1] != '/') || buffer[i] == '\0') {
                if(3 == sscanf(buffer.c_str(), "%d %d %d", &a, &b, &c)) {
                    magneticDate->Month       = a;
                    magneticDate->Day         = b;
                    magneticDate->Year        = c;
                    magneticDate->DecimalYear = 99999;
                } else if(1 == sscanf(buffer.c_str(), "%d %d %d", &a, &b, &c)) {
                    magneticDate->DecimalYear = a;
                    done                      = 1;
                }
                if(!(magneticDate->DecimalYear == a)) {
                    if(!DateToYear(magneticDate, Error_Message)) {
                        printf("%s", Error_Message);
                        strlcpy_equivalent(buffer, "", sizeof(buffer));
                        printf(
                            "\nError encountered, please re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a decimal "
                            "year\n");
                        while(NULL == fgets(buffer, sizeof(buffer), stdin)) {
                            printf("\nError encountered, please re-enter Date in MM/DD/YYYY or MM DD YYYY format, or as a "
                                   "decimal year\n");
                            if(buffer[strlen(buffer) - 1] != '\n') {
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
                while(NULL == fgets(buffer, sizeof(buffer), stdin)) {
                    printf("\nError encountered, please re-enter as MM/DD/YYYY, MM DD YYYY, or as YYYY.yyy:\n");
                    if(buffer[sizeof(buffer) - 1] != '\n') {
                        clear_input_buffer();  // Remove the left characters from stdin if the buffer is not able to read
                                               // the whole stdin
                    }
                }

                i = -1;
            }
            if(done) {
                if(magneticDate->DecimalYear > magneticModel->CoefficientFileEndDate ||
                   magneticDate->DecimalYear < magneticModel->min_year) {
                    switch(Warnings(4, magneticDate->DecimalYear, magneticModel)) {
                        case 0:
                            return 0;
                        case 1:
                            done = 0;
                            i    = -1;
                            strlcpy_equivalent(buffer, "", sizeof(buffer));
                            printf("\nPlease enter the decimal year or calendar date\n (YYYY.yyy, MM DD YYYY or "
                                   "MM/DD/YYYY):\n");
                            while(NULL == fgets(buffer, sizeof(buffer), stdin)) {
                                printf("\nPlease enter the decimal year or calendar date\n (YYYY.yyy, MM DD YYYY or "
                                       "MM/DD/YYYY):\n");
                                if(buffer[sizeof(buffer) - 1] != '\n') {
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
        free(Qstring);
        return true;
    } /*GetUserInput*/

    int GetUserGrid(CoordGeodetic *minimum, CoordGeodetic *maximum, double *step_size, double *a_step_size,
                    double *step_time, Date *StartDate, Date *EndDate, int *ElementOption, int *PrintOption,
                    char *OutputFile, Geoid *Geoid, MagneticModel *model)

    /* Prompts user to enter parameters to compute a grid - for use with the grid function
    Note: The user entries are not validated before here. The function populates the input variables & data structures.

    UPDATE : minimum Pointer to data structure with the following elements
                    double lambda; (longitude)
                    double phi; ( geodetic latitude)
                    double HeightAboveEllipsoid; (height above the ellipsoid (HaE) )
                    double HeightAboveGeoid;(height above the Geoid )

                    maximum   -same as the above -USE_GEOID
                    step_size  : double pointer : spatial step size, in decimal degrees
                    a_step_size : double pointer :  double altitude step size (km)
                    step_time : double pointer : time step size (decimal years)
                    StartDate : pointer to data structure with the following elements updates
                                            double DecimalYear;     ( decimal years )
                    EndDate :	Same as the above
    CALLS : none


     */
    {
        FILE *fileout;
        char filename[] = "GridProgramDirective.txt";
        char buffer[20];
        int strcopy_size         = 32;
        char default_outfile[16] = "GridResults.txt";

        int dummy;
        double lat_bound[2] = {LAT_BOUND_MIN, LAT_BOUND_MAX};
        double lon_bound[2] = {LON_BOUND_MIN, LON_BOUND_MAX};


        printf("Please Enter Minimum Latitude (in decimal degrees):\n");
        char var_name[20] = "latitude";
        GetMinGridInput(&minimum->phi, lat_bound, var_name);

        strlcpy_equivalent(buffer, "", sizeof(buffer));
        lat_bound[0] = minimum->phi;
        printf("Please Enter Maximum Latitude (in decimal degrees):\n");
        GetMaxGridInput(&maximum->phi, lat_bound, var_name);

        strlcpy_equivalent(buffer, "", sizeof(buffer));
        printf("Please Enter Minimum Longitude (in decimal degrees):\n");
        strlcpy_equivalent(var_name, "lontitude", sizeof(var_name));
        GetMinGridInput(&minimum->lambda, lon_bound, var_name);

        strlcpy_equivalent(buffer, "", sizeof(buffer));
        lon_bound[0] = minimum->lambda;
        printf("Please Enter Maximum Longitude (in decimal degrees):\n");
        GetMaxGridInput(&maximum->lambda, lon_bound, var_name);

        strlcpy_equivalent(buffer, "", sizeof(buffer));
        printf("Please Enter Step Size (in decimal degrees):\n");
        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", step_size) != 1) {
            *step_size = fmax(maximum->phi - minimum->phi, maximum->lambda - minimum->lambda);
            printf("Unrecognized input default %lf used\n", *step_size);

        } else {
            sscanf(buffer, "%lf", step_size);
        }


        strlcpy_equivalent(buffer, "", sizeof(buffer));

        printf("Select height (default : above MSL) \n1. Above Mean Sea Level\n2. Above WGS-84 Ellipsoid \n");
        if(NULL == fgets(buffer, sizeof(buffer), stdin)) {
            Geoid->UseGeoid = 1;
            printf("Unrecognized option, height above MSL used.");

        } else {
            sscanf(buffer, "%d", &dummy);
            if(dummy == 2) {
                Geoid->UseGeoid = 0;
            } else {
                Geoid->UseGeoid = 1;
            }
        }


        strlcpy_equivalent(buffer, "", sizeof(buffer));
        if(Geoid->UseGeoid == 1) {
            printf("Please Enter Minimum Height above MSL (in km):\n");
            if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", &minimum->HeightAboveGeoid) != 1) {
                minimum->HeightAboveGeoid = 0;
                printf("Unrecognized input default %lf used\n", minimum->HeightAboveGeoid);
            } else {
                sscanf(buffer, "%lf", &minimum->HeightAboveGeoid);
            }
            strlcpy_equivalent(buffer, "", sizeof(buffer));
            printf("Please Enter Maximum Height above MSL (in km):\n");
            GetMaxGridInputAlt(&maximum->HeightAboveGeoid, minimum->HeightAboveGeoid);
            strlcpy_equivalent(buffer, "", sizeof(buffer));

        } else {
            printf("Please Enter Minimum Height above the WGS-84 Ellipsoid (in km):\n");
            if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", &maximum->HeightAboveGeoid) != 1) {
                maximum->HeightAboveGeoid = 0;
                printf("Unrecognized input default %lf used\n", maximum->HeightAboveGeoid);


            } else {
                sscanf(buffer, "%lf", &minimum->HeightAboveGeoid);
            }

            strlcpy_equivalent(buffer, "", sizeof(buffer));
            printf("Please Enter Maximum Height above the WGS-84 Ellipsoid (in km):\n");
            GetMaxGridInputAlt(&maximum->HeightAboveGeoid, minimum->HeightAboveGeoid);
            strlcpy_equivalent(buffer, "", sizeof(buffer));
        }
        printf("Please Enter height step size (in km):\n");
        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", a_step_size) != 1) {
            *a_step_size = maximum->HeightAboveGeoid - minimum->HeightAboveGeoid;
            printf("Unrecognized input default %lf used\n", *a_step_size);

        } else {
            sscanf(buffer, "%lf", a_step_size);
        }

        strlcpy_equivalent(buffer, "", sizeof(buffer));

        double dec_year_bound[2] = {model->min_year, model->CoefficientFileEndDate};
        printf("\nPlease Enter the decimal year starting time:\n");
        GetMinGridInputDecYear(&StartDate->DecimalYear, dec_year_bound);

        strlcpy_equivalent(buffer, "", sizeof(buffer));

        dec_year_bound[0] = StartDate->DecimalYear;
        printf("Please Enter the decimal year ending time:\n");
        GetMaxGridInputDecYear(&EndDate->DecimalYear, dec_year_bound);

        strlcpy_equivalent(buffer, "", sizeof(buffer));
        printf("Please Enter the time step size:\n");
        if(NULL == fgets(buffer, sizeof(buffer), stdin) || sscanf(buffer, "%lf", step_time) != 1) {
            *step_time = EndDate->DecimalYear - StartDate->DecimalYear;
            printf("Unrecognized input, default of %lf used\n", *step_time);

        } else {
            sscanf(buffer, "%lf", step_time);
        }


        strlcpy_equivalent(buffer, "", sizeof(buffer));
        printf("Enter a geomagnetic element to print. Your options are:\n");
        printf(" 1. Declination	9.   Ddot\n 2. Inclination	10. Idot\n 3. F		11. Fdot\n 4. H		12. Hdot\n 5. X		13. "
               "Xdot\n 6. Y		14. Ydot\n 7. Z		15. Zdot\n 8. GV		16. GVdot\nFor gradients enter: 17\n");
        if(NULL == fgets(buffer, sizeof(buffer), stdin)) {
            *ElementOption = 1;
            printf("Unrecognized input, default of %d used\n", *ElementOption);
        }
        sscanf(buffer, "%d", ElementOption);


        strlcpy_equivalent(buffer, "", sizeof(buffer));
        if(*ElementOption == 17) {
            printf("Enter a gradient element to print. Your options are:\n");
            printf(" 1. dX/dphi \t2. dY/dphi \t3. dZ/dphi\n");
            printf(" 4. dX/dlambda \t5. dY/dlambda \t6. dZ/dlambda\n");
            printf(" 7. dX/dz \t8. dY/dz \t9. dZ/dz\n");
            strlcpy_equivalent(buffer, "", sizeof(buffer));
            if(NULL == fgets(buffer, 20, stdin)) {
                *ElementOption = 1;
                printf("Unrecognized input, default of %d used\n", *ElementOption);
            } else {
                sscanf(buffer, "%d", ElementOption);
            }
            strlcpy_equivalent(buffer, "", sizeof(buffer));
            *ElementOption += 16;
        }
        printf("Select output :\n");
        printf(" 1. Print to a file \n 2. Print to Screen\n");
        if(NULL == fgets(buffer, sizeof(buffer), stdin)) {
            *PrintOption = 2;
            printf("Unrecognized input, default of printing to screen\n");

        } else {
            sscanf(buffer, "%d", PrintOption);
        }


        strlcpy_equivalent(buffer, "", sizeof(buffer));
        fileout = fopen(filename, "a");
        if(*PrintOption == 1) {
            printf("Please enter output filename\nfor default ('%s') press enter:\n", default_outfile);
            if(NULL == fgets(buffer, sizeof(buffer), stdin) || strlen(buffer) <= 1) {
                fprintf(fileout, "\nResults printed in: %s\n", default_outfile);
                // The size of OutputFile is 32 which is defined in #127 in GeomagnetismHeader.h
                strcopy_size = (sizeof(default_outfile) > 32) ? 32 : sizeof(default_outfile);
                strlcpy_equivalent(OutputFile, "GridResults.txt", strcopy_size);

            } else {
                sscanf(buffer, "%s", OutputFile);
                fprintf(fileout, "\nResults printed in: %s\n", OutputFile);
            }
            /*strcpy(OutputFile, buffer);*/

            strlcpy_equivalent(buffer, "", sizeof(buffer));
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
}  // namespace wmm