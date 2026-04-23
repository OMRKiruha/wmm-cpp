//
// Created by Kiryuhin Viacheslav on 04/21/26
//

#pragma once

namespace wmm {
    // Constants
    // These error values come from the ISCWSA error model:
    // http://www.copsegrove.com/Pages/MWDGeomagneticModels.aspx

    static constexpr double INCL_ERROR_BASE        = 0.20;
    static constexpr double DECL_ERROR_OFFSET_BASE = 0.36;
    static constexpr int F_ERROR_BASE              = 130;
    static constexpr int DECL_ERROR_SLOPE_BASE     = 5000;
    static constexpr double WMM_ERROR_MULTIPLIER   = 1.21;
    static constexpr double IGRF_ERROR_MULTIPLIER  = 1.21;

    static constexpr int WMMHR_UNCERTAINTY_F           = 134;
    static constexpr int WMMHR_UNCERTAINTY_H           = 130;
    static constexpr int WMMHR_UNCERTAINTY_X           = 135;
    static constexpr int WMMHR_UNCERTAINTY_Y           = 85;
    static constexpr int WMMHR_UNCERTAINTY_Z           = 134;
    static constexpr double WMMHR_UNCERTAINTY_I        = 0.19;
    static constexpr double WMMHR_UNCERTAINTY_D_OFFSET = 0.25;
    static constexpr int WMMHR_UNCERTAINTY_D_COEF      = 5205;

    // These error values are the NCEI error model
    static constexpr int WMM_UNCERTAINTY_F           = 138;
    static constexpr int WMM_UNCERTAINTY_H           = 133;
    static constexpr int WMM_UNCERTAINTY_X           = 137;
    static constexpr int WMM_UNCERTAINTY_Y           = 89;
    static constexpr int WMM_UNCERTAINTY_Z           = 141;
    static constexpr double WMM_UNCERTAINTY_I        = 0.20;
    static constexpr double WMM_UNCERTAINTY_D_OFFSET = 0.26;
    static constexpr int WMM_UNCERTAINTY_D_COEF      = 5417;

    static constexpr int PS_MIN_LAT_DEGREE     = -55;    // Minimum Latitude for  Polar Stereographic projection in degrees
    static constexpr int PS_MAX_LAT_DEGREE     = 55;     // Maximum Latitude for Polar Stereographic projection in degrees
    static constexpr double UTM_MIN_LAT_DEGREE = -80.5;  // Minimum Latitude for UTM projection in degrees
    static constexpr double UTM_MAX_LAT_DEGREE = 84.5;   // Maximum Latitude for UTM projection in degrees

    static constexpr double GEO_POLE_TOLERANCE = 1e-5;

    static constexpr int LAT_BOUND_MIN         = -90;
    static constexpr int LAT_BOUND_MAX         = 90;
    static constexpr int LON_BOUND_MIN         = -180;
    static constexpr int LON_BOUND_MAX         = 360;
    static constexpr int ALT_BOUND_MIN         = -10;
    static constexpr int NO_ALT_MAX            = -99999;
    static constexpr int USER_GAVE_UP          = -1;
    static constexpr double DEC_YEAR_BOUND_MIN = 2024.866;
    static constexpr double DEC_YEAR_BOUND_MAX = 2030;

}  // namespace wmm