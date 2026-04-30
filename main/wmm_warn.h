//
// Created by Kiryuhin Viacheslav on 23.04.2026.
//

#pragma once

constexpr auto BOZ_WARN_TEXT_STRONG =
    "Warning: some calculated locations are in the blackout zone around the magnetic pole\n"
    "as defined by the WMM military specification \n"
    "(https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf).\n"
    "Compass accuracy is highly degraded in this region.\n";
constexpr auto BOZ_WARN_TEXT_WEAK = "Caution: some calculated locations approach the blackout zone around the magnetic\n"
                                    "pole as defined by the WMM military specification \n"
                                    "(https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf).\n"
                                    "Compass accuracy may be degraded in this region.\n";

constexpr auto WMM_MileSpec_INFO =
    "Warning: The height validity of the geomagnetic components is dependent on the geomagnetic activity level. For more "
    "information see \n"
    "(https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model)\n";
constexpr auto WMM_MileSpec_WARN =
    "Warning: WMM will not meet MilSpec at this altitude. For more information see \n"
    "(https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model)\n";