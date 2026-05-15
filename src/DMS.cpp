//
// Created by Kiryuhin Viacheslav on 30.04.2026.
//

#include "DMS.h"

#include <cmath>

namespace wmm {

    bool DMS::isValid() const {
        return std::isfinite(angle) && !std::isnan(angle);
        //               && !(m_isLatitude && m_angle > 90.0) &&               !(!m_isLatitude && m_angle > 180.0);
    }
}  // namespace wmm