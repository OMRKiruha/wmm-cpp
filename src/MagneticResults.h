//
// Created by Professional on 22.04.2026.
//

#pragma once

namespace wmm {
    struct LegendreFunction;
    struct MagneticModel;
    struct SphericalHarmonicVariables;
    struct CoordSpherical;
    struct CoordGeodetic;

    struct MagneticResults {
        void Summation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                       const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical);

        void SecVarSummation(const LegendreFunction &legendreFunction, MagneticModel &magneticModel,
                             const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical);

        void RotateMagneticVector(const CoordSpherical &CoordSpherical, const CoordGeodetic &coordGeodetic,
                                  const MagneticResults &magneticResultsSph);

        void GradYSummation(const LegendreFunction &LegendreFunction, const MagneticModel &magneticModel,
                            const SphericalHarmonicVariables &SphVariables, const CoordSpherical &coordSpherical);

        double Bx{};  // North
        double By{};  // East
        double Bz{};  // Down

    private:
        void SummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                              const CoordSpherical &coordSpherical);

        void SecVarSummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                                    const CoordSpherical &coordSpherical);
    };

}  // namespace wmm