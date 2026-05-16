//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

namespace wmm {
    struct LegendreFunction;
    struct MagneticModel;
    struct SphericalHarmonicVariables;
    struct CoordSpherical;
    struct CoordGeodetic;

    struct MagneticResults {
        void summation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                       const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical);

        void secVarSummation(const LegendreFunction &legendreFunction, MagneticModel &magneticModel,
                             const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical);

        void rotateMagneticVector(const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
                                  const MagneticResults &magneticResultsSph);

        void GradYSummation(const LegendreFunction &legendreFunction, const MagneticModel &magneticModel,
                            const SphericalHarmonicVariables &sphVariables, const CoordSpherical &coordSpherical);

        double Bx{0.0};  // North
        double By{0.0};  // East
        double Bz{0.0};  // Down

    private:
        void summationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                              const CoordSpherical &coordSpherical);

        void secVarSummationSpecial(const MagneticModel &magneticModel, const SphericalHarmonicVariables &sphVariables,
                                    const CoordSpherical &coordSpherical);
    };

}  // namespace wmm