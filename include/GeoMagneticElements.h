//
// Created by Kiryuhin Viacheslav on 22.04.2026.
//

#pragma once

namespace wmm {
    struct Ellipsoid;
    struct CoordGeodetic;
    struct CoordSpherical;
    struct MagneticModel;
    struct MagneticResults;

    struct GeoMagneticElements {
        void calculate(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
                       MagneticModel &timedMagneticModel);

        void calculate(const Ellipsoid &ellip, const CoordGeodetic &coordGeodetic, MagneticModel &timedMagneticModel);

        void gradY(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
                   const MagneticModel &timedMagneticModel, const GeoMagneticElements &geoMagneticElements);

        void calculateGridVariation(const CoordGeodetic &location);

        void WMMerrorCalc(double H_);

        [[maybe_unused]] void WMMHRerrorCalc(double H_);

        [[maybe_unused]] [[nodiscard]] GeoMagneticElements scaled(double factor) const;

        void scale(double factor);

        GeoMagneticElements operator-(const GeoMagneticElements &subtrahend) const;

        [[maybe_unused]] void errorCalc(const GeoMagneticElements &B);

        double Decl{0.0};     // 1. Angle between the magnetic field vector and true north, positive east
        double Incl{0.0};     // 2. Angle between the magnetic field vector and the horizontal plane, positive down
        double F{0.0};        // 3. Magnetic Field Strength
        double H{0.0};        // 4. Horizontal Magnetic Field Strength
        double X{0.0};        // 5. Northern component of the magnetic field vector
        double Y{0.0};        // 6. Eastern component of the magnetic field vector
        double Z{0.0};        // 7. Downward component of the magnetic field vector
        double GV{0.0};       // 8. The Grid Variation
        double Decldot{0.0};  // 9. Yearly Rate of change in declination
        double Incldot{0.0};  // 10. Yearly Rate of change in inclination
        double Fdot{0.0};     // 11. Yearly rate of change in Magnetic field strength
        double Hdot{0.0};     // 12. Yearly rate of change in horizontal field strength
        double Xdot{0.0};     // 13. Yearly rate of change in the northern component
        double Ydot{0.0};     // 14. Yearly rate of change in the eastern component
        double Zdot{0.0};     // 15. Yearly rate of change in the downward component
        double GVdot{0.0};    // 16. Yearly rate of change in grid variation

    private:
        void calculate(const MagneticResults &magneticResultsGeo);

        void calculateSecularVariation(const MagneticResults &magneticVariation);

        void calculateGradientElements(const MagneticResults &gradResults, const GeoMagneticElements &magElem);
    };

    void BaseErrors(double declCoef, double declBaseline, double inclOffset, double fOffset, double multiplier, double H,
                    double &declErr, double &inclErr, double &fErr);
}  // namespace wmm