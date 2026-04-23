//
// Created by Professional on 22.04.2026.
//

#pragma once

namespace wmm {
    struct Ellipsoid;
    struct CoordGeodetic;
    struct CoordSpherical;
    struct MagneticModel;
    struct MagneticResults;

    struct GeoMagneticElements {
        void geomag(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
                    MagneticModel &timedMagneticModel);

        void gradY(const Ellipsoid &ellip, const CoordSpherical &coordSpherical, const CoordGeodetic &coordGeodetic,
                   MagneticModel &timedMagneticModel, const GeoMagneticElements &geoMagneticElements);

        int CalculateGridVariation(const CoordGeodetic &location);

        void WMMErrorCalc(double H_);

        void WMMHRErrorCalc(double H_);

        [[nodiscard]] GeoMagneticElements scaled(double factor) const;

        void scale(double factor);

        GeoMagneticElements operator-(const GeoMagneticElements &subtrahend) const;

        void ErrorCalc(GeoMagneticElements B);

        double Decl{};     // 1. Angle between the magnetic field vector and true north, positive east
        double Incl{};     // 2. Angle between the magnetic field vector and the horizontal plane, positive down
        double F{};        // 3. Magnetic Field Strength
        double H{};        // 4. Horizontal Magnetic Field Strength
        double X{};        // 5. Northern component of the magnetic field vector
        double Y{};        // 6. Eastern component of the magnetic field vector
        double Z{};        // 7. Downward component of the magnetic field vector
        double GV{};       // 8. The Grid Variation
        double Decldot{};  // 9. Yearly Rate of change in declination
        double Incldot{};  // 10. Yearly Rate of change in inclination
        double Fdot{};     // 11. Yearly rate of change in Magnetic field strength
        double Hdot{};     // 12. Yearly rate of change in horizontal field strength
        double Xdot{};     // 13. Yearly rate of change in the northern component
        double Ydot{};     // 14. Yearly rate of change in the eastern component
        double Zdot{};     // 15. Yearly rate of change in the downward component
        double GVdot{};    // 16. Yearly rate of change in grid variation

    private:
        void calculate(const MagneticResults &magneticResultsGeo);

        void calculateSecularVariation(const MagneticResults &magneticVariation);

        void calculateGradientElements(const MagneticResults &gradResults, const GeoMagneticElements &magneticElements);
    };
}  // namespace wmm