//
// Created by dasha on 25.04.2023.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H
#include <iostream>
#include <vector>
#include "easy_image.h"

class Vector3D;
class Face;


class Figure
{
public:
    Figure();

    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
    img::Color fullAmbientReflection;
    img::Color diffuseReflection;
    img::Color specularReflection;
    double reflectionCoefficient;
    void cube();
    void tetrahedron();
    void octahedron();
    void dodecahedron();
    void generateFractal(std::vector<Figure> & fractal, const int nr_iterations, const double scale) const;
    void scaleFigure(const double scale);
    void translate(const Vector3D &vector);
//    Figure(Figure& fig);
};


#endif //ENGINE_FIGURE_H
