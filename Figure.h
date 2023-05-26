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

//const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
class Figure
{
public:
    Figure();

    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
    std::vector<int> textureNrs;
    int cubeMapNr = -1;
    std::vector<double> fullAmbientReflection;
    std::vector<double> diffuseReflection;
    std::vector<double> specularReflection;
    double reflectionCoefficient;
    void cube();
    void tetrahedron();
    void octahedron();
    void icosahedron();
    void dodecahedron();
    void sphere(int n);
    void buckyBall();
    void mengerSponge();
    void cylinderZonderVlakken(int n, double height);
    static std::vector<Figure> generateThickFigure(const Figure& to_enlarge, const double r, const int n, const int m);
    void generateFractal(std::vector<Figure> & fractal, const int nr_iterations, const double scale) const;
    void scaleFigure(const double scale);
    void translate(const Vector3D &vector);
    static Figure copyFigure(const Figure& fig);
    static Figure parseObj(const std::string &src, std::vector<Texture *> &vector, Matrix eyeTransf, bool normalOn);
    void herschaalPuntenBal(Vector3D& punt);
//    Figure(Figure& fig);
};


#endif //ENGINE_FIGURE_H
