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
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
    void cube();
    void tetrahedron();
    void octahedron();
    void dodecahedron();
};


#endif //ENGINE_FIGURE_H
