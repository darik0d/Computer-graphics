//
// Created by dasha on 25.04.2023.
//

#include "Figure.h"
#include "Face.h"
#include "vector3d.h"
#include <cmath>
#define _USE_MATH_DEFINES

Figure::Figure() {}
//Figure::Figure(Figure& fig){
//    points = std::vector<Vector3D> (fig.points);
//    faces = std::vector<Face> (fig.faces);
//    color = fig.color;
//}
void Figure::scaleFigure(const double scale){
    Matrix m;
    for(auto i:{1,2,3}){
        m(i,i) = m(i,i)*scale;
    }
    // Apply scale
    for(Vector3D & point: points) {
        point = point*m;
    }
}
void Figure::translate(const Vector3D &vector){
    Matrix to_return;
    to_return(4,1) = vector.x;
    to_return(4,2) = vector.y;
    to_return(4,3) = vector.z;
    // Apply scale
    for(Vector3D & point: points) {
        point = point*to_return;
    }
}
void Figure::generateFractal(std::vector<Figure> & fractal, const int nr_iterations, const double scale) const{
    std::vector<Figure> figures_to_fractalise;
    // Copy constructor
    Figure figure;
    figure.color = color;
    figure.faces = std::vector<Face>(faces);
    figure.points = std::vector<Vector3D>(points);
    figures_to_fractalise.push_back(figure);
    for(int n = 0; n < nr_iterations; n++) {
        for(const auto& fig:figures_to_fractalise) {
            for (int i = 0; i < fig.points.size(); i++) {
                Figure fract_fig;
                // Copy constructor
                fract_fig.color = fig.color;
                fract_fig.faces = std::vector<Face>(fig.faces);
                fract_fig.points = std::vector<Vector3D>(fig.points);
                // Scale figure
                fract_fig.scaleFigure(1 / scale);
                // Calculate vector
                Vector3D to_move = fig.points[i] - fract_fig.points[i];
                // Move
                fract_fig.translate(to_move);
                // Add figure to fractal
                fractal.push_back(fract_fig);
                // Do it for more iterations
            }
        }
        figures_to_fractalise = fractal;
        fractal.clear();
    }
    fractal = figures_to_fractalise;
}
void Figure::cube(){
    points.push_back(Vector3D::point(1,-1,-1));
    points.push_back(Vector3D::point(-1,1,-1));
    points.push_back(Vector3D::point(1,1,1));
    points.push_back(Vector3D::point(-1,-1,1));
    points.push_back(Vector3D::point(1,1,-1));
    points.push_back(Vector3D::point(-1,-1,-1));
    points.push_back(Vector3D::point(1,-1,1));
    points.push_back(Vector3D::point(-1,1,1));

    Face f1;
    f1.point_indexes.push_back(0);
    f1.point_indexes.push_back(4);
    f1.point_indexes.push_back(2);
    f1.point_indexes.push_back(6);
    Face f2;
    f2.point_indexes.push_back(4);
    f2.point_indexes.push_back(1);
    f2.point_indexes.push_back(7);
    f2.point_indexes.push_back(2);
    Face f3;
    f3.point_indexes.push_back(1);
    f3.point_indexes.push_back(5);
    f3.point_indexes.push_back(3);
    f3.point_indexes.push_back(7);
    Face f4;
    f4.point_indexes.push_back(5);
    f4.point_indexes.push_back(0);
    f4.point_indexes.push_back(6);
    f4.point_indexes.push_back(3);
    Face f5;
    f5.point_indexes.push_back(6);
    f5.point_indexes.push_back(2);
    f5.point_indexes.push_back(7);
    f5.point_indexes.push_back(3);
    Face f6;
    f6.point_indexes.push_back(0);
    f6.point_indexes.push_back(5);
    f6.point_indexes.push_back(1);
    f6.point_indexes.push_back(4);

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
}
void Figure::tetrahedron(){
    points.push_back(Vector3D::point(1,-1,-1));
    points.push_back(Vector3D::point(-1,1,-1));
    points.push_back(Vector3D::point(1,1,1));
    points.push_back(Vector3D::point(-1,-1,1));

    Face f1;
    f1.point_indexes.push_back(0);
    f1.point_indexes.push_back(1);
    f1.point_indexes.push_back(2);
    Face f2;
    f2.point_indexes.push_back(1);
    f2.point_indexes.push_back(3);
    f2.point_indexes.push_back(2);
    Face f3;
    f3.point_indexes.push_back(0);
    f3.point_indexes.push_back(3);
    f3.point_indexes.push_back(1);
    Face f4;
    f4.point_indexes.push_back(0);
    f4.point_indexes.push_back(2);
    f4.point_indexes.push_back(3);

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
}
void Figure::octahedron(){
    points.push_back(Vector3D::point(1,0,0));
    points.push_back(Vector3D::point(0,1,0));
    points.push_back(Vector3D::point(-1,0,0));
    points.push_back(Vector3D::point(0,-1,0));
    points.push_back(Vector3D::point(0,0,-1));
    points.push_back(Vector3D::point(0,0,1));

    Face f1;
    f1.point_indexes.push_back(0);
    f1.point_indexes.push_back(1);
    f1.point_indexes.push_back(5);
    Face f2;
    f2.point_indexes.push_back(1);
    f2.point_indexes.push_back(2);
    f2.point_indexes.push_back(5);
    Face f3;
    f3.point_indexes.push_back(2);
    f3.point_indexes.push_back(3);
    f3.point_indexes.push_back(5);
    Face f4;
    f4.point_indexes.push_back(3);
    f4.point_indexes.push_back(0);
    f4.point_indexes.push_back(5);
    Face f5;
    f5.point_indexes.push_back(1);
    f5.point_indexes.push_back(0);
    f5.point_indexes.push_back(4);
    Face f6;
    f6.point_indexes.push_back(2);
    f6.point_indexes.push_back(1);
    f6.point_indexes.push_back(4);
    Face f7;
    f7.point_indexes.push_back(3);
    f7.point_indexes.push_back(2);
    f7.point_indexes.push_back(4);
    Face f8;
    f8.point_indexes.push_back(0);
    f8.point_indexes.push_back(3);
    f8.point_indexes.push_back(4);

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
}
void Figure::dodecahedron() {
    std::vector<Vector3D> icoPoints;
    icoPoints.push_back(Vector3D::point(0,0, std::sqrt(5)/2));
    for(int l = 2; l < 7; l++){
        icoPoints.push_back(Vector3D::point(std::cos(2*M_PI*(l-2)/5), std::sin(2*M_PI*(l-2)/5), 0.5));
    }
    for(int l = 7; l < 12; l++){
        icoPoints.push_back(Vector3D::point(std::cos((M_PI/5)+((l-7)*2*M_PI/5)), std::sin((M_PI/5)+((l-7)*2*M_PI/5)), -0.5));
    }
    icoPoints.push_back(Vector3D::point(0,0, -std::sqrt(5)/2));
    // Here i'd like to give a point from old faces
    std::vector<std::vector<int>> icoInd = {{0,1,2}, {0,2,3}, {0,3,4},
                                            {0,4,5}, {0,5,1}, {1,6,2},
                                            {2,6,7}, {2,7,3}, {3,7,8},
                                            {3,8,4}, {4,8,9}, {4,9,5},
                                            {5,9,10}, {5,10,1}, {1,10,6},
                                            {11,7,6}, {11,8,7}, {11,9,8},
                                            {11,10,9}, {11,6,10}};
    for(auto point: icoInd){
        points.push_back(Vector3D::point((icoPoints[point[0]].x + icoPoints[point[1]].x + icoPoints[point[2]].x)/3,
                                         (icoPoints[point[0]].y + icoPoints[point[1]].y + icoPoints[point[2]].y)/3,
                                         (icoPoints[point[0]].z + icoPoints[point[1]].z + icoPoints[point[2]].z)/3));
    }
    Face f1 = Face({0,1,2,3,4});
    Face f2 = Face({0,5,6,7,1});
    Face f3 = Face({1,7,8,9,2});
    Face f4 = Face({2,9,10,11,3});
    Face f5 = Face({3,11,12,13,4});
    Face f6 = Face({4,13,14,5,0});
    Face f7 = Face({19,18,17,16,15});
    Face f8;
    f8.point_indexes.push_back(19);
    f8.point_indexes.push_back(14);
    f8.point_indexes.push_back(13);
    f8.point_indexes.push_back(12);
    f8.point_indexes.push_back(18);
    Face f9;
    f9.point_indexes.push_back(18);
    f9.point_indexes.push_back(12);
    f9.point_indexes.push_back(11);
    f9.point_indexes.push_back(10);
    f9.point_indexes.push_back(17);
    Face f10;
    f10.point_indexes.push_back(17);
    f10.point_indexes.push_back(10);
    f10.point_indexes.push_back(9);
    f10.point_indexes.push_back(8);
    f10.point_indexes.push_back(16);
    Face f11;
    f11.point_indexes.push_back(16);
    f11.point_indexes.push_back(8);
    f11.point_indexes.push_back(7);
    f11.point_indexes.push_back(6);
    f11.point_indexes.push_back(15);
    Face f12;
    f12.point_indexes.push_back(15);
    f12.point_indexes.push_back(6);
    f12.point_indexes.push_back(5);
    f12.point_indexes.push_back(14);
    f12.point_indexes.push_back(19);

    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
    faces.push_back(f9);
    faces.push_back(f10);
    faces.push_back(f11);
    faces.push_back(f12);
}
