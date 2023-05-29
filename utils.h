//
// Created by dasha on 12.05.2023.
//

#ifndef ENGINE_UTILS_H
#define ENGINE_UTILS_H

#include <vector>
#include "Point2D.h"
#include "vector3d.h"

class Figure;
class Matrix;
class Line2D;
typedef std::vector<Line2D> Lines2D;

struct utils {
    static void applyTransformation(Figure &fig, const Matrix& finalTrans);
    static Lines2D doProjection(const Figure & figuur);
    static void doProjection(const std::vector<Figure> & figures, Lines2D& to_return);
    static void getDxDyD(const Lines2D &lines, int size, double &dx, double &dy, double &d);
    static std::pair<double,double> getMinimum(const Lines2D &lines);
    static std::pair<double,double> getMaximum(const Lines2D &lines);
    static Point2D doProjection(const Vector3D &point, const double d);
    static void applyTransformation(std::vector<Figure> &figs, const Matrix &m);
    static double getCubeSizeRadius(const std::vector<Figure> & all_projected_figures);
    static void getExtrema(const std::vector<Figure> & all_projected_figures, double& x_min, double& x_max, double& y_min, double& y_max, double& z_min, double& z_max);
    static Vector3D findMiddle(Vector3D a, Vector3D b);
    static std::vector<Figure> generateThickFigure(const Figure& to_enlarge, const double r, const int n, const int m);
    static void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
    static double to_radialen(double graden);
    static double to_degrees(double graden);
    static Matrix rotateX(const double angle);
    static Matrix rotateY(const double cangle);
    static Matrix rotateZ(const double cangle);
};


#endif //ENGINE_UTILS_H
