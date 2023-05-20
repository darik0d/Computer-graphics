//
// Created by dasha on 12.05.2023.
//

#include <cmath>
#include "utils.h"
#include "vector3d.h"
#include "Figure.h"
#include "Face.h"
#include "Line2D.h"

void utils::applyTransformation(Figure &fig, const Matrix& finalTrans){
    for(Vector3D & point: fig.points) {
        point = point*finalTrans;
    }
}
Point2D utils::doProjection(const Vector3D &point, const double d){
    Point2D to_return;
    to_return.x = (d*point.x)/(-point.z);
    to_return.y = (d*point.y)/(-point.z);
    if(point.z == 0) {
        to_return.x = 0;
        to_return.y = 0;
    }
    to_return.z = point.z;
    return to_return;
}
Lines2D utils::doProjection(const Figure & figuur){
    Lines2D to_return;
    std::vector<Point2D> currentPoints;
    for(auto point: figuur.points){
        currentPoints.push_back(doProjection(point, 1));
    }
    // Now we have to work with faces. We add a new line with color of the figure (line goes from prev to current point)
    for(auto face: figuur.faces){
        int currentIndex = face.point_indexes[0];
        int nextIndex = face.point_indexes[1];
        if(face.point_indexes.size() == 2) to_return.push_back(Line2D(currentPoints[currentIndex],currentPoints[nextIndex], figuur.color, currentPoints[currentIndex].z, currentPoints[nextIndex].z));
        else {
            int i = 0;
            while(i < face.point_indexes.size()-1){
                to_return.push_back(Line2D(currentPoints[face.point_indexes[i]],currentPoints[face.point_indexes[i+1]], figuur.color, currentPoints[face.point_indexes[i]].z,currentPoints[face.point_indexes[i+1]].z));
                i++;
            }
            to_return.push_back(Line2D(currentPoints[face.point_indexes[face.point_indexes.size()-1]],currentPoints[face.point_indexes[0]], figuur.color, currentPoints[face.point_indexes[face.point_indexes.size()-1]].z,currentPoints[face.point_indexes[0]].z));
        }
    }
    return to_return;
}
std::pair<double,double> utils::getMinimum(const Lines2D &lines){
    double x = lines[0].a.x;
    double y = lines[0].a.y;
    for(auto lijn: lines){
        if (std::isinf(lijn.a.x) || std::isinf(-lijn.a.x)){
            lijn.a.x = 0;
            std::cout << "Inf detected \n";
        }
        if (std::isinf(lijn.a.y) || std::isinf(-lijn.a.y)){
            lijn.a.y = 0;
            std::cout << "Inf detected \n";
        }
        if (lijn.a.x < x) x = lijn.a.x;
        if (lijn.b.x < x) x = lijn.b.x;
        if (lijn.a.y < y) y = lijn.a.y;
        if (lijn.b.y < y) y = lijn.b.y;
    }
    std::pair<double,double> to_return = std::make_pair(x,y);
    return to_return;
}
std::pair<double,double> utils::getMaximum(const Lines2D &lines){
    double x = lines[0].a.x;
    double y = lines[0].a.y;
    for(auto lijn: lines){
        if (lijn.a.x > x) x = lijn.a.x;
        if (lijn.b.x > x) x = lijn.b.x;
        if (lijn.a.y > y) y = lijn.a.y;
        if (lijn.b.y > y) y = lijn.b.y;
    }
    std::pair<double,double> to_return = std::make_pair(x,y);
    return to_return;
}
void utils::getDxDyD(const Lines2D &lines, const int size, double &dx, double &dy, double &d){
    double x_min = getMinimum(lines).first;
    double y_min = getMinimum(lines).second;
    double x_max = getMaximum(lines).first;
    double y_max = getMaximum(lines).second;
    // Bereken x_range en y_range
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    // Bereken imagex
    int imagex = static_cast<int> (std::round(size*x_range/std::max(x_range, y_range)));
    // Bereken imagey
    int imagey = static_cast<int> (std::round(size*y_range/std::max(x_range, y_range)));
    // Bereken d
    d = 0.95*imagex/x_range;
    // Bereken DCx
    double dcx = d*(x_min+x_max)/2.0;
    // Bereken DCy
    double dcy = d*(y_min+y_max)/2.0;
    //dx
    dx = imagex/2.0 - dcx;
    //dy
    dy = imagey/2.0 - dcy;
}
void utils::applyTransformation(std::vector<Figure> &figs, const Matrix &m){
    for(auto &fig: figs){
        for(auto &point: fig.points) point = point * m;
    }
}
void utils::doProjection(const std::vector<Figure> & figures, Lines2D& to_return){
    for(const Figure& fig: figures){
        Lines2D to_insert = doProjection(fig);
        to_return.insert(to_return.end(), to_insert.begin(), to_insert.end());
    }
}
double utils::getCubeSizeRadius(const std::vector<Figure> & all_projected_figures){
    double x_min = -1;
    double y_min = -1;
    double z_min = -1;
    double x_max = 1;
    double y_max = 1;
    double z_max = 1;
    for(const Figure& figure:all_projected_figures){
        for(const Vector3D& point: figure.points){
            if(x_min > point.x) x_min = point.x;
            if(y_min > point.y) y_min = point.y;
            if(z_min > point.z) z_min = point.z;
            if(x_max < point.x) x_max = point.x;
            if(y_max < point.y) y_max = point.y;
            if(z_max < point.z) z_max = point.z;
        }
    }
    x_min = std::abs(x_min);
    y_min = std::abs(y_min);
    z_min = std::abs(z_min);
    x_max = std::abs(x_max);
    y_max = std::abs(y_max);
    z_max = std::abs(z_max);
    double c = std::max(std::max(std::max(x_min,x_max),y_min),std::max(std::max(y_max,z_min),z_max));
    //c /= 0.95;
    return c;
}