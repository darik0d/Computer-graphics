//
// Created by dasha on 31.03.2023.
//

#ifndef ENGINE_POINT2D_H
#define ENGINE_POINT2D_H
#include <iostream>
class Point2D{
public:
    Point2D(){};
    Point2D(double xy, double yy){
        x=xy;
        y=yy;
    }
    Point2D(std::pair<double,double> co){
        x = co.first;
        y = co.second;
    }
    Point2D(double co1, double co2, double co3){
        x = co1;
        y = co2;
        z = co3;
    }
    bool operator!=(Point2D p){
        if(x == p.x && y == p.y && z == p.z) return true;
        return false;
    }
    double x;
    double y;
    double z;
};
#endif //ENGINE_POINT2D_H
