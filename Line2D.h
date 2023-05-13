//
// Created by dasha on 13.05.2023.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H
#include "Point2D.h"
#include "easy_image.h" // Ik hoop dat het geen circulaire import zal geven...

class Line2D{
public:
    Point2D a;
    Point2D b;
    img::Color color;
    double z1;
    double z2;
    Line2D(){};
    Line2D(Point2D ap, Point2D bp, img::Color kleur){
        a=ap;
        b=bp;
        color = kleur;
    }
    Line2D(Point2D ap, Point2D bp, img::Color kleur, double z_1, double z_2){
        a=ap;
        b=bp;
        color = kleur;
        z1 = z_1;
        z2 = z_2;
    }

    ~Line2D() = default;
};

#endif //ENGINE_LINE2D_H
