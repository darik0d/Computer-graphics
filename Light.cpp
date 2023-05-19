#include <cmath>
#include "Light.h"

bool Light::pointIsVisible(double x_ac, double y_ac, double z_e, double eye_d, double eye_dx, double eye_dy, const Matrix &eyeTransf) const{
    double x_e = -z_e*(x_ac - eye_dx)/eye_d;
    double y_e = -z_e*(y_ac - eye_dy)/eye_d;
    Vector3D worldCo = Vector3D::point(x_e, y_e, z_e)*Matrix::inv(eyeTransf);
    Vector3D lightCo = worldCo*eye;
    double x_lac = (-d*lightCo.x)/lightCo.z + dx;
    double y_lac = (-d*lightCo.y)/lightCo.z + dy;
    double een_z_l = 1/lightCo.z;
//    if(y_lac > shadowMask.size()-1 || x_lac > shadowMask.size()-1) return false;
//    if(y_lac < 0 || x_lac < 0) return false;
    // Get all 1/z neighbours
    double een_z_a = shadowMask[std::ceil(y_lac)][std::floor(x_lac)];
    double een_z_b = shadowMask[std::ceil(y_lac)][std::ceil(x_lac)];
    double een_z_c = shadowMask[std::floor(y_lac)][std::floor(x_lac)];
    double een_z_d = shadowMask[std::floor(y_lac)][std::ceil(x_lac)];
    // Get a_x and a_y
    double a_x = x_lac - std::floor(x_lac);
    double a_y = y_lac - std::floor(y_lac);
    // een_z_e betekent 1/z van punt E
    double een_z_e = (1 - a_x)*een_z_a + a_x*een_z_b;
    double een_z_f = (1 - a_x)*een_z_c + a_x*een_z_d;

    double diepte = a_y*een_z_e + (1 - a_y)*een_z_f;
    //if(diepte == INFINITY || een_z_l == INFINITY) return false;
    if(diepte > een_z_l) return true;
    if(std::abs(diepte - een_z_l) < std::pow(10, -5)) return true;
    return false;
}

Light::~Light() {
    //std::cout << "The light is destroyed!" << std::endl;
}
