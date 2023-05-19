//
// Created by dasha on 14.05.2023.
//

#include "Texture.h"
#include <math.h>

Texture::~Texture() {
    delete image;
    //std::cout << "The texture is destroyed!" << std::endl;
}
Vector3D Texture::getUVW(const Vector3D& xyz) const{
    Vector3D c = Vector3D::cross(a, b);
    Vector3D uvw;
    Matrix to_inverse;
    to_inverse(1,1) = a.x;
    to_inverse(1,2) = a.y;
    to_inverse(1,3) = a.z;
    to_inverse(2,1) = b.x;
    to_inverse(2,2) = b.y;
    to_inverse(2,3) = b.z;
    to_inverse(3,1) = c.x;
    to_inverse(3,2) = c.y;
    to_inverse(3,3) = c.z;
    Vector3D xminp = xyz - p;
    Matrix inverted = Matrix::inv(to_inverse);
    uvw = xminp*inverted;
    // TODO: check if uvw coos are always in [0,1]
    // TODO: check if uvw is a vector. No, it's a point
    // bool testpoint = uvw.is_point();
    return uvw;
}
