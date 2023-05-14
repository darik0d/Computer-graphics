//
// Created by dasha on 14.05.2023.
//

#ifndef ENGINE_TEXTURE_H
#define ENGINE_TEXTURE_H

#include "vector3d.h"
#include "easy_image.h"

struct Texture {
    int number;
//    int width;
//    int height;
    Vector3D p, a, b;
    img::EasyImage* image;
    Vector3D getUVW(const Vector3D& xyz) const;
    virtual ~Texture();
};


#endif //ENGINE_TEXTURE_H
