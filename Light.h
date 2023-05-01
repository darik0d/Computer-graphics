//
// Created by dasha on 01.05.2023.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "easy_image.h"
#include "vector3d.h"

class Light
{
public:
    //de ambiente licht component
    img::Color ambientLight;
    //de diffuse licht component
    img::Color diffuseLight;
    //de diffuse licht component
    img::Color specularLight;
};

class InfLight: public Light
{
public:
    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;
};

class PointLight: public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
};



#endif //ENGINE_LIGHT_H
