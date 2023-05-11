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
    std::vector<double> ambientLight;
    //de diffuse licht component
    std::vector<double> diffuseLight;
    //de diffuse licht component
    std::vector<double> specularLight;
    // Infinitive?
    bool inf;
    //de richting waarin het
    //licht schijnt (of location als puntbron)
    Vector3D ldVector = Vector3D::vector(0,0,0);
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
};



#endif //ENGINE_LIGHT_H
