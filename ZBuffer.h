//
// Created by dasha on 22.03.2023.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include "vector"

class ZBuffer: public std::vector<std::vector<double>>
{
public:
//Constructor: maakt een Z-Buffer van de correcte
//grootte aan en initialiseert alle velden op +inf
    //std::vector<std::vector<double>> z_waarden;
    double min;
    double max;
    ZBuffer();

    ZBuffer(const int width, const int height);
};


#endif //ENGINE_ZBUFFER_H
