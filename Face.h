//
// Created by dasha on 25.04.2023.
//

#ifndef ENGINE_FACE_H
#define ENGINE_FACE_H
#include <iostream>
#include <vector>

class Face
{
public:
    Face(){};
    Face(std::vector<int> inds){
        point_indexes = inds;
    }
    std::vector<int> point_indexes;
    std::vector<std::vector<double> > uv;
    int map_Ka = -1;
    int map_Kd = -1;
    int map_Ks = -1;
};


#endif //ENGINE_FACE_H
