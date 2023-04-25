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
};


#endif //ENGINE_FACE_H
