//
// Created by dasha on 25.04.2023.
//

#ifndef ENGINE_FACE_H
#define ENGINE_FACE_H
#include <iostream>
#include <vector>

class Vector3D;
class Face
{
public:
    Face(){};
    Face(std::vector<int> inds){
        point_indexes = inds;
    }
    std::vector<int> point_indexes;
    std::vector<std::vector<double> > uv;
    std::vector<Vector3D* > norm;
    int map_Ka = -1;
    int map_Kd = -1;
    int map_Ks = -1;
    static Vector3D getNorm(const Vector3D &A, const Vector3D& B, const Vector3D& C, Vector3D* norma, Vector3D* normb, Vector3D* normc, const Vector3D& point);
};


#endif //ENGINE_FACE_H
