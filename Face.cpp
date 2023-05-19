//
// Created by dasha on 25.04.2023.
//

#include "Face.h"
#include "vector3d.h"

Vector3D Face::getNorm(const Vector3D &A, const Vector3D& B, const Vector3D& C, Vector3D* norma, Vector3D* normb, Vector3D* normc, const Vector3D& point) {
    Matrix to_inverse;
    to_inverse(1,1) =  A.x;
    to_inverse(1,2) =  A.y;
    to_inverse(1,3) =  A.z;
    to_inverse(2,1) =  B.x;
    to_inverse(2,2) =  B.y;
    to_inverse(2,3) =  B.z;
    to_inverse(3,1) =  C.x;
    to_inverse(3,2) =  C.y;
    to_inverse(3,3) =  C.z;
    Matrix inverted = Matrix::inv(to_inverse);
    // Get 1-p-q, p, q
    Vector3D coefs = point*inverted;
    double n_x = coefs.x*norma->x + coefs.y*normb->x + coefs.z*normc->x;
    double n_y = coefs.x*norma->y + coefs.y*normb->y + coefs.z*normc->y;
    double n_z = coefs.x*norma->z + coefs.y*normb->z + coefs.z*normc->z;
    return Vector3D::vector(n_x, n_y, n_z);
}