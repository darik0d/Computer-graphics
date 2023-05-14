//
// Created by dasha on 22.03.2023.
//

#include "ZBuffer.h"
#include "limits"

double posInf = std::numeric_limits<double>::infinity();
double negInf = -std::numeric_limits<double>::infinity();

ZBuffer::ZBuffer(const int width, const int height){
    for(int h = 0; h < height; h++){
        std::vector<double> newVec;
        for(int w = 0; w < width; w++){
            newVec.push_back(posInf);
        }
        this->push_back(newVec);
    }
    min = posInf;
    max = negInf;
}

ZBuffer::ZBuffer() {
    min = posInf;
    max = negInf;}
