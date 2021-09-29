//
// Created by Johannes Martin on 20.09.21.
//

#ifndef PARALOBSTAR_BOX_H
#define PARALOBSTAR_BOX_H

#include "global.h"
#include "Particle.h"

class Box {
public:
    double lower[global::dim] {};
    double upper[global::dim] {};

    int sonBoxAndIndex(Box &sonBox, Particle &p);
    double getLength();
    bool particleWithin(Particle &p);
};

#endif //PARALOBSTAR_BOX_H
