//
// Created by Johannes Martin on 21.09.21.
//

#ifndef PARALOBSTAR_PARTICLE_H
#define PARALOBSTAR_PARTICLE_H

#include <cmath>

#include "global.h"

class Particle {
public:
    //TODO: implement overloaded constructor
    double m { -1. }; // mass
    double x[global::dim] {}; // position
    double v[global::dim] {}; // velocity
    double F[global::dim] {}; // force
    double Fn[global::dim] {}; // force last time step
    bool toDelete { false };
    bool moved { false };

    void force(Particle &p);
    void updateX(double dt);
    void updateV(double dt);
};


#endif //PARALOBSTAR_PARTICLE_H
