//
// Created by Johannes Martin on 21.09.21.
//

#ifndef PARALOBSTAR_PARTICLE_H
#define PARALOBSTAR_PARTICLE_H

#include <cmath>
#include <boost/serialization/access.hpp>

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

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & m;
        ar & x;
        ar & v;
        ar & F;
        ar & Fn;
        ar & toDelete;
        ar & moved;
    }
};


#endif //PARALOBSTAR_PARTICLE_H
