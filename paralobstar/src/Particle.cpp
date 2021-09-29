//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/Particle.h"

void Particle::force(Particle &p){
    double dsqr = 0.;
    for (int d=0; d<global::dim; ++d){
        dsqr += pow(p.x[d] - x[d], 2.);
    }
    double f = global::G * m * p.m / (sqrt(dsqr) * dsqr);
    for (int d=0; d<global::dim; ++d){
        F[d] += f * (p.x[d] - x[d]);
    }
}

void Particle::updateX(double dt){
    double a = dt * .5 / m;
    for (int d=0; d<global::dim; ++d){
        x[d] += dt * (v[d] + a * F[d]);
        Fn[d] = F[d];
    }
}

void Particle::updateV(double dt){
    double a = dt * .5 / m;
    for (int d=0; d<global::dim; ++d){
        v[d] += a * (F[d] + Fn[d]);
    }
}