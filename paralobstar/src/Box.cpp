//
// Created by Johannes Martin on 20.09.21.
//

#include "../include/Box.h"

int Box::sonBoxAndIndex(Box &sonBox, Particle &p) {
    int sonIndex_ = 0;
    for (int d=global::dim-1; d>=0; --d){
        if (p.x[d] < .5 * (upper[d] + lower[d])){
            sonIndex_ = 2 * sonIndex_;
            sonBox.lower[d] = lower[d];
            sonBox.upper[d] = .5 * (lower[d] + upper[d]);
        } else {
            sonIndex_ = 2 * sonIndex_ + 1;
            sonBox.lower[d] = .5 * (lower[d] + upper[d]);
            sonBox.upper[d] = upper[d];
        }
    }
    return sonIndex_;
}

void Box::sonBoxByIndex(Box &sonBox, int sonIndex){
    for (int d=0; d<global::dim; ++d){
        if (sonIndex % 2){
            sonBox.lower[d] = .5 * (lower[d] + upper[d]);
            sonBox.upper[d] = upper[d];
        } else {
            sonBox.lower[d] = lower[d];
            sonBox.upper[d] = .5 * (lower[d] + upper[d]);
        }
        sonIndex /= 2;

    }
}

double Box::getLength(){
    double length_ = 0;
    for (int d=0; d<global::dim; ++d){
        if (upper[d] - lower[d] > length_){
            length_ = upper[d] - lower[d];
        }
    }
    return length_;
}

bool Box::particleWithin(Particle &p){
    for (int d=0; d<global::dim; ++d){
        if (p.x[d] < lower[d] || p.x[d] >= upper[d]){
            return false;
        }
    }
    return true;
}

double Box::smallestDistance(Particle &p){
    double distance_ = 0.;
    for (int d=0; d<global::dim; ++d){
        if (p.x[d] < lower[d]){
            distance_ += pow(lower[d] - p.x[d], 2.);
        } else if (p.x[d] > upper[d]) {
            distance_ += pow(p.x[d] - upper[d], 2.);
        } else {
            distance_ += 0.;
        }
    }
    return sqrt(distance_);
}