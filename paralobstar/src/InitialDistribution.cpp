//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/InitialDistribution.h"

InitialDistribution::InitialDistribution(const std::string &file){
    // open file
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet mass = h5file.getDataSet("/m"); // TODO: enable different masses for different particles
    HighFive::DataSet pos = h5file.getDataSet("/x");
    HighFive::DataSet vel = h5file.getDataSet("/v");

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);

    // sanity check
    if (x.size() == v.size()){ //TODO: also crosscheck mass
        numberOfParticles = x.size();
    } else {
        throw std::length_error("Length mismatch between position and velocity vectors.");
    }
}

void InitialDistribution::getParticles(Particle *&particles){
    std::vector<std::vector<double>>::iterator xit = x.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    // TODO: also iterate over mass

    int pCounter = 0;

    while (xit != x.end()){
        particles[pCounter].m = m;
        for (int d=0; d<global::dim; ++d){
            particles[pCounter].x[d] = (*xit)[d];
            particles[pCounter].v[d] = (*vit)[d];
        }
        ++xit;
        ++vit;
        ++pCounter;
    }
}