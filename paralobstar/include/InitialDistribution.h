//
// Created by Johannes Martin on 23.09.21.
//

#ifndef PARALOBSTAR_INITIALDISTRIBUTION_H
#define PARALOBSTAR_INITIALDISTRIBUTION_H

#include <highfive/H5File.hpp>

#include "Particle.h"

class InitialDistribution {
public:
    InitialDistribution(const std::string &file);

    int getNumberOfParticles() const { return numberOfParticles; };
    void getParticles(Particle *&particles);

private:
    // containers to be filled from hdf5 file
    std::vector<double> m {};
    std::vector<std::vector<double>> x {}, v {};
    int numberOfParticles { 0 };
};


#endif //PARALOBSTAR_INITIALDISTRIBUTION_H
