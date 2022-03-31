//
// Created by Johannes Martin on 23.09.21.
//

#ifndef PARALOBSTAR_INITIALDISTRIBUTION_H
#define PARALOBSTAR_INITIALDISTRIBUTION_H

#include <highfive/H5File.hpp>

#include "Particle.h"

class InitialDistribution {
public:
    InitialDistribution(const std::string &file, bool _traceMaterial=false);

    int getNumberOfParticles() const { return numberOfParticles; };
    void getAllParticles(Particle *&particles);
    void getParticles(Particle *&particles, int offset, int amount);

private:
    // containers to be filled from hdf5 file
    std::vector<double> m {};
    std::vector<std::vector<double>> x {}, v {};
    std::vector<int> materialId {};
    int numberOfParticles { 0 };
    bool traceMaterial;
};


#endif //PARALOBSTAR_INITIALDISTRIBUTION_H
