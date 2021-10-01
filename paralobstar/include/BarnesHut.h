//
// Created by Johannes Martin on 23.09.21.
//

#ifndef PARALOBSTAR_BARNESHUT_H
#define PARALOBSTAR_BARNESHUT_H

#include <highfive/H5File.hpp>

#include "Logger.h"
#include "ConfigParser.h"
#include "H5Profiler.h"
#include "InitialDistribution.h"
#include "Particle.h"
#include "DomainTree.h"

class BarnesHut {
public:
    BarnesHut(ConfigParser confP, int _myRank, int _numProcs);
    ~BarnesHut();

    void run();

private:
    // values read in from config file
    double domainSize {};
    std::string initFile {};
    bool parallel {};
    double timeStep {};
    double timeEnd {};
    double theta {};
    int h5DumpInterval {};
    int loadBalancingInterval {};

    // mpi related variables
    int myRank {};
    int numProcs {};

    // storing reference to singleton profiler instance
    H5Profiler &profiler = H5Profiler::getInstance();

    // internal variables initialized in constructor
    int N; // number of particles
    Particle *particles;
    Tree *tree;

};


#endif //PARALOBSTAR_BARNESHUT_H
