//
// Created by Johannes Martin on 06.10.21.
//

#ifndef PARALOBSTAR_SUBDOMAINTREE_H
#define PARALOBSTAR_SUBDOMAINTREE_H

#include <boost/mpi.hpp>

#include "Tree.h"

namespace mpi = boost::mpi;

class SubDomainTree : public Tree {
public:
    SubDomainTree(double domainSize, double theta, double timeStep);
    ~SubDomainTree();

    void compPseudoParticles() override;
    void compForce() override;
    void compPosition() override;
    void compVelocity() override;
    void moveParticles() override;
    int getParticleData(std::vector<double> &m,
                        std::vector<std::vector<double>> &x,
                        std::vector<std::vector<double>> &v,
                        std::vector<keytype> &k) override;

    void guessRanges() override;

private:
    mpi::communicator comm;
    int numProcs;
    int myRank;
    keytype *range;
    int numParticles;

    void guessRanges(TreeNode &t, int &pCounter, int &rangeIndex, keytype k, int lvl);

    void insertParticle(Particle &p, TreeNode &t) override;

};


#endif //PARALOBSTAR_SUBDOMAINTREE_H
