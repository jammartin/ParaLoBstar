//
// Created by Johannes Martin on 06.10.21.
//

#ifndef PARALOBSTAR_SUBDOMAINTREE_H
#define PARALOBSTAR_SUBDOMAINTREE_H

#include <map>
#include <boost/mpi.hpp>

#include "Tree.h"

namespace mpi = boost::mpi;

class SubDomainTree : public Tree {
public:
    SubDomainTree(double domainSize, double theta, double timeStep);
    ~SubDomainTree();

    void compPseudoParticles() override;
    void compForce() override;
    void dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet, HighFive::DataSet &vDataSet,
                   HighFive::DataSet &kDataSet) override;
    void guessRanges() override;
    void sendParticles() override;
    void buildCommonCoarseTree() override;
    std::vector<keytype> getRanges() override;
    int countParticles() override;

private:
    mpi::communicator comm;
    int numProcs;
    int myRank;
    keytype *range;
    int numParticles;

    static constexpr int mpiTag { 17 }; // arbitrary tag

    void insertParticle(Particle &p, TreeNode &t) override;
    void insertSubTree(Particle &p, TreeNode &t);
    void compPseudoParticles(TreeNode &t);
    void fillCommonCoarseLeavesVector(std::vector<Particle> &ccLeaves2send, TreeNode &t);
    void compCommonCoarseNodes(std::vector<Particle>::iterator &ccLeavesIt, TreeNode &t);
    void particles2sendByTheta(std::map<keytype, Particle> *&particles2send, TreeNode &t, keytype k, int lvl);
    void particles2sendByTheta(TreeNode &cc, std::map<keytype, Particle> &particles4proc, TreeNode &t,
                               double l, keytype k);
    void repair(TreeNode &t) override;
    void compForce(TreeNode &t, keytype k, int lvl);
    void compPosition(TreeNode &t) override;
    void compVelocity(TreeNode &t) override;
    void moveParticles(TreeNode &t) override;
    void guessRanges(TreeNode &t, int &pCounter, int &rangeIndex, keytype k, int lvl);
    int key2proc(keytype k);
    void fillSendVectors(std::vector<Particle> *&particles2send, TreeNode &t, keytype k, int lvl);
    int particleExchange(std::vector<Particle> *&particles2send, Particle *&particles2receive);
    void buildCommonCoarseTree(TreeNode &t, keytype k, int lvl);
    void getParticleData(std::vector<double> &m,
                         std::vector<std::vector<double>> &x,
                         std::vector<std::vector<double>> &v,
                         std::vector<keytype> &keys, TreeNode &t, keytype k, int lvl);

};


#endif //PARALOBSTAR_SUBDOMAINTREE_H
