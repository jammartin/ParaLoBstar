//
// Created by Johannes Martin on 06.10.21.
//

#ifndef PARALOBSTAR_SUBDOMAINTREE_H
#define PARALOBSTAR_SUBDOMAINTREE_H

#include <map>
#include <boost/mpi.hpp>

#include "H5Profiler.h"
#include "Tree.h"

namespace mpi = boost::mpi;

class SubDomainTree : public Tree {
public:
    SubDomainTree(double domainSize, double theta, double softening, double timeStep, bool hilbert);
    ~SubDomainTree();

    void compPseudoParticles() override;
    void compForce() override;
    void dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet, HighFive::DataSet &vDataSet,
                   HighFive::DataSet &kDataSet) override;
    void guessRanges() override;
    void newLoadDistribution() override;
    void getRanges(std::vector<keytype> &ranges) override;
    void moveParticles() override;
    int countParticles() override;
    double totalEnergy() override;
    void angularMomentum(std::vector<double> &L) override;

private:
    mpi::communicator comm;
    int numProcs;
    int myRank;
    keytype *range;

    static constexpr int mpiTag { 17 }; // arbitrary tag

    // storing reference to singleton profiler instance
    H5Profiler &profiler = H5Profiler::getInstance();

    bool hilbertFlag { false };
    static keytype Lebesgue(keytype k_, int _){ return k_; }
    static keytype Lebesgue2Hilbert(keytype lebesgue, int level);
    keytype (*getKey)(keytype, int){ &Lebesgue };
    static std::string key2str(const keytype &key);

    void insertSubTree(Particle &p, TreeNode &t);
    void fillCommonCoarseLeavesVector(std::vector<Particle> &ccLeaves2send, TreeNode &t);
    void compCommonCoarseNodes(std::vector<Particle>::iterator &ccLeavesIt, TreeNode &t);
    void particles2sendByTheta(std::map<keytype, Particle> *&particles2send, TreeNode &t, keytype k, int lvl);
    void particles2sendByTheta(TreeNode &cc, std::map<keytype, Particle> &particles4proc, TreeNode &t,
                               double l, keytype k);
    void compForce(TreeNode &t, keytype k, int lvl);
    void guessRanges(int &maxLvl, int &pCounter, int &rangeIndex, TreeNode &t, keytype k, int lvl);
    void sendParticles();
    void buildCommonCoarseTree();
    void updateRanges(int &myDistr, int &rangeIndex, int newDistribution[], TreeNode &t, keytype k, int lvl);
    int key2proc(keytype k);
    void fillSendVectors(std::vector<Particle> *&particles2send, TreeNode &t, keytype k, int lvl);
    int particleExchange(std::vector<Particle> *&particles2send, Particle *&particles2receive);
    void buildCommonCoarseTree(TreeNode &t, keytype k, int lvl);
    void clearCommonCoarseTree(TreeNode &t);
    void getParticleData(std::vector<double> &m,
                         std::vector<std::vector<double>> &x,
                         std::vector<std::vector<double>> &v,
                         std::vector<keytype> &keys, TreeNode &t, keytype k, int lvl);

};


#endif //PARALOBSTAR_SUBDOMAINTREE_H
