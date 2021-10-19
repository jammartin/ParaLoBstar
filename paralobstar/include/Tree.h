//
// Created by Johannes Martin on 21.09.21.
//

#ifndef PARALOBSTAR_TREE_H
#define PARALOBSTAR_TREE_H

#include <cmath>
#include <cstdint>
#include <vector>

#include "TreeNode.h"
#include "Logger.h"
#include "lookup.h"
#include <highfive/H5File.hpp>

typedef std::uint_fast64_t keytype;

// abstract class
class Tree {
public:
    Tree(double domainSize, double _theta, double _timeStep);
    virtual ~Tree();

    virtual void compPseudoParticles() = 0;
    virtual void compForce() = 0;
    virtual void dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet, HighFive::DataSet &vDataSet,
                           HighFive::DataSet &kDataSet) = 0;

    void insertParticle(Particle &p);
    void compPosition();
    void compVelocity();
    virtual void moveParticles();
    virtual int countParticles();
    virtual double totalEnergy();
    virtual void angularMomentum(std::vector<double> &L_tot);
    virtual void getRanges(std::vector<keytype> &ranges);
    int getNumParticles() const { return numParticles; }

    // placeholders for functions only implemented by SubDomainTree for the parallel mode
    virtual void guessRanges(){}
    virtual void sendParticles(){}
    virtual void buildCommonCoarseTree(){}
    virtual void newLoadDistribution(){}


protected:
    static constexpr keytype keyMax { std::numeric_limits<keytype>::max() };
    TreeNode root {};
    double theta;
    double timeStep;
    int numParticles;

    void forceBH(TreeNode &leaf, TreeNode &t, double l);


private:
    virtual void insertParticle(Particle &p, TreeNode &t) = 0;
    virtual void compPosition(TreeNode &t) = 0;
    virtual void compVelocity(TreeNode &t) = 0;
    virtual void moveParticles(TreeNode &t) = 0;
    virtual void repair(TreeNode &t) = 0;

    void compEnergy(double &E, TreeNode &t);
    void compAngularMomentum(std::vector<double> &L, TreeNode &t);
    void countParticles(TreeNode &t, int &N);

    void resetFlags(TreeNode &t);
    void deallocate(TreeNode &t);
};


#endif //PARALOBSTAR_TREE_H
