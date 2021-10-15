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
    virtual std::vector<keytype> getRanges();

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

    void countParticles(TreeNode &t, int &N);
    void forceBH(TreeNode &leaf, TreeNode &t, double l);


private:
    virtual void insertParticle(Particle &p, TreeNode &t) = 0;
    virtual void compPosition(TreeNode &t) = 0;
    virtual void compVelocity(TreeNode &t) = 0;
    virtual void moveParticles(TreeNode &t) = 0;
    virtual void repair(TreeNode &t) = 0;

    void resetFlags(TreeNode &t);
    void deallocate(TreeNode &t);
};


#endif //PARALOBSTAR_TREE_H
