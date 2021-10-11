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

typedef std::uint_fast64_t keytype;

// abstract class
class Tree {
public:
    Tree(double domainSize, double _theta, double _timeStep);
    virtual ~Tree();

    virtual void compPseudoParticles() = 0;
    virtual void compForce() = 0;
    virtual void compPosition() = 0;
    virtual void compVelocity() = 0;
    virtual void moveParticles() = 0;
    virtual int getParticleData(std::vector<double> &m,
                                std::vector<std::vector<double>> &x,
                                std::vector<std::vector<double>> &v,
                                std::vector<keytype> &k) = 0;

    virtual void insertParticle(Particle &p) final;
    int countParticles();
    virtual std::vector<keytype> getRanges();

    // placeholders for functions only implemented by SubDomainTree for the parallel mode
    virtual void guessRanges(){}
    virtual void sendParticles(){}
    virtual void buildCommonCoarseTree(){}


protected:
    static constexpr keytype keyMax { std::numeric_limits<keytype>::max() };
    TreeNode root {};
    double theta;
    double timeStep;

    void forceBH(TreeNode &leaf, TreeNode &t, double l);

private:
    virtual void insertParticle(Particle &p, TreeNode &t) = 0;
    void countParticles(TreeNode &t, int &N);
    void deallocate(TreeNode &t);
};


#endif //PARALOBSTAR_TREE_H
