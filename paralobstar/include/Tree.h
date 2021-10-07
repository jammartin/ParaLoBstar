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
    static constexpr keytype keyMax { std::numeric_limits<keytype>::max() };

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

    virtual void guessRanges(){}

protected:
    TreeNode root {};
    double theta;
    double timeStep;

private:
    virtual void insertParticle(Particle &p, TreeNode &t) = 0;
    void countParticles(TreeNode &t, int &N);
    void deallocate(TreeNode &t);
};


#endif //PARALOBSTAR_TREE_H
