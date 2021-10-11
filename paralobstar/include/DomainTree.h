//
// Created by Johannes Martin on 21.09.21.
//

#ifndef PARALOBSTAR_DOMAINTREE_H
#define PARALOBSTAR_DOMAINTREE_H

#include "Tree.h"

class DomainTree : public Tree {
public:
    DomainTree(double domainSize, double theta, double timeStep);

    void compPseudoParticles() override;
    void compForce() override;
    void compPosition() override;
    void compVelocity() override;
    void moveParticles() override;
    int getParticleData(std::vector<double> &m,
                        std::vector<std::vector<double>> &x,
                        std::vector<std::vector<double>> &v,
                        std::vector<keytype> &k) override;

private:
    void insertParticle(Particle &p, TreeNode &t) override;
    void compPseudoParticles(TreeNode &t);
    void compForce(TreeNode &t);
    void compPosition(TreeNode &t);
    void compVelocity(TreeNode &t);
    void resetFlags(TreeNode &t);
    void moveLeaves(TreeNode &t);
    void repair(TreeNode &t);
    void getParticleData(TreeNode &t, std::vector<double> &m,
                        std::vector<std::vector<double>> &x,
                        std::vector<std::vector<double>> &v,
                        std::vector<keytype> &k, int &N);

};


#endif //PARALOBSTAR_DOMAINTREE_H
