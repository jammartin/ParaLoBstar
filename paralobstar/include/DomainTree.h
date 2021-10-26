//
// Created by Johannes Martin on 21.09.21.
//

#ifndef PARALOBSTAR_DOMAINTREE_H
#define PARALOBSTAR_DOMAINTREE_H

#include "Tree.h"

class DomainTree : public Tree {
public:
    DomainTree(double domainSize, double theta, double softening, double timeStep);

    void compPseudoParticles() override;
    void compForce() override;
    void dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet,
                   HighFive::DataSet &vDataSet, HighFive::DataSet &kDataSet) override;

private:
    //void insertParticle(Particle &p, TreeNode &t) override;
    void compPseudoParticles(TreeNode &t);
    void compForce(TreeNode &t);
    void compPosition(TreeNode &t) override;
    void compVelocity(TreeNode &t) override;
    void moveParticles(TreeNode &t) override;
    void repair(TreeNode &t) override;
    void getParticleData(std::vector<double> &m,
                        std::vector<std::vector<double>> &x,
                        std::vector<std::vector<double>> &v,
                        std::vector<keytype> &k, TreeNode &t);

};


#endif //PARALOBSTAR_DOMAINTREE_H
