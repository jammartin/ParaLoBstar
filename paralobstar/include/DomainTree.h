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
                   HighFive::DataSet &vDataSet, HighFive::DataSet &kDataSet,
                   HighFive::DataSet &matDataSet, int traceMaterial) override;

private:
    void compForce(TreeNode &t);
    void getParticleData(std::vector<double> &m,
                        std::vector<std::vector<double>> &x,
                        std::vector<std::vector<double>> &v,
                        std::vector<keytype> &k,
                        std::vector<int> &matIds, TreeNode &t);

};


#endif //PARALOBSTAR_DOMAINTREE_H
