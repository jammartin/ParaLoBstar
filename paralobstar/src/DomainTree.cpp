//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/DomainTree.h"

DomainTree::DomainTree(double domainSize, double theta, double softening,
                       double timeStep) : Tree(domainSize, theta, softening, timeStep){}

void DomainTree::compPseudoParticles(){
    Tree::compPseudoParticles(root);
}

void DomainTree::compForce(){
    compForce(root);
}

void DomainTree::compForce(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compForce(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle){
        // actual force calculation
        for (int d=0; d<global::dim; ++d){
            t.p.F[d] = 0.;
        }
        t.p.U = 0.; // reset particle's energy
        forceBH(t, root, root.box.getLength());
    }
}

void DomainTree::dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet,
                           HighFive::DataSet &vDataSet, HighFive::DataSet &kDataSet) {
    // containers for particle data
    std::vector<double> m {};
    std::vector<std::vector<double>> x {};
    std::vector<std::vector<double>> v {};
    std::vector<keytype> k {};

    getParticleData(m, x, v, k, root);

    mDataSet.write(m);
    xDataSet.write(x);
    vDataSet.write(v);
    kDataSet.write(k);

}

void DomainTree::getParticleData(std::vector<double> &m,
                                 std::vector<std::vector<double>> &x,
                                 std::vector<std::vector<double>> &v,
                                 std::vector<keytype> &k, TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            getParticleData(m, x, v, k, *t.son[i]);
        }
    }
    if (t.isLeaf()){
        m.push_back(t.p.m);
        std::vector<double> xBuffer {};
        std::vector<double> vBuffer {};
        for(int d=0; d<global::dim; ++d){
            xBuffer.push_back(t.p.x[d]);
            vBuffer.push_back(t.p.v[d]);
        }
        x.push_back(xBuffer);
        v.push_back(vBuffer);
        k.push_back(0UL); // dummy key
    }
}