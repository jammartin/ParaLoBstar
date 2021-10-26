//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/DomainTree.h"

DomainTree::DomainTree(double domainSize, double theta, double softening,
                       double timeStep) : Tree(domainSize, theta, softening, timeStep){}

/*void DomainTree::insertParticle(Particle &p, TreeNode &t){
    Box sonBox {};
    int i = t.box.sonBoxAndIndex(sonBox, p);
    if (t.son[i] == nullptr){
        if (t.isLeaf()){ // t.p is a particle
            Particle pBuffer = t.p;
            t.son[i] = new TreeNode();
            t.son[i]->p = p;
            t.son[i]->box = sonBox;
            insertParticle(pBuffer, t);
        } else {
            t.son[i] = new TreeNode();
            t.son[i]->p = p;
            t.son[i]->box = sonBox;
        }
    } else {
        insertParticle(p, *t.son[i]);
    }
}*/

void DomainTree::compPseudoParticles(){
    compPseudoParticles(root);
}

void DomainTree::compPseudoParticles(TreeNode &t){
    for (int i=0; i<global::powdim; ++i) {
        if (t.son[i] != nullptr){
            compPseudoParticles(*t.son[i]);
        }
    }
    if (!t.isLeaf()){ // t.p is a pseudo-particle
        t.p.m = 0.;
        for (int d=0; d<global::dim; ++d){
            t.p.x[d] = 0.;
        }
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] != nullptr){
                t.p.m += t.son[i]->p.m;
                for (int d=0; d<global::dim; ++d){
                    t.p.x[d] += t.son[i]->p.m * t.son[i]->p.x[d];
                }
            }
        }
        for (int d=0; d<global::dim; ++d){
            t.p.x[d] = t.p.x[d] / t.p.m;
        }
    }
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
    if (t.isLeaf()){
        // actual force calculation
        for (int d=0; d<global::dim; ++d){
            t.p.F[d] = 0.;
        }
        t.p.U = 0.; // reset particle's energy
        forceBH(t, root, root.box.getLength());
    }
}

void DomainTree::compPosition(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compPosition(*t.son[i]);
        }
    }
    if (t.isLeaf()){
        t.p.updateX(timeStep);
    }
}

void DomainTree::compVelocity(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compVelocity(*t.son[i]);
        }
    }
    if (t.isLeaf()){
        t.p.updateV(timeStep);
    }
}

void DomainTree::moveParticles(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            moveParticles(*t.son[i]);
        }
    }
    if (t.isLeaf() && !t.p.moved){
        t.p.moved = true;
        if (!t.box.particleWithin(t.p)){
            if (root.box.particleWithin(t.p)){
                insertParticle(t.p, root);
                t.p.toDelete = true;
            } else {
                Logger(WARN) << "\t\tmoveParticles(): Particle left system. x = ("
                             << t.p.x[0] << ", " << t.p.x[1] << ", " << t.p.x[2] << ")";
                t.p.toDelete = true;
            }
        }
    }
}

void DomainTree::repair(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            repair(*t.son[i]);
        }
    }
    if (!t.isLeaf()){
        int numberOfSons = 0;
        int sonIndex;
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] != nullptr){
                if (t.son[i]->p.toDelete){
                    delete t.son[i];
                    t.son[i] = nullptr;
                } else {
                    ++numberOfSons;
                    sonIndex = i;
                }
            }
        }
        if (numberOfSons == 0){
            t.p.toDelete = true;
        } else if (numberOfSons == 1 && t.son[sonIndex]->isLeaf()){
            t.p = t.son[sonIndex]->p;
            t.type = t.son[sonIndex]->type;
            delete t.son[sonIndex];
            t.son[sonIndex] = nullptr;

        }
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