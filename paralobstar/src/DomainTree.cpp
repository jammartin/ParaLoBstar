//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/DomainTree.h"

DomainTree::DomainTree(double domainSize, double theta, double timeStep) : Tree(domainSize, theta, timeStep){}

void DomainTree::insertParticle(Particle &p, TreeNode &t){
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
}

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
        t.p.m = 0;
        for (int d=0; d<global::dim; ++d){
            t.p.x[d] = 0;
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
        forceBH(t, root, root.box.getLength());
    }
}

void DomainTree::forceBH(TreeNode &leaf, TreeNode &t, double l){
    if (&leaf != &t){
        double distance = 0.;
        for (int d=0; d<global::dim; ++d){
            distance += pow(t.p.x[d] - leaf.p.x[d], 2.);
        }
        distance = sqrt(distance);
        if (t.isLeaf() || l < theta * distance){
            leaf.p.force(t.p);
        } else {
            for (int i=0; i<global::powdim; ++i){
                if (t.son[i] != nullptr){
                    forceBH(leaf, *t.son[i], .5 * l);
                }
            }
        }
    }
}

void DomainTree::compPosition(){
    compPosition(root);
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

void DomainTree::compVelocity(){
    compVelocity(root);
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

void DomainTree::moveParticles(){
    resetFlags(root);
    moveLeaves(root);
    repair(root);
}

void DomainTree::resetFlags(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            resetFlags(*t.son[i]);
        }
    }
    t.p.moved = false;
    t.p.toDelete = false;
}

void DomainTree::moveLeaves(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            moveLeaves(*t.son[i]);
        }
    }
    if (t.isLeaf() && !t.p.moved){
        t.p.moved = true;
        if (!t.box.particleWithin(t.p)){
            if (root.box.particleWithin(t.p)){
                insertParticle(t.p, root);
                t.p.toDelete = true;
            } else {
                Logger(DEBUG) << "\t\tmoveLeaves(): Particle left system. x = ("
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
        } else if (numberOfSons == 1){
            t.p = t.son[sonIndex]->p;
            if(t.son[sonIndex]->isLeaf()){
                delete t.son[sonIndex];
                t.son[sonIndex] = nullptr;
            }
        }
    }
}

int DomainTree::getParticleData(std::vector<double> &m,
                                std::vector<std::vector<double>> &x,
                                std::vector<std::vector<double>> &v,
                                std::vector<keytype> &k){
    int N_ = 0;
    getParticleData(root, m, x, v, k, N_);
    return N_;
}

void DomainTree::getParticleData(TreeNode &t, std::vector<double> &m,
                                std::vector<std::vector<double>> &x,
                                std::vector<std::vector<double>> &v,
                                 std::vector<keytype> &k, int &N){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            getParticleData(*t.son[i], m, x, v, k, N);
        }
    }
    if (t.isLeaf()){
        ++N;
        m.push_back(t.p.m);
        std::vector<double> xBuffer {};
        std::vector<double> vBuffer {};
        for(int d=0; d<global::dim; ++d){
            xBuffer.push_back(t.p.x[d]);
            vBuffer.push_back(t.p.v[d]);
        }
        x.push_back(xBuffer);
        v.push_back(vBuffer);
        // dummy key
        k.push_back(0UL);
    }
}