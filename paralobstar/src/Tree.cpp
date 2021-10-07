//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/Tree.h"

Tree::Tree(double domainSize, double _theta, double _timeStep) : theta { _theta }, timeStep { _timeStep }{
    for (int d=0; d<global::dim; ++d){
        root.box.lower[d] = - .5 * domainSize;
        root.box.upper[d] = .5 * domainSize;
    }
}

Tree::~Tree(){
    deallocate(root);
}

void Tree::insertParticle(Particle &p){
    if (root.box.particleWithin(p)){
        if (root.isEmpty()){
            root.p = p;
        } else {
            insertParticle(p, root);
        }
    } else {
        Logger(WARN) << "insertTree(): Particle not in domain. x = ("
                     << p.x[0] << ", " << p.x[1] << ", " << p.x[2] << ")";
    }
}

int Tree::countParticles(){
    int N_ = 0;
    countParticles(root, N_);
    return N_;
}

void Tree::countParticles(TreeNode &t, int &N){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            countParticles(*t.son[i], N);
        }
    }
    if (t.isLeaf() && t.type == NodeType::particle) ++N;
}

void Tree::deallocate(TreeNode &t){
    for(int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            if (t.son[i]->isLeaf()){
                delete t.son[i];
                t.son[i] = nullptr;
            } else {
                deallocate(*t.son[i]);
            }
        }
    }
}
