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

void Tree::compPosition(){
    compPosition(root);
}

void Tree::compVelocity(){
    compVelocity(root);
}

void Tree::moveParticles(){
    resetFlags(root);
    moveParticles(root);
    repair(root);
}

void Tree::resetFlags(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            resetFlags(*t.son[i]);
        }
    }
    t.p.moved = false;
    t.p.toDelete = false;
}

void Tree::forceBH(TreeNode &leaf, TreeNode &t, double l){
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
    if (t.type == NodeType::particle && t.isLeaf()) ++N;
}

std::vector<keytype> Tree::getRanges(){
    std::vector<keytype> rangesVec_ { 0UL, keyMax };
    return rangesVec_;
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
