//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/Tree.h"

Tree::Tree(double _theta, double _timeStep) : theta { _theta }, timeStep { _timeStep }{}

Tree::~Tree(){
    deallocate(root);
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