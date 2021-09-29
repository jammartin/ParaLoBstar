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