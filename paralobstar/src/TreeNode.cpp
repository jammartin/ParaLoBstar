//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/TreeNode.h"

TreeNode::TreeNode(){
    for (int i=0; i<global::powdim; ++i){
        son[i] = nullptr;
    }
}

bool TreeNode::isLeaf(){
    for (int i=0; i<global::powdim; ++i){
        if(son[i] != nullptr) return false;
    }
    return true;
}

bool TreeNode::isCommonCoarseLeaf(){
    if (type == NodeType::commonCoarse){
        for (int i=0; i<global::powdim; ++i){
            if (son[i] != nullptr && son[i]->type == NodeType::commonCoarse) return false;
        }
        return true;
    } else {
        return false;
    }
}