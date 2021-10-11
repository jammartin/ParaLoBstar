//
// Created by Johannes Martin on 21.09.21.
//

#ifndef PARALOBSTAR_TREENODE_H
#define PARALOBSTAR_TREENODE_H

#include "global.h"
#include "Particle.h"
#include "Box.h"

enum class NodeType {
    particle, pseudoParticle, commonCoarse
};

class TreeNode {
public:
    Particle p {};
    Box box {};
    TreeNode *son[global::powdim];
    NodeType type { NodeType::particle };

    TreeNode();

    bool isEmpty() const { return p.m < 0.; }
    bool isLeaf();
    bool isCommonCoarseLeaf();

};


#endif //PARALOBSTAR_TREENODE_H
