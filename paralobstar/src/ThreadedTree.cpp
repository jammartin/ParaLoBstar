//
// Created by Johannes Martin on 03.02.22.
//

#include "../include/ThreadedTree.h"

ThreadedTree::ThreadedTree(double domainSize, double theta, double softening, double timeStep,
                           bool hilbert) : SubDomainTree(domainSize, theta, softening, timeStep, hilbert){}

void ThreadedTree::compForce(TreeNode &t, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compForce(*t.son[i], k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    std::vector<TreeNode> forceLeaves {};
    if (key2proc(getKey(k, lvl)) == myRank && t.type == NodeType::particle){
        for (int d=0; d<global::dim; ++d){
            t.p.F[d] = 0.;
        }
        t.p.U = 0.; // reset particle's potential energy
        // store particles to calculate force in parallel
        forceLeaves.push_back(t);
    }

#pragma omp parallel for
    for (std::vector<TreeNode>::iterator pIt = forceLeaves.begin(); pIt != forceLeaves.end(); ++pIt){
        forceBH(*pIt, root, root.box.getLength());
    }
}

