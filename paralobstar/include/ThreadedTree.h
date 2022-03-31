//
// Created by Johannes Martin on 03.02.22.
//

#ifndef PARALOBSTAR_THREADEDTREE_H
#define PARALOBSTAR_THREADEDTREE_H

#include <omp.h>

#include "SubDomainTree.h"

class ThreadedTree : public SubDomainTree {
public:
    ThreadedTree(double domainSize, double theta, double softening, double timeStep, bool hilbert);

private:
    void compForce(TreeNode &t, keytype k, int lvl) override;
};


#endif //PARALOBSTAR_THREADEDTREE_H
