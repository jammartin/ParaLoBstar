//
// Created by Johannes Martin on 06.10.21.
//

#include "../include/SubDomainTree.h"

SubDomainTree::SubDomainTree(double domainSize, double theta, double timeStep) : Tree(domainSize, theta, timeStep){
    numProcs = comm.size();
    myRank = comm.rank();
    range = new keytype[numProcs+1];
}

SubDomainTree::~SubDomainTree(){
    delete [] range;
}

void SubDomainTree::insertParticle(Particle &p, TreeNode &t) {
    Box sonBox{};
    int i = t.box.sonBoxAndIndex(sonBox, p);
    if (t.son[i] == nullptr) {
        if (t.isLeaf()) {
            if (t.type != NodeType::commonCoarse){
                t.type = NodeType::pseudoParticle;
            }
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
        t.son[i]->box = sonBox; // setting boxes for common coarse tree nodes
        insertParticle(p, *t.son[i]);
    }
}

void SubDomainTree::compPseudoParticles(){

}

void SubDomainTree::compForce(){

}

void SubDomainTree::compPosition(){

}

void SubDomainTree::compVelocity(){

}

void SubDomainTree::moveParticles(){

}

int SubDomainTree::getParticleData(std::vector<double> &m,
                    std::vector<std::vector<double>> &x,
                    std::vector<std::vector<double>> &v,
                    std::vector<keytype> &k){
    return 0;
}

void SubDomainTree::guessRanges(){
    numParticles = countParticles();
    Logger(DEBUG) << "\tNumber of particles on process = " << numParticles;
    range[0] = 0UL;
    range[numProcs] = keyMax;
    int pCounter { 0 };
    int rangeIndex { 1 };
    guessRanges(root, pCounter, rangeIndex, 0UL, 0);
    Logger(DEBUG) << "\tGuessed ranges on process:";
    for (int j=0; j<=numProcs; ++j){
        Logger(DEBUG) << "\t\trange[" << j << "] = " << range[j];
    }

    keytype sendRange[numProcs+1];
    for (int j=0; j<=numProcs; ++j){
        sendRange[j] = range[j]/numProcs;
    }

    mpi::all_reduce(comm, sendRange, numProcs+1, range, std::plus<keytype>());

    Logger(DEBUG) << "\tAveraged guessed ranges:";
    for (int j=0; j<=numProcs; ++j){
        Logger(DEBUG) << "\t\trange[" << j << "] = " << range[j];
    }
}

void SubDomainTree::guessRanges(TreeNode &t, int &pCounter, int &rangeIndex, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && rangeIndex < numProcs){
            guessRanges(*t.son[i], pCounter, rangeIndex,
                        k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    if (t.isLeaf() && t.type == NodeType::particle){
        ++pCounter;
        if (pCounter > round((double)numParticles/(double)numProcs)*rangeIndex){
            range[rangeIndex] = k;
            ++rangeIndex;
        }
    }
}