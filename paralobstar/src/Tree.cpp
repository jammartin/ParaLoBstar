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
        if ((t.isLeaf() || l < theta * distance) && !t.isEmpty()){ // skip empty domain list nodes
            leaf.p.force(t.p);
            leaf.p.U += -.5 * global::G*leaf.p.m*t.p.m/distance; // track gravitational energy on the fly
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
    numParticles = 0;
    countParticles(root, numParticles);
    return numParticles;
}

void Tree::countParticles(TreeNode &t, int &N){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            countParticles(*t.son[i], N);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()) ++N;
}

// compForce() has to be called before
double Tree::totalEnergy(){
    double E_ = 0.;
    compEnergy(E_, root);
    return E_;
}

void Tree::compEnergy(double &E, TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compEnergy(E, *t.son[i]);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        double vSqr = 0.;
        for (int d=0; d<global::dim; ++d){
            vSqr += pow(t.p.v[d], 2.);
        }
        E += t.p.U + .5 * t.p.m * vSqr; // adding potential energy and kinetic energy for each particle
    }
}

void Tree::angularMomentum(std::vector<double> &L){
    compAngularMomentum(L, root);
}

void Tree::compAngularMomentum(std::vector<double> &L, TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compAngularMomentum(L, *t.son[i]);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        // Note: only works for three dimensions
        for (int d=0; d<global::dim; ++d){
            int i = lookup::VectorProduct[d][0];
            int j = lookup::VectorProduct[d][1];
            L[d] += t.p.m * (t.p.x[i] * t.p.v[j] - t.p.x[j] * t.p.v[i]);
        }
    }
}

void Tree::getRanges(std::vector<keytype> &ranges){
    ranges.assign({ 0UL, keyMax }); // dummy ranges
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
