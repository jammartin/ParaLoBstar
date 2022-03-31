//
// Created by Johannes Martin on 21.09.21.
//

#include "../include/Tree.h"

Tree::Tree(double domainSize, double _theta, double _softening,
           double _timeStep) : theta { _theta }, softening { _softening }, timeStep { _timeStep }{
    for (int d=0; d<global::dim; ++d){
        root.box.lower[d] = - .5 * domainSize;
        root.box.upper[d] = .5 * domainSize;
    }
}

Tree::~Tree(){
    deallocate(root);
    Logger(DEBUG) << "Tree destroyed.";
}

void Tree::insertParticle(Particle &p){
    if (root.box.particleWithin(p)){
        if (root.isEmpty()){
            root.p = p;
        } else {
            insertParticle(p, root);
        }
    } else {
        Logger(WARN) << "insertParticle(): Particle not in domain. x = ("
                     << p.x[0] << ", " << p.x[1] << ", " << p.x[2] << ")";
    }
}

void Tree::insertParticle(Particle &p, TreeNode &t) {
    Box sonBox {};
    int i = t.box.sonBoxAndIndex(sonBox, p);

    if (t.son[i] == nullptr) {
        if (t.type == NodeType::particle){
            t.type = NodeType::pseudoParticle;
            Particle pBuffer = t.p;
            t.son[i] = new TreeNode();
            t.son[i]->p = p;
            t.son[i]->box = sonBox;
            insertParticle(pBuffer, t);
            // needed if particle has been flagged for deletion in moveParticles()
            // repair() will take care of deletion if necessary (numberOfSons == 1)
            t.p.toDelete = false;
        } else {
            t.son[i] = new TreeNode();
            t.son[i]->p = p;
            t.son[i]->box = sonBox;
        }
    } else {
        insertParticle(p, *t.son[i]);
    }
}

void Tree::compPseudoParticles(TreeNode &t){
    for (int i=0; i<global::powdim; ++i) {
        if (t.son[i] != nullptr){
            compPseudoParticles(*t.son[i]);
        }
    }
    if (t.type == NodeType::pseudoParticle || t.isCommonCoarseLeaf()){
        t.p.m = 0.;
        for (int d=0; d<global::dim; ++d){
            t.p.x[d] = 0.;
        }
        if (!t.isLeaf()){ // skip calculations for empty common coarse leaf nodes
            for (int i=0; i<global::powdim; ++i){
                if (t.son[i] != nullptr){
                    t.p.m += t.son[i]->p.m;
                    for (int d=0; d<global::dim; ++d){
                        t.p.x[d] += t.son[i]->p.m * t.son[i]->p.x[d];
                    }
                }
            }
            for (int d=0; d<global::dim; ++d){
                t.p.x[d] = t.p.x[d] / t.p.m;
            }
        }
    }
}

void Tree::compPosition(){
    compPosition(root);
}

void Tree::compPosition(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compPosition(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle){
        t.p.updateX(timeStep);
    }
}

void Tree::compVelocity(){
    compVelocity(root);
}

void Tree::compVelocity(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compVelocity(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle){
        t.p.updateV(timeStep);
    }
}

void Tree::moveParticles(){
    resetFlags(root);
    moveParticles(root);
    repair(root);
}

void Tree::moveParticles(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            moveParticles(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle && !t.p.moved){
        t.p.moved = true;
        if (!t.box.particleWithin(t.p)){
            if (root.box.particleWithin(t.p)){
                insertParticle(t.p, root);
                t.p.toDelete = true;
            } else {
                Logger(WARN) << "\t\tmoveParticles(): Particle left system. x = ("
                             << t.p.x[0] << ", " << t.p.x[1] << ", " << t.p.x[2] << ")";
                t.p.toDelete = true;
            }
        }
    }
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

void Tree::repair(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            repair(*t.son[i]);
        }
    }
    if (!t.isLeaf()){
        int numberOfSons = 0;
        int sonIndex;
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] != nullptr){
                if (t.son[i]->p.toDelete){
                    delete t.son[i];
                    t.son[i] = nullptr;
                } else {
                    ++numberOfSons;
                    sonIndex = i;
                }
            }
        }
        if (t.type != NodeType::commonCoarse){
            if (numberOfSons == 0){
                t.p.toDelete = true;
            } else if (numberOfSons == 1 && t.son[sonIndex]->isLeaf()){
                t.p = t.son[sonIndex]->p;
                t.type = t.son[sonIndex]->type;
                delete t.son[sonIndex];
                t.son[sonIndex] = nullptr;
            }
        }
    }
}

void Tree::forceBH(TreeNode &leaf, TreeNode &t, double l){
    if (&leaf != &t){
        double distance = 0.;
        for (int d=0; d<global::dim; ++d){
            distance += pow(t.p.x[d] - leaf.p.x[d], 2.);
        }
        distance = sqrt(distance);
        if ((t.isLeaf() || l < theta * distance) && !t.isEmpty()){ // skip empty domain list nodes
            leaf.p.force(t.p, softening);
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
    if (t.type == NodeType::particle) ++N;
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
    if (t.type == NodeType::particle){
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
    if (t.type == NodeType::particle){
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

void Tree::getCenterOfMass(std::vector<double> &com){
    com.assign(root.p.x, root.p.x+global::dim);
}

void Tree::deallocate(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            deallocate(*t.son[i]);
        }
    }
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && t.son[i]->isLeaf()){
            delete t.son[i];
            t.son[i] = nullptr;
        }
    }
}
