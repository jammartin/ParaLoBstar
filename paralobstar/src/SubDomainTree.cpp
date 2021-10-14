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
    Box sonBox {};
    int i = t.box.sonBoxAndIndex(sonBox, p);

    if (t.son[i] == nullptr) {
        if (t.type == NodeType::particle && t.isLeaf()) {
            t.type = NodeType::pseudoParticle;
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
        insertParticle(p, *t.son[i]);
    }
}

void SubDomainTree::insertSubTree(Particle &p, TreeNode &t) {
    Box sonBox {};
    int i = t.box.sonBoxAndIndex(sonBox, p);

    if (t.son[i] == nullptr) {
        t.son[i] = new TreeNode();
        t.son[i]->p = p;
        t.son[i]->p.toDelete = true;
        t.son[i]->box = sonBox;
    } else {
        insertSubTree(p, *t.son[i]);
    }
}

void SubDomainTree::compPseudoParticles(){
    compPseudoParticles(root);

    std::vector<Particle> ccLeaves2send {};
    std::vector<Particle> ccLeaves2receive {};

    fillCommonCoarseLeavesVector(ccLeaves2send, root);
    const size_t ccLeavesLength { ccLeaves2send.size() };

    Logger(DEBUG) << "\tUpdating " << ccLeavesLength << " common coarse leaves ...";

    mpi::all_gather(comm, ccLeaves2send.data(), ccLeavesLength, ccLeaves2receive);

    std::vector<Particle> ccLeaves2insert { ccLeavesLength };
    // setting relevant particle values to 0
    for (std::vector<Particle>::iterator ccIt=ccLeaves2insert.begin(); ccIt != ccLeaves2insert.end(); ++ccIt){
        ccIt->m = 0.;
        for (int d=0; d<global::dim; ++d){
            ccIt->x[d] = 0.;
        }
    }
    // concatenating common coarse leaves from all processes
    for (int proc=0; proc<numProcs; ++proc){
        for (int ccIndex=0; ccIndex<ccLeaves2insert.size(); ++ccIndex){
            ccLeaves2insert[ccIndex].m += ccLeaves2receive[ccIndex+proc*ccLeavesLength].m;
            for (int d=0; d<global::dim; ++d){
                ccLeaves2insert[ccIndex].x[d] += ccLeaves2receive[ccIndex+proc*ccLeavesLength].x[d];
            }
        }
    }

    Logger(DEBUG) << "\t... done. Computing common coarse nodes ...";

    std::vector<Particle>::iterator ccLeavesIt = ccLeaves2insert.begin();
    compCommonCoarseNodes(ccLeavesIt, root);

    Logger(DEBUG) << "\t... done.";


}

void SubDomainTree::compPseudoParticles(TreeNode &t){
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

void SubDomainTree::fillCommonCoarseLeavesVector(std::vector<Particle> &ccLeaves2send, TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            fillCommonCoarseLeavesVector(ccLeaves2send, *t.son[i]);
        }
    }
    if (t.isCommonCoarseLeaf()){
        ccLeaves2send.push_back(t.p);
    }
}

void SubDomainTree::compCommonCoarseNodes(std::vector<Particle>::iterator &ccLeavesIt, TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && t.son[i]->type == NodeType::commonCoarse){
            compCommonCoarseNodes(ccLeavesIt, *t.son[i]);
        }
    }
    if (t.isCommonCoarseLeaf()){
        t.p = *ccLeavesIt;
        ++ccLeavesIt;
    } else { // update inner common coarse node like pseudo particles
        t.p.m = 0.;
        for (int d=0; d<global::dim; ++d){
            t.p.x[d] = 0.;
        }
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] != nullptr){
                t.p.m += t.son[i]->p.m;
                for (int d=0; d<global::dim; ++d){
                    t.p.x[d] += t.son[i]->p.m * t.son[i]->p.x[d];
                }
            }
        }
        if (t.p.m > 0.){
            for (int d=0; d<global::dim; ++d){
                t.p.x[d] = t.p.x[d] / t.p.m;
            }
        }
    }
}

void SubDomainTree::compForce(){
    auto particles2send = new std::map<keytype, Particle>[numProcs];

    Logger(DEBUG) << "\tCollecting particles to send to processes ...";
    particles2sendByTheta(particles2send, root, 0UL, 0);

    // copying map entries into vector
    auto particles4procs = new std::vector<Particle>[numProcs];

    for (int proc=0; proc<numProcs; ++proc){
        Logger(DEBUG) << "\t... sending " << particles2send[proc].size() << " particles to process " << proc;
        for (std::map<keytype, Particle>::iterator pIt = particles2send[proc].begin();
                pIt != particles2send[proc].end(); ++pIt){
            particles4procs[proc].push_back(pIt->second);
        }
    }

    Particle *particles2receive;
    int totalReceiveLength = particleExchange(particles4procs, particles2receive);

    delete [] particles4procs;

    for (int i=0; i<totalReceiveLength; ++i){
        insertSubTree(particles2receive[i], root);
    }

    compForce(root, 0UL, 0);
    repair(root);

    delete[] particles2receive;
}

void SubDomainTree::compForce(TreeNode &t, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compForce(*t.son[i], k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    if (key2proc(k) == myRank && t.type != NodeType::commonCoarse && t.isLeaf()){
        // actual force calculation
        for (int d=0; d<global::dim; ++d){
            t.p.F[d] = 0.;
        }
        forceBH(t, root, root.box.getLength());
    }
}

void SubDomainTree::particles2sendByTheta(std::map<keytype, Particle> *&particles2send, TreeNode &t,
                                          keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && t.son[i]->type==NodeType::commonCoarse){
            particles2sendByTheta(particles2send, *t.son[i],
                                  k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    int proc = key2proc(k);
    if (proc != myRank && lvl > 0){ // exclude root
        particles2sendByTheta(t, particles2send[proc], root, root.box.getLength(), 1UL);
    }
}

void SubDomainTree::particles2sendByTheta(TreeNode &cc, std::map<keytype, Particle> &particles4proc, TreeNode &t,
                                          double l, keytype k) {
    if (t.type != NodeType::commonCoarse) {
        particles4proc[k] = t.p; // map is filled with unique keys and therefore overwrites duplicates
    }
    double distance = cc.box.smallestDistance(t.p);
    if (l >= theta * distance) {
        for (int i = 0; i < global::powdim; ++i) {
            if (t.son[i] != nullptr){
                // the key chosen below ensures uniques and ordering from root to leaves
                particles2sendByTheta(cc, particles4proc, *t.son[i],
                                      .5 * l, (k << global::dim) | (keytype)i);
            }
        }
    }
}


void SubDomainTree::compPosition(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compPosition(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        t.p.updateX(timeStep);
    }
}

void SubDomainTree::compVelocity(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compVelocity(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        t.p.updateV(timeStep);
    }
}

void SubDomainTree::moveParticles(TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            moveParticles(*t.son[i]);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf() && !t.p.moved){
        t.p.moved = true;
        if (!t.box.particleWithin(t.p)){
            if (root.box.particleWithin(t.p)){
                insertParticle(t.p, root);
                t.p.toDelete = true;
            } else {
                Logger(DEBUG) << "\t\tmoveLeaves(): Particle left system. x = ("
                              << t.p.x[0] << ", " << t.p.x[1] << ", " << t.p.x[2] << ")";
                t.p.toDelete = true;
            }
        }
    }
}

void SubDomainTree::repair(TreeNode &t){
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

void SubDomainTree::guessRanges(){
    numParticles = 0;
    Tree::countParticles(root, numParticles);
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

    Logger(INFO) << "Averaged guessed ranges:";
    for (int j=0; j<=numProcs; ++j){
        Logger(INFO) << "\trange[" << j << "] = " << range[j];
    }
}

void SubDomainTree::guessRanges(TreeNode &t, int &pCounter, int &rangeIndex, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && rangeIndex < numProcs){ // range[numProcs] is set manually to keyMax
            guessRanges(*t.son[i], pCounter, rangeIndex,
                        k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        ++pCounter;
        if (pCounter > round((double)numParticles/(double)numProcs)*rangeIndex){
            range[rangeIndex] = k;
            ++rangeIndex;
        }
    }
}

void SubDomainTree::sendParticles(){
    auto particles2send = new std::vector<Particle>[numProcs];

    fillSendVectors(particles2send, root, 0UL, 0);

    for (int proc=0; proc < numProcs; ++proc){
        Logger(DEBUG) << "\tSending " << particles2send[proc].size() << " particles to process " << proc;
    }

    repair(root);

    Particle *particles2receive;
    int totalReceiveLength = particleExchange(particles2send, particles2receive);

    delete [] particles2send;

    for (int i=0; i<totalReceiveLength; ++i){
        insertParticle(particles2receive[i], root);
    }

    delete [] particles2receive;
}

void SubDomainTree::fillSendVectors(std::vector<Particle> *&particles2send, TreeNode &t, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            fillSendVectors(particles2send, *t.son[i], k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    int particleProc;
    if (t.type == NodeType::particle && t.isLeaf() && (particleProc = key2proc(k)) != myRank){
        particles2send[particleProc].push_back(t.p);
        t.p.toDelete = true;
    }
}

int SubDomainTree::particleExchange(std::vector<Particle> *&particles2send, Particle *&particles2receive){

    auto receiveLengths = new int[numProcs];
    auto sendLengths = new int[numProcs];

    for (int proc=0; proc<numProcs; ++proc){
        if (proc != myRank){
            sendLengths[proc] = particles2send[proc].size();
        } else {
            sendLengths[proc] = -1;
            receiveLengths[proc] = -1;
        }
    }

    std::vector<mpi::request> reqLengths;
    std::vector<mpi::status> statLengths;

    for (int proc=0; proc<numProcs; ++proc){
        if (proc != myRank){
            reqLengths.push_back(comm.isend(proc, mpiTag, &sendLengths[proc], 1));
            statLengths.push_back(comm.recv(proc, mpiTag, &receiveLengths[proc], 1));
        }
    }

    mpi::wait_all(reqLengths.begin(), reqLengths.end());

    int totalReceiveLength_ = 0;
    for (int proc=0; proc<numProcs; ++proc){
        if (proc != myRank){
            totalReceiveLength_ += receiveLengths[proc];
        }
    }

    Logger(INFO) << "particleExchange(): Receiving " << totalReceiveLength_ << " particles in total ...";

    particles2receive = new Particle[totalReceiveLength_];

    std::vector<mpi::request> reqParticles;
    std::vector<mpi::status> statParticles;

    int offset = 0;
    for (int proc=0; proc<numProcs; ++proc){
        if (proc != myRank){
            reqParticles.push_back(comm.isend(proc, mpiTag, particles2send[proc].data(), sendLengths[proc]));
            statParticles.push_back(comm.recv(proc, mpiTag, particles2receive + offset, receiveLengths[proc]));
            offset += receiveLengths[proc];
        }
    }

    mpi::wait_all(reqParticles.begin(), reqParticles.end());

    Logger(INFO) << "                 ... done.";

    delete [] receiveLengths;
    delete [] sendLengths;

    return totalReceiveLength_;
}


void SubDomainTree::buildCommonCoarseTree() {
    buildCommonCoarseTree(root, 0UL, 0);
}

void SubDomainTree::buildCommonCoarseTree(TreeNode &t, keytype k, int lvl) {
    t.type = NodeType::commonCoarse;
    int proc1 = key2proc(k);
    int proc2 = key2proc(k | keyMax >> (global::dim*lvl+1));
    if (proc1 != proc2){
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] == nullptr){
                t.son[i] = new TreeNode();
                t.box.sonBoxByIndex(t.son[i]->box, i);
            } else if (t.son[i]->type == NodeType::particle && t.son[i]->isLeaf()){
                t.son[i]->type = NodeType::commonCoarse; // has to be set before calling insertParticle
                insertParticle(t.son[i]->p, *t.son[i]); // move particle one level deeper
                t.son[i]->p = Particle(); // reset common coarse node's particle to default
            }
            buildCommonCoarseTree(*t.son[i], k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
}

int SubDomainTree::key2proc(keytype k){
    for (int j=0; j<numProcs; ++j){
        if (k >= range[j] && k < range[j+1]) return j;
    }
    Logger(ERROR) << "key2proc(): Key " << k << " cannot be assigned to any process.";
    return -1;
}

int SubDomainTree::countParticles(){
    numParticles = 0;
    Tree::countParticles(root, numParticles);
    int N_ = 0;
    mpi::all_reduce(comm, numParticles, N_, std::plus<keytype>());
    return N_;
}

void SubDomainTree::dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet,
                              HighFive::DataSet &vDataSet, HighFive::DataSet &kDataSet){
    // containers for particle data
    std::vector<double> m {};
    std::vector<std::vector<double>> x {};
    std::vector<std::vector<double>> v {};
    std::vector<keytype> k {};

    getParticleData(m, x, v, k, root, 0UL, 0);

    std::vector<int> procsNumParticles {};
    mpi::all_gather(comm, numParticles, procsNumParticles);

    std::size_t offset = 0;
    for(int proc=0; proc<myRank; ++proc){
        offset += procsNumParticles[proc];
    }
    Logger(DEBUG) << "\tWriting " << numParticles << " particles @" << offset;

    // assuming countParticles() has been called immediately before
    mDataSet.select({offset}, {std::size_t(numParticles)}).write(m);
    xDataSet.select({offset , 0},
                    {std::size_t(numParticles), std::size_t(global::dim)}).write(x);
    vDataSet.select({offset , 0},
                    {std::size_t(numParticles), std::size_t(global::dim)}).write(v);
    kDataSet.select({offset}, {std::size_t(numParticles)}).write(k);
}

void SubDomainTree::getParticleData(std::vector<double> &m,
                                    std::vector<std::vector<double>> &x,
                                    std::vector<std::vector<double>> &v,
                                    std::vector<keytype> &keys, TreeNode &t, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            getParticleData(m, x, v, keys,
                            *t.son[i], k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        m.push_back(t.p.m);
        std::vector<double> xBuffer {};
        std::vector<double> vBuffer {};
        for(int d=0; d<global::dim; ++d){
            xBuffer.push_back(t.p.x[d]);
            vBuffer.push_back(t.p.v[d]);
        }
        x.push_back(xBuffer);
        v.push_back(vBuffer);
        keys.push_back(k);
    }
}

std::vector<keytype> SubDomainTree::getRanges(){
    std::vector<keytype> rangesVec_ {};
    rangesVec_.assign(range, range+numProcs+1);
    return rangesVec_;
}