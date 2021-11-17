//
// Created by Johannes Martin on 06.10.21.
//

#include "../include/SubDomainTree.h"

SubDomainTree::SubDomainTree(double domainSize, double theta, double softening, double timeStep,
                             bool hilbert) : Tree(domainSize, theta, softening, timeStep){
    numProcs = comm.size();
    myRank = comm.rank();
    range = new keytype[numProcs+1];
    hilbertFlag = hilbert;
    if (hilbertFlag) getKey = &Lebesgue2Hilbert;
}

SubDomainTree::~SubDomainTree(){
    delete [] range;
}

keytype SubDomainTree::Lebesgue2Hilbert(keytype lebesgue, int level){
    keytype hilbert = 0UL;
    int dir = 0;
    for (int lvl=global::maxTreeLvl; lvl>0; --lvl){
        int cell = (lebesgue >> ((lvl-1)*global::dim)) & (keytype)((1<<global::dim)-1);
        hilbert = hilbert << global::dim;
        if (lvl>global::maxTreeLvl-level){
            hilbert += lookup::HilbertTable[dir][cell];
        }
        dir = lookup::DirTable[dir][cell];
    }
    return hilbert;
}

std::string SubDomainTree::key2str(const keytype &key){
    int levels[global::maxTreeLvl];
    for (int lvl=0; lvl<global::maxTreeLvl; ++lvl) {
        levels[lvl] = (key >> global::dim*lvl) & (keytype)7;
    }
    std::string str_ = "#|";
    for (int rLvl=global::maxTreeLvl-1; rLvl>=0; --rLvl) {
        str_ += std::to_string(levels[rLvl]);
        str_ += "|";
    }
    return str_;
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
    Tree::compPseudoParticles(root);

    profiler.time(ProfilerIds::timeCommonCoarse);
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
    for (int ccIndex=0; ccIndex<ccLeaves2insert.size(); ++ccIndex){
        // each common coarse leave can only have children on one process
        for (int proc=0; proc<numProcs; ++proc){
            if (ccLeaves2receive[ccIndex+proc*ccLeavesLength].m > 0.){
                ccLeaves2insert[ccIndex].m = ccLeaves2receive[ccIndex+proc*ccLeavesLength].m;
                for (int d=0; d<global::dim; ++d){
                    ccLeaves2insert[ccIndex].x[d] = ccLeaves2receive[ccIndex+proc*ccLeavesLength].x[d];
                }
            }
        }
    }

    Logger(DEBUG) << "\t... done. Computing common coarse nodes ...";

    std::vector<Particle>::iterator ccLeavesIt = ccLeaves2insert.begin();
    compCommonCoarseNodes(ccLeavesIt, root);

    profiler.time2file(ProfilerIds::timeCommonCoarse);
    Logger(DEBUG) << "\t... done.";
}

void SubDomainTree::fillCommonCoarseLeavesVector(std::vector<Particle> &ccLeaves2send, TreeNode &t){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && t.son[i]->type == NodeType::commonCoarse){
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

    profiler.time(ProfilerIds::timeForceExchange);
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

    delete [] particles2send;

    Particle *particles2receive;
    int totalReceiveLength = particleExchange(particles4procs, particles2receive);
    profiler.value2file<int>(ProfilerIds::forceRcv, totalReceiveLength);

    delete [] particles4procs;

    for (int i=0; i<totalReceiveLength; ++i){
        insertSubTree(particles2receive[i], root);
    }

    delete[] particles2receive;
    profiler.time2file(ProfilerIds::timeForceExchange);

    compForce(root, 0UL, 0);
    repair(root);
}

void SubDomainTree::compForce(TreeNode &t, keytype k, int lvl){
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr){
            compForce(*t.son[i], k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
        }
    }
    if (key2proc(getKey(k, lvl)) == myRank && t.type == NodeType::particle && t.isLeaf()){
        // actual force calculation
        for (int d=0; d<global::dim; ++d){
            t.p.F[d] = 0.;
        }
        t.p.U = 0.; // reset particle's potential energy
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
    int proc = key2proc(getKey(k, lvl));
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

void SubDomainTree::moveParticles(){
    Tree::moveParticles();
    profiler.disableWrite();
    sendParticles();
    profiler.enableWrite();
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
    profiler.value2file<int>(ProfilerIds::lbRcv, totalReceiveLength);

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
    if (t.type == NodeType::particle && t.isLeaf() && (particleProc = key2proc(getKey(k, lvl))) != myRank){
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

    Logger(DEBUG) << "particleExchange(): Receiving " << totalReceiveLength_ << " particles in total ...";

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

    Logger(DEBUG) << "                 ... done.";

    delete [] receiveLengths;
    delete [] sendLengths;

    return totalReceiveLength_;
}


void SubDomainTree::buildCommonCoarseTree() {
    buildCommonCoarseTree(root, 0UL, 0);
}

void SubDomainTree::buildCommonCoarseTree(TreeNode &t, keytype k, int lvl) {
    t.type = NodeType::commonCoarse;
    keytype key = getKey(k, lvl);
    int proc1 = key2proc(key);
    int proc2 = key2proc(key | keyMax >> (global::dim*lvl+1));
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

void SubDomainTree::clearCommonCoarseTree(TreeNode &t){
    t.type = NodeType::pseudoParticle;
    for (int i=0; i<global::powdim; ++i){
        if (t.son[i] != nullptr && t.son[i]->type == NodeType::commonCoarse){
            clearCommonCoarseTree(*t.son[i]);
            if (t.son[i]->isLeaf() && t.son[i]->type == NodeType::pseudoParticle){
                delete t.son[i];
                t.son[i] = nullptr;
            }
        }
    }
}

void SubDomainTree::guessRanges(){
    numParticles = Tree::countParticles();
    Logger(DEBUG) << "\tNumber of particles on process = " << numParticles;
    int pCounter = 0;
    int rangeIndex = 1;
    int maxLvl = 0;
    guessRanges(maxLvl, pCounter, rangeIndex, root, 0UL, 0);

    Logger(DEBUG) << "\tGuessed ranges on process:";
    for (int j=1; j<numProcs; ++j){
        Logger(DEBUG) << "\t\trange[" << j << "] = " << key2str(range[j]);
    }

    keytype sendRange[numProcs+1];
    for (int j=1; j<numProcs; ++j){
        sendRange[j] = range[j]/numProcs;
    }
    mpi::all_reduce(comm, sendRange, numProcs+1, range, std::plus<keytype>());

    int maxLevel = 0;
    mpi::all_reduce(comm, maxLvl, maxLevel, mpi::maximum<keytype>());
    range[0] = 0UL;
    for (int j=1; j<numProcs; ++j){
        // cutting off trailing bits of ranges to avoid the creation of unnecessary common coarse nodes
        range[j] = range[j] & (keyMax << ((global::maxTreeLvl-maxLevel)*global::dim));
    }
    range[numProcs] = keyMax;

    Logger(INFO) << "Averaged guessed ranges:";
    for (int j=0; j<=numProcs; ++j){
        Logger(INFO) << "\trange[" << j << "] = " << key2str(range[j]);
    }

    sendParticles();
}

void SubDomainTree::guessRanges(int &maxLvl, int &pCounter, int &rangeIndex, TreeNode &t, keytype k, int lvl){
    if (hilbertFlag){
        std::map<keytype, int> keyMap;
        for (int i = 0; i < global::powdim; ++i) {
            // sorting implicitly in ascending Hilbert key order utilizing an ordered map
            keytype hilbert = getKey(k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
            keyMap[hilbert] = i;
        }
        // actual recursion in correct order
        for (std::map<keytype, int>::iterator kIt = keyMap.begin(); kIt != keyMap.end(); ++kIt){
            if (t.son[kIt->second] != nullptr && rangeIndex < numProcs){
                guessRanges(maxLvl, pCounter, rangeIndex, *t.son[kIt->second],
                            k | ((keytype)kIt->second << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
            }
        }
    } else {
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] != nullptr && rangeIndex < numProcs){ // range[numProcs] is set manually to keyMax
                guessRanges(maxLvl, pCounter, rangeIndex, *t.son[i],
                            k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
            }
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        ++pCounter;
        if (pCounter > round((double)numParticles/(double)numProcs)*rangeIndex){
            maxLvl = lvl > maxLvl ? lvl : maxLvl;
            range[rangeIndex] = getKey(k, lvl);
            ++rangeIndex;
        }
    }
}

void SubDomainTree::newLoadDistribution(){

    int N = countParticles();
    int particlesOnProc[numProcs];
    mpi::all_gather(comm, numParticles, particlesOnProc);

    // create old and new particle distributions
    int oldDistribution[numProcs+1], newDistribution[numProcs+1];
    oldDistribution[0] = 0;
    for (int proc=0; proc<numProcs; ++proc){
        oldDistribution[proc+1] = oldDistribution[proc] + particlesOnProc[proc];
    }
    for (int j=0; j<=numProcs; ++j){
        newDistribution[j] = (j * N)/numProcs; // new load distribution
        range[j] = 0UL; // reset ranges to zero for all_reduce later on
    }
    int rangeIndex = 0;
    int myDistr = oldDistribution[myRank];
    while (myDistr > newDistribution[rangeIndex]) ++rangeIndex;
    updateRanges(myDistr, rangeIndex, newDistribution, root, 0UL, 0);

    range[0] = 0UL;
    range[numProcs] = keyMax;

    keytype sendRange[numProcs+1];
    std::copy(range, range+numProcs+1, sendRange);

    mpi::all_reduce(comm, sendRange, numProcs+1, range, mpi::maximum<keytype>());

    Logger(INFO) << "New ranges:";
    for (int j=0; j<=numProcs; ++j){
        Logger(INFO) << "\trange[" << j << "] = " << key2str(range[j]);
    }

    sendParticles();
    clearCommonCoarseTree(root);
    buildCommonCoarseTree();

}

void SubDomainTree::updateRanges(int &myDistr, int &rangeIndex, int newDistribution[], TreeNode &t, keytype k, int lvl){
    if (hilbertFlag){
        std::map<keytype, int> keyMap;
        for (int i = 0; i < global::powdim; ++i) {
            // sorting implicitly in ascending Hilbert key order utilizing an ordered map
            keytype hilbert = getKey(k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
            keyMap[hilbert] = i;
        }
        // actual recursion in correct order
        for (std::map<keytype, int>::iterator kIt = keyMap.begin(); kIt != keyMap.end(); ++kIt){
            if (t.son[kIt->second] != nullptr){
                updateRanges(myDistr, rangeIndex, newDistribution,*t.son[kIt->second],
                            k | ((keytype)kIt->second << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
            }
        }
    } else {
        for (int i=0; i<global::powdim; ++i){
            if (t.son[i] != nullptr){
                updateRanges(myDistr, rangeIndex, newDistribution, *t.son[i],
                            k | ((keytype)i << (global::dim*(global::maxTreeLvl-lvl-1))), lvl+1);
            }
        }
    }
    if (t.type == NodeType::particle && t.isLeaf()){
        while (myDistr >= newDistribution[rangeIndex]){
            range[rangeIndex] = getKey(k, lvl);
            ++rangeIndex;
        }
        ++myDistr;
    }
}

int SubDomainTree::key2proc(keytype key){
    for (int j=0; j<numProcs; ++j){
        if (key >= range[j] && key < range[j+1]) return j;
    }
    Logger(ERROR) << "key2proc(): Key " << key2str(key) << " cannot be assigned to any process.";
    return -1;
}

int SubDomainTree::countParticles(){
    Tree::countParticles(); // sets numParticles
    int N_ = 0;
    mpi::all_reduce(comm, numParticles, N_, std::plus<int>());
    return N_;
}

double SubDomainTree::totalEnergy(){
    double E = Tree::totalEnergy();
    double E_tot_ = 0;
    mpi::all_reduce(comm, E, E_tot_, std::plus<double>());
    return E_tot_;
}

void SubDomainTree::angularMomentum(std::vector<double> &L_tot){
    std::vector<double> L(global::dim, 0.);
    Tree::angularMomentum(L);
    mpi::all_reduce(comm, L.data(), global::dim, L_tot.data(), std::plus<double>());
}

void SubDomainTree::dump2file(HighFive::DataSet &mDataSet, HighFive::DataSet &xDataSet,
                              HighFive::DataSet &vDataSet, HighFive::DataSet &kDataSet){
    // containers for particle data
    std::vector<double> m {};
    std::vector<std::vector<double>> x {};
    std::vector<std::vector<double>> v {};
    std::vector<keytype> k {};

    getParticleData(m, x, v, k, root, 0UL, 0);

    std::vector<int> particlesOnProc {};
    mpi::all_gather(comm, numParticles, particlesOnProc);

    std::size_t offset = 0;
    for(int proc=0; proc<myRank; ++proc){
        offset += particlesOnProc[proc];
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
        keys.push_back(getKey(k, lvl));
    }
}

void SubDomainTree::getRanges(std::vector<keytype> &ranges){
    ranges.assign(range, range+numProcs+1);
}