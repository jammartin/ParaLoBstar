//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/BarnesHut.h"

BarnesHut::BarnesHut(ConfigParser confP) : domainSize { confP.getVal<double>("domainSize") },
                                           initFile { confP.getVal<std::string>("initFile") },
                                           outDir { confP.getVal<std::string>("outDir") },
                                           parallel { confP.getVal<bool>("parallel") },
                                           timeStep { confP.getVal<double>("timeStep") },
                                           timeEnd { confP.getVal<double>("timeEnd") },
                                           h5DumpInterval { confP.getVal<int>("h5DumpInterval") },
                                           loadBalancingInterval { confP.getVal<int>("loadBalancingInterval") }
{
    numProcs = comm.size();
    myRank = comm.rank();

    Logger(INFO) << "Reading in initial distribution ...";
    InitialDistribution initDist { initFile };

    N = initDist.getNumberOfParticles();
    Logger(INFO) << "... done. Total number of particles N = " << N;

    steps = (int)round(timeEnd/timeStep);
    Logger(INFO) << "Number of time steps = " << steps;

    if (steps % loadBalancingInterval != 0){
        Logger(ERROR) << "Number of time steps (" << steps
                      << ") must have load balancing interval (" << loadBalancingInterval << ") as divisor. Aborting.";
        throw std::invalid_argument("Number of time steps must have load balancing interval as divisor.");
    }

    // Initialize profiler data sets for both serial and parallel mode
    profiler.createValueDataSet<int>(ProfilerIds::N, steps+1);
    profiler.createTimeDataSet(ProfilerIds::timePos, steps);
    profiler.createTimeDataSet(ProfilerIds::timeMove, steps);
    profiler.createTimeDataSet(ProfilerIds::timePseudo, steps);
    profiler.createTimeDataSet(ProfilerIds::timeForce, steps);
    profiler.createTimeDataSet(ProfilerIds::timeVel, steps);
    // Initialize only parallel mode data sets
    if (parallel){
        profiler.createTimeDataSet(ProfilerIds::timeCommonCoarse, steps);
        profiler.createTimeDataSet(ProfilerIds::timeForceExchange, steps);
        profiler.createValueDataSet<int>(ProfilerIds::forceRcv, steps);
        profiler.createTimeDataSet(ProfilerIds::timeLb, steps/loadBalancingInterval);
        profiler.createValueDataSet<int>(ProfilerIds::lbRcv, steps/loadBalancingInterval);
    }

    int particlesPerProc = N/numProcs;

    if (N % numProcs != 0){
        Logger(ERROR) << "Number of particles (" << N
                      << ") must have number of processes (" << numProcs << ") as divisor. Aborting.";
        throw std::domain_error("Total amount of particles must have number of processes as divisor.");
    }

    particles = new Particle[particlesPerProc]; // allocate memory

    if (parallel){
        if (numProcs == 1){
            Logger(WARN) << "Parallel execution on one process is not intended. Use serial mode instead.";
        }
        initDist.getParticles(particles, myRank*particlesPerProc, particlesPerProc);
        tree = new SubDomainTree(domainSize, confP.getVal<double>("theta"),
                                 confP.getVal<double>("softening"), timeStep, confP.getVal<bool>("hilbert"));

    } else if (numProcs == 1){
        initDist.getAllParticles(particles);
        tree = new DomainTree(domainSize, confP.getVal<double>("theta"), confP.getVal<double>("softening"), timeStep);

    } else {
        Logger(ERROR) << "Serial execution with more than one process is not possible. Aborting.";
        throw std::invalid_argument("Serial execution must not use more than one process.");
    }

    Logger(INFO) << "Inserting particles into tree ...";

    for (int i=0; i<particlesPerProc; ++i){
        tree->insertParticle(particles[i]);
    }
    delete [] particles;

    profiler.disableWrite();
    if (parallel){
        Logger(INFO) << "...done. Creating load distribution via space-filling curves ...";
        tree->guessRanges(); // guessing some ranges
        tree->newLoadDistribution();
    }

    Logger(INFO) << "... done. Computing pseudo-particles ...";
    tree->compPseudoParticles();

    Logger(INFO) << "... done. Computing forces ...";
    tree->compForce();

    profiler.enableWrite();
    N = tree->countParticles();
    profiler.value2file(ProfilerIds::N, tree->getNumParticles());

    Logger(INFO) << "... done.";
    Logger(WARN) << "Initialization finished! N = " << N;
}

BarnesHut::~BarnesHut(){
    delete tree;
}

void BarnesHut::run(){
    double t = 0.;
    int step = 0;
    while (step <= steps){
        if (step % h5DumpInterval == 0){
            Logger(INFO) << "Dumping particles to h5 ...";
            N = tree->countParticles();
            Logger(DEBUG) << "\t... writing " << N << " particles to file ...";

            std::vector<size_t> dataSpaceDims(2);
            dataSpaceDims[0] = std::size_t(N); // number of particles
            dataSpaceDims[1] = global::dim;

            std::stringstream stepss;
            stepss << std::setw(6) << std::setfill('0') << step;
            HighFive::File h5File { outDir + "/ts" + stepss.str() + ".h5",
                                   HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
                                   HighFive::MPIOFileDriver(comm, MPI_INFO_NULL) };

            // creating data sets for global information on all processes
            std::vector<keytype> ranges {};
            tree->getRanges(ranges);
            HighFive::DataSet rangesDataSet = h5File.createDataSet<keytype>("/ranges", HighFive::DataSpace::From(ranges));

            HighFive::DataSet tDataSet = h5File.createDataSet<double>("/t", HighFive::DataSpace::From(t));

            std::vector<double> com {};
            tree->getCenterOfMass(com);
            HighFive::DataSet comDataSet = h5File.createDataSet<double>("/COM", HighFive::DataSpace::From(com));

            double E_tot = tree->totalEnergy();
            HighFive::DataSet energyDataSet = h5File.createDataSet<double>("/E_tot", HighFive::DataSpace::From(E_tot));
            Logger(DEBUG) << "Total energy E_tot = " << E_tot;

            std::vector<double> L_tot(global::dim, 0.);
            tree->angularMomentum(L_tot);
            Logger(DEBUG) << "Angular Momentum L_tot = ("
                            << L_tot[0] << ", " << L_tot[1] << ", " << L_tot[2] << ")";
            HighFive::DataSet angMomDataSet = h5File.createDataSet<double>("/L_tot", HighFive::DataSpace::From(L_tot));

            if (myRank == 0){
                rangesDataSet.write(ranges);
                tDataSet.write(t);
                comDataSet.write(com);
                energyDataSet.write(E_tot);
                angMomDataSet.write(L_tot);
            }

            // each process writes its particle data
            HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
            HighFive::DataSet xDataSet = h5File.createDataSet<double>("/x", HighFive::DataSpace(dataSpaceDims));
            HighFive::DataSet vDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
            HighFive::DataSet kDataSet = h5File.createDataSet<keytype>("/key", HighFive::DataSpace(N));

            tree->dump2file(mDataSet, xDataSet, vDataSet, kDataSet);
            Logger(INFO) << "... done.";
        }

        t += timeStep;
        ++step;
        if (step > steps){
            Logger(WARN) << "Finished! N = " << N;
            break;
        }

        Logger(INFO) << "Timestep t = " << t << " ...";
        Logger(DEBUG) << "\tComputing positions and updating tree ...";
        profiler.time(ProfilerIds::timePos);
        tree->compPosition();
        profiler.time2file(ProfilerIds::timePos);

        profiler.time(ProfilerIds::timeMove);
        tree->moveParticles();
        profiler.time2file(ProfilerIds::timeMove);

        if (parallel && step % loadBalancingInterval == 0){
            Logger(INFO) << "\t... load balancing ...";
            profiler.setStep(step/loadBalancingInterval-1);
            profiler.time(ProfilerIds::timeLb);
            tree->newLoadDistribution();
            profiler.time2file(ProfilerIds::timeLb);
            profiler.setStep(step-1);
        }

        profiler.time(ProfilerIds::timePseudo);
        tree->compPseudoParticles();
        profiler.time2file(ProfilerIds::timePseudo);

        Logger(DEBUG) << "\t... computing forces ...";
        profiler.time(ProfilerIds::timeForce);
        tree->compForce();
        profiler.time2file(ProfilerIds::timeForce);

        Logger(DEBUG) << "\t... computing velocities ...";
        profiler.time(ProfilerIds::timeVel);
        tree->compVelocity();
        profiler.time2file(ProfilerIds::timeVel);
        Logger(INFO) << "... done.";

        profiler.setStep(step);
        N = tree->countParticles();
        Logger(INFO) << "Total amount of particles in simulation domain N = " << N;
        profiler.value2file(ProfilerIds::N, tree->getNumParticles());
    }
}