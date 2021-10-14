//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/BarnesHut.h"

BarnesHut::BarnesHut(ConfigParser confP) : domainSize { confP.getVal<double>("domainSize") },
                                           initFile { confP.getVal<std::string>("initFile") },
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

    int steps = (int)round(timeEnd/timeStep);
    Logger(INFO) << "Number of time steps = " << steps;

    // Initialize profiler data sets for both serial and parallel mode
    profiler.createValueDataSet<int>(ProfilerIds::N, steps);
    profiler.createTimeDataSet(ProfilerIds::timePos, steps);
    profiler.createTimeDataSet(ProfilerIds::timeMove, steps);
    profiler.createTimeDataSet(ProfilerIds::timePseudo, steps);
    profiler.createTimeDataSet(ProfilerIds::timeForce, steps);
    profiler.createTimeDataSet(ProfilerIds::timeVel, steps);

    int particlesPerProc = N/numProcs;

    if (N % numProcs != 0){
        Logger(ERROR) << "Number of particles (" << N
                      << ") must have number of processes (" << numProcs << ") as divisor. Aborting.";
        throw std::domain_error("Total amount of particles must have number of processes as divisor.");
    }

    particles = new Particle[particlesPerProc]; // allocate memory

    if (parallel && numProcs > 1){
        initDist.getParticles(particles, myRank*particlesPerProc, particlesPerProc);
        tree = new SubDomainTree(domainSize, confP.getVal<double>("theta"), timeStep);

    } else if (!parallel && numProcs == 1){
        initDist.getAllParticles(particles);
        tree = new DomainTree(domainSize,confP.getVal<double>("theta"), timeStep);

    } else {
        Logger(ERROR) << "Serial execution with more than one process is not possible. Aborting.";
        throw std::invalid_argument("Serial execution must not use more than one process.");
    }

    Logger(INFO) << "Inserting particles into tree ...";

    for (int i=0; i<particlesPerProc; ++i){
        tree->insertParticle(particles[i]);
    }
    delete [] particles;

    if (parallel){
        tree->guessRanges();
        Logger(INFO) << "...done. Sending particles ...";
        tree->sendParticles();
        Logger(INFO) << "...done. Building common coarse tree ...";
        tree->buildCommonCoarseTree();
    }

    Logger(INFO) << "... done. Computing pseudo-particles ...";
    tree->compPseudoParticles();

    Logger(INFO) << "... done. Computing forces ...";
    tree->compForce();

    Logger(INFO) << "... done. Initialization finished!";
}

BarnesHut::~BarnesHut(){
    delete tree;
}

void BarnesHut::run(){
    double t = 0.;
    int step = 0;
    while (t <= timeEnd){
        if (step % h5DumpInterval == 0){
            Logger(INFO) << "Dumping particles to h5 ...";
            N = tree->countParticles();
            Logger(DEBUG) << "\t... writing " << N << " particles to file ...";

            std::vector<size_t> dataSpaceDims(2);
            dataSpaceDims[0] = std::size_t(N); // number of particles
            dataSpaceDims[1] = global::dim;

            std::stringstream stepss;
            stepss << std::setw(6) << std::setfill('0') << step;
            HighFive::File h5File {"output/ts" + stepss.str() + ".h5",
                                   HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
                                   HighFive::MPIOFileDriver(comm, MPI_INFO_NULL) };

            HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
            HighFive::DataSet xDataSet = h5File.createDataSet<double>("/x", HighFive::DataSpace(dataSpaceDims));
            HighFive::DataSet vDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
            HighFive::DataSet kDataSet = h5File.createDataSet<keytype>("/key", HighFive::DataSpace(N));

            std::vector<keytype> ranges = tree->getRanges();
            HighFive::DataSet rangesDataSet = h5File.createDataSet<keytype>("/ranges", HighFive::DataSpace::From(ranges));
            // actually writing ranges only once
            if (myRank == 0){
                rangesDataSet.write(ranges);
            }

            tree->dump2file(mDataSet, xDataSet, vDataSet, kDataSet);
            Logger(INFO) << "... done.";
        }
        if (t >= timeEnd){
            Logger(INFO) << "Finished!";
            break;
        }

        profiler.setStep(step);
        t += timeStep;
        ++step;

        Logger(INFO) << "Timestep t = " << t << " ...";
        Logger(DEBUG) << "\tComputing positions and updating tree ...";
        profiler.time(ProfilerIds::timePos);
        tree->compPosition();
        profiler.time2file(ProfilerIds::timePos);

        profiler.time(ProfilerIds::timeMove);
        tree->moveParticles();
        profiler.time2file(ProfilerIds::timeMove);

        if (parallel){
            tree->sendParticles();
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

        profiler.value2file(ProfilerIds::N, tree->countParticles());
    }
}