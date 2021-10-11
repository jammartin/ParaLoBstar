//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/BarnesHut.h"

BarnesHut::BarnesHut(ConfigParser confP,
                     int _myRank, int _numProcs) : domainSize { confP.getVal<double>("domainSize") },
                                                   initFile { confP.getVal<std::string>("initFile") },
                                                   parallel { confP.getVal<bool>("parallel") },
                                                   timeStep { confP.getVal<double>("timeStep") },
                                                   timeEnd { confP.getVal<double>("timeEnd") },
                                                   h5DumpInterval { confP.getVal<int>("h5DumpInterval") },
                                                   loadBalancingInterval { confP.getVal<int>("loadBalancingInterval") },
                                                   myRank { _myRank }, numProcs { _numProcs }
{
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

    if (parallel){
        initDist.getParticles(particles, myRank*particlesPerProc, particlesPerProc);
        tree = new SubDomainTree(domainSize, confP.getVal<double>("theta"), timeStep);

    } else if (numProcs == 1){
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
    delete [] particles;
    delete tree;
}

void BarnesHut::run(){
    double t = 0.;
    int step = 0;
    while (t <= timeEnd){
        if (step % h5DumpInterval == 0){
            Logger(INFO) << "Dumping particles to h5 ...";
            // containers for particle data
            std::vector<double> m {};
            std::vector<std::vector<double>> x {};
            std::vector<std::vector<double>> v {};
            std::vector<keytype> k {};
            std::vector<keytype> ranges = tree->getRanges();

            Logger(DEBUG) << "\t... collecting particle data ...";
            N = tree->getParticleData(m, x, v, k);
            Logger(DEBUG) << "\t... done. Writing " << N << " particles to file ...";

            std::stringstream stepss;
            stepss << std::setw(6) << std::setfill('0') << step;
            HighFive::File h5File {"output/ts" + stepss.str() + ".h5",
                                   HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
                                   HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL) };

            HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m",  HighFive::DataSpace::From(m));
            HighFive::DataSet xDataSet = h5File.createDataSet<double>("/x",  HighFive::DataSpace::From(x));
            HighFive::DataSet vDataSet = h5File.createDataSet<double>("/v",  HighFive::DataSpace::From(v));
            HighFive::DataSet kDataSet = h5File.createDataSet<keytype>("/key", HighFive::DataSpace::From(k));
            HighFive::DataSet rangesDataSet = h5File.createDataSet<keytype>("/ranges", HighFive::DataSpace::From(ranges));

            mDataSet.write(m);
            xDataSet.write(x);
            vDataSet.write(v);
            kDataSet.write(k);
            rangesDataSet.write(ranges);
            Logger(INFO) << "... done.";
        }
        if (t == timeEnd){
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