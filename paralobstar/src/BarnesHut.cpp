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
    particles = new Particle[N]; // allocate memory
    initDist.getParticles(particles);

    Logger(INFO) << "... done. Number of particles N = " << N;

    int steps = (int)round(timeEnd/timeStep);

    if (parallel){
        //TODO: Implement
    } else if (numProcs == 1){
        // initialize profiler data sets
        profiler.createValueDataSet<int>(ProfilerIds::N, steps);
        profiler.createTimeDataSet(ProfilerIds::timePos, steps);
        profiler.createTimeDataSet(ProfilerIds::timeMove, steps);
        profiler.createTimeDataSet(ProfilerIds::timePseudo, steps);
        profiler.createTimeDataSet(ProfilerIds::timeForce, steps);
        profiler.createTimeDataSet(ProfilerIds::timeVel, steps);

        tree = new DomainTree(domainSize,confP.getVal<double>("theta"), timeStep);
    } else {
        Logger(ERROR) << "Serial execution with more than one process is not possible. Aborting.";
        exit(0);
    }

    Logger(INFO) << "Insertion particles into tree ...";

    for (int i=0; i<N; ++i){
        tree->insertParticle(particles[i]);
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
            std::vector<keytype> ranges { 0UL, DomainTree::keyMax };

            Logger(DEBUG) << "\t... collecting particle data ...";
            N = tree->getParticleData(m, x, v, k);
            Logger(DEBUG) << "\t... done. Writing " << N << " particles to file ...";

            std::stringstream stepss;
            stepss << std::setw(6) << std::setfill('0') << step;
            HighFive::File h5File {"output/" + stepss.str() + ".h5",
                                   HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate };

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