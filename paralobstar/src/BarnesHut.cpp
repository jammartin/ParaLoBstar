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
    Logger(INFO) << "Reading in initial distribution ...";
    InitialDistribution initDist { initFile };

    N = initDist.getNumberOfParticles();
    particles = new Particle[N]; // allocate memory
    initDist.getParticles(particles);

    Logger(INFO) << "... done. Number of particles N = " << N;

    // TODO: implement parallel case
    //tree = parallel ? SubDomainTree() : DomainTree();

    tree = new DomainTree(domainSize,confP.getVal<double>("theta"), timeStep);

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
            Logger(INFO) << "\tDumping particles to h5 ...";
            // containers for particle data
            std::vector<double> m {};
            std::vector<std::vector<double>> x {};
            std::vector<std::vector<double>> v {};
            std::vector<keytype> k {};
            std::vector<keytype> ranges { 0UL, DomainTree::keyMax };
            N = tree->getParticleData(m, x, v, k);

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
            Logger(INFO) << "\t... done.";
        }
        if (t == timeEnd){
            break;
        }
        t += timeStep;
        ++step;
        Logger(INFO) << "Timestep t = " << t << " ...";
        Logger(DEBUG) << "\tComputing positions and updating tree ...";
        tree->compPosition();
        tree->moveParticles();
        Logger(DEBUG) << "\t... computing pseudo-particles and forces ...";
        tree->compPseudoParticles();
        tree->compForce();
        Logger(DEBUG) << "\t... computing velocities ...";
        tree->compVelocity();
        Logger(INFO) << "... done.";
    }
}