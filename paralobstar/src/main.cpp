//
// Created by Johannes Martin on 17.09.21.
//

#include <mpi.h>

// header only libraries
#include <cxxopts.hpp>

#include "../include/Logger.h"
#include "../include/ConfigParser.h"
#include "../include/H5Profiler.h"
#include "../include/BarnesHut.h"

structlog LOGCFG = {};

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    cxxopts::Options cmdLineOptions { "paralobstar",
                                      "Run parallel load balanced N-body simulations via the Barnes-Hut method." };
    cmdLineOptions.add_options()
            ("c,config", "Path to config file", cxxopts::value<std::string>()->default_value("config.info"))
            ("p,profiling", "Path to h5 profiling file", cxxopts::value<std::string>()->default_value("profiling.h5"))
            ("v,verbose", "More printouts for debugging")
            ("h,help", "Show this help");

    auto cmdLineOpts = cmdLineOptions.parse(argc, argv);

    ConfigParser configParser { cmdLineOpts["config"].as<std::string>() };

    // initialize Logger
    LOGCFG.headers = true;
    LOGCFG.level = cmdLineOpts.count("verbose") ? DEBUG : INFO;
    LOGCFG.myRank = myRank;
    LOGCFG.outputRank = configParser.getVal<int>("outputRank");

    // create singleton instance
    H5Profiler &profiler = H5Profiler::getInstance(cmdLineOpts["profiling"].as<std::string>(), myRank, numProcs);

    BarnesHut algorithm { configParser, myRank, numProcs };

    algorithm.run();

    MPI_Finalize();
    return 0;
}