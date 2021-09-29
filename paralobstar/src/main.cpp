//
// Created by Johannes Martin on 17.09.21.
//

// header only libraries
#include <cxxopts.hpp>

#include "../include/Logger.h"
#include "../include/ConfigParser.h"
#include "../include/BarnesHut.h"

structlog LOGCFG = {};

int main(int argc, char *argv[]){

    cxxopts::Options cmdLineOptions { "paralobstar",
                                      "Run parallel load balanced N-body simulations via the Barnes-Hut method." };
    cmdLineOptions.add_options()
            ("c,config", "Path to config file", cxxopts::value<std::string>()->default_value("config.info"))
            ("v,verbose", "More printouts for debugging")
            ("h,help", "Show this help");

    auto cmdLineOpts = cmdLineOptions.parse(argc, argv);

    LOGCFG.headers = true;
    LOGCFG.level = cmdLineOpts.count("verbose") ? DEBUG : INFO;

    ConfigParser configParser { cmdLineOpts["config"].as<std::string>() };

    BarnesHut algorithm { configParser };

    algorithm.run();

    return 0;
}