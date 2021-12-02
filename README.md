# ParaLoBstar

![paralobstarIcon](media/paralobstar.png)

**Para**lell **Lo**ad **B**alanced Barnes-Hut algorithm for gravitational N-body simulations involving a lot of **star**s. 

## Examples

![mpN1000160lb1np140Hi_th0_5](media/mpN1000160lb1np140Hi_th0_5.gif)

The example above shows 4 [Plummer](https://en.wikipedia.org/wiki/Plummer_model) spheres consisting of 1000160 particles in total interacting with each other gravitationally computed on 140 cores. The left shows the x-y-plane, the right the x-z-plane and the different colors indicate the computation domain of the different processes. For dynamic load-balancing Hilbert space-filling curves are used.

### Scaling


## Installation

ParaLoBstar uses the header-only libraries [HighFive](https://github.com/BlueBrain/HighFive) for [HDF5](https://www.hdfgroup.org/solutions/hdf5/) file-I/O and [cxxopts](https://github.com/jarro2783/cxxopts) for command line argument handling. These dependencies can be installed with the `install` target included in the `Makefile`. 

```
$ cd paralobstar
$ make install
```

### Installation on BinAC

The target platform for ParaLoBstar is the [BinAC](https://wiki.bwhpc.de/e/Category:BwForCluster_BinAC) cluster operated by [bwHPC](https://www.bwhpc.de/index.html).

As [Boost.MPI](https://www.boost.org/doc/libs/1_77_0/doc/html/mpi.html) is used for distributed memory communication we install [boost](https://www.boost.org/users/history/version_1_77_0.html) locally before we can `make` ParaLoBstar as described in the previous section. 

#### Installing boost locally
Download the source code into some directory

```
$ mkdir src
$ cd src
$ wget https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.gz
$ tar -xvf boost_1_77_0.tar.gz
``` 
and load the mpi module needed to compile.

```
$ module load mpi/openmpi/4.1-gnu-9.2
$ which mpic++
/opt/bwhpc/common/mpi/openmpi/4.1.1-gnu-9.2/bin/mpic++
``` 
Now we compile the boost library including mpi locally

```
$ cd boost_1_77_0
$ ./bootstrap.sh --prefix=$HOME/local --without-libraries=
```
by adding `using mpi : /opt/bwhpc/common/mpi/openmpi/4.1.1-gnu-9.2/bin/mpic++ ;` to the file `project-config.jam` and finally building the library.

```
$ ./b2 install | tee install.log
```

### Building the source
After succesful istallation ParaLoBstar can be compiled selecting one of the following targets:

```
$ make
$ make debug
$ make binac
```
**Note:** Before using `make [TARGET]` on BinAC you might want to `source binacLoadModules.sh` or `module load mpi/openmpi/4.1-gnu-9.2` to load the MPI compiler `mpic++`.

## Running simulations

```
Run parallel load balanced N-body simulations via the Barnes-Hut method.
Usage:
  paralobstar [OPTION...]

  -c, --config arg     Path to config file (default: config.info)
  -p, --profiling arg  Path to h5 profiling file (default: profiling.h5)
  -v, --verbose        More printouts for debugging
  -s, --silent         Suppress normal printouts
  -h, --help           Show this help
```

### Configuration file

The simulation parameters are read from a configuration file, which must be provided. The file has the following contents:

```
domainSize 10.          // box side length of the simulation domain
initFile plN4096.h5     // h5 file containing the initial particle distribution
outDir output           // output directory for particle distribution h5 files
parallel true           // switch for serial or parallel mode
timeStep .025           // fixed time step
timeEnd 25.             // end of the simulation
theta 1.                // clumping parameter
softening 0.032         // softening parameter
hilbert true            // switch for Lebesgue or Hilbert curves used for domain decomposition
h5DumpInterval 1        // output particle distribution each h5DumpInterval steps
loadBalancingInterval 5 // obtain a new load distribution each loadBalancingInterval steps
outputRank -1           // standard output only by process with rank outputRank (-1: all ranks)
``` 
An appropriate initial distribution file can be generated via the [ParticleDistributor](https://github.com/MichaelSt98/ParticleDistributor).

## Tools

### h5renderer

### jobCreatorBinAC
As setting up job scripts on [BinAC](https://wiki.bwhpc.de/e/Category:BwForCluster_BinAC) can be a pain in the ass for a large amount of simulations, it can be automated by `tools/jobCreatorBinAC`. 

1. Update the configuration template file `templates/config.info` but do *not* change the capital written values e.g. `INITFILE` as they are placeholders to be automatically updated
2. Update the job file `templates/job.sh` but leave all `_XYZ_` values e.g. `_JOBNAME_`  as is for the same reason.
3. Edit the following part of `createJobs.sh`

```
# BEGIN CONFIGURATION
SIM_DIR="../../simulations"
modes=("Se" "Le" "Hi")
lbIntervals=("10" "50" "100")
numProcs=("10" "28" "56")
# END CONFIGURATION
```
Finally run `./createJobs.sh [INIT_FILE_STEM]` and a directory structure is created in `$SIM_DIR`. Finally you have to copy your initial distribution file called `INIT_FILE_STEM.h5` into `$SIM_DIR` and run the autogenerated script `submit.sh` to launch your jobs all at once. 

## Sources

The algorithm is a C++ implementation of the approach described in 
> M. Griebel, S. Knapek, and G. Zumbusch. **Numerical Simulation in Molecular Dynamics**. pp 313–370. Springer 2007. ISBN: 9783540680949.



## Acknowledgements
- Utitlities used for logging and configuration taken from [CppUtils](https://github.com/MichaelSt98/CppUtils)
- Pre-developement on Branch [MolecularDynamics in NNS](https://github.com/MichaelSt98/NNS/tree/MolecularDynamics)
- ParaLoBstar's icon has been created from free content [lobster image](https://pixabay.com/images/id-2027717/) and [star eyes image](https://pixabay.com/images/id-303363/)