# ParaLoBstar

![paralobstarIcon](icons/paralobstar.png)

**Para**lell **Lo**ad **B**alanced Barnes-Hut algorithm for gravitational N-body simulations involving a lot of **star**s. 

## Installation

ParaLoBstar uses the header-only libraries [HighFive](https://github.com/BlueBrain/HighFive) for [HDF5](https://www.hdfgroup.org/solutions/hdf5/) file-I/O and [cxxopts](https://github.com/jarro2783/cxxopts) for command line argument handling. These dependencies can be installed with the `install` target included in the `Makefile`. 

```
$ cd paralobstar
$ make install
```

### Installation on BinAC

The target platform for ParaLoBstar is the [BinAC](https://wiki.bwhpc.de/e/Category:BwForCluster_BinAC) cluster operated by [bwHPC](https://www.bwhpc.de/index.html).

As [Boost.MPI](https://www.boost.org/doc/libs/1_77_0/doc/html/mpi.html) is used for distributed memory communication we install [boost](https://www.boost.org/users/history/version_1_77_0.html) locally before we can `make` ParaLoBstar as described in the above section. 

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
domainSize 50.
initFile plummerN10000seed3777336863.h5
outDir output
parallel true
timeStep .001
timeEnd 20.
theta .6
softening 0.01
hilbert true
h5DumpInterval 50
loadBalancingInterval 10
outputRank -1
``` 
An appropriate initial distribution file can be generated via the [ParticleDistributor](https://github.com/MichaelSt98/ParticleDistributor).

## Acknowledgements
- [CppUtils](https://github.com/MichaelSt98/CppUtils)
- [NNS](https://github.com/MichaelSt98/NNS)
- [lobster image](https://pixabay.com/images/id-2027717/)
- [star eyes image](https://pixabay.com/images/id-303363/)