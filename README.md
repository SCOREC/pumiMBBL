
## layout

- src - source codes (will be installed with `make
install`)
- test - example code that exercises APIs
- CMakeLists.txt - cmake files that specify how to build the library and test

## setup

```
git clone https://github.com/SCOREC/pumiMBBL.git
```

## build

The following assumes that a valid C and C++ compiler, and `cmake`, are in your PATH.

`CMAKE_INSTALL_PREFIX` is the path where the library, headers, and test binary
are installed.

`kk` is the path where kokkos is installed (either GPU or OpenMP based installation).

Set path as `export kk=/path/to/kokkos/install`

`kk_compiler` is the path to kokkos compiler
Set path as `export kk_compiler=/path/to/nvcc_wrapper` for GPU build.

Ignore this for OpenMP build

Load necessary modules:
```
module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core
module load \
gcc/10.1.0 \
openmpi \
cmake/3.20.0  \
cuda/11.4
```

Building with GPU
```
mkdir build-GPU
cd build-GPU
export CMAKE_PREFIX_PATH=$kk/lib64/cmake/Kokkos:$CMAKE_PREFIX_PATH
cmake ../pumiMBBL -DCMAKE_CXX_COMPILER=$kk_compiler -DCMAKE_INSTALL_PREFIX=$PWD/install # on GPU
make -j 8
make install
ctest
```

Building with OpenMP
```
mkdir build-omp
cd build-omp
export CMAKE_PREFIX_PATH=$kk/lib64/cmake/Kokkos:$CMAKE_PREFIX_PATH
cmake ../pumiMBBL -DCMAKE_INSTALL_PREFIX=$PWD/install # on OpenMP
make -j 8
make install
ctest
```

## Building with [spack](https://github.com/spack/spack) package manager
```
# Configuring spack environment
# One time step -- Create spack environment <env_name> as
# spack env create <env_name>

spack env create pumimbbl_serial # for serial build
spack env create pumimbbl_omp # for openmp build
spack env create pumimbbl_cuda # for cuda build

# Activate environment
spack env activate <env_name>

# One time step -- Add necessary dependency for build
spack add kokkos # for serial build
spack add kokkos +openmp # for openmp build
spack add kokkos +cuda +cuda_lambda +wrapper cuda_arch=XX # for cuda build, ensure correct value for XX

# One time step -- Install package
spack install

# Building pumiMBBL with spack

# For serial pumiMBBL build
spack activate pumimbbl_serial
mkdir build-serial
cd build-serial
cmake ../pumiMBBL -DCMAKE_INSTALL_PREFIX=$PWD/install
make -j install
ctest

# For openMP pumiMBBL build
spack activate pumimbbl_omp
mkdir build-omp
cd build-omp
cmake ../pumiMBBL -DCMAKE_INSTALL_PREFIX=$PWD/install
make -j install
ctest

# For cuda pumiMBBL build
spack activate pumimbbl_cuda
mkdir build-cuda
cd build-cuda
cmake ../pumiMBBL -DCMAKE_INSTALL_PREFIX=$PWD/install -DCMAKE_CXX_COMPILER=nvcc_wrapper
make -j install
ctest
```

## Documentation

Building Doxygen Documentation
```
cd pumiMBBL/
doxygen
```
HTML file `index.html` inside folder `html` will contain the generated code documentation

Building Library Reference Documentation
```
cd pumiMBBL/doc
pdflatex pumiMBBL-GPU.tex
```
PDF file `pumiMBBL-GPU.pdf` containing implemented mesh concepts and documentation will be generated

## pseudo particle push tests

```
./install/bin/pumiMBBL1D_Demo # for 1D particle tracking test
./install/bin/pumiMBBL2D_Demo # for 2D particle tracking test
```
