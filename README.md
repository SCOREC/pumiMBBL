
## layout

- src - source codes (will be installed with `make
install`)
- test - example code that exercises APIs
- CMakeLists.txt - cmake files that specify how to build the library and test

## setup

```
git clone https://github.com/SCOREC/pumiMBBL.git
mkdir build-GPU
```

## build

The following assumes that a valid C and C++ compiler, and `cmake`, are in your PATH.

`CMAKE_INSTALL_PREFIX` is the path where the library, headers, and test binary
are installed.

`kk` is the path where kokkos is installed
Set path as `export kk=/path/to/kokkos/install`

Load necessary modules:
```
module load gcc/7.3.0-bt47fwr mpich/3.3-diz4f6i cmake/3.20.0 cuda/10.2
```

```
cd build-GPU
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$kk
cmake ../pumiMBBL -DCMAKE_CXX_COMPILER=path/to/nvcc_wrapper -DCMAKE_INSTALL_PREFIX=$PWD/install
make -j 8
make install
```

## test

```
./install/bin/pumiMBBL2D_Demo
```

