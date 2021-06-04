
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
# export kk=/path/to/kokkos/install

```
cd build-GPU
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$kk
cmake ../pumiMBBL -DCMAKE_CXX_COMPILER=/lore/vittav/Kokkos/kokkos/bin/nvcc_wrapper -DCMAKE_INSTALL_PREFIX=$PWD/install
make -j 8
make install
```

## test

```
./install/bin/pumiMBBL2D_Demo
```

