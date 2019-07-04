Derived from: https://cmake.org/examples/

## layout

- src - source codes (will be installed with `make
install`)
- test - example code that exercises APIs
- CMakeLists.txt - cmake files that specify how to build the library and test

## setup

```
git clone https://github.com/SCOREC/pumiMBBL.git
mkdir buildpumiMBBL
```

## build

The following assumes that a valid C and C++ compiler, and `cmake`, are in your PATH.

`CMAKE_INSTALL_PREFIX` is the path where the library, headers, and test binary
are installed.

```
cd buildpumiMBBL
cmake ../pumiMBBL -DCMAKE_INSTALL_PREFIX=$PWD/install
make 
make install
```

## test

```
./install/bin/pumiMBBL_Demo
```

