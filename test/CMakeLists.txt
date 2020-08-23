add_executable (pumiMBBL1D_Demo dummy_pumiMBBL1D_main.c)
add_executable (pumiMBBL2D_Demo dummy_pumiMBBL2D_main.c)

target_link_libraries (pumiMBBL1D_Demo LINK_PUBLIC pumiMBBL LINK_PUBLIC m)
target_link_libraries (pumiMBBL2D_Demo LINK_PUBLIC pumiMBBL LINK_PUBLIC m)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

install(TARGETS pumiMBBL1D_Demo DESTINATION bin)
install(TARGETS pumiMBBL2D_Demo DESTINATION bin)
