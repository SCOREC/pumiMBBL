add_executable (pumiMBBL1D_Demo pumiMBBL_1D_test.cpp)
add_executable (pumiMBBL2D_Demo pumiMBBL_2D_test.cpp)

target_link_libraries (pumiMBBL2D_Demo LINK_PUBLIC pumiMBBL LINK_PUBLIC m)
target_link_libraries (pumiMBBL1D_Demo LINK_PUBLIC pumiMBBL LINK_PUBLIC m)

install(TARGETS pumiMBBL2D_Demo DESTINATION bin)
install(TARGETS pumiMBBL1D_Demo DESTINATION bin)
