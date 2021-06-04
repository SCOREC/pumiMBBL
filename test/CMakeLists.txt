add_executable (pumiMBBL2D_Demo pumiMBBL_test.cpp)

target_link_libraries (pumiMBBL2D_Demo LINK_PUBLIC pumiMBBL LINK_PUBLIC m)

install(TARGETS pumiMBBL2D_Demo DESTINATION bin)
