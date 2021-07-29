add_library (pumiMBBL INTERFACE pumiMBBLGPU.hpp pumiMBBL_initiate.hpp pumiMBBL_utils.hpp pumiMBBL_routines.hpp)

target_include_directories(pumiMBBL INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pumiMBBL INTERFACE Kokkos::kokkos)

install(TARGETS pumiMBBL DESTINATION lib)
install(FILES pumiMBBLGPU.hpp pumiMBBL_utils.hpp pumiMBBL_initiate.hpp pumiMBBL_routines.hpp DESTINATION include)
