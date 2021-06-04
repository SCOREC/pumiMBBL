add_library (pumiMBBL INTERFACE pumiMBBLGPU.hpp pumi_utils.hpp)

target_include_directories(pumiMBBL INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pumiMBBL INTERFACE Kokkos::kokkos)

install(TARGETS pumiMBBL DESTINATION lib)
install(FILES pumiMBBLGPU.hpp pumi_utils.hpp DESTINATION include)
