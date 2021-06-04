add_library (pumiMBBL INTERFACE pumiMBBLGPU.hpp pumi_utils.hpp)

target_include_directories(pumiMBBL PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(pumiMBBL Kokkos::kokkos)

install(TARGETS pumiMBBL DESTINATION lib)
install(FILES pumiMBBLGPU.hpp pumi_utils.hpp DESTINATION include)
