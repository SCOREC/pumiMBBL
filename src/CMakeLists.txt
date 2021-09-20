set(HEADERS
  pumiMBBLGPU.hpp
  pumiMBBLGPU_impl.hpp
  pumiMBBL_utils.hpp
  pumiMBBL_routines.hpp
  pumiMBBL_initiate.hpp
)

set(SOURCES
  pumiMBBL_initiate.cpp
  pumiMBBL_routines.cpp
)

add_library (pumiMBBL-core ${SOURCES})

target_include_directories(pumiMBBL-core INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(pumiMBBL-core Kokkos::kokkos)

pumiMBBL_export_lib(pumiMBBL-core "${HEADERS}")

add_library (pumiMBBL INTERFACE)
target_link_libraries(pumiMBBL INTERFACE ${pumiMBBL_EXPORTED_TARGETS})
export_target(pumiMBBL)
