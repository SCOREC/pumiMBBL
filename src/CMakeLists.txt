add_library (pumiMBBL pumi_initiate.c pumi_routines.c)

target_include_directories(pumiMBBL PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS pumiMBBL DESTINATION lib)
install(FILES pumiMBBL.h pumi_initiate.h pumi_routines.h DESTINATION include)
