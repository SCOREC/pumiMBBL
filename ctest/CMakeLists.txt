add_executable (particle_search_1d particle_search_1d.cpp)

target_link_libraries (particle_search_1d pumiMBBL m)

add_test(Particle_Search_1d particle_search_1d)
