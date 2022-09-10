
add_executable (particle_search_1d particle_search_1d.cpp)
add_executable (particle_search_2d particle_search_2d.cpp)

target_link_libraries (particle_search_1d pumiMBBL m)
target_link_libraries (particle_search_2d pumiMBBL m)

add_test(Particle_Search_1d particle_search_1d)
add_test(Particle_Search_2d particle_search_2d)

install(TARGETS particle_search_2d DESTINATION bin)
