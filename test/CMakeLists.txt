add_executable (pumiMBBL_Demo dummy_pumi_main.c)

target_link_libraries (pumiMBBL_Demo LINK_PUBLIC pumiMBBL LINK_PUBLIC m)
install(TARGETS pumiMBBL_Demo DESTINATION bin)
