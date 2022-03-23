# This is p2pt/tests/CMakeLists.txt

#--------------------------------------------------------------

set(p2pt_test_all_sources
  test-driver.cxx
  test-p2pt.cxx
)

aux_source_directory(Templates p2pt_test_all_sources)
add_executable(p2pt_test_all ${p2pt_test_all_sources})

target_link_libraries(p2pt_test_all p2pt ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}vbl ${VXL_LIB_PREFIX}testlib)
add_test( NAME test-p2pt COMMAND p2pt_test_all test_p2pt )

add_executable(p2pt_test_include test-include.cxx)
target_link_libraries(p2pt_test_include p2pt)
