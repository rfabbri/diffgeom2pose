# optional, header-only library can be used without this library
# library is useful when the function will be used in many places
# to avoid having to include p2pt many times and slowdown the compiler
set(p2pt_sources
  p2pt.hxx p2pt.h 
)
aux_source_directory(Templates p2pt_sources)

# this library can be used if you want to precompile instantiations
# but header-only will also work, see README.md

include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_library( p2pt ${p2pt_sources} )
set_target_properties(p2pt PROPERTIES LINKER_LANGUAGE CXX)

#add_library( p2pt ${p2pt_sources}  "poly.h")
#set_target_properties(p2pt PROPERTIES LINKER_LANGUAGE CXX)
