# CMakeList.txt : CMake project for diffgeom2pose, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (cmd "diffgeom2pose.cxx" )
target_link_libraries(cmd p2pt)

# TODO: Add tests and install targets if needed.
