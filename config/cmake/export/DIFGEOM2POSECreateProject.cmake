# This DIFFGEOM2POSECreateProject.cmake file handles the creation of files needed by
# other client projects that use DIFFGEOM2POSE.  Nothing is built by this
# CMakeLists.txt file.  This CMakeLists.txt file must be processed by
# CMake after all the other CMakeLists.txt files in the DIFFGEOM2POSE tree,
# which is why the add_subdirectory(config/cmake/export) command is at the end
# of the top level CMakeLists.txt file.

# Save library dependencies.
#set(DIFFGEOM2POSE_CMAKE_DOXYGEN_DIR  ${DIFFGEOM2POSE_ROOT_SOURCE_DIR}/config/cmake/doxygen)

get_property(DIFFGEOM2POSETargets_MODULES GLOBAL PROPERTY DIFFGEOM2POSETargets_MODULES)

set(DIFFGEOM2POSE_CONFIG_CMAKE_DIR "share/DIFFGEOM2POSE/cmake")
if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
   set(INTERFACE_LINK_OPTION "")
else()
   set(INTERFACE_LINK_OPTION "EXPORT_LINK_INTERFACE_LIBRARIES")
endif()

if(DIFFGEOM2POSETargets_MODULES)
  export(TARGETS
    ${DIFFGEOM2POSETargets_MODULES}
    APPEND
    FILE "${CMAKE_CURRENT_BINARY_DIR}/DIFFGEOM2POSETargets.cmake"
    ${INTERFACE_LINK_OPTION}
  )
  install(EXPORT ${DIFFGEOM2POSE_INSTALL_EXPORT_NAME} DESTINATION ${DIFFGEOM2POSE_CONFIG_CMAKE_DIR}
          COMPONENT Development)
endif()

# Create the DIFFGEOM2POSEConfig.cmake file for the build tree.
configure_file(${DIFFGEOM2POSE_CMAKE_DIR}/DIFFGEOM2POSEConfig.cmake.in
               ${PROJECT_BINARY_DIR}/DIFFGEOM2POSEConfig.cmake @ONLY)

configure_file(${DIFFGEOM2POSE_CMAKE_DIR}/DIFFGEOM2POSEConfig_export.cmake.in
               ${PROJECT_BINARY_DIR}/config/cmake/export/DIFFGEOM2POSEConfig.cmake
               @ONLY)

install(FILES
  ${PROJECT_BINARY_DIR}/config/cmake/export/DIFFGEOM2POSEConfig.cmake
  ${DIFFGEOM2POSE_CMAKE_DIR}/DIFFGEOM2POSEStandardOptions.cmake
  ${DIFFGEOM2POSE_CMAKE_DIR}/UseDIFFGEOM2POSE.cmake
  DESTINATION ${DIFFGEOM2POSE_CONFIG_CMAKE_DIR}
)
