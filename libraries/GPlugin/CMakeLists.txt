add_library(GPLUGIN SHARED
  GPluginLoader.cxx
  GPluginManager.cxx
  GPluginManifest.cxx
  GPluginRegistry.cxx
)

target_include_directories(GPLUGIN PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(GPLUGIN PUBLIC ROOT::Core)
target_compile_features(GPLUGIN PUBLIC cxx_std_17)

set_target_properties(GPLUGIN PROPERTIES
  EXPORT_NAME Plugin
)

install(
  TARGETS GPLUGIN
  EXPORT GrootPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
