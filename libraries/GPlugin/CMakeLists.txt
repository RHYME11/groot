add_library(GPLUGIN SHARED
  GPluginLoader.cxx
  GPluginManager.cxx
  GPluginManifest.cxx
  GPluginRegistry.cxx
  GPluginSessionRegistry.cxx
)

target_include_directories(GPLUGIN PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(GPLUGIN PUBLIC ROOT::Core)
target_compile_features(GPLUGIN PUBLIC cxx_std_17)
