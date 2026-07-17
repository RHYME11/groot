add_library(GPluginMock SHARED MockPlugin.cxx)
target_include_directories(GPluginMock PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_library(GPluginWrongApiMock SHARED MockPlugin.cxx)
target_include_directories(GPluginWrongApiMock PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_definitions(GPluginWrongApiMock PRIVATE MOCK_API_VERSION=999)

add_library(GPluginMissingSymbolsMock SHARED MissingSymbolsPlugin.cxx)

add_executable(GPluginManagerTest GPluginManagerTest.cxx)
target_link_libraries(GPluginManagerTest PRIVATE GPLUGIN ROOT::Core)

add_test(
  NAME GPluginManagerTest
  COMMAND GPluginManagerTest
    $<TARGET_FILE:GPluginMock>
    $<TARGET_FILE:GPluginWrongApiMock>
    $<TARGET_FILE:GPluginMissingSymbolsMock>
    ${CMAKE_CURRENT_BINARY_DIR}/runtime
)
