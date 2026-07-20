add_executable(GCanvasPluginInteractionTest
  GCanvasPluginInteractionTest.cxx
)

target_link_libraries(GCanvasPluginInteractionTest PRIVATE
  GROOTLIB
  GPLUGIN
  ROOT::Core
  ROOT::Graf
  ROOT::Gpad
  ROOT::Hist
)

add_test(NAME GCanvasPluginInteraction COMMAND GCanvasPluginInteractionTest)
