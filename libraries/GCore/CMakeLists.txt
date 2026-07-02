
#this will compile the current library.
add_library(GCORE SHARED
  GHTTPConnection.cxx
  GObjectManager.cxx
  GOptions.cxx
)

target_include_directories(GCORE PUBLIC
  ${CMAKE_SOURCE_DIR}/include 
  ${ROOT_INCLUDE_DIRS}
) 

root_generate_dictionary(G__GCORE
  GObjectManager.h
  MODULE GCORE
  LINKDEF Linkdef.h)
target_link_libraries(GCORE PUBLIC ${ROOT_LIBRARIES} GROOTLIB)
