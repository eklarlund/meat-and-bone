cmake_minimum_required(VERSION 2.8.12)
project(libmeat_and_bone)

set(LIBIGL_LOCATION "${PROJECT_SOURCE_DIR}/../libigl/"
    CACHE FILEPATH "Location of libIGL")


### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   OFF)

### libIGL options: choose your dependencies (by default everything is OFF, in this example we need the viewer) ###
option(LIBIGL_WITH_BBW              "Use BBW"            ON)
find_package(CGAL QUIET)
option(LIBIGL_WITH_CGAL             "Use CGAL"           "${CGAL_FOUND}")
option(LIBIGL_WITH_COMISO           "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK             "Use CORK"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         ON)
option(LIBIGL_WITH_LIM              "Use LIM"            ON)
find_package(MATLAB QUIET)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         "${MATLAB_FOUND}")
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          "${MOSEK_FOUND}")
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         ON)
option(LIBIGL_WITH_PNG              "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       ON)
option(LIBIGL_WITH_XML              "Use XML"            ON)
### End   to be tested ----

### libIGL options: decide if you want to use the functionalities that depends on cgal
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, cgal strange build system fails without it
  find_package(CGAL REQUIRED)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

### Adding libIGL: choose the path to your local copy libIGL ###
### This is going to compile everything you requested ###
add_subdirectory("${LIBIGL_LOCATION}/shared/cmake" "libigl")


### Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w") # disable all warnings (not ideal but...)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # disable all warnings (not ideal but...)
endif()



include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})


########

# SOURCES

set(LIBMEAT_AND_BONE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB mysources "${CMAKE_CURRENT_SOURCE_DIR}/meat_and_bone/*.cpp")

# TARGETS

add_library(libmeat_and_bone_lib STATIC ${mysources})
add_executable(libmeat_and_bone_bin main.cpp)

# TARGETS - HOW

# Library
target_include_directories(libmeat_and_bone_lib PUBLIC
  ${LIBMEAT_AND_BONE_INCLUDE_DIRS}
  ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(libmeat_and_bone_lib PUBLIC 
  ${LIBIGL_DEFINITIONS})
target_link_libraries(libmeat_and_bone_lib PUBLIC
  ${LIBIGL_LIBRARIES}
  ${LIBIGL_VIEWER_EXTRA_LIBRARIES})


# Executable
target_include_directories(libmeat_and_bone_bin PRIVATE 
  ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(libmeat_and_bone_bin PRIVATE 
  ${LIBIGL_DEFINITIONS})
target_link_libraries(libmeat_and_bone_bin PRIVATE
  libmeat_and_bone_lib
  ${LIBIGL_LIBRARIES} 
  ${LIBIGL_VIEWER_EXTRA_LIBRARIES})
