###################################################################
# The following code is modified from: https://github.com/libigl/libigl-example-project/blob/master/CMakeLists.txt
# Per 3/14/2017
#This code may need to be updated as libigl evolves

cmake_minimum_required(VERSION 3.5)
project(libmeat_and_bone)

if (NOT DEFINED )
  if ("^${CMAKE_SOURCE_DIR}$" STREQUAL "^${PROJECT_SOURCE_DIR}$")
    set (MEAT_AND_BONE_IS_SUBPROJECT FALSE)
  else ()
    set (MEAT_AND_BONE_IS_SUBPROJECT TRUE)
  endif ()
endif ()

if (NOT MEAT_AND_BONE_IS_SUBPROJECT)
# Assume libigl is installed as sibling directory or given in variable LIBIGL_LOCATION
  set(LIBIGL_LOCATION "${PROJECT_SOURCE_DIR}/../libigl/"
      CACHE FILEPATH "Location of libIGL")

  ### BEGIN FROM LIBIGL-EXAMPLE-PROJECT with few changes

  # Compilation flags: adapt to your needs 
  if(MSVC)
    # Enable parallel compilation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") 
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
  else()
    # Libigl requires a modern C++ compiler that supports c++11
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()

  # Must use libigl as static
  option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)

  # add a customizable menu bar
  option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)

  # libigl options: choose your dependencies (by default everything is OFF except opengl) 
  option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
  option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
  option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)

  if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
    find_package(CGAL REQUIRED COMPONENTS Core)
    set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
    include(${CGAL_USE_FILE})
  endif()

  add_subdirectory("${LIBIGL_LOCATION}/shared/cmake" "libigl")
    
else (NOT MEAT_AND_BONE_IS_SUBPROJECT)

  message("meat-and-bone: using existing libigl")
  
endif (NOT MEAT_AND_BONE_IS_SUBPROJECT)

# libigl information 
message("meat-and-bone: libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("meat-and-bone: libigl libraries: ${LIBIGL_LIBRARIES}")
message("meat-and-bone: libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("meat-and-bone: libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("meat-and-bone: libigl definitions: ${LIBIGL_DEFINITIONS}")
message("meat-and-bone: LIBIGL_USE_STATIC_LIBRARY: ${LIBIGL_USE_STATIC_LIBRARY}")
message("meat-and-bone: CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")


# SOURCES

set(LIBMEAT_AND_BONE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB mysources "${CMAKE_CURRENT_SOURCE_DIR}/meat_and_bone/*.cpp")

# TARGETS
# We assume static builds

add_library(meat_and_bone_lib STATIC ${mysources})
add_executable(meat_and_bone_bin main.cpp)
add_executable(meat_and_bone_test test.cpp)

# TARGETS - HOW

# Library
target_include_directories(meat_and_bone_lib PUBLIC
  ${LIBMEAT_AND_BONE_INCLUDE_DIRS}
  ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(meat_and_bone_lib PUBLIC 
  ${LIBIGL_DEFINITIONS})
target_link_libraries(meat_and_bone_lib PUBLIC
  ${LIBIGL_LIBRARIES}
  ${LIBIGL_EXTRA_LIBRARIES}
  ${LIBIGL_VIEWER_EXTRA_LIBRARIES}
  ${LIBIGL_CGAL_EXTRA_LIBRARIES})

# Executable
target_include_directories(meat_and_bone_bin PRIVATE 
  ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(meat_and_bone_bin PRIVATE 
  ${LIBIGL_DEFINITIONS})
target_link_libraries(meat_and_bone_bin PRIVATE
  meat_and_bone_lib
  ${LIBIGL_LIBRARIES} 
  ${LIBIGL_VIEWER_EXTRA_LIBRARIES})
  
target_include_directories(meat_and_bone_test PRIVATE 
  ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(meat_and_bone_test PRIVATE 
  ${LIBIGL_DEFINITIONS})
target_link_libraries(meat_and_bone_test PRIVATE
  meat_and_bone_lib
  ${LIBIGL_LIBRARIES} 
  ${LIBIGL_VIEWER_EXTRA_LIBRARIES})



  # Tests

  add_test(NAME unit_tests COMMAND meat_and_bone_test)

  configure_file(testdata/palate_example/laser_meat.obj testdata/palate_example/laser_meat.obj)
  configure_file(testdata/palate_example/ct_bone.obj testdata/palate_example/ct_bone.obj)
  configure_file(testdata/test/ref.obj testdata/test/ref.obj)
  configure_file(testdata/test/ref.obj testdata/test/ref.obj)

  add_test(NAME test_palate COMMAND meat_and_bone_bin testdata/palate_example/laser_meat.obj testdata/palate_example/ct_bone.obj result.obj)