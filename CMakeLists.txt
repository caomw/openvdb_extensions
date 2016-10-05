cmake_minimum_required (VERSION 2.8.3)

project("OpenVDB")

if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ext/glfw/src")
  message(FATAL_ERROR "The OpenVDB dependency repositories (GLFW, etc.) are missing! "
    "You probably did not clone the project with --recursive. It is possible to recover "
    "by calling \"git submodule update --init --recursive\"")
endif()

option(OPENVDB_BUILD_VIEWER "Build OpenVDB viewer application?" ON)
option(OPENVDB_BUILD_SHARED "Build OpenVDB as a shared library?" ON)

# Required libraries for linking against OpenVDB (all targets)
set(OPENVDB_EXTRA_INCLUDE_DIRS "")
set(OPENVDB_EXTRA_LIBRARY_DIRS "")
set(OPENVDB_EXTRA_LIBRARIES "")

if(APPLE AND OPENVDB_BUILD_SHARED)
  set(CMAKE_MACOSX_RPATH ON)
endif()

include(CMakeParseArguments)

# Set a default build configuration (Release)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# This build file 
include(ExternalProject)
set(CMAKE_EXTERNAL_ARGUMENTS -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_GENERATOR=${CMAKE_GENERATOR} -Wno-dev
  -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/ext_build/dist)

# Parse arguments for nested build process
if(NOT CMAKE_EXTERNAL_ARGUMENTS)
  set(CMAKE_EXTERNAL_ARGUMENTS_LIST "")
else()
  string(REPLACE "$" ";" CMAKE_EXTERNAL_ARGUMENTS "${CMAKE_EXTERNAL_ARGUMENTS}")
  set(CMAKE_EXTERNAL_ARGUMENTS_LIST ${CMAKE_EXTERNAL_ARGUMENTS})
endif()


if (MSVC)
  # Windows-specific build flags

  # Don't advertise awkward and non-standard "secure" C++ functions
  add_definitions (/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (/D "NOMINMAX")
  include_directories(ext/glew/include)
  set(OPENVDB_EXTRA_LIBRARIES opengl32 glew)

  set(OPENVDB_EXTRA_LIBRARIES ${OPENVDB_EXTRA_LIBRARIES}
    $<$<CONFIG:Debug>:zlibstaticd>
    $<$<CONFIG:RelWithDebInfo>:zlibstatic>
    $<$<CONFIG:Release>:zlibstatic> 
    $<$<CONFIG:MinSizeRel>:zlibstatic>
  )

  # Statically link against the C++ runtime library, also apply these settings to nested projects
  set(CompilerFlags
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
  foreach(CompilerFlag ${CompilerFlags})
    string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
    string(REPLACE "/MDd" "/MTd" ${CompilerFlag} "${${CompilerFlag}}")
    set(CMAKE_EXTERNAL_ARGUMENTS ${CMAKE_EXTERNAL_ARGUMENTS} -D${CompilerFlag}:INTERNAL=${${CompilerFlag}})
  endforeach()
elseif(APPLE)
  # OSX-specific build flags

  # Find standard libraries needed for OpenGL/GLFW
  find_library(cocoa_library Cocoa)
  find_library(opengl_library OpenGL)
  find_library(corevideo_library CoreVideo)
  find_library(iokit_library IOKit)
  set(OPENVDB_EXTRA_LIBRARIES ${cocoa_library} ${opengl_library} ${corevideo_library} ${iokit_library} z)

  # Compile in C++11 mode
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")

  # Try to auto-detect a suitable OSX SDK
  execute_process(COMMAND bash -c "xcodebuild -version -sdk | grep MacOSX | grep Path | head -n 1 | cut -f 2 -d ' '" OUTPUT_VARIABLE CMAKE_OSX_SYSROOT)
  string(REGEX REPLACE "(\r?\n)+$" "" CMAKE_OSX_SYSROOT "${CMAKE_OSX_SYSROOT}")
  string(REGEX REPLACE "^.*X([0-9.]*).sdk$" "\\1" CMAKE_OSX_DEPLOYMENT_TARGET "${CMAKE_OSX_SYSROOT}")

  # Also use the same SDK in nested projects
  set(CMAKE_EXTERNAL_ARGUMENTS ${CMAKE_EXTERNAL_ARGUMENTS} -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET} -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT})
elseif("${CMAKE_SYSTEM}" MATCHES "Linux")
  # Linux-specific build flags
  set(OPENVDB_EXTRA_LIBRARIES dl GL GLU Xxf86vm Xrandr Xinerama Xcursor Xi X11 pthread z)

  # Compile in C++11 mode
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Disable some clang warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-private-field")
endif()

if(WIN32)
  # Build zlib (only needed on Windows, Linux/OSX already have it)
  ExternalProject_Add(zlib
    URL             "${CMAKE_CURRENT_SOURCE_DIR}/ext/zlib"
    PREFIX          "ext_build"
    CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS}
  )
endif(WIN32)

# Build Intel Thread Building Blocks
ExternalProject_Add(tbb_vdb
  URL             "${CMAKE_CURRENT_SOURCE_DIR}/ext/tbb"
  PREFIX          "ext_build"
  CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS} -DBUILD_SHARED_LIBS=OFF -DNAMESPACE_VERSIONING=NO
)

# Build IlmBase, a support library needed by OpenEXR
ExternalProject_Add(IlmBase
  URL             "${CMAKE_CURRENT_SOURCE_DIR}/ext/openexr/IlmBase"
  PREFIX          "ext_build"
  CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS} -DILMBASE_BUILD_SHARED_LIBS=OFF -DBUILD_SHARED_LIBS=OFF -DNAMESPACE_VERSIONING=NO
)

# Build OpenEXR
ExternalProject_Add(OpenEXR
  URL             "${CMAKE_CURRENT_SOURCE_DIR}/ext/openexr/OpenEXR"
  PREFIX          "ext_build"
  CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS} -DOPENEXR_BUILD_SHARED_LIBS=OFF -DILMBASE_PACKAGE_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/ext_build/dist -DBUILD_SHARED_LIBS=OFF -DNAMESPACE_VERSIONING=NO
)

# Build GLFW
ExternalProject_Add(glfw
  URL             "${CMAKE_CURRENT_SOURCE_DIR}/ext/glfw"
  PREFIX          "ext_build"
  CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS} -DGLFW_BUILD_TESTS=0 -DGLFW_BUILD_EXAMPLES=0 -DGLFW_BUILD_DOCS=0
)

# Build Alloy
ExternalProject_Add(alloy
  URL             "${CMAKE_CURRENT_SOURCE_DIR}/ext/alloy"
  PREFIX          "ext_build"
  CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS} -DALLOY_BUILD_EXAMPLE=0
)

# Make sure everything is built in the right oder
add_dependencies(OpenEXR IlmBase)
if(WIN32)
  add_dependencies(OpenEXR zlib)
endif(WIN32)

# Pass external arguments to nested build processes
string(REGEX REPLACE ";" "$" CMAKE_EXTERNAL_ARGUMENTS_STR "${CMAKE_EXTERNAL_ARGUMENTS}")
list(APPEND CMAKE_EXTERNAL_ARGUMENTS -D CMAKE_EXTERNAL_ARGUMENTS:STRING=${CMAKE_EXTERNAL_ARGUMENTS_STR})

# Boost
find_package(Boost REQUIRED COMPONENTS system iostreams thread)
list(APPEND OPENVDB_EXTRA_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
list(APPEND OPENVDB_EXTRA_LIBRARY_DIRS ${Boost_LIBRARY_DIRS})

# Include & library directories
list(APPEND OPENVDB_EXTRA_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ext_build/dist/include)
list(APPEND OPENVDB_EXTRA_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ext_build/dist/include/OpenEXR)
list(APPEND OPENVDB_EXTRA_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ext_build/dist/lib)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/openvdb)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/zlib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/tbb)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/tbb/include)

set(OPENVDB_EXTRA_LIBRARIES IlmImf-2_2 IlmThread-2_2 Iex-2_2 IexMath-2_2 Imath-2_2 Half tbb glfw3 ${OPENVDB_EXTRA_LIBRARIES} ${Boost_LIBRARIES})

include_directories(${OPENVDB_EXTRA_INCLUDE_DIRS})
link_directories(${OPENVDB_EXTRA_LIBRARY_DIRS})

#--------------------------------------------------------------------
# Define defintions
add_definitions(-DOPENVDB_TOOLS_RAYTRACER_USE_EXR)
if(NOT BUILD_SHARED_LIBS)
    add_definitions(-DOPENVDB_LIB)
    add_definitions(-DOPENVDB_STATIC)
endif()

#--------------------------------------------------------------------
# Include directories


# openvdb header files
file(GLOB_RECURSE openvdb_io_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/io/*.h)
file(GLOB_RECURSE openvdb_math_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/math/*.h)
file(GLOB_RECURSE openvdb_metadata_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/metadata/*.h)
file(GLOB_RECURSE openvdb_tools_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/tools/*.h)
file(GLOB_RECURSE openvdb_tree_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/tree/*.h)
file(GLOB_RECURSE openvdb_util_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/util/*.h)

set(openvdb_headers 
  openvdb/openvdb/Exceptions.h
  openvdb/openvdb/Grid.h
  openvdb/openvdb/Metadata.h
  openvdb/openvdb/openvdb.h
  openvdb/openvdb/Platform.h
  openvdb/openvdb/PlatformConfig.h
  openvdb/openvdb/Types.h
  openvdb/openvdb/version.h
)

list(APPEND openvdb_headers 
  ${openvdb_io_headers} 
  ${openvdb_math_headers} 
  ${openvdb_metadata_headers} 
  ${openvdb_tools_headers}
  ${openvdb_tree_headers}
  ${openvdb_util_headers}
)

# openvdb source files
file(GLOB_RECURSE openvdb_io_sources ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/io/*.cc)
file(GLOB_RECURSE openvdb_math_sources ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/math/*.cc)
file(GLOB_RECURSE openvdb_metadata_sources ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/metadata/*.cc)
file(GLOB_RECURSE openvdb_util_sources ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/util/*.cc)

set(openvdb_sources 
  openvdb/openvdb/Grid.cc
  openvdb/openvdb/openvdb.cc
  openvdb/openvdb/Platform.cc
)

list(APPEND openvdb_sources
  ${openvdb_io_sources} 
  ${openvdb_math_sources} 
  ${openvdb_metadata_sources}
  ${openvdb_util_sources}
)

# set library type
if(OPENVDB_BUILD_SHARED)
  set(OPENVDB_LIBRARY_TYPE "SHARED")
else()
  set(OPENVDB_LIBRARY_TYPE "STATIC")
endif()

# openvdb library
add_library(openvdb ${OPENVDB_LIBRARY_TYPE} ${openvdb_sources} ${openvdb_headers})
target_compile_definitions(openvdb PRIVATE OPENVDB_PRIVATE)
#set_target_properties(openvdb PROPERTIES COMPILE_DEFINITIONS "OPENVDB_PRIVATE")
target_link_libraries(openvdb ${OPENVDB_EXTRA_LIBRARIES})

# openvdb viewer library
file(GLOB_RECURSE openvdb_viewer_headers ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/viewer/*.h)
file(GLOB_RECURSE openvdb_viewer_sources ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/viewer/*.cc)
add_library(openvdb_viewer STATIC ${openvdb_viewer_sources} ${openvdb_viewer_headers})
target_compile_definitions(openvdb_viewer PRIVATE OPENVDB_USE_GLFW_3 GL_GLEXT_PROTOTYPES)

# command line utilities
add_executable(openvdb_print openvdb/openvdb/cmd/openvdb_print/main.cc)
target_link_libraries(openvdb_print openvdb)
# add_executable(openvdb_render openvdb/openvdb/cmd/openvdb_render/main.cc)
# target_link_libraries(openvdb_render openvdb)
add_executable(openvdb_view openvdb/openvdb/cmd/openvdb_view/main.cc)
target_link_libraries(openvdb_view openvdb_viewer openvdb)
add_custom_command(TARGET openvdb_viewer PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/viewer
                   ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb/openvdb_viewer
)

# dependencies to ensure the correct build order
add_dependencies(openvdb OpenEXR)
add_dependencies(openvdb tbb_vdb)
add_dependencies(openvdb alloy)
add_dependencies(openvdb glfw)
add_dependencies(openvdb_print openvdb)
add_dependencies(openvdb_view openvdb)
# add_dependencies(openvdb_render openvdb)

if(OPENVDB_BUILD_SHARED)
  if(NOT MSVC)
    install(TARGETS openvdb LIBRARY DESTINATION lib)
  endif()
else()
  install(TARGETS openvdb ARCHIVE DESTINATION lib)
endif()
install(TARGETS openvdb_print openvdb_view RUNTIME DESTINATION bin)
install(DIRECTORY openvdb/openvdb DESTINATION include)

get_directory_property(OPENVDB_HAS_PARENT PARENT_DIRECTORY)
if(OPENVDB_HAS_PARENT)
  # This project is included from somewhere else. Export OPENVDB_EXTRA_XXX variables
  set(OPENVDB_EXTRA_INCLUDE_DIRS ${OPENVDB_EXTRA_INCLUDE_DIRS} PARENT_SCOPE)
  set(OPENVDB_EXTRA_LIBRARY_DIRS ${OPENVDB_EXTRA_LIBRARY_DIRS} PARENT_SCOPE)
  set(OPENVDB_EXTRA_LIBRARIES ${OPENVDB_EXTRA_LIBRARIES} PARENT_SCOPE)
endif()
