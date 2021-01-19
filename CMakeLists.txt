cmake_minimum_required(VERSION 3.1)
project(Assignment1)

#-------------------------------------------------------------------------------
# Build options
#-------------------------------------------------------------------------------
option(BUILD_DEBUG     "Build with debug settings"    OFF)
option(BUILD_DOCS      "Build documentation"          OFF)

#-------------------------------------------------------------------------------
# Set target
#-------------------------------------------------------------------------------
set(APPLICATION_SOURCE
    src/texture.cpp
    src/triangulation.cpp
    src/svgparser.cpp
    src/transforms.cpp
    src/rasterizer.cpp
    src/drawrend.cpp
    src/svg.cpp
    src/main.cpp
    src/rasterizer.h
    src/drawrend.h
    src/svg.h
    src/svgparser.h
    src/texture.h
    src/transforms.h
    src/triangulation.h
)

# dirent on Windows
if (WIN32)
    list(APPEND APPLICATION_SOURCE src/dirent/dirent.c)
endif()

add_executable(draw ${APPLICATION_SOURCE})

target_link_libraries(draw PUBLIC CGL)

#-------------------------------------------------------------------------------
# Find dependencies
#-------------------------------------------------------------------------------
add_subdirectory(CGL)
target_include_directories(draw PUBLIC ${CGL_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGL_CXX_FLAGS}")

# dirent on Windows
if (WIN32)
    target_include_directories(draw PUBLIC ${CMAKE_SOURCE_DIR}/src/dirent)
endif()

#-------------------------------------------------------------------------------
# Add subdirectories
#-------------------------------------------------------------------------------

# build documentation
if(BUILD_DOCS)
  find_package(DOXYGEN)
  if(DOXYGEN_FOUND AND BUILD_DOCS)
    add_subdirectory(docs)
  endif()
endif()

# Install settings
set(CMAKE_INSTALL_PREFIX "${Assignment1_SOURCE_DIR}/")
