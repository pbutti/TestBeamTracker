cmake_minimum_required(VERSION 3.5)

project(TestBeamTracker LANGUAGES CXX)

set(ROOT_COMPONENTS Core GenVector Hist Tree TreePlayer)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# required packages

# Enable C++11 threading support, with a preference for pthread
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
#find_package(Boost 1.64 REQUIRED COMPONENTS program_options)
find_package(ROOT 6.10 REQUIRED COMPONENTS ${ROOT_COMPONENTS} Geom)

# ACTS
#add_subdirectory(external/acts-core)
find_package(Acts REQUIRED)
message(STATUS "ACTS FOUND::INCLUDE DIR=" ${ACTS_INCLUDE_DIR})

#Eigen
find_package(Eigen3 3.2.9 REQUIRED)
message(STATUS "EIGEN::" ${EIGEN3_INCLUDE_DIRS})

#boost
find_package(Boost 1.62 REQUIRED COMPONENTS program_options unit_test_framework)

file(GLOB_RECURSE TestBeamTracker_files "src/*.cpp" "src/*.cxx" "include/TestBeamTracker/*.hpp" "include/TestBeamTracker/*.h")
add_library(TestBeamTracker SHARED ${TestBeamTracker_files})
target_include_directories (TestBeamTracker PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(TestBeamTracker PUBLIC ${ROOT_INCLUDE_DIRS} ${ACTS_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
target_link_libraries(TestBeamTracker PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(TestBeamTracker PUBLIC ${TestBeamTrackerClassesDict})
target_link_directories(TestBeamTracker PUBLIC ${ACTS_LIBRARY_DIR} ${ROOT_LIBRARY_DIR})
install(TARGETS TestBeamTracker LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_subdirectory(util)

