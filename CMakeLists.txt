cmake_minimum_required(VERSION 3.0.0)
project(MDNCPP VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include(CTest)
enable_testing()
include(ExternalProject)

set(CMAKE_STAGING_PREFIX ${PROJECT_SOURCE_DIR}/build/stage)

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    message("Defining debug")
    add_definitions(-DMDNCPP_BUILD_DEBUG)
endif()

add_definitions(-DMDNCPP_VERSION="0.1.1")


set(CMAKE_CXX_FLAGS_Release "-O3")
set(CMAKE_CXX_FLAGS_Debug "-Wall")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ADIOS2 REQUIRED)

if(NOT TARGET MPI::MPI_CXX)
    add_library(MPI::MPI_CXX IMPORTED INTERFACE)

    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_FLAGS})
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
endif()

if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()

ExternalProject_Add(spdlog
    PREFIX spdlog
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib/spdlog
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_STAGING_PREFIX}
    # -DSPDLOG_MASTER_PROJECT=OFF
    # -DSPDLOG_BUILD_EXAMPLE=OFF
    # -DSPDLOG_BUILD_EXAMPLE_HO=OFF
    # -DSPDLOG_BUILD_TESTS=OFF
    # -DSPDLOG_BUILD_TESTS_HO=OFF
    # -DSPDLOG_BUILD_SHARED=OFF
)

ExternalProject_Add(argparse
    PREFIX argparse
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib/argparse
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_STAGING_PREFIX}
    -DBUILD_SHARED_LIBS=OFF
    -DARGPARSE_BUILD_SAMPLES=OFF
    -DARGPARSE_BUILD_TESTS=OFF
)

# ExternalProject_Add(adios2
#     PREFIX adios2
#     SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib/adios2
#     CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#     -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
#     -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
#     -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
#     -DCMAKE_INSTALL_PREFIX=${CMAKE_STAGING_PREFIX}
#     -DBUILD_SHARED_LIBS=OFF
#     -DADIOS2_BUILD_EXAMPLES=OFF
#     -DBUILD_TESTING=OFF
#     -DADIOS2_USE_MPI=ON
#     # -DMPI_CXX_FOUND=${MPI_CXX_FOUND}
#     # -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
#     # -DMPI_CXX_COMPILE_FLAGS=${MPI_CXX_COMPILE_FLAGS}
#     # -DMPI_CXX_INCLUDE_PATH=${MPI_CXX_INCLUDE_PATH}
#     # -DMPI_CXX_LINK_FLAGS=${MPI_CXX_LINK_FLAGS}
#     # -DMPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}
#     # -DMPI_C_FOUND=${MPI_C_FOUND}
#     # -DMPI_C_COMPILER=${MPI_C_COMPILER}
#     # -DMPI_C_COMPILE_FLAGS=${MPI_C_COMPILE_FLAGS}
#     # -DMPI_C_INCLUDE_PATH=${MPI_C_INCLUDE_PATH}
#     # -DMPI_C_LINK_FLAGS=${MPI_C_LINK_FLAGS}
#     # -DMPI_C_LIBRARIES=${MPI_C_LIBRARIES}
# )

ExternalProject_Add(json
    PREFIX json
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib/json
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_STAGING_PREFIX}
    -DBUILD_SHARED_LIBS=OFF
    -DJSON_BuildTests=OFF
)

# message("CMAKE_BINARY_DIR is ${CMAKE_BINARY_DIR}")

# configure_file(${PROJECT_SOURCE_DIR}/src/config.hpp.in ${CMAKE_BINARY_DIR}/generated/config.hpp @ONLY)
# include_directories(${CMAKE_BINARY_DIR}/generated/)
configure_file(${PROJECT_SOURCE_DIR}/src/config.hpp.in ${PROJECT_SOURCE_DIR}/src/config.hpp @ONLY)
include_directories(${PROJECT_SOURCE_DIR}/lib/spdlog/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/argparse/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/json/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/yas/include)
include_directories(${ADIOS2_DIR})


add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/MDNCPP.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE adios2::cxx11_mpi MPI::MPI_C)

add_dependencies(${PROJECT_NAME} spdlog)
add_dependencies(${PROJECT_NAME} argparse)
add_dependencies(${PROJECT_NAME} json)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
