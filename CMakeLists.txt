cmake_minimum_required(VERSION 3.0.0)
project(MDNCPP VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)

include(CTest)
enable_testing()
include(ExternalProject)

set(CMAKE_STAGING_PREFIX ${PROJECT_SOURCE_DIR}/build/stage)
set(STAGING_DIR ${PROJECT_SOURCE_DIR}/build/stage)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

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
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_INSTALL_PREFIX=${STAGING_DIR}
    -DSPDLOG_BUILD_SHARED=OFF
)

ExternalProject_Add(argparse
    PREFIX argparse
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib/argparse
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_INSTALL_PREFIX=${STAGING_DIR}
    -DBUILD_SHARED_LIBS=OFF
)

ExternalProject_Add(adios2
    PREFIX adios2
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib/adios2
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_INSTALL_PREFIX=${STAGING_DIR}
    -DBUILD_SHARED_LIBS=OFF
    -DADIOS2_BUILD_EXAMPLES=OFF
    -DBUILD_TESTING=OFF
)

include_directories(${PROJECT_SOURCE_DIR}/lib/spdlog/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/argparse/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/adios2/bindings)

add_executable(MDNCPP ${PROJECT_SOURCE_DIR}/src/MDNCPP.cpp)

message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
target_link_libraries(MDNCPP PUBLIC MPI::MPI_CXX)
target_link_libraries(MDNCPP PUBLIC OpenMP::OpenMP_CXX)

add_dependencies(${PROJECT_NAME} spdlog)
add_dependencies(${PROJECT_NAME} argparse)
add_dependencies(${PROJECT_NAME} adios2)
# target_link_libraries(MDNCPP STATIC ${STAGING_DIR}/lib/spdlog.a)
# target_link_libraries(MDNCPP PUBLIC ${STAGING_DIR}/include)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
