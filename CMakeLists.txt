cmake_minimum_required(VERSION 3.18.0)
project(MDNCPP VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include(CTest)
enable_testing()

add_definitions(-DMDNCPP_VERSION="0.1.2")

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    message("Defining debug")
    add_definitions(-DMDNCPP_CONF_BUILD_DEBUG)
endif()

if(DEFINED MDNCPP_ENABLE_PROFILING)
    if(MDNCPP_ENABLE_PROFILING STREQUAL ON)
        message("Profiling: on")
        add_definitions(-DMDNCPP_CONF_ENABLE_PROFILING)
    else()
        message("Profiling: off")
    endif()
else()
    message("Profiling: off")
endif()

if(DEFINED MDNCPP_CALC_ENTHROPY)
    if(MDNCPP_CALC_ENTHROPY STREQUAL ON)
        message("Enthropy calculation: on")
        add_definitions(-DMDNCPP_CONF_CALC_ENTHROPY)
    else()
        message("Enthropy calculation: off")
    endif()
else()
    message("Enthropy calculation: off")
endif()

if(DEFINED MDNCPP_KE_PE_PRESENT)
    if(MDNCPP_KE_PE_PRESENT STREQUAL ON)
        message("Ke, pe pesent: yes")
        add_definitions(-DMDNCPP_CONF_KE_PE_PRESENT)
    else()
        message("Ke, pe pesent: no")
    endif()
else()
    message("Ke, pe pesent: no")
endif()

if(DEFINED MDNCPP_TRACE_OUT)
    if(MDNCPP_TRACE_OUT STREQUAL ON)
        message("Trace out: on")
        add_definitions(-DMDNCPP_CONF_TRACE_OUT)
    else()
        message("Trace out: off")
    endif()
else()
    message("Trace out: off")
endif()

set(CMAKE_STAGING_PREFIX ${PROJECT_SOURCE_DIR}/build/stage)
set(CMAKE_CXX_FLAGS_Release "-O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_Debug "-O3 -Wall -Wextra -Wpedantic")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(fmt REQUIRED)

if(spdlog_FOUND)
    message("spdlog found")
else()
    message("spdlog NOT found")
endif()

include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

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
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()

include_directories(${PROJECT_SOURCE_DIR}/src/CSV/CSV)

configure_file(${PROJECT_SOURCE_DIR}/src/config.hpp.in ${PROJECT_SOURCE_DIR}/src/config.hpp @ONLY)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/MDNCPP.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE adios2::cxx11_mpi MPI::MPI_C)
target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog)
target_link_libraries(${PROJECT_NAME} PRIVATE argparse)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
