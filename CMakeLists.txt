cmake_minimum_required(VERSION 3.0.0)
project(MDNCPP VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

# find_package(spdlog REQUIRED)
# if(NOT spdlog_FOUND)
# 	message(SEND_ERROR "Failed to find spdlog")
# 	return()
# else()
#     include_directories(${spdlog_INCLUDE_DIR})
# endif()

add_library(spdlog STATIC lib/spdlog/include/spdlog/spdlog.h)
# set(TARGET_scr src/MDNCPP.cpp)
# set(TARGET src/MDNCPP)

add_executable(MDNCPP src/MDNCPP.cpp)

# target_link_libraries(main ${spdlog_LIBRARY})
# target_link_libraries(MDNCPP spdlog::spdlog)
target_link_libraries(MDNCPP spdlog)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
