cmake_minimum_required(VERSION 3.16)

project(lbm-sim)

# generate compile_commands for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ensure c++ version
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(LbmLib "src/lbm.cpp")

target_include_directories(LbmLib PUBLIC "include/")

add_executable(LBMSim "src/main.cpp")

target_link_libraries(LBMSim PRIVATE LbmLib)
