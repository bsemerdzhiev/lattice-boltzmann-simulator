cmake_minimum_required(VERSION 3.16)

project(lbm-sim)

# generate compile_commands for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ensure c++ version
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(LbmLib "src/lbm.cpp" "src/benchmark_simulation.cpp"
                   "src/render_simulation.cpp" "src/cell.cpp")

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)

add_library(glad STATIC "frameworks/glad/src/glad.c")
target_include_directories(glad PUBLIC "frameworks/glad/include")

target_include_directories(LbmLib PUBLIC "include/")
target_link_libraries(LbmLib PRIVATE glfw glad OpenGL::GL)

add_executable(LBMSim "src/main.cpp")

target_link_libraries(LBMSim PRIVATE LbmLib)
