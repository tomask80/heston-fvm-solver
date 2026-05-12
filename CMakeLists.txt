cmake_minimum_required(VERSION 3.14)
project(HestonModel VERSION 1.0 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    Eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
)
FetchContent_MakeAvailable(Eigen)

# Pridanie hlavičkových súborov z tvojej include zložky
include_directories(include)


add_executable(HestonSolver src/main.cpp src/heston.cpp)

target_link_libraries(HestonSolver PRIVATE Eigen3::Eigen)