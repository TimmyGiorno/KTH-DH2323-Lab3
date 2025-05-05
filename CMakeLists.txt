cmake_minimum_required(VERSION 3.15)
project(DH2323Lab2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O2 -Wall")

find_package(SDL2 CONFIG REQUIRED)

include_directories(
        ${CMAKE_SOURCE_DIR}/glm
        ${CMAKE_SOURCE_DIR}/SDL2Auxiliary
)

add_executable(DH2323SkeletonSDL2
        skeletonSDL2.cpp
        ${CMAKE_SOURCE_DIR}/SDL2Auxiliary/SDL2Auxiliary.cpp
)

target_link_libraries(DH2323SkeletonSDL2
        PRIVATE
        SDL2::SDL2main
        SDL2::SDL2
)