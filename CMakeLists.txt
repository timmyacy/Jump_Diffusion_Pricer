cmake_minimum_required(VERSION 3.15)
project(JumpDiffusionPricer)

set(CMAKE_CXX_STANDARD 17)
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/models
    ${CMAKE_SOURCE_DIR}/pricers
    ${CMAKE_SOURCE_DIR}/fft
)
add_executable(JumpDiffusionPricer
    main.cpp
    models/black_scholes.cpp
    pricers/carr_madan.cpp
    pricers/mc.cpp
    fft/fft.cpp
    models/variance_gamma.cpp
)

target_include_directories(JumpDiffusionPricer PRIVATE ${CMAKE_SOURCE_DIR})
