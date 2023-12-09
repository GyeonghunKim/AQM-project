cmake_minimum_required(VERSION  3.6)
project(aqm_final_cversion)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
add_executable(aqm_final_cversion
#        mc.cpp
        main.cpp
)
