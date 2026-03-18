cmake_minimum_required(VERSION 3.18)
project(wmm-cpp)

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red "${Esc}[31m")
set(Green "${Esc}[32m")
set(Yellow "${Esc}[33m")
set(Blue "${Esc}[34m")

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE SOURCES "include/*" "src/*")
#message("\n${Yellow}==> Project: ${PROJECT_NAME} (${CMAKE_CURRENT_LIST_DIR}) ${ColourReset}")
#message("\n${Yellow}==> Sources: ${SOURCES} ${ColourReset}")

add_library(${PROJECT_NAME} STATIC ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC "include")