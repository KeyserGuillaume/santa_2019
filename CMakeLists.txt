cmake_minimum_required(VERSION 3.10)
project(Santa2019)

# set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11")
set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -std=c++14 -DIL_STD -m64 -Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG           "${CMAKE_CXX_FLAGS_DEBUG} -Og -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE         "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3 -flto")
set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(LINKER_OPTIONS                  -flto -Wl,--no-as-needed)

if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
    set(CMAKE_CXX_FLAGS_DEBUG       "${CMAKE_CXX_FLAGS_DEBUG} -save-temps=obj")
endif()

find_package(Cplex)
include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})

set(SOURCE_FILES
        src/main.cpp src/Flow.cpp src/Flow.h src/constants.h src/tools.h src/tools.cpp src/LocalSearch.cpp src/LocalSearch.h src/Assignment.cpp src/Assignment.h src/Presets.cpp src/Presets.h src/day_sub_problem.cpp src/day_sub_problem.h src/day_cost_lower_bound.cpp src/day_cost_lower_bound.h src/BranchAndBound.cpp src/BranchAndBound.h src/dwlb.cpp src/dwlb.h)

add_executable(         santa ${SOURCE_FILES})
target_link_libraries(  santa PRIVATE ${CPLEX_LIBRARIES})
target_link_libraries(  santa PRIVATE dl)
