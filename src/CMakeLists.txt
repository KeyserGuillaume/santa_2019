cmake_minimum_required(VERSION 3.10)
project(Santa2019)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11")

add_executable(santa main.cpp Flow.cpp Flow.h constants.h tools.h tools.cpp LocalSearch.cpp LocalSearch.h Assignment.cpp Assignment.h Presets.cpp Presets.h day_sub_problem.cpp day_sub_problem.h day_cost_lower_bound.cpp day_cost_lower_bound.h BranchAndBound.cpp BranchAndBound.h)