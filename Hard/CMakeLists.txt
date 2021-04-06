cmake_minimum_required (VERSION 2.6)
# Include the base directories so that every module has it included
include_directories(eigen)
include_directories(statistic)

# Add the subdirectories to the project
add_subdirectory(statistic)
