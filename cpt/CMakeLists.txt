cmake_minimum_required(VERSION "3.15.0")

project("pt" VERSION 1.0)

# C++ standard
add_library("${PROJECT_NAME}_compiler_flags" INTERFACE)

#Selecting C++17 to enable using the `reduce` function
#target_compile_features("${PROJECT_NAME}_compiler_flags" INTERFACE cxx_std_14)
target_compile_features("${PROJECT_NAME}_compiler_flags" INTERFACE cxx_std_17)

# add compiler warning flags just when building this project via
# the BUILD_INTERFACE genex
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
#The following generator expression was introduced in CMake 3.15
target_compile_options("${PROJECT_NAME}_compiler_flags" INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused>>"
  "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

# set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

# set building shared libraries option.
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

# add libraries
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/lib_Utils")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/lib_Grid")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/lib_Particle_Tracking")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/lib_Event_Execution_Time")

# add executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} 
    PUBLIC "${PROJECT_NAME}_compiler_flags" lib_Particle_Tracking lib_Utils lib_Event_Execution_Time
    )

# install libs
set(installable_libs lib_Grid lib_Particle_Tracking lib_Utils lib_Event_Execution_Time "${PROJECT_NAME}_compiler_flags")
#install(TARGETS ${installable_libs} ARCHIVE DESTINATION lib)
install(TARGETS ${installable_libs} DESTINATION lib)

# add the install targets
#install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

# setup installer
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)
