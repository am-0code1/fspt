add_library(lib_Particle_Tracking 
Particle_Tracking.cpp
Particle.cpp
)

target_link_libraries(lib_Particle_Tracking 
PUBLIC "${PROJECT_NAME}_compiler_flags" lib_Grid lib_Utils
)

# Any target linking to this target may need to include its source dir
target_include_directories(lib_Particle_Tracking INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

#shared lib
if (BUILD_SHARED_LIBS)
  target_compile_definitions(lib_Particle_Tracking PRIVATE "THIS_IS_SHARED_LIB")
endif()


# install include headers
install(FILES 
Particle_Tracking.h
Particle.h
DESTINATION include)
