add_library(lib_Utils 
Utils.cpp
)

target_link_libraries(lib_Utils
PUBLIC "${PROJECT_NAME}_compiler_flags"
)

# Any target linking to this target may need to include its source dir
target_include_directories(lib_Utils INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

#shared lib
if (BUILD_SHARED_LIBS)
  target_compile_definitions(lib_Utils PRIVATE "THIS_IS_SHARED_LIB")
endif()

# install include headers
install(FILES 
Utils.h
DESTINATION include)
