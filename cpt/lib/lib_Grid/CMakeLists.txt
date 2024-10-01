add_library(lib_Grid 
Grid.cpp
Cell.cpp 
Face.cpp
Node.cpp
Zone.cpp
)

target_link_libraries(lib_Grid 
PUBLIC "${PROJECT_NAME}_compiler_flags" lib_Utils lib_Event_Execution_Time
)

# Any target linking to this target may need to include its source dir
target_include_directories(lib_Grid INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

#shared lib
if (BUILD_SHARED_LIBS)
  target_compile_definitions(lib_Grid PRIVATE "THIS_IS_SHARED_LIB")
endif()

#OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lib_Grid PUBLIC OpenMP::OpenMP_CXX)
endif()


# install include headers
install(FILES 
Bin.h 
Cell.h
Face.h
Grid.h
Node.h
Props.h
Zone.h
DESTINATION include)
