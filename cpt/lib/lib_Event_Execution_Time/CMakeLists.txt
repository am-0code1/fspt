add_library(lib_Event_Execution_Time 
Event_Execution_Time.cpp
)

target_link_libraries(lib_Event_Execution_Time 
PUBLIC "${PROJECT_NAME}_compiler_flags"
)

# Any target linking to this target may need to include its source dir
target_include_directories(lib_Event_Execution_Time INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

#shared lib
if (BUILD_SHARED_LIBS)
  target_compile_definitions(lib_Event_Execution_Time PRIVATE "THIS_IS_SHARED_LIB")
endif()

# install include headers
install(FILES 
Event_Execution_Time.h 
DESTINATION include)
