cmake_minimum_required(VERSION 3.13.0)
project(SPH)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED) 
find_package(GLUT REQUIRED) 
find_package(glm REQUIRED)

add_executable(sph "")
target_include_directories(sph PRIVATE ${GLM_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${GLUT_INCLUDE_DIR})
target_link_libraries(sph PRIVATE GLUT::GLUT GLEW::GLEW OpenGL::GL)

add_subdirectory(src)
add_subdirectory(resources)
