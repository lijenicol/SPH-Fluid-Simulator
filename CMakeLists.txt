cmake_minimum_required(VERSION 3.25.0)
project(SPH LANGUAGES CXX CUDA)
set(CMAKE_CUDA_STANDARD 17)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED) 
find_package(GLUT REQUIRED) 
find_package(glm REQUIRED)

add_executable(sph "")
target_include_directories(sph PRIVATE ${GLM_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${GLUT_INCLUDE_DIR} src)
target_link_libraries(sph PRIVATE GLUT::GLUT GLEW::GLEW OpenGL::GL)

add_subdirectory(src)
add_subdirectory(resources)
