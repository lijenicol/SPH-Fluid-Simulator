target_sources(sph 
    PRIVATE 
        Camera.cpp
        Geometry.cpp
        neighborTable.cpp
        Shader.cpp
        SPHSystem.cpp
        sph.cpp
        Tester.cpp
)

add_subdirectory(imgui)
add_subdirectory(kernels)
