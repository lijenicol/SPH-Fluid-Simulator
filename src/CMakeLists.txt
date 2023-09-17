target_sources(sph 
    PRIVATE 
        Camera.cpp
        Geometry.cpp
        Shader.cpp
        SPHSystem.cpp
        sph.cpp
        Tester.cpp
)

add_subdirectory(imgui)