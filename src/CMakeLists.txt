target_sources(sph 
    PRIVATE 
        Camera.cpp
        Camera.h
        Core.h
        Geometry.cpp
        Geometry.h
        Particle.cpp
        Particle.h
        Shader.cpp
        Shader.h
        SPHSystem.cpp
        SPHSystem.h
        Tester.cpp
        Tester.h
)

add_subdirectory(imgui)