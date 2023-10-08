#pragma once
#include <glm/glm.hpp>

struct Particle
{
	glm::vec3 position, velocity, acceleration, force;
    uint32_t id;
    float density;
    float pressure;
    uint16_t hash;
};
