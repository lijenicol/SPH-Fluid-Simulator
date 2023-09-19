#pragma once
#include <glm/glm.hpp>

struct Particle
{
	glm::vec3 position, velocity, acceleration, force;

    int id;
    /// ID of next particle (-1 for no next)
	int next;

	float density;
	float pressure;
};
