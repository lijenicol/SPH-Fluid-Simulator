#pragma once
#include <glm/glm.hpp>

struct Particle
{
	glm::vec3 position, velocity, acceleration, force;
	Particle* next;
	float density;
	float pressure;
};
