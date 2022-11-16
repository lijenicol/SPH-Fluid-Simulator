#pragma once
#include "GL/glew.h"
#include <glm/glm.hpp>

class Particle
{
private:
	//for creating ids
	static int pCount;
public:
	// Attributes of particle
	float mass, size, elasticity;
	glm::vec3 position, velocity, acceleration;
	glm::vec3 force;

	// For linked list
	Particle* next;

	float density;
	float pressure;
	float id;

	Particle(float mass, float size, glm::vec3 position, glm::vec3 velocity);
	~Particle();
};

