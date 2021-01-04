#include "Particle.h"
#include <iostream>
#include <glm/gtx/string_cast.hpp>

int Particle::pCount = 0;

Particle::Particle(float mass, float size, glm::vec3 position, glm::vec3 velocity) {
	this->mass = mass;
	this->size = size;
	this->position = position;
	this->velocity = velocity;

	force = glm::vec3(0);
	acceleration = glm::vec3(0);

	density = 0.0f;
	pressure = 0.0f;
	id = pCount++;
}

Particle::~Particle() {}