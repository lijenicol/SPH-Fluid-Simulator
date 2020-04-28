#include "Particle.h"
#include <iostream>
#include <glm/gtx/string_cast.hpp>

int Particle::pCount = 0;

Particle::Particle(float mass, float size, float elasticity, float lifeTime, glm::vec3 position, glm::vec3 velocity, bool fixed) {
	this->mass = mass;
	this->size = size;
	this->elasticity = elasticity;
	this->lifeTime = lifeTime;
	this->position = position;
	this->velocity = velocity;
	this->fixed = fixed;

	currentLife = 0;
	force = glm::vec3(0);
	wallForce = glm::vec3(0);
	acceleration = glm::vec3(0);

	density = 0.0f;
	pressure = 0.0f;
	id = pCount++;
	
	//initialize model
	//model.MakeBox(glm::vec3(-1), glm::vec3(1));
}

Particle::~Particle() {}

void Particle::update(float deltaTime) {
	//zero out wall force
	wallForce = glm::vec3(0);

	if (!fixed) {
		//update the position of the particle
		calculatePosition(deltaTime);
	}

	//zero out the force to be added next frame
	force = glm::vec3(0);
	acceleration = glm::vec3(0);

	//apply wall force to be used next frame
	//applyForce(wallForce);
}

//void Particle::draw(glm::mat4 viewProjMtx, GLuint shader) {
//	//generate model matrix
//	glm::mat4 M = glm::translate(glm::mat4(1), position) * glm::scale(glm::mat4(1), glm::vec3(size));
//
//	//draw the model 
//	model.Draw(M, viewProjMtx, shader);
//}

void Particle::calculateAcceleration() {
	acceleration += force / mass;
}

void Particle::calculateVelocity(float deltaTime) {
	velocity += acceleration * deltaTime;
}

void Particle::calculatePosition(float deltaTime) {	
	//calculate acceleration and velocity
	calculateAcceleration();
	calculateVelocity(deltaTime);

	//clamp velocity to a max velocity governed by the CFL law
	float maxVelocity = 0.4 * size / deltaTime;
	if (glm::length(velocity) > maxVelocity) {
		//svelocity = velocity / glm::length(velocity) * maxVelocity;
	}

	position += velocity * deltaTime;

	// Handle collisions with box
	float boxWidth = 1.0f;
	if (position.y < size) {
		position.y = -position.y + 2 * size + 0.001f;
		velocity.y = -elasticity * velocity.y;
	}

	if (position.x < size - boxWidth) {
		position.x = -position.x + 2 * (size - boxWidth) + 0.001f;
		velocity.x = -elasticity * velocity.x;
	}

	if (position.x > -size + boxWidth) {
		position.x = -position.x + 2 * -(size - boxWidth) - 0.001f;
		velocity.x = -elasticity * velocity.x;
	}

	if (position.z < size - boxWidth) {
		position.z = -position.z + 2 * (size - boxWidth) + 0.001f;
		velocity.z = -elasticity * velocity.z;
	}

	if (position.z > -size + boxWidth) {
		position.z = -position.z + 2 * -(size - boxWidth) - 0.001f;
		velocity.z = -elasticity * velocity.z;
	}
}

void Particle::applyForce(glm::vec3 force) {
	this->force += force;
}

void Particle::applyAcceleration(glm::vec3 acceleration) {
	this->acceleration += acceleration;
}

void Particle::setFixed(bool val) { fixed = val; }

glm::vec3 Particle::getPosition() { return position; }
glm::vec3 Particle::getVelocity() { return velocity; }
glm::vec3 Particle::getForce() { return force; }
float Particle::getMass() { return mass; }