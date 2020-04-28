#pragma once
#include "GL/glew.h"
#include <glm/glm.hpp>
#include "Model.h"

class Particle
{
private:
	//particle properties
	float lifeTime, currentLife;
	float mass, size, elasticity;
	glm::vec3 position, velocity, acceleration;
	glm::vec3 force;
	glm::vec3 wallForce;
	
	//interactive properties
	bool fixed;

	//method for calculating new position of a particle
	void calculatePosition(float deltaTime);

	//calculates the acceleration of the particle using Newtonian mechanics
	void calculateAcceleration();

	//calculates the velocity of the particle
	void calculateVelocity(float deltaTime);
	
	//display particle currently as a box
	//Model model;

	//for creating ids
	static int pCount;
public:
	void update(float deltaTime);
	//void draw(glm::mat4 viewProjMtx, GLuint shader);

	//applies force to particle
	void applyForce(glm::vec3 force);

	//applies acceleration to particle
	void applyAcceleration(glm::vec3 acceleration);

	glm::vec3 getForce();
	glm::vec3 getPosition();
	glm::vec3 getVelocity();
	float getMass();

	void setFixed(bool val);

	float density;
	float pressure;
	float id;

	Particle(float mass, float size, float elasticity, float lifeTime, glm::vec3 position, glm::vec3 velocity, bool fixed);
	~Particle();
};

