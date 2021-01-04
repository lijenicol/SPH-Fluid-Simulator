#pragma once
#include <vector>
#include "GL/glew.h"
#include <glm/glm.hpp>
#include "Particle.h"
#include "Geometry.h"

class SPHSystem
{
private:
	//particle data
	unsigned int numParticles;
	std::vector<Particle*> particles;
	std::vector<std::vector<Particle*>> neighbouringParticles;
	bool started;

	//kernel/fluid constants
	float POLY6, SPIKY_GRAD, SPIKY_LAP, GAS_CONSTANT, MASS, H2, SELF_DENS;

	//fluid properties
	float restDensity;
	float viscosity, h, g, tension;

	//initializes the particles that will be used
	void initParticles();

	//searches the neighbours of all particles
	void searchNeighbours();

	// Updates position of all particles
	// SPH should handle this instead of the particle class
	void updateParticles(float deltaTime);

	// Sphere geometry for rendering
	Geometry* sphere;
	glm::mat4 sphereScale;

public:
	SPHSystem(unsigned int numParticles, float mass, float restDensity, float gasConst, float viscosity, float h, float g, float tension);
	~SPHSystem();

	//updates the SPH system
	void update(float deltaTime);

	//draws the SPH system & its particles
	void draw(glm::mat4 viewProjMtx, GLuint shader);

	void print();

	void reset();
	void startSimulation();
};

