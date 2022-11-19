#pragma once

#include <vector>
#include <thread>

#include "GL/glew.h"
#include <glm/glm.hpp>

#include "Particle.h"
#include "Geometry.h"

#define THREAD_COUNT 8

class SPHSystem
{
private:
	//particle data
	unsigned int numParticles;
	std::vector<std::vector<Particle*>> neighbouringParticles;
	bool started;

	//initializes the particles that will be used
	void initParticles();

	// Creates hash table for particles in infinite domain
	void buildTable();

	// Sphere geometry for rendering
	Geometry* sphere;
	glm::mat4 sphereScale;
	glm::mat4* sphereModelMtxs;
	GLuint vbo;

	// Threads and thread blocks
	std::thread threads[THREAD_COUNT];
	int blockBoundaries[THREAD_COUNT + 1];
	int tableBlockBoundaries[THREAD_COUNT + 1];

public:
	SPHSystem(unsigned int numParticles, float mass, float restDensity, float gasConst, float viscosity, float h, float g, float tension);
	~SPHSystem();

	//kernel/fluid constants
	float POLY6, SPIKY_GRAD, SPIKY_LAP, GAS_CONSTANT, MASS, H2, SELF_DENS;

	//fluid properties
	float restDensity;
	float viscosity, h, g, tension;

	std::vector<Particle*> particles;
	Particle** particleTable;
	glm::ivec3 getCell(Particle *p) const;

	// std::mutex mtx;

	//updates the SPH system
	void update(float deltaTime);

	//draws the SPH system & its particles
	void draw(const glm::mat4& viewProjMtx, GLuint shader);

	void print();

	void reset();
	void startSimulation();
};

