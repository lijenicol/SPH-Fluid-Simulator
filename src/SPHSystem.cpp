#include "SPHSystem.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <mutex>

#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>

#define PI 3.14159265f
#define TABLE_SIZE 1000000

// This will lock across all instances of SPHSystem's,
// however since we only really have one instance, this
// should be okay for now  
std::mutex mtx;

/**
 * Hashes the position of a cell, giving where
 * the cell goes in the particle table 
 */
uint getHash(const glm::ivec3& cell) {
	return (
		(uint)(cell.x * 73856093) 
	  ^ (uint)(cell.y * 19349663) 
	  ^ (uint)(cell.z * 83492791)
	) % TABLE_SIZE;
}

SPHSystem::SPHSystem(unsigned int numParticles, float mass, float restDensity, float gasConst, float viscosity, float h, float g, float tension) {
	this->numParticles = numParticles;
	this->restDensity = restDensity;
	this->viscosity = viscosity;
	this->h = h;
	this->g = g;
	this->tension = tension;

	POLY6 = 315.0f / (64.0f * PI * pow(h, 9));
	SPIKY_GRAD = -45.0f / (PI * pow(h, 6));
	SPIKY_LAP = 45.0f / (PI * pow(h, 6));
	MASS = mass;
	GAS_CONSTANT = gasConst;
	H2 = h * h;
	SELF_DENS = MASS * POLY6 * pow(h, 6);

	//setup densities & volume
	int cbNumParticles = numParticles * numParticles * numParticles;
	neighbouringParticles.resize(cbNumParticles);
	particles.resize(cbNumParticles);

	//initialize particles
	initParticles();

	// Load in sphere geometry and allocate matrice space
	sphere = new Geometry("resources/lowsphere.obj");
	sphereScale = glm::scale(glm::vec3(h/2.f));
	sphereModelMtxs = new glm::mat4[cbNumParticles];
	
	// Generate VBO for sphere model matrices
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(glm::mat4), &sphereModelMtxs[0], GL_DYNAMIC_DRAW);

	// Setup instance VAO
	glBindVertexArray(sphere->vao);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), 0);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)sizeof(glm::vec4));
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(sizeof(glm::vec4)*2));
	glEnableVertexAttribArray(5);
	glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(sizeof(glm::vec4)*3));

	glVertexAttribDivisor(2,1);
	glVertexAttribDivisor(3,1);
	glVertexAttribDivisor(4,1);
	glVertexAttribDivisor(5,1);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	//start init
	started = false;

	// Allocate table memory
	particleTable = (Particle **)malloc(sizeof(Particle *) * TABLE_SIZE);

	// Init block boundaries (for funcs that loop through particles)
	blockBoundaries[0] = 0;
	int blockSize = particles.size() / THREAD_COUNT;
	for (int i = 1; i < THREAD_COUNT; i++) {
		blockBoundaries[i] = i * blockSize;
	}
	blockBoundaries[THREAD_COUNT] = particles.size();

	// Init table block boundaries (for table clearing func)
	tableBlockBoundaries[0] = 0;
	blockSize = TABLE_SIZE / THREAD_COUNT;
	for (int i = 1; i < THREAD_COUNT; i++) {
		tableBlockBoundaries[i] = i * blockSize;
	}
	tableBlockBoundaries[THREAD_COUNT] = TABLE_SIZE;
}

SPHSystem::~SPHSystem() {
	// free table
	free(particleTable);
	free(sphereModelMtxs);

	//delete particles
	particles.clear();
	particles.shrink_to_fit();

	//delete neighbouring particles
	neighbouringParticles.clear();
	neighbouringParticles.shrink_to_fit();
}

void SPHSystem::initParticles() {
	std::srand(1024);
	float particleSeperation = h + 0.01f;
	for (int i = 0; i < numParticles; i++) {
		for (int j = 0; j < numParticles; j++) {
			for (int k = 0; k < numParticles; k++) {
				// dam like particle positions
				float ranX = (float(rand()) / float((RAND_MAX)) * 0.5f - 1) * h / 10;
				float ranY = (float(rand()) / float((RAND_MAX)) * 0.5f - 1) * h / 10;
				float ranZ = (float(rand()) / float((RAND_MAX)) * 0.5f - 1) * h / 10;
				glm::vec3 nParticlePos = glm::vec3(i * particleSeperation + ranX - 1.5f, j * particleSeperation + ranY + h + 0.1f, k * particleSeperation + ranZ - 1.5f);

				//create new particle
				Particle* nParticle = new Particle(MASS, h,	nParticlePos, glm::vec3(0));

				//append particle
				particles[i + (j + numParticles * k) * numParticles] = nParticle;
			}
		}
	}
}

/**
 * Parallel computation function for calculating density
 * and pressures of particles in the given SPH System.
 */
void parallelDensityAndPressures(const SPHSystem& sphSystem, int start, int end) {
	float massPoly6Product = sphSystem.MASS * sphSystem.POLY6;
	
	for (int i = start; i < end; i++) {
		float pDensity = 0;
		Particle* pi = sphSystem.particles[i];
		glm::ivec3 cell = sphSystem.getCell(pi);

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					glm::ivec3 near_cell = cell + glm::ivec3(x, y, z);
					uint index = getHash(near_cell);
					Particle* pj = sphSystem.particleTable[index];
					
					// Iterate through cell linked list
					while (pj != NULL) {
						float dist2 = glm::length2(pj->position - pi->position);
						if (dist2 < sphSystem.H2 && pi != pj) {
							pDensity += massPoly6Product * glm::pow(sphSystem.H2 - dist2, 3);
						}
						pj = pj->next;
					}
				}
			}
		}
		
		// Include self density (as itself isn't included in neighbour)
		pi->density = pDensity + sphSystem.SELF_DENS;

		// Calculate pressure
		float pPressure = sphSystem.GAS_CONSTANT * (pi->density - sphSystem.restDensity);
		pi->pressure = pPressure;
	}
}

/**
 * Parallel computation function for calculating forces
 * of particles in the given SPH System.
 */
void parallelForces(const SPHSystem& sphSystem, int start, int end) {
	for (int i = start; i < end; i++) {
		Particle* pi = sphSystem.particles[i];
		pi->force = glm::vec3(0);
		glm::ivec3 cell = sphSystem.getCell(pi);

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					glm::ivec3 near_cell = cell + glm::ivec3(x, y, z);
					uint index = getHash(near_cell);
					Particle* pj = sphSystem.particleTable[index];

					// Iterate through cell linked list
					while (pj != NULL) {
						float dist2 = glm::length2(pj->position - pi->position);
						if (dist2 < sphSystem.H2 && pi != pj) {
							//unit direction and length
							float dist = sqrt(dist2);
							glm::vec3 dir = glm::normalize(pj->position - pi->position);

							//apply pressure force
							glm::vec3 pressureForce = -dir * sphSystem.MASS * (pi->pressure + pj->pressure) / (2 * pj->density) * sphSystem.SPIKY_GRAD;
							pressureForce *= std::pow(sphSystem.h - dist, 2);
							pi->force += pressureForce;

							//apply viscosity force
							glm::vec3 velocityDif = pj->velocity - pi->velocity;
							glm::vec3 viscoForce = sphSystem.viscosity * sphSystem.MASS * (velocityDif / pj->density) * sphSystem.SPIKY_LAP * (sphSystem.h - dist);
							pi->force += viscoForce;
						}
						pj = pj->next;
					}
				}
			}
		}
	}
}

/**
 * Parallel computation function moving positions
 * of particles in the given SPH System.
 */
void parallelUpdateParticlePositions(const SPHSystem& sphSystem, float deltaTime, int start, int end) {
	for (int i = start; i < end; i++) {
		Particle *p = sphSystem.particles[i];

		//calculate acceleration and velocity
		glm::vec3 acceleration = p->force / p->density + glm::vec3(0, sphSystem.g, 0);
		p->velocity += acceleration * deltaTime;
		
		// Update position
		p->position += p->velocity * deltaTime;

		// Handle collisions with box
		float boxWidth = 3.f;
		float elasticity = 0.5f;
		if (p->position.y < p->size) {
			p->position.y = -p->position.y + 2 * p->size + 0.0001f;
			p->velocity.y = -p->velocity.y * elasticity;
		}

		if (p->position.x < p->size - boxWidth) {
			p->position.x = -p->position.x + 2 * (p->size - boxWidth) + 0.0001f;
			p->velocity.x = -p->velocity.x * elasticity;
		}

		if (p->position.x > -p->size + boxWidth) {
			p->position.x = -p->position.x + 2 * -(p->size - boxWidth) - 0.0001f;
			p->velocity.x = -p->velocity.x * elasticity;
		}

		if (p->position.z < p->size - boxWidth) {
			p->position.z = -p->position.z + 2 * (p->size - boxWidth) + 0.0001f;
			p->velocity.z = -p->velocity.z * elasticity;
		}

		if (p->position.z > -p->size + boxWidth) {
			p->position.z = -p->position.z + 2 * -(p->size - boxWidth) - 0.0001f;
			p->velocity.z = -p->velocity.z * elasticity;
		}
	}
}

void SPHSystem::update(float deltaTime) {
	if (!started) return;
	
	// To increase system stability, a fixed deltaTime is set
	deltaTime = 0.003f;

	// Build particle hash table
	buildTable();

	// Calculate densities and pressures of particles
	for (int i = 0; i < THREAD_COUNT; i++) {
		threads[i] = std::thread(parallelDensityAndPressures, std::ref(*this), blockBoundaries[i], blockBoundaries[i + 1]);
	}
	for (std::thread& thread : threads) {
		thread.join();
	}

	// Calclulate forces of particles
	for (int i = 0; i < THREAD_COUNT; i++) {
		threads[i] = std::thread(parallelForces, std::ref(*this), blockBoundaries[i], blockBoundaries[i + 1]);
	}
	for (std::thread& thread : threads) {
		thread.join();
	}

	// Update positions of all particles
	for (int i = 0; i < THREAD_COUNT; i++) {
		threads[i] = std::thread(parallelUpdateParticlePositions, std::ref(*this), deltaTime, blockBoundaries[i], blockBoundaries[i + 1]);
	}
	for (std::thread& thread : threads) {
		thread.join();
	}
}

void SPHSystem::draw(const glm::mat4& viewProjMtx, GLuint shader) {
	// Calculate model matrices for each particle
	for (int i = 0; i < particles.size(); i++) {
		glm::mat4 translate = glm::translate(particles[i]->position);
		sphereModelMtxs[i] = translate * sphereScale;
	}

	// Send matrix data to GPU
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	void* data = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	memcpy(data, sphereModelMtxs, sizeof(glm::mat4) * particles.size());
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Draw instanced particles
	glUseProgram(shader);
	glUniformMatrix4fv(glGetUniformLocation(shader, "viewProjMtx"), 1, false, (float*)&viewProjMtx);
	glBindVertexArray(sphere->vao);
	glDrawElementsInstanced(GL_TRIANGLES, sphere->indices.size(), GL_UNSIGNED_INT, 0, particles.size());
	glBindVertexArray(0);
	glUseProgram(0);
}

/**
 * Parallel helper for clearing table 
 */
void tableClearHelper(SPHSystem& sphSystem, int start, int end) {
	for (int i = start; i < end; i++) {
		sphSystem.particleTable[i] = NULL;
	}
}

/**
 * Parallel helper for building table 
 */
void buildTableHelper(SPHSystem& sphSystem, int start, int end) {
	for (int i = start; i < end; i++) {
		Particle* pi = sphSystem.particles[i];

		// Calculate hash index using hashing formula
		uint index = getHash(sphSystem.getCell(pi));

		// Setup linked list if need be
		mtx.lock();
		if (sphSystem.particleTable[index] == NULL) {
			pi->next = NULL;
			sphSystem.particleTable[index] = pi;
		}
		else {
			pi->next = sphSystem.particleTable[index];
			sphSystem.particleTable[index] = pi;
		}
		mtx.unlock();
	}
}

void SPHSystem::buildTable() {
	// Parallel empty the table
	for (int i = 0; i < THREAD_COUNT; i++) {
		threads[i] = std::thread(tableClearHelper, std::ref(*this), tableBlockBoundaries[i], tableBlockBoundaries[i + 1]);
	}
	for (std::thread& thread : threads) {
		thread.join();
	}

	// Parallel build table
	for (int i = 0; i < THREAD_COUNT; i++) {
		threads[i] = std::thread(buildTableHelper, std::ref(*this), blockBoundaries[i], blockBoundaries[i + 1]);
	}
	for (std::thread& thread : threads) {
		thread.join();
	}
}

glm::ivec3 SPHSystem::getCell(Particle* p) const {
	return glm::ivec3(p->position.x / h, p->position.y / h, p->position.z / h);
}

void SPHSystem::print() {
	std::cout << "PARTICLE DENSITY AND PRESSURES: " << std::endl;
	for (int i = 0; i < particles.size(); i++) {
		std::cout << "PARTICLE: " << particles[i]->id << " " << std::endl;
		std::cout << "DENSITY: " << particles[i]->density << " " << std::endl;
		std::cout << "PRESSURE: " << particles[i]->pressure << " " << std::endl;
		std::cout << "NEIGHBOURS [" << neighbouringParticles[i].size() << "]: ";
		for (int j = 0; j < neighbouringParticles[i].size(); j++) {
			std::cout << neighbouringParticles[i][j]->id << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
}

void SPHSystem::reset() {
	initParticles();
	started = false;
}

void SPHSystem::startSimulation() {
	started = true;
}