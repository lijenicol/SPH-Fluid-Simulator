#include "SPHSystem.h"
#include <iostream>
#include <cstdlib>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#define PI 3.14159265f

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

	// Load in sphere geometry for rendering particles
	sphere = new Geometry("sphere.obj");
	sphereScale = glm::scale(glm::vec3(h/2.f));

	//start init
	started = false;
}

SPHSystem::~SPHSystem() {
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
				glm::vec3 nParticlePos = glm::vec3(i * particleSeperation + ranX - 0.25f, j * particleSeperation + ranY + h + 0.5f, k * particleSeperation + ranZ - 0.25f);

				//create new particle
				Particle* nParticle = new Particle(MASS, h,	nParticlePos, glm::vec3(0));

				//append particle
				particles[i + (j + numParticles * k) * numParticles] = nParticle;
			}
		}
	}
}

// Main update loop
// First finds neighbours of all particles, then calculates pressures
// Calculates pressure force and viscosity force
// Lastly updates particle positions
void SPHSystem::update(float deltaTime) {
	if (!started) return;
	deltaTime = 0.001f;
		
	//find the neighbours of all particles
	searchNeighbours();

	//CALCULATE DENSITY AND PRESSURES
	for (int i = 0; i < particles.size(); i++) {
		float pDensity = 0;
		Particle* pi = particles[i];

		//DENSITY USING ALL PARTICLES
		for (int j = 0; j < neighbouringParticles[i].size(); j++) {
			Particle* pj = neighbouringParticles[i][j];
			float dist2 = glm::length2(pj->position - pi->position);
			pDensity += MASS * POLY6 * glm::pow(H2 - dist2, 3);
		}
		
		// Include self density (as itself isn't included in neighbour)
		pi->density = pDensity + SELF_DENS;

		// Calculate pressure
		float pPressure = GAS_CONSTANT * (pi->density - restDensity);
		pi->pressure = pPressure;
	}

	//CALCULATE FORCES
	for (int i = 0; i < particles.size(); i++) {
		Particle* pi = particles[i];
		pi->force = glm::vec3(0);

		for (int j = 0; j < neighbouringParticles[i].size(); j++) {
			Particle* pj = neighbouringParticles[i][j];

			//unit direction and length
			float dist = glm::length(pj->position - pi->position);
			glm::vec3 dir = glm::normalize(pj->position - pi->position);

			//apply pressure force
			glm::vec3 pressureForce = -dir * MASS * (pi->pressure + pj->pressure) / (2 * pj->density) * SPIKY_GRAD;
			pressureForce *= std::pow(h - dist, 2);
			pi->force += pressureForce;

			//apply viscosity force
			glm::vec3 velocityDif = pj->velocity - pi->velocity;
			glm::vec3 viscoForce = viscosity * MASS * (velocityDif / pj->density) * SPIKY_LAP * (h - dist);
			pi->force += viscoForce;
		}
	}

	//update particle positions
	updateParticles(deltaTime);
}

void SPHSystem::updateParticles(float deltaTime) {
	for (int i = 0; i < particles.size(); i++) {
		Particle *p = particles[i];

		//calculate acceleration and velocity
		glm::vec3 acceleration = p->force / p->density + glm::vec3(0,g,0);
		p->velocity += acceleration * deltaTime;
		
		// Update position
		p->position += p->velocity * deltaTime;

		// Handle collisions with box
		float boxWidth = 0.5f;
		float elasticity = 0.15f;
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

void SPHSystem::draw(glm::mat4 viewProjMtx, GLuint shader) {
	//draw each of the particles
	for (int i = 0; i < particles.size(); i++) {
		//calculate transformation matrix for sphere
		glm::mat4 translate = glm::translate(particles[i]->position);
		sphere->draw(translate * sphereScale, viewProjMtx, shader);
	}
}

void SPHSystem::searchNeighbours() {
	for (int i = 0; i < particles.size(); i++) {
		std::vector<Particle*> nbParticles;

		for (int j = 0; j < particles.size(); j++) {
			//dont check the particle against itself
			if (i == j) continue;

			//find distance
			glm::vec3 s1 = particles[i]->position;
			glm::vec3 s2 = particles[j]->position;
			float distance = glm::length(s2 - s1);
			
			//append if distance is less than h
			if (distance <= h) {
				nbParticles.push_back(particles[j]);
			}
		}
		neighbouringParticles[i] = nbParticles;
	}
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