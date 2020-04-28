#include "SPHSystem.h"
#include <iostream>
#include <cstdlib>
#include <glm/gtx/string_cast.hpp>

SPHSystem::SPHSystem(unsigned int numParticles, float mass, float restDensity, float viscosity, float h, float g, float tension) {
	this->numParticles = numParticles;
	this->restDensity = restDensity;
	this->viscosity = viscosity;
	this->h = h;
	this->g = g;
	this->tension = tension;

	POLY6 = 315.0f / (64.0f * glm::pi<float>() * pow(h, 9));
	SPIKY_GRAD = -45.0f / (glm::pi<float>() * pow(h, 6));
	SPIKY_LAP = 45.0f / (glm::pi<float>() * pow(h, 6));
	MASS = mass;
	GAS_CONSTANT = 2000.f;

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
	float particleSeperation = h;
	for (int i = 0; i < numParticles; i++) {
		for (int j = 0; j < numParticles; j++) {
			for (int k = 0; k < numParticles; k++) {
				// dam like particle positions
				float ranX = (float(rand()) / float((RAND_MAX)) * 2.f - 1) * h / 10;
				float ranY = (float(rand()) / float((RAND_MAX)) * 2.f - 1) * h / 10;
				float ranZ = (float(rand()) / float((RAND_MAX)) * 2.f - 1) * h / 10;
				glm::vec3 nParticlePos = glm::vec3(i * h + ranX, j * h + ranY + h, k * h + ranZ);

				//create new particle
				Particle* nParticle = new Particle(MASS, h / 2.f, 0.04f, 0.0f,
					nParticlePos, glm::vec3(0), false);

				//append particle
				particles[i + (j + numParticles * k) * numParticles] = nParticle;
			}
		}
	}
}

void SPHSystem::update(float deltaTime) {
	if (!started) return;

	//find the neighbours of all particles
	searchNeighbours();

	//EDIT THIS OUT AFTERWARDS
	deltaTime = 0.003f;

	//CALCULATE DENSITY AND PRESSURES
	for (int i = 0; i < particles.size(); i++) {
		float pDensity = 0;
		Particle* pi = particles[i];

		//DENSITY USING ALL PARTICLES
		for (int j = 0; j < neighbouringParticles[i].size(); j++) {
			Particle* pj = neighbouringParticles[i][j];
			float dist = glm::length(pj->getPosition() - pi->getPosition());
			pDensity += pj->getMass() * POLY6 * glm::pow((h * h) - (dist * dist), 3);
		}
		pi->density = pDensity + pi->getMass() * POLY6 * glm::pow(h, 6);

		//calculate pressure
		float pPressure = GAS_CONSTANT * (glm::pow((pi->density / restDensity), 7) - 1.f);
		pi->pressure = pPressure; 
	}

	//CALCULATE FORCES
	for (int i = 0; i < particles.size(); i++) {
		Particle* pi = particles[i];
		float vi = pi->getMass() / pi->density;	// Volume of particle i
		for (int j = 0; j < neighbouringParticles[i].size(); j++) {
			Particle* pj = neighbouringParticles[i][j];

			//unit direction and length
			float dist = glm::length(pj->getPosition() - pi->getPosition());
			glm::vec3 dir = glm::normalize(pj->getPosition() - pi->getPosition());

			//calculate volumes of pi and pj
			float vj = pj->getMass() / pj->density;

			//apply pressure force
			glm::vec3 pressureForce = -dir * (vi * vj * ((pi->pressure + pj->pressure)) / 2.0f) * SPIKY_GRAD * glm::pow(h - dist, 2);
			pi->applyForce(pressureForce);

			//apply viscosity force
			glm::vec3 velocityDif = pj->getVelocity() - pi->getVelocity();
			glm::vec3 viscoForce = dir * (vi * vj * viscosity * velocityDif) * SPIKY_LAP * (h - dist);
			pi->applyForce(viscoForce);
		}

		//apply gravity force
		pi->applyForce(glm::vec3(0, g * pi->getMass(), 0));
	}

	//update particle positions
	for (int i = 0; i < particles.size(); i++) {
		particles[i]->update(deltaTime);
	}
}

void SPHSystem::draw(glm::mat4 viewProjMtx, GLuint shader) {
	//draw each of the particles
	for (int i = 0; i < particles.size(); i++) {
		//calculate transformation matrix for sphere
		glm::mat4 translate = glm::translate(particles[i]->getPosition());
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
			glm::vec3 s1 = particles[i]->getPosition();
			glm::vec3 s2 = particles[j]->getPosition();
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