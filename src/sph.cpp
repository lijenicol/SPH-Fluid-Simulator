//
// Created by lijenicol on 16/09/23.
//

#include <glm/gtx/norm.hpp>

#include <neighborTable.h>
#include <sph.h>
#include <kernels/sphGPU.h>

/// Parallel computation function for calculating density
/// and pressures of particles in the given SPH System.
void parallelDensityAndPressures(
    Particle *particles, int start, int end, const int *particleTable,
    const SPHSettings &settings)
{
	float massPoly6Product = settings.mass * settings.poly6;

	for (int i = start; i < end; i++) {
		float pDensity = 0;
		Particle* pi = &particles[i];
		glm::ivec3 cell = getCell(pi, settings.h);

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					glm::ivec3 near_cell = cell + glm::ivec3(x, y, z);
					uint index = getHash(near_cell);
                    int pjIndex = particleTable[index];

					// Iterate through cell linked list
					while (pjIndex != -1) {
                        Particle *pj = &particles[pjIndex];
						float dist2 = glm::length2(pj->position - pi->position);
						if (dist2 < settings.h2 && pi != pj) {
							pDensity += massPoly6Product
                                * glm::pow(settings.h2 - dist2, 3);
						}
						pjIndex = pj->next;
					}
				}
			}
		}

		// Include self density (as itself isn't included in neighbour)
		pi->density = pDensity + settings.selfDens;

		// Calculate pressure
		float pPressure
            = settings.gasConstant * (pi->density - settings.restDensity);
		pi->pressure = pPressure;
	}
}

/// Parallel computation function for calculating forces
/// of particles in the given SPH System.
void parallelForces(
    Particle *particles, int start, int end, const int *particleTable,
    const SPHSettings &settings)
{
	for (int i = start; i < end; i++) {
		Particle* pi = &particles[i];
		pi->force = glm::vec3(0);
		glm::ivec3 cell = getCell(pi, settings.h);

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					glm::ivec3 near_cell = cell + glm::ivec3(x, y, z);
					uint index = getHash(near_cell);
                    int pjIndex = particleTable[index];

					// Iterate through cell linked list
					while (pjIndex != -1) {
                        Particle *pj = &particles[pjIndex];
						float dist2 = glm::length2(pj->position - pi->position);
						if (dist2 < settings.h2 && pi != pj) {
							//unit direction and length
							float dist = sqrt(dist2);
							glm::vec3 dir = glm::normalize(pj->position - pi->position);

							//apply pressure force
							glm::vec3 pressureForce = -dir * settings.mass * (pi->pressure + pj->pressure) / (2 * pj->density) * settings.spikyGrad;
							pressureForce *= std::pow(settings.h - dist, 2);
							pi->force += pressureForce;

							//apply viscosity force
							glm::vec3 velocityDif = pj->velocity - pi->velocity;
							glm::vec3 viscoForce = settings.viscosity * settings.mass * (velocityDif / pj->density) * settings.spikyLap * (settings.h - dist);
							pi->force += viscoForce;
						}
                        pjIndex = pj->next;
					}
				}
			}
		}
	}
}

/// Parallel computation function moving positions
/// of particles in the given SPH System.
void parallelUpdateParticlePositions(
    Particle *particles, int start, int end, glm::mat4 *particleTransforms,
    const SPHSettings &settings, const float &deltaTime)
{
    glm::mat4 sphereScale = glm::scale(glm::vec3(settings.h / 2.f));
    float boxWidth = 3.f;
    float elasticity = 0.5f;

	for (size_t i = start; i < end; i++) {
		Particle *p = &particles[i];

		//calculate acceleration and velocity
		glm::vec3 acceleration = p->force / p->density + glm::vec3(0, settings.g, 0);
		p->velocity += acceleration * deltaTime;

		// Update position
		p->position += p->velocity * deltaTime;

		// Handle collisions with box
		if (p->position.y < settings.h) {
			p->position.y = -p->position.y + 2 * settings.h + 0.0001f;
			p->velocity.y = -p->velocity.y * elasticity;
		}

		if (p->position.x < settings.h - boxWidth) {
			p->position.x = -p->position.x + 2 * (settings.h - boxWidth) + 0.0001f;
			p->velocity.x = -p->velocity.x * elasticity;
		}

		if (p->position.x > -settings.h + boxWidth) {
			p->position.x = -p->position.x + 2 * -(settings.h - boxWidth) - 0.0001f;
			p->velocity.x = -p->velocity.x * elasticity;
		}

		if (p->position.z < settings.h - boxWidth) {
			p->position.z = -p->position.z + 2 * (settings.h - boxWidth) + 0.0001f;
			p->velocity.z = -p->velocity.z * elasticity;
		}

		if (p->position.z > -settings.h + boxWidth) {
			p->position.z = -p->position.z + 2 * -(settings.h - boxWidth) - 0.0001f;
			p->velocity.z = -p->velocity.z * elasticity;
		}

        particleTransforms[i] = glm::translate(p->position) * sphereScale;
	}
}

/// CPU update particles implementation
void updateParticlesCPU(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime)
{
    const size_t threadCount = std::thread::hardware_concurrency();
    std::thread threads[threadCount];

    size_t blockBoundaries[threadCount + 1];
    blockBoundaries[0] = 0;
    size_t blockSize = particleCount / threadCount;
    for (size_t i = 1; i < threadCount; i++) {
        blockBoundaries[i] = i * blockSize;
    }
    blockBoundaries[threadCount] = particleCount;

    int *particleTable
        = createNeighborTable(particles, particleCount, settings);

    // Calculate densities and pressures
    for (int i = 0; i < threadCount; i++) {
        threads[i] = std::thread(
            parallelDensityAndPressures, particles, blockBoundaries[i],
            blockBoundaries[i + 1], particleTable, settings);
    }
    for (std::thread& thread : threads) {
        thread.join();
    }

    // Calculate forces
    for (int i = 0; i < threadCount; i++) {
        threads[i] = std::thread(
            parallelForces, particles, blockBoundaries[i],
            blockBoundaries[i + 1], particleTable, settings);
    }
    for (std::thread& thread : threads) {
        thread.join();
    }

    // Update particle positions
    for (int i = 0; i < threadCount; i++) {
        threads[i] = std::thread(
            parallelUpdateParticlePositions, particles, blockBoundaries[i],
            blockBoundaries[i + 1], particleTransforms, settings, deltaTime);
    }
    for (std::thread& thread : threads) {
        thread.join();
    }

    delete(particleTable);
}

void updateParticles(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime, const bool onGPU)
{
    if (onGPU) {
        updateParticlesGPU(
            particles, particleTransforms, particleCount, settings, deltaTime);
    }
    else {
        updateParticlesCPU(
            particles, particleTransforms, particleCount, settings, deltaTime);
    }
}