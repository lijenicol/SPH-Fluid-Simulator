//
// Created by lijenicol on 16/09/23.
//

#include <algorithm>
#include <numeric>

#include <glm/gtx/norm.hpp>

#include <neighborTable.h>
#include <sph.h>
#include <kernels/sphGPU.h>
#include <timer.h>


/// Calculates and stores particle hashes.
void parallelCalculateHashes(
    Particle *particles, size_t start, size_t end, const SPHSettings &settings)
{
    for (size_t i = start; i < end; i++) {
        Particle *particle = &particles[i];
        particle->hash = getHash(getCell(particle, settings.h));
    }
}

/// Parallel computation function for calculating density
/// and pressures of particles in the given SPH System.
void parallelDensityAndPressures(
    Particle *particles, const size_t particleCount, const size_t start,
    const size_t end, const uint32_t *particleTable,
    const SPHSettings &settings)
{
	float massPoly6Product = settings.mass * settings.poly6;

	for (size_t piIndex = start; piIndex < end; piIndex++) {
		float pDensity = 0;
		Particle* pi = &particles[piIndex];
		glm::ivec3 cell = getCell(pi, settings.h);

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					uint16_t cellHash = getHash(cell + glm::ivec3(x, y, z));
                    uint32_t pjIndex = particleTable[cellHash];
                    if (pjIndex == NO_PARTICLE) {
                        continue;
                    }
					while (pjIndex < particleCount) {
                        if (pjIndex == piIndex) {
                            pjIndex++;
                            continue;
                        }
                        Particle *pj = &particles[pjIndex];
                        if (pj->hash != cellHash) {
                            break;
                        }
						float dist2 = glm::length2(pj->position - pi->position);
						if (dist2 < settings.h2) {
							pDensity += massPoly6Product
                                * glm::pow(settings.h2 - dist2, 3);
						}
                        pjIndex++;
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
    Particle *particles, const size_t particleCount, const size_t start,
    const size_t end, const uint32_t *particleTable,
    const SPHSettings &settings)
{
	for (size_t piIndex = start; piIndex < end; piIndex++) {
		Particle* pi = &particles[piIndex];
		pi->force = glm::vec3(0);
		glm::ivec3 cell = getCell(pi, settings.h);

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
                    uint16_t cellHash = getHash(cell + glm::ivec3(x, y, z));
                    uint32_t pjIndex = particleTable[cellHash];
                    if (pjIndex == NO_PARTICLE) {
                        continue;
                    }
                    while (pjIndex < particleCount) {
                        if (pjIndex == piIndex) {
                            pjIndex++;
                            continue;
                        }
                        Particle *pj = &particles[pjIndex];
                        if (pj->hash != cellHash) {
                            break;
                        }
						float dist2 = glm::length2(pj->position - pi->position);
						if (dist2 < settings.h2) {
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
                        pjIndex++;
					}
				}
			}
		}
	}
}

/// Parallel computation function moving positions
/// of particles in the given SPH System.
void parallelUpdateParticlePositions(
    Particle *particles, const size_t particleCount, const size_t start,
    const size_t end, glm::mat4 *particleTransforms,
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

        particleTransforms[p->id]
            = glm::translate(p->position) * sphereScale;
	}
}

/// Sort particles by the particle's hash
void sortParticles(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t &particleCount)
{
    std::vector<size_t> indices(particleCount);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(
        indices.begin(), indices.end(),
        [&](size_t i, size_t j) {
            return particles[i].hash < particles[j].hash;
        }
    );

    std::vector<bool> done(particleCount);
    for (size_t i = 0; i < particleCount; ++i) {
        if (done[i]) {
            continue;
        }
        done[i] = true;
        size_t prev_j = i;
        size_t j = indices[i];
        while (i != j) {
            std::swap(particles[prev_j], particles[j]);
            std::swap(particleTransforms[prev_j], particleTransforms[j]);
            done[j] = true;
            prev_j = j;
            j = indices[j];
        }
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

    // Calculate hashes
    {
//        Timer timer("hashes");
        for (int i = 0; i < threadCount; i++) {
            threads[i] = std::thread(
                parallelCalculateHashes, particles, blockBoundaries[i],
                blockBoundaries[i + 1], settings);
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
    }

    // Sort particles and transforms
    {
//        Timer timer("sort");
        sortParticles(particles, particleTransforms, particleCount);
    }

    uint32_t *particleTable
        = createNeighborTable(particles, particleCount);

    // Calculate densities and pressures
    {
//        Timer timer("densities");
        for (int i = 0; i < threadCount; i++) {
            threads[i] = std::thread(
                parallelDensityAndPressures, particles, particleCount,
                blockBoundaries[i], blockBoundaries[i + 1], particleTable,
                settings);
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
    }

    // Calculate forces
    {
//        Timer timer("forces");
        for (int i = 0; i < threadCount; i++) {
            threads[i] = std::thread(
                parallelForces, particles, particleCount, blockBoundaries[i],
                blockBoundaries[i + 1], particleTable, settings);
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
    }

    // Update particle positions
    {
//        Timer timer("positions");
        for (int i = 0; i < threadCount; i++) {
            threads[i] = std::thread(
                parallelUpdateParticlePositions, particles, particleCount,
                blockBoundaries[i], blockBoundaries[i + 1], particleTransforms,
                settings, deltaTime);
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
    }

    delete(particleTable);
}

void updateParticles(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime, const bool onGPU)
{
    if (onGPU) {
//        updateParticlesGPU(
//            particles, particleTransforms, particleCount, settings, deltaTime);
    }
    else {
        updateParticlesCPU(
            particles, particleTransforms, particleCount, settings, deltaTime);
    }
}