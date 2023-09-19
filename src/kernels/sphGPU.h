#ifndef SPH_SPHGPU_H
#define SPH_SPHGPU_H

#include <Particle.h>
#include <SPHSystem.h>

/// Update attrs of particles in place.
void updateParticlesGPU(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime);

#endif // SPH_SPHGPU_H