//
// Created by lijenicol on 16/09/23.
//

#ifndef SPH_SPH_H
#define SPH_SPH_H

#include "Particle.h"
#include "SPHSystem.h"

#define TABLE_SIZE 1000000

/// Update attrs of particles in place.
void updateParticles(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime, const bool onGPU);

#endif //SPH_SPH_H
