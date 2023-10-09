/// Methods for interacting with the particle neighbor hash table

#ifndef SPH_NEIGHBORTABLE_H
#define SPH_NEIGHBORTABLE_H

#include <Particle.h>
#include <SPHSystem.h>

const uint32_t TABLE_SIZE = 262144;
const uint32_t NO_PARTICLE = 0xFFFFFFFF;

/// Returns a hash of the cell position
uint32_t getHash(const glm::ivec3 &cell);

/// Get the cell that the particle is in.
glm::ivec3 getCell(Particle *p, float h);

/// Creates the particle neighbor hash table.
/// It is the caller's responsibility to free the table.
uint32_t* createNeighborTable(
    Particle *sortedParticles, const size_t &particleCount);

#endif //SPH_NEIGHBORTABLE_H
