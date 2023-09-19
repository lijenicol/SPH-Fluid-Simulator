/// Methods for interacting with the particle neighbor hash table

#ifndef SPH_NEIGHBORTABLE_H
#define SPH_NEIGHBORTABLE_H

#include <Particle.h>
#include <SPHSystem.h>

#define TABLE_SIZE 1000000

/// Returns a hash of the cell position
uint getHash(const glm::ivec3 &cell);

/// Get the cell that the particle is in.
glm::ivec3 getCell(Particle *p, float h);

/// Creates the particle neighbor hash table.
/// It is the caller's responsibility to free the table.
int* createNeighborTable(
    Particle *particles, const size_t &particleCount,
    const SPHSettings &settings);

#endif //SPH_NEIGHBORTABLE_H
