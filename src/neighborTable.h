/// Methods for interacting with the particle neighbor hash table

#ifndef SPH_NEIGHBORTABLE_H
#define SPH_NEIGHBORTABLE_H

#include <Particle.h>
#include <SPHSystem.h>

// Note: Since hashes are stored as uint16, this is the max table size.
#define TABLE_SIZE 65535

// Note: Hashes are in range 0x0000 - 0xFFFF, so we allocate the next
// bit 0x10000 to indicate no particle has this hash.
const uint32_t NO_PARTICLE = 0x10000;

/// Returns a hash of the cell position
uint16_t getHash(const glm::ivec3 &cell);

/// Get the cell that the particle is in.
glm::ivec3 getCell(Particle *p, float h);

/// Creates the particle neighbor hash table.
/// It is the caller's responsibility to free the table.
uint32_t* createNeighborTable(
    Particle *sortedParticles, const size_t &particleCount);

#endif //SPH_NEIGHBORTABLE_H
