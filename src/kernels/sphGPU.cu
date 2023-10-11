#include <glm/gtx/norm.hpp>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include <neighborTable.h>
#include <kernels/sphGPU.h>

static const float BOX_WIDTH = 10.f;
static const float BOX_COLLISION_ELASTICITY = 1.f;
static const float BOX_COLLISION_OFFSET = 0.00001;

static const uint16_t MAX_NEIGHBORS = 32;

/// Returns a hash of the cell position
__device__ uint32_t getHashDevice(const glm::ivec3 &cell)
{
    return (
        (uint)(cell.x * 73856093)
        ^ (uint)(cell.y * 19349663)
        ^ (uint)(cell.z * 83492791)
    ) % TABLE_SIZE;
}

/// Get the cell that the particle is in.
__device__ glm::ivec3 getCellDevice(Particle *p, float h)
{
    return {p->position.x / h, p->position.y / h, p->position.z / h};
}

/// Calculate particle hashes
__global__ void calculateHashesKernel(
    Particle *particles, const size_t particleCount, float h)
{
    size_t particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (particleIndex >= particleCount) {
        return;
    }
    Particle *particle = &particles[particleIndex];
    particle->hash = getHashDevice(getCellDevice(particle, h));
}

/// Comparison struct for sorting particles by hash.
struct HashComp
{
    __host__ __device__ bool operator()(
            const Particle& p1, const Particle& p2) {
        return p1.hash < p2.hash;
    }
};

/// Constructs the hash to particle map, storing the result
/// in `hashToParticleMap`.
__global__ void constructHashToParticleMap(
    Particle *sortedParticles, const size_t particleCount,
    uint32_t *hashToParticleMap)
{
    size_t particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (particleIndex >= particleCount) {
        return;
    }

    uint32_t prevHash = particleIndex == 0
        ? NO_PARTICLE : sortedParticles[particleIndex-1].hash;
    uint32_t currentHash = sortedParticles[particleIndex].hash;
    if (currentHash != prevHash) {
        hashToParticleMap[currentHash] = particleIndex;
    }
}

/// Calculates neighbors and stores the result in `neighborList`
__global__ void initNeighborList(
    Particle *particles, const size_t particleCount,
    const uint32_t *hashToParticleMap, const SPHSettings settings,
    Particle **neighborList)
{
    size_t piIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (piIndex > particleCount) {
        return;
    }
    Particle *pi = &particles[piIndex];
    glm::ivec3 cell = getCellDevice(pi, settings.h);
    size_t neighborCount = 0;
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            for (int z = -1; z <= 1; z++) {
                uint16_t cellHash
                    = getHashDevice(cell + glm::ivec3(x, y, z));
                uint32_t pjIndex = hashToParticleMap[cellHash];
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
                    if (dist2 < settings.h2 && pi != pj) {
                        neighborList[piIndex * MAX_NEIGHBORS + neighborCount]
                            = pj;
                        neighborCount++;
                    }
                    pjIndex++;
                }
            }
        }
    }

    // Fill the rest of the table with null pointers
    for (size_t i = neighborCount; i < MAX_NEIGHBORS; i++) {
        neighborList[piIndex * MAX_NEIGHBORS + i] = nullptr;
    }
}

/// Kernel computation function for calculating density
/// and pressures of particles in the given SPH System.
__global__ void calculateDensitiesAndPressuresKernel(
    Particle *particles, const size_t particleCount,
    Particle **neighborList, const SPHSettings settings)
{
    size_t piIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (piIndex > particleCount) {
        return;
    }
    Particle *pi = &particles[piIndex];
    pi->density = settings.selfDens;
    glm::ivec3 cell = getCellDevice(pi, settings.h);

    size_t neighborOffset = piIndex * MAX_NEIGHBORS;
    const Particle *pj = neighborList[neighborOffset];
    while (pj) {
        float dist2 = glm::length2(pj->position - pi->position);
        pi->density += settings.massPoly6Product
                    * glm::pow(settings.h2 - dist2, 3);
        neighborOffset++;
        pj = neighborList[neighborOffset];
    }

    pi->pressure = settings.gasConstant * (pi->density - settings.restDensity);
}

/// Parallel computation function for calculating forces
/// of particles in the given SPH System.
__global__ void calculateForcesKernel(
    Particle *particles, const size_t particleCount,
    Particle **neighborList, const SPHSettings settings)
{
    size_t piIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (piIndex > particleCount) {
        return;
    }
    Particle *pi = &particles[piIndex];
    pi->force = glm::vec3(0);
    glm::ivec3 cell = getCellDevice(pi, settings.h);

    size_t neighborOffset = piIndex * MAX_NEIGHBORS;
    const Particle *pj = neighborList[neighborOffset];
    while (pj) {
        // Unit direction and length
        float dist = glm::length(pj->position - pi->position);
        glm::vec3 dir = glm::normalize(pj->position - pi->position);

        // Apply pressure force
        glm::vec3 pressureForce
            = -dir * settings.mass
              * (pi->pressure + pj->pressure) / (2 * pj->density)
              * settings.spikyGrad;
        pressureForce *= std::pow(settings.h - dist, 2);
        pi->force += pressureForce;

        // Apply viscosity force
        glm::vec3 velocityDif = pj->velocity - pi->velocity;
        glm::vec3 viscoForce
            = settings.viscosity * settings.mass
              * (velocityDif / pj->density) * settings.spikyLap
              * (settings.h - dist);
        pi->force += viscoForce;

        neighborOffset++;
        pj = neighborList[neighborOffset];
    }
}

/// Parallel computation function moving positions
/// of particles in the given SPH System.
__global__ void updateParticlePositionsKernel(
    Particle *particles, const size_t particleCount,
    glm::mat4 *particleTransforms, const SPHSettings settings,
    const float deltaTime)
{
    size_t pIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (pIndex > particleCount) {
        return;
    }
    Particle *p = &particles[pIndex];

    //calculate acceleration and velocity
    glm::vec3 acceleration
        = p->force / p->density + glm::vec3(0, settings.g, 0);
    p->velocity += acceleration * deltaTime;

    // Update position
    p->position += p->velocity * deltaTime;

    // Handle collisions with box
    if (p->position.y < settings.h) {
        p->position.y = -p->position.y + 2 * settings.h
                        + BOX_COLLISION_OFFSET;
        p->velocity.y = -p->velocity.y * BOX_COLLISION_ELASTICITY;
    }

    if (p->position.x < settings.h - BOX_WIDTH) {
        p->position.x = -p->position.x + 2 * (settings.h - BOX_WIDTH)
                        + BOX_COLLISION_OFFSET;
        p->velocity.x = -p->velocity.x * BOX_COLLISION_ELASTICITY;
    }

    if (p->position.x > -settings.h + BOX_WIDTH) {
        p->position.x = -p->position.x + 2 * -(settings.h - BOX_WIDTH)
                        - BOX_COLLISION_OFFSET;
        p->velocity.x = -p->velocity.x * BOX_COLLISION_ELASTICITY;
    }

    if (p->position.z < settings.h - BOX_WIDTH) {
        p->position.z = -p->position.z + 2 * (settings.h - BOX_WIDTH)
                        + BOX_COLLISION_OFFSET;
        p->velocity.z = -p->velocity.z * BOX_COLLISION_ELASTICITY;
    }

    if (p->position.z > -settings.h + BOX_WIDTH) {
        p->position.z = -p->position.z + 2 * -(settings.h - BOX_WIDTH)
                        - BOX_COLLISION_OFFSET;
        p->velocity.z = -p->velocity.z * BOX_COLLISION_ELASTICITY;
    }

    particleTransforms[pIndex]
        = glm::translate(p->position) * settings.sphereScale;
}

void updateParticlesGPU(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime)
{
    size_t threadsPerBlock = 512;
    size_t gridSize = particleCount / threadsPerBlock + 1;

    // Start by copying particles to GPU.
    thrust::device_vector<Particle> dParticleVector(
        particles, particles + particleCount);
    Particle *dParticles = thrust::raw_pointer_cast(dParticleVector.data());

    // Set up hashes for all particles
    calculateHashesKernel<<<gridSize, threadsPerBlock>>>(
        dParticles, particleCount, settings.h);
    cudaDeviceSynchronize();

    // Sort particles by hash.
    thrust::sort(
        dParticleVector.begin(), dParticleVector.end(), HashComp());

    // Set up the hash->particle map.
    thrust::device_vector<uint32_t> dHashToParticleMapVector(
        TABLE_SIZE, NO_PARTICLE);
    uint32_t *dHashToParticleMap
        = thrust::raw_pointer_cast(dHashToParticleMapVector.data());
    constructHashToParticleMap<<<gridSize, threadsPerBlock>>>(
        dParticles, particleCount, dHashToParticleMap);
    cudaDeviceSynchronize();

    // Calculate neighbors for every particle.
    Particle **dNeighborList;
    cudaMalloc((void**)&dNeighborList,
               sizeof(Particle *) * particleCount * MAX_NEIGHBORS);
    initNeighborList<<<gridSize, threadsPerBlock>>>(
        dParticles, particleCount, dHashToParticleMap, settings,
        dNeighborList);
    cudaDeviceSynchronize();

    // Calculate densities and pressures.
    calculateDensitiesAndPressuresKernel<<<gridSize, threadsPerBlock>>>(
        dParticles, particleCount, dNeighborList, settings);
    cudaDeviceSynchronize();

    // Calculate forces.
    calculateForcesKernel<<<gridSize, threadsPerBlock>>>(
        dParticles, particleCount, dNeighborList, settings);
    cudaDeviceSynchronize();

    cudaFree(dNeighborList);

    // Update positions and transforms.
    glm::mat4 *dParticleTransforms;
    size_t transformsSize = sizeof(glm::mat4) * particleCount;
    cudaMalloc((void**)&dParticleTransforms, transformsSize);
    updateParticlePositionsKernel<<<gridSize, threadsPerBlock>>>(
        dParticles, particleCount, dParticleTransforms, settings, deltaTime);
    cudaDeviceSynchronize();

    cudaMemcpy(particles, dParticles, sizeof(Particle) * particleCount,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(particleTransforms, dParticleTransforms, transformsSize,
               cudaMemcpyDeviceToHost);
    cudaFree(dParticleTransforms);
}