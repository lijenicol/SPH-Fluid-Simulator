#include <glm/gtx/norm.hpp>

#include <kernels/sphGPU.h>
#include <neighborTable.h>
#include <sph.h>
#include <timer.h>

/// Returns a hash of the cell position
__device__ uint16_t getHashDevice(const glm::ivec3 &cell)
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

///// Calculate particle hashes
//__global__ void calculateHashesKernel(
//    Particle *particles, const size_t particleCount, float h)
//{
//    size_t particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
//    if (particleIndex >= particleCount) {
//        return;
//    }
//    Particle *particle = &particles[particleIndex];
//    particles->hash = getHashDevice(getCellDevice(particle, h));
//}
//
//struct HashComp
//{
//    __host__ __device__ bool operator()(
//            const Particle& p1, const Particle& p2) {
//        return p1.hash < p2.hash;
//    }
//};

///// Sort particles by hash
//__global__ void sortParticles(Particle *particles, const size_t particleCount)
//{
//    thrust::sort(particles, particles + particleCount, HashComp());
//}
//
///// Constructs the neighbor table and stores the result in `createdTable`.
//__global__ void constructNeighborTable(
//    Particle *sortedParticles, const size_t particleCount,
//    uint32_t **createdTable)
//{
//    size_t globalThreadIdx = blockIdx.x * blockDim.x + threadIdx.x;
//    if (globalThreadIdx >= 1) {
//        // This method should only be run on one thread. This is here
//        // for safety.
//        return;
//    }
//
//    uint32_t *particleTable
//        = (uint32_t *)malloc(sizeof(uint32_t) * TABLE_SIZE);
//    for (size_t i = 0; i < TABLE_SIZE; ++i) {
//        particleTable[i] = NO_PARTICLE;
//    }
//
//    uint32_t prevHash = NO_PARTICLE;
//    for (size_t i = 0; i < particleCount; ++i) {
//        uint16_t currentHash = sortedParticles[i].hash;
//        if (currentHash != prevHash) {
//            particleTable[currentHash] = i;
//            prevHash = currentHash;
//        }
//    }
//
//    createdTable = &particleTable;
//}

/// Kernel computation function for calculating density
/// and pressures of particles in the given SPH System.
__global__ void calculateDensitiesAndPressuresKernel(
    Particle *particles, const size_t particleCount,
    const uint32_t *particleTable, const SPHSettings settings)
{
    size_t piIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (piIndex > particleCount) {
        return;
    }
    Particle *pi = &particles[piIndex];

    // TODO: Try make use of CUDA constants
    float massPoly6Product = settings.mass * settings.poly6;
    glm::ivec3 cell = getCellDevice(pi, settings.h);

    float pDensity = 0;

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            for (int z = -1; z <= 1; z++) {
                uint16_t cellHash
                    = getHashDevice(cell + glm::ivec3(x, y, z));
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
                    if (dist2 < settings.h2 && pi != pj) {
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

/// Parallel computation function for calculating forces
/// of particles in the given SPH System.
__global__ void calculateForcesKernel(
    Particle *particles, const size_t particleCount,
    const uint32_t *particleTable, const SPHSettings settings)
{
    size_t piIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (piIndex > particleCount) {
        return;
    }
    Particle *pi = &particles[piIndex];

    // Another constant
    glm::ivec3 cell = getCellDevice(pi, settings.h);

    pi->force = glm::vec3(0);

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            for (int z = -1; z <= 1; z++) {
                uint16_t cellHash
                    = getHashDevice(cell + glm::ivec3(x, y, z));
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
                    pjIndex++;
                }
            }
        }
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

    // TODO: These should be constants somewhere else.
    glm::mat4 sphereScale = glm::scale(glm::vec3(settings.h / 2.f));
    float boxWidth = 3.f;
    float elasticity = 0.5f;

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

    particleTransforms[pIndex] = glm::translate(p->position) * sphereScale;
}

void updateParticlesGPU(
    Particle *particles, glm::mat4 *particleTransforms,
    const size_t particleCount, const SPHSettings &settings,
    float deltaTime)
{
//    std::cout << "Running on GPU" << std::endl;

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
        Timer timer("hashes");
        for (int i = 0; i < threadCount; i++) {
            threads[i] = std::thread(
                parallelCalculateHashes, particles, blockBoundaries[i],
                blockBoundaries[i + 1], settings);
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
    }

    // Sort particles
    {
        Timer timer("sort");
        sortParticles(particles, particleCount);
    }

    // Copy particles
    Particle *dParticles;
    size_t particlesSize = sizeof(Particle) * particleCount;
    cudaMalloc((void**)&dParticles, particlesSize);
    cudaMemcpy(dParticles, particles, particlesSize, cudaMemcpyHostToDevice);

    // Copy particle transforms
    glm::mat4 *dParticleTransforms;
    size_t transformsSize = sizeof(glm::mat4) * particleCount;
    cudaMalloc((void**)&dParticleTransforms, transformsSize);
    cudaMemcpy(
        dParticleTransforms, particleTransforms, transformsSize,
        cudaMemcpyHostToDevice);

    // Create and copy particle table
    uint32_t *dParticleTable;
    {
        Timer timer("tableCreation");
        uint32_t *particleTable
            = createNeighborTable(particles, particleCount);
        size_t tableSize = sizeof(uint32_t) * TABLE_SIZE;
        cudaMalloc((void**)&dParticleTable, tableSize);
        cudaMemcpy(
            dParticleTable, particleTable, tableSize, cudaMemcpyHostToDevice);
        // Can free the host table since it won't be of use.
        free(particleTable);
    }

//    const size_t blockSize = 512;
//    size_t gridSize = blockSize / particleCount + 1;
//    calculateHashesKernel<<<gridSize, blockSize>>>(
//        dParticles, particleCount, settings.h);
//    cudaDeviceSynchronize();
//
//    uint32_t *dParticleTable;
//    constructNeighborTable<<<1, 1>>>(
//        particles, particleCount, &dParticleTable);
//    cudaDeviceSynchronize();

    size_t threadsPerBlock = 512;
    size_t gridSize = particleCount / threadsPerBlock + 1;

    {
        Timer timer("densities");
        calculateDensitiesAndPressuresKernel<<<gridSize, threadsPerBlock>>>(
            dParticles, particleCount, dParticleTable, settings);
        cudaDeviceSynchronize();
    }
//    // Check for errors
//    cudaError_t error = cudaGetLastError();
//    if (error != cudaSuccess) {
//        printf("CUDA error: %s\n", cudaGetErrorString(error));
//    }

    {
        Timer timer("forces");
        calculateForcesKernel<<<gridSize, threadsPerBlock>>>(
            dParticles, particleCount, dParticleTable, settings);
        cudaDeviceSynchronize();
    }

    {
        Timer timer("positions");
        updateParticlePositionsKernel<<<gridSize, threadsPerBlock>>>(
            dParticles, particleCount, dParticleTransforms, settings, deltaTime);
        cudaDeviceSynchronize();
    }

    cudaMemcpy(particles, dParticles, particlesSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(
        particleTransforms, dParticleTransforms, transformsSize,
        cudaMemcpyDeviceToHost);

    // Free allocated memory
    cudaFree(dParticles);
    cudaFree(dParticleTransforms);
    cudaFree(dParticleTable);
}