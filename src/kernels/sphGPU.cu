//#include <glm/gtx/norm.hpp>
//
//#include <kernels/sphGPU.h>
//#include <neighborTable.h>
//
///// Main kernel for running SPH and updating particle positions
//__global__ void updateParticlesKernel(
//    Particle *particles, glm::mat4 *particleTransforms,
//    const size_t particleCount, const SPHSettings &settings,
//    float deltaTime)
//{
//    size_t blockWidth = particleCount / blockDim.x + 1;
//    size_t blockOffset = threadIdx.x * blockWidth;
//    size_t maybeEnd = blockWidth + blockOffset;
//    size_t end = maybeEnd < particleCount ? maybeEnd : particleCount;
//    for (size_t i = blockOffset; i < end; i++) {
//        particleTransforms[i] *= glm::translate(glm::vec3(0, deltaTime * 10, 0));
//    }
//}
//
///// Returns a hash of the cell position
//__device__ uint getHashDevice(const glm::ivec3 &cell)
//{
//    return (
//        (uint)(cell.x * 73856093)
//        ^ (uint)(cell.y * 19349663)
//        ^ (uint)(cell.z * 83492791)
//    ) % TABLE_SIZE;
//}
//
///// Get the cell that the particle is in.
//__device__ glm::ivec3 getCellDevice(Particle *p, float h)
//{
//    return {p->position.x / h, p->position.y / h, p->position.z / h};
//}
//
///// Kernel computation function for calculating density
///// and pressures of particles in the given SPH System.
//__global__ void calculateDensitiesAndPressuresKernel(
//    Particle *particles, const size_t particleCount, const int *particleTable,
//    const SPHSettings settings)
//{
//    // TODO: Try make use of CUDA constants
//	float massPoly6Product = settings.mass * settings.poly6;
//
//    // TODO: Try make use of CUDA constants
//    size_t blockWidth = particleCount / (blockDim.x * gridDim.x) + 1;
//    size_t start = (blockIdx.x * blockDim.x + threadIdx.x) * blockWidth;
//    size_t maybeEnd = start + blockWidth;
//    size_t end = maybeEnd < particleCount ? maybeEnd : particleCount;
//
//    for (size_t i = start; i < end; i++) {
//		float pDensity = 0;
//        Particle* pi = &particles[i];
//        glm::ivec3 cell = getCellDevice(pi, settings.h);
//
//		for (int x = -1; x <= 1; x++) {
//			for (int y = -1; y <= 1; y++) {
//				for (int z = -1; z <= 1; z++) {
//					glm::ivec3 near_cell = cell + glm::ivec3(x, y, z);
//					uint index = getHashDevice(near_cell);
//
//                    // Iterate through cell linked list
//                    int pjIndex = particleTable[index];
//                    while (pjIndex != -1) {
//                        Particle *pj = &particles[pjIndex];
//                        float dist2 = glm::length2(pj->position - pi->position);
//                        if (dist2 < settings.h2 && pi != pj) {
//                            pDensity += massPoly6Product
//                                * glm::pow(settings.h2 - dist2, 3);
//                        }
//                        pjIndex = pj->next;
//                    }
//				}
//			}
//		}
//
//		// Include self density (as itself isn't included in neighbour)
//		pi->density = pDensity + settings.selfDens;
//
//		// Calculate pressure
//		float pPressure
//            = settings.gasConstant * (pi->density - settings.restDensity);
//		pi->pressure = pPressure;
//	}
//}
//
///// Parallel computation function for calculating forces
///// of particles in the given SPH System.
//__global__ void calculateForcesKernel(
//    Particle *particles, const size_t particleCount, const int *particleTable,
//    const SPHSettings settings)
//{
//	// TODO: Try make use of CUDA constants
//    size_t blockWidth = particleCount / (blockDim.x * gridDim.x) + 1;
//    size_t start = (blockIdx.x * blockDim.x + threadIdx.x) * blockWidth;
//    size_t maybeEnd = start + blockWidth;
//    size_t end = maybeEnd < particleCount ? maybeEnd : particleCount;
//
//    for (size_t i = start; i < end; i++) {
//		Particle* pi = &particles[i];
//		pi->force = glm::vec3(0);
//		glm::ivec3 cell = getCellDevice(pi, settings.h);
//
//		for (int x = -1; x <= 1; x++) {
//			for (int y = -1; y <= 1; y++) {
//				for (int z = -1; z <= 1; z++) {
//					glm::ivec3 near_cell = cell + glm::ivec3(x, y, z);
//					uint index = getHashDevice(near_cell);
//                    int pjIndex = particleTable[index];
//
//					// Iterate through cell linked list
//					while (pjIndex != -1) {
//                        Particle *pj = &particles[pjIndex];
//						float dist2 = glm::length2(pj->position - pi->position);
//						if (dist2 < settings.h2 && pi != pj) {
//							//unit direction and length
//							float dist = sqrt(dist2);
//							glm::vec3 dir = glm::normalize(pj->position - pi->position);
//
//							//apply pressure force
//							glm::vec3 pressureForce = -dir * settings.mass * (pi->pressure + pj->pressure) / (2 * pj->density) * settings.spikyGrad;
//							pressureForce *= std::pow(settings.h - dist, 2);
//							pi->force += pressureForce;
//
//							//apply viscosity force
//							glm::vec3 velocityDif = pj->velocity - pi->velocity;
//							glm::vec3 viscoForce = settings.viscosity * settings.mass * (velocityDif / pj->density) * settings.spikyLap * (settings.h - dist);
//							pi->force += viscoForce;
//						}
//                        pjIndex = pj->next;
//					}
//				}
//			}
//		}
//	}
//}
//
///// Parallel computation function moving positions
///// of particles in the given SPH System.
//__global__ void updateParticlePositionsKernel(
//    Particle *particles, const size_t particleCount,
//    glm::mat4 *particleTransforms, const SPHSettings settings,
//    const float deltaTime)
//{
//    glm::mat4 sphereScale = glm::scale(glm::vec3(settings.h / 2.f));
//    float boxWidth = 3.f;
//    float elasticity = 0.5f;
//
//    // TODO: Try make use of CUDA constants
//    size_t blockWidth = particleCount / (blockDim.x * gridDim.x) + 1;
//    size_t start = (blockIdx.x * blockDim.x + threadIdx.x) * blockWidth;
//    size_t maybeEnd = start + blockWidth;
//    size_t end = maybeEnd < particleCount ? maybeEnd : particleCount;
//
//    for (size_t i = start; i < end; i++) {
//		Particle *p = &particles[i];
//
//		//calculate acceleration and velocity
//		glm::vec3 acceleration = p->force / p->density + glm::vec3(0, settings.g, 0);
//		p->velocity += acceleration * deltaTime;
//
//		// Update position
//		p->position += p->velocity * deltaTime;
//
//		// Handle collisions with box
//		if (p->position.y < settings.h) {
//			p->position.y = -p->position.y + 2 * settings.h + 0.0001f;
//			p->velocity.y = -p->velocity.y * elasticity;
//		}
//
//		if (p->position.x < settings.h - boxWidth) {
//			p->position.x = -p->position.x + 2 * (settings.h - boxWidth) + 0.0001f;
//			p->velocity.x = -p->velocity.x * elasticity;
//		}
//
//		if (p->position.x > -settings.h + boxWidth) {
//			p->position.x = -p->position.x + 2 * -(settings.h - boxWidth) - 0.0001f;
//			p->velocity.x = -p->velocity.x * elasticity;
//		}
//
//		if (p->position.z < settings.h - boxWidth) {
//			p->position.z = -p->position.z + 2 * (settings.h - boxWidth) + 0.0001f;
//			p->velocity.z = -p->velocity.z * elasticity;
//		}
//
//		if (p->position.z > -settings.h + boxWidth) {
//			p->position.z = -p->position.z + 2 * -(settings.h - boxWidth) - 0.0001f;
//			p->velocity.z = -p->velocity.z * elasticity;
//		}
//
//        particleTransforms[i] = glm::translate(p->position) * sphereScale;
//	}
//}
//
//void updateParticlesGPU(
//    Particle *particles, glm::mat4 *particleTransforms,
//    const size_t particleCount, const SPHSettings &settings,
//    float deltaTime)
//{
//    std::cout << "Running on GPU" << std::endl;
//    // Create the neighbor table first as it adjusts some of the
//    // particle properties
//    int *neighborTable
//        = createNeighborTable(particles, particleCount, settings);
//
//    // Copy particles
//    Particle *dParticles;
//    size_t particlesSize = sizeof(Particle) * particleCount;
//    cudaMalloc((void**)&dParticles, particlesSize);
//    cudaMemcpy(dParticles, particles, particlesSize, cudaMemcpyHostToDevice);
//
//    // Copy particle transforms
//    glm::mat4 *dParticleTransforms;
//    size_t transformsSize = sizeof(glm::mat4) * particleCount;
//    cudaMalloc((void**)&dParticleTransforms, transformsSize);
//    cudaMemcpy(
//        dParticleTransforms, particleTransforms, transformsSize,
//        cudaMemcpyHostToDevice);
//
//    // Copy neighbor table
//    int *dParticleTable;
//    size_t particleTableSize = sizeof(int) * TABLE_SIZE;
//    cudaMalloc((void**)&dParticleTable, particleTableSize);
//    cudaMemcpy(
//        dParticleTable, neighborTable, particleTableSize,
//        cudaMemcpyHostToDevice);
//    free(neighborTable);
//
//    calculateDensitiesAndPressuresKernel<<<8, 512>>>(
//        dParticles, particleCount, dParticleTable, settings);
//    cudaDeviceSynchronize();
////    // Check for errors
////    cudaError_t error = cudaGetLastError();
////    if (error != cudaSuccess) {
////        printf("CUDA error: %s\n", cudaGetErrorString(error));
////    }
//
//    calculateForcesKernel<<<8, 512>>>(
//        dParticles, particleCount, dParticleTable, settings);
//    cudaDeviceSynchronize();
//
//    updateParticlePositionsKernel<<<8, 512>>>(
//        dParticles, particleCount, dParticleTransforms, settings, deltaTime);
//    cudaDeviceSynchronize();
//
//    cudaMemcpy(particles, dParticles, particlesSize, cudaMemcpyDeviceToHost);
//    cudaMemcpy(
//        particleTransforms, dParticleTransforms, transformsSize,
//        cudaMemcpyDeviceToHost);
//
//    // Free allocated memory
//    cudaFree(dParticles);
//    cudaFree(dParticleTransforms);
//    cudaFree(dParticleTable);
//}