#include <cstring>

#include "SPHSystem.h"
#include "sph.h"

#define PI 3.14159265f

SPHSettings::SPHSettings(
    float mass, float restDensity, float gasConst, float viscosity, float h,
    float g, float tension)
    : mass(mass)
    , restDensity(restDensity)
    , gasConstant(gasConst)
    , viscosity(viscosity)
    , h(h)
    , g(g)
    , tension(tension)
{
    poly6 = 315.0f / (64.0f * PI * pow(h, 9));
    spikyGrad = -45.0f / (PI * pow(h, 6));
    spikyLap = 45.0f / (PI * pow(h, 6));
    h2 = h * h;
    selfDens = mass * poly6 * pow(h, 6);
    massPoly6Product = mass * poly6;
    sphereScale = glm::scale(glm::vec3(h / 2.f));
}

SPHSystem::SPHSystem(
    size_t particleCubeWidth, const SPHSettings &settings,
    const bool &runOnGPU)
    : particleCubeWidth(particleCubeWidth)
    , settings(settings)
    , runOnGPU(runOnGPU)
{
    particleCount = particleCubeWidth * particleCubeWidth * particleCubeWidth;
    particles = (Particle*)malloc(sizeof(Particle) * particleCount);

    // Load in sphere geometry and allocate matrice space
    sphere = new Geometry("resources/lowsphere.obj");
    sphereModelMtxs = new glm::mat4[particleCount];

    initParticles();

	// Generate VBO for sphere model matrices
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, particleCount * sizeof(glm::mat4), &sphereModelMtxs[0], GL_DYNAMIC_DRAW);

	// Setup instance VAO
	glBindVertexArray(sphere->vao);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), 0);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)sizeof(glm::vec4));
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(sizeof(glm::vec4)*2));
	glEnableVertexAttribArray(5);
	glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(sizeof(glm::vec4)*3));

	glVertexAttribDivisor(2,1);
	glVertexAttribDivisor(3,1);
	glVertexAttribDivisor(4,1);
	glVertexAttribDivisor(5,1);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	//start init
	started = false;
}

SPHSystem::~SPHSystem() {
	free(sphereModelMtxs);
    delete[](particles);
}

void SPHSystem::initParticles()
{
	std::srand(1024);
	float particleSeperation = settings.h + 0.01f;
	for (int i = 0; i < particleCubeWidth; i++) {
		for (int j = 0; j < particleCubeWidth; j++) {
			for (int k = 0; k < particleCubeWidth; k++) {
				float ranX
                    = (float(rand()) / float((RAND_MAX)) * 0.5f - 1)
                    * settings.h / 10;
				float ranY
                    = (float(rand()) / float((RAND_MAX)) * 0.5f - 1)
                    * settings.h / 10;
				float ranZ
                    = (float(rand()) / float((RAND_MAX)) * 0.5f - 1)
                    * settings.h / 10;
				glm::vec3 nParticlePos = glm::vec3(
                    i * particleSeperation + ranX - 1.5f,
                    j * particleSeperation + ranY + settings.h + 0.1f,
                    k * particleSeperation + ranZ - 1.5f);

                size_t particleIndex
                    = i + (j + particleCubeWidth * k) * particleCubeWidth;
                Particle* particle = &particles[particleIndex];
                particle->position = nParticlePos;
                particle->velocity = glm::vec3(0);

                sphereModelMtxs[particleIndex]
                    = glm::translate(particle->position) * settings.sphereScale;
			}
		}
	}
}

void SPHSystem::update(float deltaTime) {
	if (!started) return;
	// To increase system stability, a fixed deltaTime is set
	deltaTime = 0.003f;
    updateParticles(
        particles, sphereModelMtxs, particleCount, settings, deltaTime,
        runOnGPU);
}

void SPHSystem::draw(const glm::mat4& viewProjMtx, GLuint shader) {
	// Send matrix data to GPU
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	void* data = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	memcpy(data, sphereModelMtxs, sizeof(glm::mat4) * particleCount);
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Draw instanced particles
	glUseProgram(shader);
	glUniformMatrix4fv(glGetUniformLocation(shader, "viewProjMtx"), 1, false, (float*)&viewProjMtx);
	glBindVertexArray(sphere->vao);
	glDrawElementsInstanced(GL_TRIANGLES, sphere->indices.size(), GL_UNSIGNED_INT, 0, particleCount);
	glBindVertexArray(0);
	glUseProgram(0);
}

void SPHSystem::reset() {
	initParticles();
	started = false;
}

void SPHSystem::startSimulation() {
	started = true;
}