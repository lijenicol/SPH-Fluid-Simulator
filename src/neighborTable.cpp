#include <mutex>

#include <neighborTable.h>

/// Global mutex for CPU methods
std::mutex mtx;

/// Returns a hash of the cell position
uint getHash(const glm::ivec3 &cell)
{
    return (
        (uint)(cell.x * 73856093)
        ^ (uint)(cell.y * 19349663)
        ^ (uint)(cell.z * 83492791)
    ) % TABLE_SIZE;
}

/// Get the cell that the particle is in.
glm::ivec3 getCell(Particle *p, float h)
{
    return {p->position.x / h, p->position.y / h, p->position.z / h};
}

/// Populate the particle table
void populateTable(
    Particle *particles, int start, int end, int *particleTable,
    const SPHSettings &settings)
{
    for (int i = start; i < end; i++) {
        Particle* pi = &particles[i];
        uint index = getHash(getCell(pi, settings.h));
        mtx.lock();
        if (particleTable[index] == -1) {
            pi->next = -1;
            particleTable[index] = pi->id;
        }
        else {
            pi->next = particleTable[index];
            particleTable[index] = pi->id;
        }
        mtx.unlock();
    }
}

/// Initialize the table with empty values
void initTable(
    int *particleTable, int start, int end)
{
    for (int i = start; i < end; i++) {
        particleTable[i] = -1;
    }
}

int* createNeighborTable(
    Particle *particles, const size_t &particleCount,
    const SPHSettings &settings)
{
    int *particleTable
        = (int *)malloc(sizeof(int) * TABLE_SIZE);

    const size_t threadCount = std::thread::hardware_concurrency();
    std::thread threads[threadCount];

    size_t blockBoundaries[threadCount + 1];
    size_t tableBoundaries[threadCount + 1];
    blockBoundaries[0] = 0;
    tableBoundaries[0] = 0;
    size_t blockSize = particleCount / threadCount;
    size_t tableBlockSize = TABLE_SIZE / threadCount;
    for (size_t i = 1; i < threadCount; i++) {
        blockBoundaries[i] = i * blockSize;
        tableBoundaries[i] = i * tableBlockSize;
    }
    blockBoundaries[threadCount] = particleCount;
    tableBoundaries[threadCount] = TABLE_SIZE;

    // Init neighbor table
    for (int i = 0; i < threadCount; i++) {
        threads[i] = std::thread(
            initTable, particleTable, tableBoundaries[i],
            tableBoundaries[i + 1]);
    }
    for (std::thread& thread : threads) {
        thread.join();
    }

    // Construct neighbor table
    for (int i = 0; i < threadCount; i++) {
        threads[i] = std::thread(
            populateTable, particles, blockBoundaries[i],
            blockBoundaries[i + 1], particleTable, settings);
    }
    for (std::thread& thread : threads) {
        thread.join();
    }

    return particleTable;
}
