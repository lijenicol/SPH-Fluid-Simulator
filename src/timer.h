#ifndef SPH_TABLE_H
#define SPH_TABLE_H

#include <chrono>
#include <string>

/// \class Timer
///
/// RAII style class which prints out elapsed time between object
/// creation and object deletion.
class Timer
{
public:
    Timer(const std::string &name);
    ~Timer();

private:
    std::string name;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};

#endif // SPH_TABLE_H
