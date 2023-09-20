#include <timer.h>

#include <iostream>

Timer::Timer(const std::string &name)
    : name(name)
{
    start = std::chrono::high_resolution_clock::now();
}

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();

    auto s = std::chrono::time_point_cast<std::chrono::microseconds>(start)
        .time_since_epoch()
        .count();
    auto e = std::chrono::time_point_cast<std::chrono::microseconds>(end)
        .time_since_epoch()
        .count();

    auto duration = (e - s) * 0.001;
    std::cerr << "<" << name << "> clock finished: " << duration << std::endl;
}
