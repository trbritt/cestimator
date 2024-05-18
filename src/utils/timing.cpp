#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include "utils.hpp"

void Utils::Timer::start() {
    this->epoch = Utils::Timer::Clock::now();
}

void Utils::Timer::stop() {
    this->telapsed = Utils::Timer::Clock::now() - this->epoch;
}

void Utils::Timer::print() {
    size_t num_spaces = std::max(0, static_cast<int>(30-this->name.length()));
    std::string padding(num_spaces, ' ');

    auto t =  std::chrono::duration_cast<std::chrono::microseconds>(this->telapsed).count();
    std::cout  << this->name+padding << "\t:\t" << t << " μs" << std::endl;
}
