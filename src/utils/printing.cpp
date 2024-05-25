#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <ctime>
#include <cctype>
#include <iterator>
#include "utils.hpp"

std::string getFormattedTime(const std::tm&time_struct, const std::string&format) {
   // Create a stringstream to capture the formatted output
   std::stringstream ss;

   // Use std::put_time to format the time and write to the stringstream
   ss << std::put_time(&time_struct, format.c_str());

   // Extract the formatted string from the stringstream
   return ss.str();
}

void Cestimator::Utils::banner(std::string fname) {
   std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
   std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
   std::tm     local_time = *std::localtime(&now_time_t);

   std::string header = "";
   header += Cestimator::Utils::colors::HEADER + "###################################################################\n";
   header += "# Cestimator v0.1\n";
   header += "# Running on " + fname + "\n";
   header += "# Analysing all columns\n";
   header += "# Starting on " + getFormattedTime(local_time, "%Y-%m-%d %H:%M:%S") + "\n";
   header += "###################################################################" + Cestimator::Utils::colors::ENDC + "\n";

   std::cout << header << std::endl;
}

void Cestimator::Utils::goodbye(std::vector <Cestimator::Utils::Timer> timers) {
   // print timing information
   std::cout << std::endl << Cestimator::Utils::colors::UNDERLINE << "Timings\n" << Cestimator::Utils::colors::ENDC << std::endl;

   for (int i = 0; i < timers.size(); i++) {
      timers[i].print();
   }
   Cestimator::Utils::goodbye();
}

void Cestimator::Utils::goodbye() {
   // Get the current time
   std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

   // Convert the time to a time_t representation (seconds since epoch)
   std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

   // Convert time_t to a local time representation (struct tm)
   std::tm local_time = *std::localtime(&now_time_t);
   std::cout << Cestimator::Utils::colors::OKGREEN << "###################################################################" << std::endl;
   std::cout << "# Cleanly exiting; ended on " << std::put_time(&local_time, "%Y-%m-%d %H:%M:%S") << " ... Goodbye." << std::endl;
   std::cout << "###################################################################" << Cestimator::Utils::colors::ENDC << std::endl;
}
