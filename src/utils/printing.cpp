#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <ctime>
#include <cctype>
#include <iterator>
#include "utils.hpp"

std::string getFormattedTime(const std::tm &time_struct, const std::string &format) {
  // Create a stringstream to capture the formatted output
  std::stringstream ss;

  // Use std::put_time to format the time and write to the stringstream
  ss << std::put_time(&time_struct, format.c_str());

  // Extract the formatted string from the stringstream
  return ss.str();
}

//pretty print vector
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{";
    if(vec.size()!=0)
    {
        std::copy(vec.begin(), vec.end()-1, std::ostream_iterator<T>(os, " "));
        os<<vec.back();
    }
    os<<"}";
    return os;
}

std::string Utils::word_wrap(std::string text, unsigned per_line)
{
    unsigned line_begin = 0;

    while (line_begin < text.size())
    {
        const unsigned ideal_end = line_begin + per_line ;
        unsigned line_end = ideal_end < text.size() ? ideal_end : text.size()-1;

        if (line_end == text.size() - 1)
            ++line_end;
        else if (std::isspace(text[line_end]))
        {
            text[line_end] = '\n';
            ++line_end;
        }
        else    // backtrack
        {
            unsigned end = line_end;
            while ( end > line_begin && !std::isspace(text[end]))
                --end;

            if (end != line_begin)                  
            {                                       
                line_end = end;                     
                text[line_end++] = '\n';            
            }                                       
            else                                    
                text.insert(line_end++, 1, '\n');
        }

        line_begin = line_end;
    }

    return text;
}

void Utils::banner(std::string fname, int dim1, int dim2) {
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm local_time = *std::localtime(&now_time_t);

    std::string header = "";
    header += Utils::Colors::HEADER + "###################################################################\n";
    header += "# Cestimator v0.1\n";
    header += "# Running on " + fname + "\n";
    header += "# Analysing columns" + std::to_string(dim1) + " and " + std::to_string(dim2) + "\n";
    header += "# Starting on " + getFormattedTime(local_time, "%Y-%m-%d %H:%M:%S") + "\n";
    header += "###################################################################" + Utils::Colors::ENDC + "\n";

    std::cout << header << std::endl;

}

void Utils::goodbye(std::vector<Utils::Timer> timers){
    // print timing information
    std::cout << std::endl << "-------------------------------------------------------------------" << std::endl;
    for (int i = 0; i < timers.size(); i++)
    {
        timers[i].print();
    }
    std::cout << "-------------------------------------------------------------------" << std::endl;
    // Get the current time
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

    // Convert the time to a time_t representation (seconds since epoch)
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

    // Convert time_t to a local time representation (struct tm)
    std::tm local_time = *std::localtime(&now_time_t);

    std::cout << Utils::Colors::OKGREEN << "###################################################################" << std::endl;
    std::cout << "# Cleanly exiting; ended on " << std::put_time(&local_time, "%Y-%m-%d %H:%M:%S") << " ... Goodbye." << std::endl;
    std::cout << "###################################################################" << Utils::Colors::ENDC << std::endl;

}

