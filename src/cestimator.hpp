/*
    author: Tristan Britt
    email: hello@tbritt.xyz
    
    file: cestimator.hpp
    license: gpl_v3

    this is the header for all the things that makes the main
    routine work, including the referrals to the estimator
    implementations

*/

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <regex>
#include <vector>

#include <Eigen/Dense>

#include "utils/utils.hpp"
#include "estimators/estimators.hpp"

using namespace Eigen;

template<typename M>
M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    /*this is a very slow way to read a file. The problem
      is that the builtin regexing to find \n is quite slow,
      but this normally isnt a problem for small ish files.
      For maximum optimization, you'll require either std::sregex_token_iterator to get content in between the 
      matches, but honestly would recommend a port to boost::regex for files > 10s GB for hyperfine data
      see https://github.com/trbritt/epw_gkk_processor/blob/main/process_g.cpp
    */
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ' ')) {
            // std::cout << std::stod(cell) << " " ;
            values.push_back(std::stod(cell));
        }
        // std::cout << std::endl;
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}

#define assertm(exp, msg) assert(((void)msg, exp));
