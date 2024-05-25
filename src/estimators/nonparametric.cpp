/*
 *  author: Tristan Britt
 *  email: hello@tbritt.xyz
 *
 *  file: nonparametric.cpp
 *  license: gpl_v3
 *
 *  this is the main code for the nonparametric estimator
 *
 */

#include "estimators.hpp"

int Cestimator::non_parametric::run() noexcept{
   //non parametric estimates done in construction
   mu    = mu_no_par;
   sigma = sigma_no_par;
   return 0;
}
