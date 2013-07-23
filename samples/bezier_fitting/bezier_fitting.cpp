/*
 Copyright 2012 Allan Yoshio Hasegawa

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 -----------------------------------------------------------------------------
 */

/*

Bezier Fitting Sample

*/


#include <cstdio>
#include <random>
#include <array>
#include <tuple>
#include <cmath>

#include "pdebc/SequentialDE.hpp"

constexpr std::array<double,2> POINT {{23.345, -10.009}};

constexpr int POPULATION_SIZE {8};
constexpr int POPULATION_DIM {2};
using POPULATION_TYPE = double;
using ERROR_TYPE = double;

constexpr double DOMAIN_LIMITS = 128;


int main(int argc, char *argv[]) {
  using pdebc::SequentialDE;
  using pdebc::BaseDE;
  using namespace std;

  random_device rd;
  mt19937 emt(rd());
  uniform_real_distribution<POPULATION_TYPE> ud(-DOMAIN_LIMITS, +DOMAIN_LIMITS);
  auto rand_domain = bind(ud, emt);

  uniform_int_distribution<uint32_t> up(0,POPULATION_SIZE);
  auto rand_pop = bind(up, emt);

  auto calc_error = [](const std::array<POPULATION_TYPE, POPULATION_DIM>& arr) -> ERROR_TYPE {
    double e = std::sqrt(
      std::pow(arr[0]-POINT[0],2) +
      std::pow(arr[1]-POINT[1],2)
      );
    //printf("ce (%g,%g): %g\n", arr[0], arr[1], e);
    return e;
  };

  auto error_evaluation = [](const ERROR_TYPE& a, const ERROR_TYPE& b) {
    return a < b;
  };
  

  SequentialDE<POPULATION_TYPE,POPULATION_DIM,POPULATION_SIZE,ERROR_TYPE> de {
    0.5, 0.8,
    std::move(rand_domain), //std::function<POP_TYPE()>&& callback_population_generator
    std::move(rand_pop), //std::function<uint32_t()>&& callback_population_picker
    std::move(calc_error), //std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error
    std::move(error_evaluation) //std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation)
  };

  
  for (int i = 0; i < 20; i++) {
    de.solveOneGeneration();
  }

  auto bc_error = get<0>(de.getBestCandidate());
  auto bc_point = get<1>(de.getBestCandidate());
  printf("POINT: (%g,%g)\n", POINT[0],POINT[1]);
  printf("Best candidate point: (%g,%g)\n", bc_point[0], bc_point[1]);
  printf("Best Candidate error: %g\n", std::sqrt(bc_error));
}