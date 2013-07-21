#include <cstdio>
#include <random>
#include <array>

#include "SequentialDE.hpp"

constexpr int POPULATION_SIZE {512};
constexpr int POPULATION_DIM {2};
using POPULATION_TYPE = double;
using ERROR_TYPE = double;

constexpr double DOMAIN_LIMITS = 128;

int main(int argc, char *argv[]) {
  using pdebc::SequentialDE;
  using namespace std;

  printf("Hello World\n");

  random_device rd;
  mt19937 emt(rd());
  uniform_real_distribution<POPULATION_TYPE> ud(-DOMAIN_LIMITS, +DOMAIN_LIMITS);
  auto rand_domain = bind(ud, emt);

  uniform_int_distribution<uint32_t> up(0,POPULATION_SIZE);
  auto rand_pop = bind(up, emt);


  auto calc_error = [](const std::array<POPULATION_TYPE, POPULATION_DIM>& arr) -> ERROR_TYPE {
    return 1;
  };

  auto error_evaluation = [](const ERROR_TYPE& a, const ERROR_TYPE& b) {
    return a < b;
  };

  SequentialDE<POPULATION_TYPE,2,512,ERROR_TYPE> de {
    std::move(rand_domain), //std::function<POP_TYPE()>&& callback_population_generator
    std::move(rand_pop), //std::function<uint32_t()>&& callback_population_picker
    std::move(calc_error), //std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error
    std::move(error_evaluation) //std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation)
  };

  de.solveOneGeneration();
  printf("%g\n", get<0>(de.getBestCandidate()));
}
