
#ifndef SEQUENTIALDE_HPP_
#define SEQUENTIALDE_HPP_

#include <array>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <functional>
#include <random>

//#include "BaseDE.hpp"

namespace pdebc {

template <class T, unsigned I, unsigned J>
using Matrix = std::array<std::array<T, J>, I>;

template <class POP_TYPE, int POP_DIM, int POP_SIZE, class ERROR_TYPE>
struct SequentialDE {
	double kCR_;
	double kF_;

	Matrix<POP_TYPE, POP_SIZE, POP_DIM> population_;

	std::function<POP_TYPE()> callback_population_generator_;
	std::function<uint32_t()> callback_population_picker_;
	std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)> callback_calc_error_;
	std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)> callback_error_evaluation_;

	SequentialDE(
		std::function<POP_TYPE()>&& callback_population_generator,
		std::function<uint32_t()>&& callback_population_picker,
		std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error,
		std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation) :
			callback_population_generator_{callback_population_generator},
			callback_population_picker_{callback_population_picker},
			callback_calc_error_{callback_calc_error},
			callback_error_evaluation_{callback_error_evaluation} {

		// Initialize random_cr_
		using namespace std;
		random_device rd;
  		mt19937 emt(rd());
  		uniform_real_distribution<double> ud(0.0, 1.0);
  		random_cr_ = bind(ud, emt);

  		// Initialize random_trials_
  		mt19937 emt2(rd());
  		uniform_int_distribution<uint32_t> ui2(0, POP_SIZE-1);
  		random_trials_ = bind(ui2, emt2);

  		// Initialize random_j_
  		mt19937 emt3(rd());
  		uniform_int_distribution<uint32_t> ui3(0.0, POP_DIM-1);
  		random_j_ = bind(ui3, emt3);

		generatePopulation();
		calcGenerationError();
	}

	~SequentialDE() {

	}

	void solveOneGeneration() {
		for (uint32_t i = 0; i < POP_SIZE; i++) {
			mutation(i);
			select(i);
		}
	}

	std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() {
		
		auto e = std::min_element(pop_errors_.begin(), pop_errors_.end(),
			callback_error_evaluation_);
		auto min = std::distance(pop_errors_.begin(), e);

		std::array<POP_TYPE,POP_DIM> r;
		for (int d = 0; d < POP_DIM; d++) {
			r[d] = population_[min][d];
		}

		return {pop_errors_[min],r};
	}


private:

	std::function<double()> random_cr_;
	std::function<uint32_t()> random_trials_;
	std::function<uint32_t()> random_j_;

	Matrix<POP_TYPE, 3, POP_DIM> pop_trials_; // Used in "mutation"
	std::array<POP_TYPE, POP_DIM> pop_candidate_;
	std::array<ERROR_TYPE, POP_SIZE> pop_errors_;

	void generatePopulation() {
		for (uint32_t i = 0; i < POP_SIZE; ++i) {
			for (int d = 0; d < POP_DIM; ++d) {
				population_[i][d] = callback_population_generator_();
			}
		}
	}

	void calcGenerationError() {
		for (uint32_t i = 0; i < POP_SIZE; ++i) {
			for (int d = 0; d < POP_DIM; ++d) {
				pop_candidate_[d] = population_[i][d];
			}
			pop_errors_[i] = callback_calc_error_(pop_candidate_);
		}
	}

	void mutation(const uint32_t actual_index) {
		int j = random_j_();

		const uint32_t it0 = random_trials_();
		const uint32_t it1 = random_trials_();
		const uint32_t it2 = random_trials_();

		for (int d = 0; d < POP_DIM; d++) {
			pop_trials_[0][d] = population_[it0][d];
			pop_trials_[1][d] = population_[it1][d];
			pop_trials_[2][d] = population_[it2][d];
		}

		pop_candidate_[j] = pop_trials_[0][j] + kF_ * (pop_trials_[1][j] - pop_trials_[2][j]);
		j = (j + 1) % POP_DIM;

		for (int k = 1; k < POP_DIM; ++k) {
			if (random_cr_() <= kCR_) {
				pop_candidate_[j] = pop_trials_[0][j] + kF_ * (pop_trials_[1][j] - pop_trials_[2][j]);
		    } else {
		      	pop_candidate_[j] = population_[actual_index][j];
		    }
		    j = (j + 1) % POP_DIM;
	    }
	}

	void select(const uint32_t actual_index) {
		ERROR_TYPE error_new = callback_calc_error_(pop_candidate_);

		if (callback_error_evaluation_(error_new, pop_errors_[actual_index])) {
			for (int d = 0; d < POP_DIM; d++) {
				population_[actual_index][d] = pop_candidate_[d];
			}
			pop_errors_[actual_index] = error_new;
		}
	}
};

} // end namespace pdebc

#endif /* SEQUENTIALDE_HPP_ */