
#ifndef SEQUENTIALDE_HPP_
#define SEQUENTIALDE_HPP_

#include <array>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <functional>

//#include "BaseDE.hpp"

namespace pdebc {

template <class T, unsigned I, unsigned J>
using Matrix = std::array<std::array<T, J>, I>;

template <class POP_TYPE = double, int POP_DIM = 2, int POP_SIZE = 512, class ERROR_TYPE = double>
struct SequentialDE {
	double kCR_;
	double kF_;

	Matrix<POP_TYPE, POP_SIZE, POP_DIM> population_;

	std::function<POP_TYPE()> callback_population_generator_;
	std::function<uint32_t(const uint32_t)> callback_population_picker;
	std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)> callback_calc_error_;
	std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)> callback_error_evaluation_;

	SequentialDE(
		std::function<POP_TYPE()>&& callback_population_generator,
		std::function<uint32_t()>&& callback_population_picker,
		std::function<ERROR_TYPE(...)>&& callback_calc_error,
		std::function<bool(ERROR_TYPE,ERROR_TYPE)>&& callback_error_evaluation) {

		generatePopulation();
		calcGenerationError();
	}

	~SequentialDE() {

	}

	void solveOneGeneration() {
		ERROR_TYPE error_best_candidate = callback_calc_error_();
		for (uint32_t i = 0; i < POP_SIZE; i++) {
			mutation(i);
			ERROR_TYPE error_new = callback_calc_error_(pop_candidate_);
			if (callback_error_evaluation_(pop_errors_[i], error_new)) {
				for (int d = 0; d < POP_DIM; d++) {
					population_[i][d] = pop_candidate_[d];
				}
				pop_errors_[i] = error_new;
			}
		}
	}

	std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() {
		
		ERROR_TYPE e = std::min_element(pop_errors_.begin(), pop_errors_.end(),
			callback_error_evaluation_);
		auto min = std::distance(pop_errors_.begin(), e);

		std::array<POP_TYPE,POP_DIM> r;
		for (int d = 0; d < POP_DIM; d++) {
			r[d] = population_[d][min];
		}

		return {e,r};
	}


private:

	Matrix<POP_TYPE, 3, POP_DIM> pop_trials_; // Used in "mutation"
	std::array<POP_TYPE, POP_DIM> pop_candidate_;
	std::array<ERROR_TYPE, POP_SIZE> pop_errors_;

	void generatePopulation() {
		for (int d = 0; d < POP_DIM; ++d) {
			for (uint32_t i = 0; i < POP_SIZE; ++i) {
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
		int j = 0;//random_j_[r_rs_];
		pop_candidate_[j] = pop_trials_[0][j] + kF_ * (pop_trials_[1][j] - pop_trials_[2][j]);
		j = (j + 1) % POP_DIM;

		for (int k = 1; k < POP_DIM; ++k) {
			if (0/*random_cr_[r_rs_][0]*/ <= kCR_) {
				pop_candidate_[j] = pop_trials_[0][k] + kF_ * (pop_trials_[1][k] - pop_trials_[2][k]);
		    } else {
		      	pop_candidate_[j] = population_[actual_index][k];
		    }
		    j = (j + 1) % POP_DIM;
	    }
	}

	void select() {

	}
};

} // end namespace pdebc

#endif /* SEQUENTIALDE_HPP_ */