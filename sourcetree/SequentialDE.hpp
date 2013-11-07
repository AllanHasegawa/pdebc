
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

#ifndef SEQUENTIALDE_HPP_
#define SEQUENTIALDE_HPP_

#include <array>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <functional>
#include <random>

#include "BaseDE.hpp"

namespace pdebc {

//! Sequential implementation of the Differential Evolution algorithm.
/*!
	\tparam POP_TYPE Population data type (usually 'double')
	\tparam POP_DIM Population dimensions (usually 2D or 3D)
	\tparam ERROR_TYPE Error type (usually 'double')
*/
template <class POP_TYPE, int POP_DIM, class ERROR_TYPE>
struct SequentialDE : public BaseDE<POP_TYPE, POP_DIM, ERROR_TYPE> {

	const uint32_t kPopSize_; ///< Population size;
	std::vector<std::array<POP_TYPE,POP_DIM>> population_; ///< Entire population.

	/*!
		\param POP_SIZE Population size.

		\param CR Mutation rate. Determines the chances
			of a mutation happening. This value must be 
			between [0,1].
		\param F Mutation weight. Determines how much
			the mutation impacts each trials. This value
			should be between [0,1].
		\param callback_population_generator Function used to generate each
			entity of the population. It must return a POP_TYPE type and use no
			parameters.
		\param callback_calc_error Function used to calculate the error with a single
			member of the population. It must return a ERROR_TYPE type and takes an
			array containg a single population entity as input parameter.
		\param callback_error_evaluation Fuction used to compare two ERROR_TYPE. It
			must return a bool. In case of true, the population from the first ERROR_TYPE
			will be picked as best candidate. Try to figure out what happens in case of false xD.
	*/
	SequentialDE(const uint32_t POP_SIZE, const double CR, const double F,
		const std::function<POP_TYPE()>&& callback_population_generator,
		const std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error,
		const std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation) :
			kPopSize_{POP_SIZE},
			BaseDE<POP_TYPE, POP_DIM, ERROR_TYPE>(
				CR, F,
				std::move(callback_population_generator),
				std::move(callback_calc_error),
				std::move(callback_error_evaluation)) {

		population_.resize(kPopSize_);
		pop_errors_.resize(kPopSize_);
		
		// Initialize random_cr_
		using namespace std;
		random_device rd;
  		mt19937 emt(rd());
  		uniform_real_distribution<double> ud(0.0, 1.0);
  		random_cr_ = bind(ud, emt);

  		// Initialize random_trials_
  		mt19937 emt2(rd());
  		uniform_int_distribution<uint32_t> ui2(0, kPopSize_-1);
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
		for (uint32_t i = 0; i < kPopSize_; i++) {
			mutation(i);
			select(i);
		}
	}

	void solveNGenerations(const uint32_t N) {
		for (uint32_t g = 0; g < N; ++g) {
			solveOneGeneration();
		}
	}

	/*!
		This operation has an O(N) complexity, where N is the population size.
	*/
	std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() {
		
		auto e = std::min_element(pop_errors_.begin(), pop_errors_.end(),
			this->callback_error_evaluation_);
		auto min = std::distance(pop_errors_.begin(), e);

		std::array<POP_TYPE,POP_DIM> r = population_[min];

		return std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>>{pop_errors_[min],r};
	}


private:
	std::function<double()> random_cr_;
	std::function<uint32_t()> random_trials_;
	std::function<uint32_t()> random_j_;

	template <class T, unsigned I, unsigned J>
	using Matrix = std::array<std::array<T, J>, I>;

	Matrix<POP_TYPE, 3, POP_DIM> pop_trials_; // Used in "mutation"
	std::array<POP_TYPE, POP_DIM> pop_candidate_;
	std::vector<ERROR_TYPE> pop_errors_;

	void generatePopulation() {
		for (uint32_t i = 0; i < kPopSize_; ++i) {
			for (int d = 0; d < POP_DIM; ++d) {
				population_[i][d] = this->callback_population_generator_();
			}
		}
	}

	void calcGenerationError() {
		for (uint32_t i = 0; i < kPopSize_; ++i) {
			pop_candidate_ = population_[i];
			pop_errors_[i] = this->callback_calc_error_(pop_candidate_);
		}
	}

	void mutation(const uint32_t actual_index) {
		int j = random_j_();

		const uint32_t it0 = random_trials_();
		uint32_t it1 = random_trials_();
		while (it1 == it0) {
			it1 = random_trials_();
		}
		uint32_t it2 = random_trials_();
		while (it2 == it1 || it2 == it0) {
			it2 = random_trials_();
		}

		for (int d = 0; d < POP_DIM; d++) {
			pop_trials_[0][d] = population_[it0][d];
			pop_trials_[1][d] = population_[it1][d];
			pop_trials_[2][d] = population_[it2][d];
		}

		pop_candidate_[j] = pop_trials_[0][j] + this->kF_ * (pop_trials_[1][j] - pop_trials_[2][j]);
		j = (j + 1) % POP_DIM;

		for (int k = 1; k < POP_DIM; ++k) {
			if (random_cr_() <= this->kCR_) {
				pop_candidate_[j] = pop_trials_[0][j] + this->kF_ * (pop_trials_[1][j] - pop_trials_[2][j]);
		    } else {
		      	pop_candidate_[j] = population_[actual_index][j];
		    }
		    j = (j + 1) % POP_DIM;
	    }
	}
	

	void select(const uint32_t actual_index) {
		ERROR_TYPE error_new = this->callback_calc_error_(pop_candidate_);

		if (this->callback_error_evaluation_(error_new, pop_errors_[actual_index])) {
			for (int d = 0; d < POP_DIM; d++) {
				population_[actual_index][d] = pop_candidate_[d];
			}
			pop_errors_[actual_index] = error_new;
		}
	}
};

} // end namespace pdebc

#endif /* SEQUENTIALDE_HPP_ */
