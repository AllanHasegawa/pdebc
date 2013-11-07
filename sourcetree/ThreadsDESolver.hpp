
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

#ifndef THREADSDESOLVER_H_
#define THREADSDESOLVER_H_

#include <thread>
#include <random>
#include <mutex>
#include <condition_variable>
#include <algorithm>
 
#include "BaseDE.hpp"

/// \cond DEV
namespace pdebc {

enum class WorkType {
	SOLVE_GENERATION,
	GET_BEST_CANDIDATE
};


//! ThreadsDE internal class.
/*!
	This class is used by ThreadsDE privately, so, Doxygen will ignore it :3
*/
template <class POP_TYPE, int POP_DIM, class ERROR_TYPE>
struct ThreadsDESolver {

	const int kID_;
	const uint32_t kPopSize_;
	BaseDE<POP_TYPE,POP_DIM,ERROR_TYPE>* base_de_;

	std::vector<std::array<POP_TYPE,POP_DIM>> population_;

	ThreadsDESolver(const int id, const uint32_t POP_SIZE,
		BaseDE<POP_TYPE,POP_DIM,ERROR_TYPE>* base_de)
		: kID_{id}, kPopSize_{POP_SIZE},
			base_de_{base_de}, finish_{false}, pending_work_{false},
			work_ready_{false} {

		population_.resize(kPopSize_);
		pop_errors_.resize(kPopSize_);
		
		using MyThreadsDESolver =
			pdebc::ThreadsDESolver<POP_TYPE,POP_DIM,ERROR_TYPE>;
		thread_ = std::thread(&MyThreadsDESolver::run,this);
	}
	~ThreadsDESolver() {
		std::unique_lock<std::mutex> lock(mutex_);
		pending_work_ = true;
		finish_ = true;
		cond_.notify_one();
		lock.unlock();
		thread_.join();
	}

	void solveOneGeneration() {
		using namespace std;
		unique_lock<mutex> lock(mutex_);
		pending_work_ = true;
		work_ready_ = false;
		work_type_ = WorkType::SOLVE_GENERATION;
		cond_.notify_one();
	}

	void solveBestCandidate() {
		using namespace std;
		unique_lock<mutex> lock(mutex_);
		pending_work_ = true;
		work_ready_ = false;
		work_type_ = WorkType::GET_BEST_CANDIDATE;
		cond_.notify_one();
	}

	std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() const {
		return best_candidate_;
	}

	void waitWork() {
		using namespace std;
		unique_lock<mutex> lock(work_ready_lock_);
		work_ready_cond_.wait(lock, [this]() {return this->work_ready_;});
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

	std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> best_candidate_;

	// Threads Flow Control
	std::thread thread_;
	std::mutex mutex_;
	std::condition_variable cond_;
	bool pending_work_;
	bool finish_;
	bool work_ready_;
	WorkType work_type_;
	std::mutex work_ready_lock_;
	std::condition_variable work_ready_cond_;

	void run() {

		using namespace std;
		{ // this scope will be called only once
			// Initialize random_cr_
			mt19937 emt(random_device{}());
			uniform_real_distribution<double> ud(0.0, 1.0);
			random_cr_ = bind(ud, emt);

			// Initialize random_trials_
			mt19937 emt2(random_device{}());
			uniform_int_distribution<uint32_t> ui2(0, kPopSize_-1);
			random_trials_ = bind(ui2, emt2);

			// Initialize random_j_
			mt19937 emt3(random_device{}());
			uniform_int_distribution<uint32_t> ui3(0, POP_DIM-1);
			random_j_ = bind(ui3, emt3);

			generatePopulation();
			calcGenerationError();
		}

		while (!finish_) {
			// non-busy wait for more work
			unique_lock<mutex> lock(mutex_);
			cond_.wait(lock, [this]() {
				return this->pending_work_;
			});
			lock.unlock();

			if (finish_) {
				continue;
			}


			if (work_type_ == WorkType::SOLVE_GENERATION) {
				for (uint32_t i = 0; i < kPopSize_; ++i) {
					mutation(i);
					select(i);
				}
			} else if (work_type_ == WorkType::GET_BEST_CANDIDATE) {
				auto e = std::min_element(pop_errors_.begin(),
					pop_errors_.end(),
					base_de_->callback_error_evaluation_);

				auto min = std::distance(pop_errors_.begin(), e);
				std::array<POP_TYPE,POP_DIM> r = population_[min];

				best_candidate_ = make_tuple(pop_errors_[min],r);
			}


			// Lets tell everyone we are DONE! <sigh>
			///   Need a full cycle to call it done!
			pending_work_ = false;
			lock_guard<mutex> work_read_lock(work_ready_lock_);
			work_ready_ = true;
			work_ready_cond_.notify_one();
		}
	}

	void generatePopulation() {
		for (uint32_t i = 0; i < kPopSize_; ++i) {
			for (int d = 0; d < POP_DIM; ++d) {
				population_[i][d] =
					base_de_->callback_population_generator_();
			}
		}
	}

	void calcGenerationError() {
		for (uint32_t i = 0; i < kPopSize_; ++i) {
			for (int d = 0; d < POP_DIM; ++d) {
				pop_candidate_[d] = population_[i][d];
			}
			pop_errors_[i] = 
				base_de_->callback_calc_error_(pop_candidate_);
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

		pop_trials_[0] = population_[it0];
		pop_trials_[1] = population_[it1];
		pop_trials_[2] = population_[it2];

		pop_candidate_[j] = pop_trials_[0][j] + base_de_->kF_
			* (pop_trials_[1][j] - pop_trials_[2][j]);
		j = (j + 1) % POP_DIM;

		for (int k = 1; k < POP_DIM; ++k) {
			if (random_cr_() <= base_de_->kCR_) {
				pop_candidate_[j] = pop_trials_[0][j]
					+ base_de_->kF_ * (pop_trials_[1][j]
					- pop_trials_[2][j]);
		    } else {
		      	pop_candidate_[j] = population_[actual_index][j];
		    }
		    j = (j + 1) % POP_DIM;
	    }
	}
	

	void select(const uint32_t actual_index) {
		ERROR_TYPE error_new = base_de_->callback_calc_error_(
				pop_candidate_);

		if (base_de_->callback_error_evaluation_(
				error_new, pop_errors_[actual_index])) {
			for (int d = 0; d < POP_DIM; d++) {
				population_[actual_index][d] =
					pop_candidate_[d];
			}
			pop_errors_[actual_index] = error_new;
		}
	}
};

};
/// \endcond

 #endif /* THREADSDESOLVER_H_ */
