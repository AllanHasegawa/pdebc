
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

#ifndef BASEDE_H_
#define BASEDE_H_

#include <functional>
#include <tuple>

namespace pdebc {


template <class POP_TYPE, int POP_DIM, int POP_SIZE, class ERROR_TYPE>
struct BaseDE {

	const double kCR_;
	const double kF_;

	const std::function<POP_TYPE()> callback_population_generator_;
	const std::function<uint32_t()> callback_population_picker_;
	const std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)> callback_calc_error_;
	const std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)> callback_error_evaluation_;

	BaseDE(const double CR, const double F,
		const std::function<POP_TYPE()>&& callback_population_generator,
		const std::function<uint32_t()>&& callback_population_picker,
		const std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error,
		const std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation) :
			kCR_{CR}, kF_{F},
			callback_population_generator_{callback_population_generator},
			callback_population_picker_{callback_population_picker},
			callback_calc_error_{callback_calc_error},
			callback_error_evaluation_{callback_error_evaluation} {

	}

	virtual void solveOneGeneration() = 0;
	virtual void solveNGenerations(const uint32_t N) = 0;
	virtual std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() const = 0;

protected:
	~BaseDE() {

	}
};

} // end namespace pdebc

#endif /* BASEDE_H_ */