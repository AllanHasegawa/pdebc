
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

//! pdebc namespace
/*!
	Every class in the pdebc library belongs in the namespace pdebc.
*/
namespace pdebc {

//! Abstract/base class for every Differential Evolution class.
/*!
	BaseDE offers a generic interface for any DE class.
	It's main methods are made public, and are overriden by
	its childrens.
	It's use is optional. The usage of a children class directly is allowed.

	\tparam POP_TYPE Population data type (usually 'double')
	\tparam POP_DIM Population dimensions (usually 2D or 3D)
	\tparam ERROR_TYPE Error type (usually 'double')
*/
template <class POP_TYPE, int POP_DIM, class ERROR_TYPE>
struct BaseDE {

	const double kCR_; ///< Mutation rate.
	const double kF_; ///< Mutation weight.

	const std::function<POP_TYPE()>
		callback_population_generator_; ///< Callback for the population generator function.
	const std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>
		callback_calc_error_; ///< Callback for the error calculator function.
	const std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>
		callback_error_evaluation_; ///< Callback for the error evaluator function.

	//! BaseDE constructor
	/*!
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
	BaseDE(const double CR, const double F,
		const std::function<POP_TYPE()>&& callback_population_generator,
		const std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error,
		const std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation) :
			kCR_{CR}, kF_{F},
			callback_population_generator_{callback_population_generator},
			callback_calc_error_{callback_calc_error},
			callback_error_evaluation_{callback_error_evaluation} {

	}

	//! It solves one generation.
	/*!
		This is a blocking method.
	*/
	virtual void solveOneGeneration() = 0;
	//! It solves `N` generations.
	/*!
		\param N Number of generations to solve.
	*/
	virtual void solveNGenerations(const uint32_t N) = 0;
	//! It gets the best candidate.
	/*!
		Will go through the entire population looking for the best candidate.
		The choice is based on the results of the BaseDE::callback_calc_error_ function.
	*/
	virtual std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() = 0;

protected:
	~BaseDE() {

	}
};

} // end namespace pdebc

#endif /* BASEDE_H_ */