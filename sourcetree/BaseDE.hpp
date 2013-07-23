

#ifndef BASEDE_H_
#define BASEDE_H_

#include <functional>

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

	~BaseDE() {

	}
	
	virtual void solveOneGeneration() = 0;
	virtual void solveNGenerations(const uint32_t N) = 0;
	virtual std::tuple<ERROR_TYPE,std::array<POP_TYPE,POP_DIM>> getBestCandidate() const = 0;
};

} // end namespace pdebc

#endif /* BASEDE_H_ */