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

#ifndef BCDESOLVER_H_
#define BCDESOLVER_H_

#include <array>
#include <vector>
#include <memory>

#include "BezierCurve.h"
#include "Vec2.h"

class BCDESolverST;
class BCDESolver {
public:
	const double kDE_F_;
	const double kDE_CR_;
	const int kPopulation_;
	const int kNProcess_;

	BezierCurve& bezier_curve_;
	const std::vector<Vec2>& data_points_;
	const std::vector<double>& parameterization_values_;


	int generation_;
	/*
	 * [NUMBER PROCESS][DIM]
	 */
	std::vector<Vec2> lowest_error_index_;

	// this is used by the migration
	// each thread calc a random number ;)
	// Callgrind said it was slow =X
  std::vector<int> random_population_interval_;

	BCDESolver(const std::vector<double>& parameterization_values,
			const std::vector<Vec2>& data_points, const int n_process,
			const double de_f, const double de_cr, const int population,
			BezierCurve& bezier_curve);
	virtual ~BCDESolver();

	void SolveOneGeneration();


private:

	std::vector<std::shared_ptr<BCDESolverST>> solvers_;

	void Initialize();
	void Migration(const int control_point);
};

#endif /* BCDESOLVER_H_ */
