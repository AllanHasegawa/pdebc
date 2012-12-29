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
#include "BezierCurve.h"
#include "Vec2.h"

class BCDESolver {
public:
	const double kDE_F_;
	const double kDE_CR_;
	const int kPopulation_;
	const int kNProcess_;

	BezierCurve& bezier_curve_;
	const std::vector<Vec2>& data_points_;
	const std::vector<double>& parameterization_values_;

	/*
	 *  [control point][population][dimension]
	 *  Here, "CP-2" is because we ignore the first and last CP
	 */
	std::vector<std::vector<Vec2>> parameters_;
	std::vector<std::vector<Vec2>> population_errors_;

	int generation_;
	/*
	 * [NUMBER PROCESS][DIM]
	 */
	std::vector<Vec2> lowest_error_index_;

	BCDESolver(const std::vector<double>& parameterization_values,
			const std::vector<Vec2>& data_points, const int n_process,
			const double de_f, const double de_cr, const int population,
			BezierCurve& bezier_curve);
	virtual ~BCDESolver();

	void SolveOneGeneration();


private:

	void Initialize();
	void Migration(const int control_point);
};

#endif /* BCDESOLVER_H_ */
