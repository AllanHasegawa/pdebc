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

#include "BCDESolver.h"

#include <vector>
#include <random>
#include <functional>

#include "Vec2.h"
#include "BCDESolverST.h"

BCDESolver::BCDESolver(const std::vector<double>& parameterization_values,
		const std::vector<Vec2>& data_points, const int n_process,
		const double de_f, const double de_cr, const int population,
		BezierCurve& bezier_curve) :
		bezier_curve_(bezier_curve), data_points_(data_points), parameterization_values_(
				parameterization_values), n_process_(n_process), kDE_F_(de_f), kDE_CR_(
				de_cr), kPopulation_(population), generation_(0) {

	parameters_.resize(bezier_curve.kNumberControlPoints_ - 2);
	for (auto i = parameters_.begin(), e = parameters_.end(); i != e; i++) {
		i->resize(kPopulation_);
	}
	lowest_error_index_.resize(n_process);
	Initialize();
}

BCDESolver::~BCDESolver() {

}

void BCDESolver::SolveOneGeneration() {
	using namespace std;

	vector<BCDESolverST> solvers;

	//for (int i = 0; i < n_process_; i++) {
		BCDESolverST solver(*this, 1, 2);
		//solvers.push_back(solver);
		solver.Start();
	//}

	//for (int i = 0; i < n_process_; i++) {
		//solvers[i].Join();
	//}
		solver.Join();
}

void BCDESolver::Initialize() {
	const int bcn = bezier_curve_.kNumberControlPoints_;

	const int w = 1024;
	const int h = 1024;

	std::uniform_real_distribution<double> distribution(0, 1);
	std::mt19937 engine; // Mersenne twister MT19937
	auto generator = std::bind(distribution, engine);

	// (bc-2): We ignore the first and last control points :3
	// random population for each control point...
	for (int i = 0; i < bcn - 2; i++) {
		for (int j = 0; j < kPopulation_; j++) {
			parameters_[i][j].x = -(w) + generator() * (w * 2);
			parameters_[i][j].y = -(h) + generator() * (h * 2);
		}
	}
}

void BCDESolver::Migration(const int control_point) {
	const double phi = 1;

	const int populationInterval = kPopulation_ / n_process_;

	std::uniform_int_distribution<int> distribution(0, populationInterval);
	std::mt19937 engine; // Mersenne twister MT19937
	auto generator = std::bind(distribution, engine);

	for (int k = 0; k < n_process_; k++) {
		if (0 <= phi) {
			int destinyIndex = k + 1;
			if (destinyIndex >= n_process_) {
				destinyIndex = 0;
			}
			int pI = generator();
			int pSDestiny = populationInterval * destinyIndex;
			parameters_[control_point][pSDestiny + pI].x =
					parameters_[control_point][lowest_error_index_[k].x].x;
			parameters_[control_point][pSDestiny + pI].y =
					parameters_[control_point][lowest_error_index_[k].y].y;
		}
	}
}
