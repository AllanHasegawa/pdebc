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

/*

Bezier Fitting Sample

*/


#include <cstdio>
#include <random>
#include <array>
#include <tuple>
#include <cmath>
#include <vector>

#include "pdebc/SequentialDE.hpp"

#include "BezierCurve.hpp"
/*
Aux function to calculate the parameterization values
of the Bezier Curve using the Chord Length method
*/
std::vector<double> calcChordLength(const std::vector<Vec2d>& data_points) {
	using namespace std;
	const int DP = data_points.size();
	vector<double> chord_length(DP);
	chord_length[0] = 0;
	chord_length[DP - 1] = 1;

	double td = 0;

	for (int i = 1; i < DP; i++) {
		const double vdx = data_points[i][0] - data_points[i - 1][0];
		const double vdy = data_points[i][1] - data_points[i - 1][1];
		const double d = sqrt(vdx * vdx + vdy * vdy);
		td += d;
	}

	double p_td = 0;
	for (int i = 1; i < DP; i++) {
		const double vdx = data_points[i][0] - data_points[i - 1][0];
		const double vdy = data_points[i][1] - data_points[i - 1][1];
		const double d = sqrt(vdx * vdx + vdy * vdy);
		p_td += d;
		chord_length[i] = p_td / td;
	}
	return chord_length;
}


constexpr int POPULATION_SIZE {128};
constexpr int POPULATION_DIM {2};
using POPULATION_TYPE = double;
using ERROR_TYPE = double;
constexpr double DOMAIN_LIMITS = 128;


int main(int argc, char *argv[]) {
	using SequentialDE =
		pdebc::SequentialDE<POPULATION_TYPE,POPULATION_DIM,POPULATION_SIZE,ERROR_TYPE>;
	using namespace std;

	/* populate the BezierCurve information */
	// 'data_points' is a vector of tuples (double,Vec2d)
	//		the double value is the parameterization values
	//		the Vec2d is it's 2D position
	// data points will have a form like a wave
	auto data_points_2dpos = vector<Vec2d>
		{{{-10,0}}, {{0,10}}, {{10,0}}, {{20,-10}}, {{30,0}}};

	auto chord_length = calcChordLength(data_points_2dpos);

	vector<std::tuple<double,Vec2d>> data_points{};
	for (int i = 0; i < chord_length.size(); i++) {
		data_points.push_back(std::tuple<double,Vec2d>{chord_length[i],data_points_2dpos[i]});
	}

	/* lets create the BezierCurve control points */
	// i will put all the points at position (0,0)
	// 		the algorith will find (fit) the proper position
	//		important now is to pass the size of the vector
	BezierCurve bezier_curve{data_points, move(vector<Vec2d>(4))};

	/* first bezier curve control point == first data point */
	/* second  "      "     "       "   == second "     "   */
	bezier_curve.control_points_[0] = data_points_2dpos[0];
	bezier_curve.control_points_[bezier_curve.control_points_.size()-1]
		= data_points_2dpos[data_points_2dpos.size()-1];

	/* lets create the callback functions */
	// generate the population at random
	random_device rd;
	mt19937 emt(rd());
	uniform_real_distribution<POPULATION_TYPE> ud(-DOMAIN_LIMITS, +DOMAIN_LIMITS);
	auto rand_domain = bind(ud, emt);

	// pick the candidates at random
	uniform_int_distribution<uint32_t> up(0,POPULATION_SIZE);
	auto rand_pop = bind(up, emt);

	// error evaluations
	auto error_evaluation =
		[](const ERROR_TYPE& a, const ERROR_TYPE& b) {
			return a < b;
	};

	/* lets create one "DE" algorithm for each control point */
	vector<SequentialDE> des;
	for (int i = 0; i < 2; i++) {
		// Since I will use the "OptimizationCache"
		// We will have to update is any time we change
		// control points
		bezier_curve.updateVariableCPForOptimizationCache(i+1);
		// each DE will have a unique error calculation function
		auto calc_error =
			[&bezier_curve,i](const array<POPULATION_TYPE, POPULATION_DIM>& arr) -> ERROR_TYPE {
				return bezier_curve.calcErrorWithOptimizationCache(arr);
		};
		

		des.push_back(SequentialDE{
			0.5, 0.8,
			std::move(rand_domain), //std::function<POP_TYPE()>&& callback_population_generator
			std::move(rand_pop), //std::function<uint32_t()>&& callback_population_picker
			std::move(calc_error), //std::function<ERROR_TYPE(const std::array<POP_TYPE,POP_DIM>&)>&& callback_calc_error
			std::move(error_evaluation) //std::function<bool(const ERROR_TYPE&,const ERROR_TYPE&)>&& callback_error_evaluation)
		});
	}


	for (int i = 0; i < 200; i++) {
		printf("%s\n", string(40,'*').c_str());
		printf("Generation %d:\n", i);
		for (int j = 0; j < 2; ++j) {
			bezier_curve.updateVariableCPForOptimizationCache(j+1);
			auto& d = des[j];
			d.solveOneGeneration();
			auto bc_error = get<0>(d.getBestCandidate());
			auto bc_point = get<1>(d.getBestCandidate());
			printf("Best candidate middle control-point: (%g,%g)\n", bc_point[0], bc_point[1]);
			printf("Best Candidate error: %g\n", std::sqrt(bc_error));
			bezier_curve.control_points_[j+1] = bc_point;
		}
	}
}