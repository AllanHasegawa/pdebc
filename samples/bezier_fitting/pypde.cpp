#include "pypde.hpp"

#include <vector>
#include <memory>
#include <tuple>

#include "pdebc/ThreadsDE.hpp"



std::vector<double> calcChordLengthSwig(const std::vector<Vec2d>& data_points) {
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

pypde::pypde(const int n_processes, const int population_size,
		const int bezier_control_points,
		std::vector<Vec2> data_points) {
	using namespace std;
	vector<Vec2d> data_points_2dpos;
	for_each(data_points.begin(),data_points.end(),
		[&data_points_2dpos](const Vec2& v){
			data_points_2dpos.push_back({{v.x,v.y}});
	});

	auto chord_length = calcChordLengthSwig(data_points_2dpos);

	vector<std::tuple<double,Vec2d>> dp{};
	for (int i = 0; i < chord_length.size(); i++) {
		dp.push_back(tuple<double,Vec2d>{chord_length[i],data_points_2dpos[i]});
	}

	bezier_curve_ = new BezierCurve(dp,
			move(vector<Vec2d>(bezier_control_points))
		);
	/* first bezier curve control point == first data point */
	/* second  "      "     "       "   == second "     "   */
	bezier_curve_->control_points_[0] = data_points_2dpos[0];
	bezier_curve_->control_points_[bezier_curve_->control_points_.size()-1]
		= data_points_2dpos[data_points_2dpos.size()-1];

	/* -2 because we dont try to fit the first and last control point */
	for (int i = 0; i < bezier_control_points-2; i++) {
		bezier_curve_->updateVariableCPForOptimizationCache(i+1);

		// population generator
		auto t1 = chrono::high_resolution_clock::now().time_since_epoch();
		mt19937 emt(chrono::duration_cast<chrono::nanoseconds>(t1).count());
		uniform_real_distribution<POPULATION_TYPE> ud(-DOMAIN_LIMITS, +DOMAIN_LIMITS);
		auto rand_domain = bind(ud, emt);

		// error evaluations
		auto error_evaluation =
			[](const ERROR_TYPE& a, const ERROR_TYPE& b) {
				return a < b;
		};
		
		// error calculation
		auto calc_error =
			[this,i](const array<POPULATION_TYPE, POPULATION_DIM>& arr) -> ERROR_TYPE {
				return this->bezier_curve_->calcErrorWithOptimizationCache(arr);
		};
		
		
		des_.push_back(make_shared<PYPDE_ThreadsDE>(
			n_processes, 1, population_size, 0.5, 0.8,
			std::move(rand_domain),
			std::move(calc_error),
			std::move(error_evaluation)
		));
	}
}

pypde::~pypde() {
	delete bezier_curve_;
}


void pypde::solveOneGeneration()  {
	using namespace std;
	for (int j = 0; j < des_.size(); ++j) {
		bezier_curve_->updateVariableCPForOptimizationCache(j+1);
		auto& d = des_[j];
		d->solveOneGeneration();
		//auto bc_error = get<0>(d.getBestCandidate());
		auto bc_point = get<1>(d->getBestCandidate());
		//printf("Best candidate middle control-point: (%g,%g)\n", bc_point[0], bc_point[1]);
		//printf("Best Candidate error: %g\n", std::sqrt(bc_error));
		bezier_curve_->control_points_[j+1] = bc_point;
	}
}

double pypde::getBestCandidateError(int i) {
	return std::get<0>(des_[i]->getBestCandidate());
}

std::vector<double> pypde::getBestCandidateCP(int i) {
	std::vector<double> v(2);
	auto p = std::get<1>(des_[i]->getBestCandidate());
	v[0] = p[0];
	v[1] = p[1];
	return v;
}