

#include "BezierCurve.hpp"

#include <array>
#include <tuple>
#include <vector>
#include <cmath>
#include <functional>
#include <cstdio>


BezierCurve::BezierCurve(
	const std::vector<std::tuple<double,Vec2d>> data_points,
	const std::vector<Vec2d> control_points) :
		kNumberControlPoints_{static_cast<uint32_t>(control_points.size())},
		data_points_{data_points},
		control_points_{control_points} {
	
	initializeOptimizationCache();

}

BezierCurve::~BezierCurve() {
}

void BezierCurve::getCurveInT(const double parameterization_value, Vec2d& out) const {
	int n = kNumberControlPoints_ - 1;
	double Bx = 0;
	double By = 0;
	for (int i = 0; i < kNumberControlPoints_; i++) {
		const double p1 = pow(parameterization_value, i);
		const double p2 = pow(1 - parameterization_value, n - i);
		const double b = binomial_cache_[n][i] * p1 * p2;
		const Vec2d& v = control_points_[i];
		Bx += b * v[0];
		By += b * v[1];
	}
	out[0] = Bx;
	out[1] = By;
}

double BezierCurve::calcError() const {
	using namespace std;
	Vec2d temp_curve_p;
	double error{0.0};

	for (auto& p : data_points_) {
		double ll = get<0>(p);
		getCurveInT(ll, temp_curve_p);
		const double dx = get<1>(p)[0] - temp_curve_p[0];
		error += dx * dx;
		const double dy = get<1>(p)[1] - temp_curve_p[1];
		error += dy * dy;
	}
	return error;
}

void BezierCurve::initializeOptimizationCache() {
	using namespace std;
	// First I optimize the p1 and p2 variables
	const int dp_s = data_points_.size();
	vector<vector<double>> p1_p2_caching;
	p1_p2_caching.resize(dp_s);

	for (auto i = p1_p2_caching.begin(); i != p1_p2_caching.end(); i++) {
		i->resize(kNumberControlPoints_);
	}

	for (int i = 0; i < dp_s; i++) {
		for (int j = 0; j < kNumberControlPoints_; j++) {
			const double pv = get<0>(data_points_[i]);
			const double p1 = pow(pv, j);
			const double p2 = pow(1 - pv,
				kNumberControlPoints_ - 1 - j);
			p1_p2_caching[i][j] = p1 * p2;
		}
	}

	// Second I calc the binomial cache
	// 	There's no need to save it in memory...
	function<long(int)> factorial;
	factorial = [&factorial](int f) -> long {
		if (f > 0) {
			return f*factorial(f-1);
		} else {
			return (long)1;
		}
	};
	for (int n = 0; n < kMaxControlPoints; n++) {
		for (int i = 0; i < kMaxControlPoints; i++) {
			binomial_cache_[n][i] = factorial(n)
			/ static_cast<double>(factorial(i) * factorial(n - i));
		}
	}
	auto calcBinomial = [&factorial,this](const int i) -> double {
		return factorial(this->kNumberControlPoints_-1)
			/ static_cast<double>(factorial(i)
				* factorial(this->kNumberControlPoints_-1 - i));
	};

	// Then I can cache the part of the bezier curve equation
	// for all data points and
	// control points
	b_caching_.resize(dp_s);
	for (auto i = b_caching_.begin(); i != b_caching_.end(); i++) {
		i->resize(kNumberControlPoints_);
	}
	for (int p = 0; p < dp_s; p++) {
		for (int i = 0; i < kNumberControlPoints_; i++) {
			b_caching_[p][i] = calcBinomial(i) * p1_p2_caching[p][i];
		}
	}
}

void BezierCurve::updateVariableCPForOptimizationCache(
	const int variable_control_point) {
  // Then I cache the control points that will remain const
  variable_control_point_ = variable_control_point;
  const int np = data_points_.size();
  const_control_point_.resize(np);
  for (int p = 0; p < np; p++) {
    double Bx = 0;
    double By = 0;
    for (int i = 0; i < kNumberControlPoints_; i++) {
      if (i == variable_control_point) {
        continue;
      }
      const double b = b_caching_[p][i];
      const Vec2d& v = control_points_[i];
      Bx += b * v[0];
      By += b * v[1];
    }
    const_control_point_[p][0] = Bx;
    const_control_point_[p][1] = By;
  }
}

double BezierCurve::calcErrorWithOptimizationCache(const Vec2d& candidate_cp) {
	using namespace std;
	const int DP = data_points_.size();

	double ex = 0;
	double ey = 0;
	Vec2d temp_curve;
	for (int k = 1; k < DP - 1; k++) {
		getCurveInTWithOptimizationCache(k, candidate_cp, temp_curve);

		const double dx = get<1>(data_points_[k])[0] - temp_curve[0];
		const double dy = get<1>(data_points_[k])[1] - temp_curve[1];
		ex += dx * dx;
		ey += dy * dy;
	}
	return ex + ey;
}

void BezierCurve::getCurveInTWithOptimizationCache(const int para_index,
	const Vec2d& candidate_cp, Vec2d& out) {
  //const Vec2d& v = control_points_[variable_control_point_];
  out[0] = candidate_cp[0] * b_caching_[para_index][variable_control_point_]
      + const_control_point_[para_index][0];
  out[1] = candidate_cp[1] * b_caching_[para_index][variable_control_point_]
      + const_control_point_[para_index][1];
}