

#include "BezierCurve.hpp"

#include <array>
#include <tuple>
#include <vector>
#include <cmath>
#include <functional>
#include <cstdio>

std::array<std::array<double, BezierCurve::kMaxControlPoints>, BezierCurve::kMaxControlPoints> BezierCurve::binomial_cache_;


BezierCurve::BezierCurve(
	const std::vector<std::tuple<double,Vec2d>>& data_points,
	const std::vector<Vec2d> control_points
	) :
	kNumberControlPoints_{static_cast<uint32_t>(control_points.size())},
	data_points_{data_points},
	control_points_{control_points} {
	calcBinomialCache();
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
		getCurveInT(get<0>(p), temp_curve_p);
		const double dx = get<1>(p)[0] - temp_curve_p[0];
		error += dx * dx;
		const double dy = get<1>(p)[1] - temp_curve_p[1];
		error += dy * dy;
	}
	return error;
}

void BezierCurve::calcBinomialCache() {
	std::function<long(int)> factorial;
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
}