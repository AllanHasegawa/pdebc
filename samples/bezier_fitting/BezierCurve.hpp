

#ifndef BEZIERCURVE_HPP_
#define BEZIERCURVE_HPP_

#include <array>
#include <tuple>
#include <vector>

using Vec2d = std::array<double,2>;

struct BezierCurve {

	static constexpr int kMaxControlPoints{20};
	static std::array<std::array<double, kMaxControlPoints>, kMaxControlPoints> binomial_cache_;

	const uint32_t kNumberControlPoints_;
	const std::vector<std::tuple<double,Vec2d>>& data_points_;
	std::vector<Vec2d> control_points_;


	BezierCurve(
		const std::vector<std::tuple<double,Vec2d>>& data_points,
		const std::vector<Vec2d> control_points
		);

	~BezierCurve();

	void getCurveInT(const double parameterization_value, Vec2d& out) const;
	double calcError() const;

private:
	void calcBinomialCache();
};

#endif /* BEZIERCURVE_HPP_ */