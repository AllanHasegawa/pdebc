

#ifndef BEZIERCURVE_HPP_
#define BEZIERCURVE_HPP_

#include <array>
#include <tuple>
#include <vector>

using Vec2d = std::array<double,2>;

struct BezierCurve {

	static constexpr int kMaxControlPoints{20};
	std::array<std::array<double, kMaxControlPoints>, kMaxControlPoints> binomial_cache_;

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


	/* Optimization Cache */
	void updateVariableCPForOptimizationCache(const int variable_control_point);
	void getCurveInTWithOptimizationCache(const int para_index,
		const Vec2d& candidate_cp, Vec2d& out);
	double calcErrorWithOptimizationCache(const Vec2d& candidate_cp);

private:
	/* Optimization Cache */
	uint32_t variable_control_point_;
	std::vector<std::vector<double>> b_caching_;
	std::vector<Vec2d> const_control_point_;


	void initializeOptimizationCache();
};

#endif /* BEZIERCURVE_HPP_ */