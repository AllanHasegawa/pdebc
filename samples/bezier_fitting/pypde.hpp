
#include <vector>
#include <memory>

#include "pdebc/ThreadsDE.hpp"
#include "BezierCurve.hpp"

struct Vec2 {
	double x;
	double y;

	Vec2() : x{}, y{} {

	}
	Vec2(double x, double y) : x{x}, y{y} {
	}
	~Vec2() {

	}
};

#define POPULATION_SIZE 128
#define POPULATION_DIM 2
#define POPULATION_TYPE double
#define ERROR_TYPE double
#define DOMAIN_LIMITS static_cast<double>(64)

using PYPDE_ThreadsDE = pdebc::ThreadsDE<POPULATION_TYPE,POPULATION_DIM,ERROR_TYPE>;

struct pypde {
	BezierCurve* bezier_curve_;
	std::vector<std::shared_ptr<PYPDE_ThreadsDE>> des_;

	pypde(const int n_processes, const int population_size,
		const int bezier_control_points,
		std::vector<Vec2> data_points);

	~pypde();

	void solveOneGeneration();

	double getBestCandidateError(int i);

	std::vector<double> getBestCandidateCP(int i);
};