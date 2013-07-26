
#include <vector>

#include "pdebc/SequentialDE.hpp"
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
#define DOMAIN_LIMITS static_cast<double>(128)

#define SequentialDE pdebc::SequentialDE<POPULATION_TYPE,POPULATION_DIM,POPULATION_SIZE,ERROR_TYPE>

struct pypde {
	BezierCurve* bezier_curve_;
	std::vector<SequentialDE> des_;

	pypde(const int bezier_control_points,
		std::vector<Vec2> data_points);

	~pypde();

	void solveOneGeneration();

	double getBestCandidateError(int i);

	std::vector<double> getBestCandidateCP(int i);
};