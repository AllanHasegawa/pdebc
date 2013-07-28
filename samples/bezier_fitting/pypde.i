%module pypde
%{
#include "pypde.hpp"
%}

%include "std_vector.i"
namespace std {
   %template(vectord) vector<double>;
   %template(vectorv) vector<Vec2>;
}

//%include "pypde.hpp"
struct Vec2 {
	double x;
	double y;

	Vec2();
	Vec2(double, double);
	~Vec2();
};

struct pypde {
	pypde(const int n_processes, const int population_size,
		const int bezier_control_points,
		std::vector<Vec2> data_points);
	~pypde();
	void solveOneGeneration();
	double getBestCandidateError(int i);
	std::vector<double> getBestCandidateCP(int i);
};