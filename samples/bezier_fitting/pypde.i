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
};

class pypde {
	pypde(const int bezier_control_points,
		std::vector<Vec2> data_points);
	~pypde();
	void solveOneGeneration();
	double getBestCandidateError(int i);
	std::vector<double> getBestCandidateCP(int i);
};