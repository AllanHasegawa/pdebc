#include <stdio.h>

#include "Globals.h"
#include "BezierCurve.h"

int main(void) {
	using namespace std;

	Globals::CalcBinomial();
	printf("Hello World\n");
	BezierCurve<10, 2> bc;
	std::array<double, 2> point;
	bc.GetCurveInT(0.5, point);

	printf("%f\n", Globals::binomial_cache_[11][10]);

	return 0;
}
