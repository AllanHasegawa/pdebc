#include <stdio.h>

#include "external/tinyxml2/tinyxml2.h"

#include <vector>
#include "Globals.h"
#include "BezierCurve.h"
#include "SVGParser.h"

int main(void) {
	using namespace std;

	Globals::CalcBinomial();
	printf("Hello World\n");
	BezierCurve<10, 2> bc;
	std::array<double, 2> point;
	bc.GetCurveInT(0.5, point);

	printf("%f\n", Globals::binomial_cache_[11][10]);

	printf("Loading XML\n");
	vector<array<double,2>> control_points;
	SVGParser::GetControlPoints("c6.svg", control_points);
	printf("N CPs: %d\n", (int)control_points.size());

	return 0;
}
