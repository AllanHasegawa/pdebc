#include <stdio.h>

#include "external/tinyxml2/tinyxml2.h"

#include <vector>
#include "Globals.h"
#include "BezierCurve.h"
#include "SVGParser.h"
#include "BCDESolver.h"
#include "Vec2.h"

int main(void) {
	using namespace std;

	Globals::CalcBinomial();
	printf("%f\n", Globals::binomial_cache_[11][10]);

	printf("Loading XML\n");
	vector<Vec2> data_points;
	SVGParser::GetDataPoints("c6.svg", data_points);
	printf("N CPs: %d\n", (int) data_points.size());

	vector<double> chord_length;
	Globals::CalcChordLength(data_points, chord_length);

	BezierCurve bc(5);
	BCDESolver de(chord_length, data_points, 1, 0.5, 0.8, 128, bc);


	return 0;
}
