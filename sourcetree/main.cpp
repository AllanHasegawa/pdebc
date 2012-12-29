#include <stdio.h>

#include "external/tinyxml2/tinyxml2.h"

#include <random>
#include <vector>
#include "Globals.h"
#include "BezierCurve.h"
#include "SVGParser.h"
#include "BCDESolver.h"
#include "Vec2.h"

int main(void) {
  using namespace std;

  Globals::CalcBinomial();
  //printf("%f\n", Globals::binomial_cache_[11][10]);

  //printf("Loading XML\n");
  vector<Vec2> data_points;
  SVGParser::GetDataPoints("g_p_1", data_points);
  //printf("N DPs: %d\n", (int) data_points.size());

  vector<double> chord_length;
  Globals::CalcChordLength(data_points, chord_length);

  /*for (auto i = chord_length.begin(); i != chord_length.end(); i++) {
   printf("%lf\n", *i);
   }
   for (auto i = data_points.begin(); i != data_points.end(); i++) {
   printf("DP; %f,%f\n", i->x, i->y);
   }*/
  BezierCurve bc(9);
  BCDESolver de(chord_length, data_points, 1, 0.8, 0.5, 128, bc);

  while (de.generation_ < 500) {
    de.SolveOneGeneration();
    Vec2 error;
    bc.CalcError(chord_length, data_points, error);
    //printf("%lf,%lf\n", error.x, error.y);
  }

  //bc.PrintControlPoints();

  string save_content = bc.SaveAsSVGPoints(128);

  printf("%s\n", save_content.c_str());
  return 0;
}
