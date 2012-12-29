#include <stdio.h>

#include "external/tinyxml2/tinyxml2.h"

#include <random>
#include <vector>
#include "Globals.h"
#include "BezierCurve.h"
#include "SVGParser.h"
#include "BCDESolver.h"
#include "Vec2.h"

void PrintUsage();

int main(int argc, char *argv[]) {
  using namespace std;

  Globals::CalcBinomial();

  bool verify_input[7];
  for (int i = 0; i < 7; i++) {
    verify_input[i] = false;
  }
  char* data_points_file;
  int n_processes;
  double de_f;
  double de_cr;
  int n_generations;
  int n_control_points;
  int n_population;

  if (argc < 2) {
    PrintUsage();
    return 0;
  }

  for (int i = 1; i < argc; i++) {
    if (i + 1 < argc) {
      if (strcmp(argv[i], "-d") == 0) {
        data_points_file = argv[i + 1];
        verify_input[0] = true;
        i++;
      } else if (strcmp(argv[i], "-b") == 0) {
        n_control_points = atoi(argv[i + 1]);
        i++;
        verify_input[1] = true;
      } else if (strcmp(argv[i], "-p") == 0) {
        n_processes = atoi(argv[i + 1]);
        i++;
        verify_input[2] = true;
      } else if (strcmp(argv[i], "-g") == 0) {
        n_generations = atoi(argv[i + 1]);
        i++;
        verify_input[3] = true;
      } else if (strcmp(argv[i], "-n") == 0) {
        n_population = atoi(argv[i + 1]);
        i++;
        verify_input[4] = true;
      } else if (strcmp(argv[i], "-f") == 0) {
        de_f = atof(argv[i + 1]);
        i++;
        verify_input[5] = true;
      } else if (strcmp(argv[i], "-c") == 0) {
        de_cr = atof(argv[i + 1]);
        i++;
        verify_input[6] = true;
      }
    }
  }

  for (int i = 0; i < 7; i++) {
    if (!verify_input[i]) {
      PrintUsage();
      return 0;
    }
  }

  vector<Vec2> data_points;
  SVGParser::GetDataPoints(data_points_file, data_points);

  vector<double> chord_length;
  Globals::CalcChordLength(data_points, chord_length);

  BezierCurve bc(n_control_points);
  BCDESolver de(chord_length, data_points, n_processes, de_f, de_cr,
                n_population, bc);

  while (de.generation_ < n_generations) {
    de.SolveOneGeneration();
  }

  string save_content = bc.SaveAsSVGPoints(128);

  printf("%s\n", save_content.c_str());
  return 0;
}

void PrintUsage() {
  printf(
      "pdebc -d <data_points_file> -b <number control points> -p <number_processes> -g <number generations> -n <number population> -f <DE_F> -c <DE_CR>\n");
}
