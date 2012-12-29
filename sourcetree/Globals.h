/*
 Copyright 2012 Allan Yoshio Hasegawa

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 -----------------------------------------------------------------------------
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_
#include <functional>
#include <array>
#include <math.h>
#include <vector>
#include "Vec2.h"

struct Globals {
  static constexpr int kMaxControlPoints = 20;
  static std::array<std::array<double, kMaxControlPoints>, kMaxControlPoints> binomial_cache_;

  static void CalcBinomial() {
    std::function<long(int)> factorial;
    factorial = [&factorial](int f) -> double {
      if (f > 0) {
        return f*factorial(f-1);
      } else {
        return (long)1;
      }
    };
    for (int n = 0; n < kMaxControlPoints; n++) {
      for (int i = 0; i < kMaxControlPoints; i++) {
        binomial_cache_[n][i] = factorial(n)
            / (factorial(i) * factorial(n - i));
      }
    }
  }

  static void CalcChordLength(const std::vector<Vec2>& data_points,
                              std::vector<double>& chord_length) {

    using namespace std;
    const int DP = data_points.size();
    chord_length.resize(DP);
    chord_length[0] = 0;
    chord_length[DP - 1] = 1;

    double td = 0;

    for (int i = 1; i < DP; i++) {
      const double vdx = data_points[i].x - data_points[i - 1].x;
      const double vdy = data_points[i].y - data_points[i - 1].y;
      const double d = sqrt(vdx * vdx + vdy * vdy);
      td += d;
    }

    double p_td = 0;
    for (int i = 1; i < DP; i++) {
      const double vdx = data_points[i].x - data_points[i - 1].x;
      const double vdy = data_points[i].y - data_points[i - 1].y;
      const double d = sqrt(vdx * vdx + vdy * vdy);
      p_td += d;
      chord_length[i] = p_td / td;
    }
  }

 private:
  Globals() {
  }
};

#endif /* GLOBALS_H_ */
