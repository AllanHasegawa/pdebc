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

#include "BezierCurve.h"
#include <stdio.h>
#include <string>

#include "Globals.h"

BezierCurve::BezierCurve(const uint32_t n_control_points)
    : kNumberControlPoints_(n_control_points),
      control_points_(n_control_points) {

}

BezierCurve::~BezierCurve() {
}

void BezierCurve::GetCurveInT(const double parameterization_value, Vec2& out) {
  int n = kNumberControlPoints_ - 1;
  double Bx = 0;
  double By = 0;
  for (int i = 0; i < kNumberControlPoints_; i++) {
    double p1 = 1;
    for (int j = 0; j < i; j++) {
      p1 *= parameterization_value;
    }
    double p2 = 1;
    for (int j = 0; j < n - i; j++) {
      p2 *= (1 - parameterization_value);
    }
    const double b = Globals::binomial_cache_[n][i] * p1 * p2;
    Bx += b * control_points_[i].x;
    By += b * control_points_[i].y;
  }
  out.x = Bx;
  out.y = By;
}

void BezierCurve::CalcError(const std::vector<double>& parameterization_values,
                            const std::vector<Vec2>& data_points, Vec2& error) {

  const int DP = data_points.size();
  error.x = error.y = 0;

  for (int k = 1; k < DP - 1; k++) {
    GetCurveInT(parameterization_values[k], temp_curve_p_);
    const double dx = data_points[k].x - temp_curve_p_.x;
    error.x += dx * dx;
    const double dy = data_points[k].y - temp_curve_p_.y;
    error.y += dy * dy;
  }
}

void BezierCurve::PrintControlPoints() {
  printf("%d\n", (int) control_points_.size());
  for (auto i = control_points_.begin(); i != control_points_.end(); i++) {
    printf("%f %f\n", i->x, i->y);
  }
}

std::string BezierCurve::SaveAsSVGPoints(const uint32_t interpolation) {
  using namespace std;
  string str = "<g stroke=\"black\" stroke-width=\"1\" fill=\"none\">\n";
  Vec2 p;
  str += "<path id=\"path_bc\" d=\"";

  GetCurveInT(0, p);

  str += "M " + to_string(p.x) + " " + to_string(p.y) + " ";

  for (double i = 1.0 / interpolation; i <= 1.0; i += 1.0 / interpolation) {
    GetCurveInT(i, p);
    str += "L " + to_string(p.x) + " " + to_string(p.y) + " ";
  }
  str += "\" />\n";
  str += "</g>";
  return str;
}
