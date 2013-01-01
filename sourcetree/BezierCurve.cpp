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
#include <math.h>

#include "Globals.h"

std::vector<Vec2> BezierCurve::const_control_point_ = std::vector<Vec2>();
std::vector<std::vector<double>> BezierCurve::b_caching_ = std::vector<
    std::vector<double>>();
int BezierCurve::variable_control_point_ = 1;

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
    const double p1 = pow(parameterization_value, i);
    const double p2 = pow(1 - parameterization_value, n - i);
    /*
     double p1 = 1;
     for (int j = 0; j < i; j++) {
     p1 *= parameterization_value;
     }
     double p2 = 1;
     for (int j = 0, je = n - i; j < je; j++) {
     p2 *= (1 - parameterization_value);
     }*/
    const double b = Globals::binomial_cache_[n][i] * p1 * p2;
    const Vec2& v = control_points_[i];
    Bx += b * v.x;
    By += b * v.y;
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

void BezierCurve::UpdateVariableCPForMadOptimizationCaching(
    const std::vector<double>& parameterization_values,
    const int variable_control_point) {
  // Then I cache the control points that will remain const
  variable_control_point_ = variable_control_point;
  const int np = parameterization_values.size();
  const_control_point_.resize(np);
  for (int p = 0; p < np; p++) {
    double Bx = 0;
    double By = 0;
    for (int i = 0; i < kNumberControlPoints_; i++) {
      if (i == variable_control_point) {
        continue;
      }
      const double b = b_caching_[p][i];
      const Vec2& v = control_points_[i];
      Bx += b * v.x;
      By += b * v.y;
    }
    const_control_point_[p].x = Bx;
    const_control_point_[p].y = By;
  }
}
void BezierCurve::SetMadOptimizationCaching(
    const std::vector<Vec2>& data_points,
    const std::vector<double>& parameterization_values) {
  local_data_points_.resize(data_points.size());
  for (int i = 0; i < data_points.size(); i++) {
    local_data_points_[i].x = data_points[i].x;
    local_data_points_[i].y = data_points[i].y;
  }

  // First I optimize the p1 and p2 variables
  std::vector<std::vector<double>> p1_p2_caching;
  p1_p2_caching.resize(parameterization_values.size());
  for (auto i = p1_p2_caching.begin(); i != p1_p2_caching.end(); i++) {
    i->resize(kNumberControlPoints_);
  }

  for (int i = 0; i < parameterization_values.size(); i++) {
    for (int j = 0; j < kNumberControlPoints_; j++) {
      const double p1 = pow(parameterization_values[i], j);
      const double p2 = pow(1 - parameterization_values[i],
                            kNumberControlPoints_ - 1 - j);
      p1_p2_caching[i][j] = p1 * p2;
    }
  }

  //then the binomial ;)
  // calc local binomial cache
  const int n = kNumberControlPoints_ - 1;
  std::function<long(int)> factorial;
  factorial = [&factorial](int f) -> double {
    if (f > 0) {
      return f*factorial(f-1);
    } else {
      return (long)1;
    }
  };
  for (int i = 0; i < Globals::kMaxControlPoints; i++) {
    local_binomial_[i] = factorial(n) / (factorial(i) * factorial(n - i));
  }

  b_caching_.resize(parameterization_values.size());
  for (auto i = b_caching_.begin(); i != b_caching_.end(); i++) {
    i->resize(kNumberControlPoints_);
  }
  for (int p = 0; p < parameterization_values.size(); p++) {
    for (int i = 0; i < kNumberControlPoints_; i++) {
      b_caching_[p][i] = local_binomial_[i] * p1_p2_caching[p][i];
    }
  }
}

void BezierCurve::CalcErrorWithMadOptimizationCaching(Vec2& error) {
  const int DP = local_data_points_.size();

  double ex = 0;
  double ey = 0;
  for (int k = 1; k < DP - 1; k++) {
    GetCurveInTWithMadOptimizationCaching(k, temp_curve_p_);

    const double dx = local_data_points_[k].x - temp_curve_p_.x;
    const double dy = local_data_points_[k].y - temp_curve_p_.y;
    ex += dx * dx;
    ey += dy * dy;
  }
  error.x = ex;
  error.y = ey;
}

void BezierCurve::GetCurveInTWithMadOptimizationCaching(const int para_index,
                                                        Vec2& out) {
  const Vec2& v = control_points_[variable_control_point_];
  out.x = v.x * b_caching_[para_index][variable_control_point_]
      + const_control_point_[para_index].x;
  out.y = v.y * b_caching_[para_index][variable_control_point_]
      + const_control_point_[para_index].y;
  /*
   double Bx = 0;
   double By = 0;
   for (int i = 0; i < kNumberControlPoints_; i++) {
   //const double p1 = p1_caching_[para_index][i];
   //const double p2 = p2_caching_[para_index][i];
   const double b = b_caching_[para_index][i];//local_binomial_[i] * p1_p2_caching_[para_index][i];
   const Vec2& v = control_points_[i];
   Bx += b * v.x;
   By += b * v.y;
   }
   out.x = Bx;
   out.y = By;*/
}

