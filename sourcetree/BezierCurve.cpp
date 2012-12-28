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

BezierCurve::BezierCurve(const uint32_t n_control_points) :
		kNumberControlPoints_(n_control_points), control_points_(
				n_control_points) {

}

BezierCurve::~BezierCurve() {
}

void BezierCurve::GetCurveInT(const double parameterization_value, Vec2& out) {
	out.x = parameterization_value;
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
