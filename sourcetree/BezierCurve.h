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

#ifndef BEZIERCURVE_H_
#define BEZIERCURVE_H_

#include <cstdint>
#include <vector>
#include <string>
#include "Vec2.h"

class BezierCurve {
public:
	const uint32_t kNumberControlPoints_;
	std::vector<Vec2> control_points_;

	BezierCurve(const uint32_t n_control_points);
	virtual ~BezierCurve();

	void GetCurveInT(const Vec2& parameterization_value, Vec2& out);
	void CalcError(const std::vector<Vec2>& parameterization_values,
			const std::vector<Vec2>& data_points, Vec2& error);
	void PrintControlPoints();
	std::string SaveAsSVGPoints(const uint32_t interpolation);

private:
	Vec2 temp_curve_p_;

};

#endif /* BEZIERCURVE_H_ */
