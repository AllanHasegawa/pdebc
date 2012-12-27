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
#include <array>

template<int CP, int DIM = 2>
class BezierCurve {
public:
	std::array<std::array<double, DIM>, CP> control_points_;

	BezierCurve();
	virtual ~BezierCurve();

	void GetCurveInT(const double parameterization_value,
			std::array<double, DIM>& out);

private:
	// [CP][Dimension]
	std::array<double, DIM> temp_curve_p_;

	template<int DP>
	void CalcError(const std::array<double, DP>& parameterization_values,
			const std::array<std::array<double, DIM>, DP>& data_points,
			std::array<double, DIM>& error);
};

template<int CP, int DIM>
BezierCurve<CP, DIM>::BezierCurve() {

}

template<int CP, int DIM>
BezierCurve<CP, DIM>::~BezierCurve() {
}

template<int CP, int DIM>
void BezierCurve<CP, DIM>::GetCurveInT(const double parameterization_value,
		std::array<double, DIM>& out) {
	out[0] = parameterization_value;
}

template<int CP, int DIM> template<int DP>
void BezierCurve<CP, DIM>::CalcError(
		const std::array<double, DP>& parameterization_values,
		const std::array<std::array<double, DIM>, DP>& data_points,
		std::array<double, DIM>& error) {

	error.fill(0);

	for (int k = 1; k < DP - 1; k++) {
		GetCurveInT(parameterization_values[k], temp_curve_p_);
		for (int i = 0; i < DIM; i++) {
			const double d = data_points[k][i] - temp_curve_p_[i];
			error[i] += d * d;
		}
	}
}

#endif /* BEZIERCURVE_H_ */
