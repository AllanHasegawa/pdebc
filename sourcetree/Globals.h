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

private:
	Globals() {}
};

std::array<std::array<double, Globals::kMaxControlPoints>,
		Globals::kMaxControlPoints> Globals::binomial_cache_ = std::array<
		std::array<double, Globals::kMaxControlPoints>,
		Globals::kMaxControlPoints>();

#endif /* GLOBALS_H_ */
