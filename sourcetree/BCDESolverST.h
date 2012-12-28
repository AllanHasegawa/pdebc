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

#ifndef BCDESOLVERST_H_
#define BCDESOLVERST_H_

#include <thread>

#include "BCDESolver.h"
#include "BezierCurve.h"

class BCDESolverST {
public:
	BCDESolverST(const BCDESolver& solver, const int control_point,
			const int process);
	virtual ~BCDESolverST();


	void Start();
	void Join();

private:
	std::thread thread_;
	const BCDESolver& solver_;
	const int kCP_;
	const int kProcess_;

	BezierCurve bezier_curve_;

	void Run();
	void Mutate(const int actual_index);
	void Select(const Vec2& trial, const Vec2& error_before, Vec2& error_new);
};

#endif /* BCDESOLVERST_H_ */
