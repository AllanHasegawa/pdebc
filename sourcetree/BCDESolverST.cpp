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

#include "BCDESolverST.h"

#include <thread>
#include <stdio.h>


#include "BezierCurve.h"
#include "Vec2.h"

BCDESolverST::BCDESolverST(const BCDESolver& solver, const int control_point,
		const int process) :
		solver_(solver), bezier_curve_(
				solver.bezier_curve_.kNumberControlPoints_), kCP_(
				control_point), kProcess_(process) {
	printf("%d\n", bezier_curve_.control_points_.size());
	BezierCurve& bc = solver.bezier_curve_;
	for (int i = 0; i < bc.kNumberControlPoints_; i++) {
		bezier_curve_.control_points_[i].x = bc.control_points_[i].x;
		bezier_curve_.control_points_[i].y = bc.control_points_[i].y;
	}
}


BCDESolverST::~BCDESolverST() {
	// TODO Auto-generated destructor stub
}

void BCDESolverST::Start() {
	printf("Starting Thread\n");
	thread_ = std::thread(&BCDESolverST::Run, this);
}

void BCDESolverST::Join() {
	thread_.join();
}


void BCDESolverST::Run() {
	printf("Running %d\n", 3);
}
