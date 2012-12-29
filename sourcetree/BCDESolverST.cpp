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
#include <vector>
#include <array>

#include "BCDESolver.h"
#include "BezierCurve.h"
#include "Vec2.h"

BCDESolverST::BCDESolverST(const int control_point, const int process,
                           BCDESolver& solver)
    : solver_(solver),
      bezier_curve_(solver.bezier_curve_.kNumberControlPoints_),
      kCP_(control_point),
      kProcess_(process),
      dt_real_0_1_(0, 1),
      dt_int_0_2_(0, 1) {

  std::random_device rd;
  engine_.seed(rd());
  const int pop_interval = (solver_.kPopulation_ / solver_.kNProcess_);
  const int pop_start = kProcess_ * pop_interval;
  const int pop_end = pop_start + pop_interval;

  dt_int_pop_interval_ = std::uniform_int_distribution<int>(pop_start,
                                                            pop_end - 1);

  BezierCurve& bc = solver.bezier_curve_;
  for (int i = 0; i < bc.kNumberControlPoints_; i++) {
    bezier_curve_.control_points_[i].x = bc.control_points_[i].x;
    bezier_curve_.control_points_[i].y = bc.control_points_[i].y;
  }
}

BCDESolverST::~BCDESolverST() {
}

void BCDESolverST::Start() {
  thread_ = std::thread(&BCDESolverST::Run, this);
}

void BCDESolverST::Join() {
  thread_.join();
}

void BCDESolverST::Run() {
  Vec2 error_before;
  Vec2 lowest_error_index;
  lowest_error_index.x = 0;
  lowest_error_index.y = 0;

  bezier_curve_.CalcError(solver_.parameterization_values_,
                          solver_.data_points_, error_before);

  const int pop_interval = (solver_.kPopulation_ / solver_.kNProcess_);
  const int pop_start = kProcess_ * pop_interval;
  const int pop_end = pop_start + pop_interval;

  Vec2 trials;

  for (int k = pop_start; k < pop_end; k++) {
    Mutate(k, trials);
    Vec2 error_now;
    Select(trials, error_before, error_now);

    if (error_now.x <= error_before.x) {
      solver_.parameters_[kCP_][k].x = trials.x;
      error_before.x = error_now.x;
    }
    if (error_now.y <= error_before.y) {
      solver_.parameters_[kCP_][k].y = trials.y;
      error_before.y = error_now.y;
    }

    solver_.population_errors_[kCP_][k].x = error_before.x;
    solver_.population_errors_[kCP_][k].y = error_before.y;

    if (error_before.x
        < solver_.population_errors_[kCP_][lowest_error_index.x].x) {
      lowest_error_index.x = k;
    }
    if (error_before.y
        < solver_.population_errors_[kCP_][lowest_error_index.y].y) {
      lowest_error_index.y = k;
    }
  }
  solver_.lowest_error_index_[kProcess_].x = lowest_error_index.x;
  solver_.lowest_error_index_[kProcess_].y = lowest_error_index.y;
}

void BCDESolverST::Mutate(const int actual_index, Vec2& trials) {
  const int D = 2;
  // final int NP = BCDESolver.N_POPULATION;

  const int pop_interval = (solver_.kPopulation_ / solver_.kNProcess_);
  const int pop_start = kProcess_ * pop_interval;
  const int pop_end = pop_start + pop_interval;

  int r1 = dt_int_pop_interval_(engine_);
  while (r1 == actual_index) {
    r1 = dt_int_pop_interval_(engine_);
  }
  int r2 = dt_int_pop_interval_(engine_);
  while (r2 == r1 || r2 == actual_index) {
    r2 = dt_int_pop_interval_(engine_);
  }
  int r3 = dt_int_pop_interval_(engine_);
  while (r3 == r2 || r3 == r1 || r3 == actual_index) {
    r3 = dt_int_pop_interval_(engine_);
  }

  int j = dt_int_0_2_(engine_);
  std::vector<Vec2>& p = solver_.parameters_[kCP_];

  for (int k = 1; k <= D; k++) {
    if (j == 0) {
      if (dt_real_0_1_(engine_) <= solver_.kDE_CR_ || k == D) {
        trials.x = p[r1].x + solver_.kDE_F_ * (p[r2].x - p[r3].x);

      } else {
        trials.x = p[actual_index].x;
      }
    } else if (j == 1) {
      if (dt_real_0_1_(engine_) <= solver_.kDE_CR_ || k == D) {
        trials.y = p[r1].y + solver_.kDE_F_ * (p[r2].y - p[r3].y);

      } else {
        trials.y = p[actual_index].y;
      }
    }
    j = (j + 1) % D;
  }
}

void BCDESolverST::Select(const Vec2& trial, const Vec2& error_before,
                          Vec2& error_new) {
  double tx = bezier_curve_.control_points_[kCP_ + 1].x;
  double ty = bezier_curve_.control_points_[kCP_ + 1].y;

  bezier_curve_.control_points_[kCP_ + 1].x = trial.x;
  bezier_curve_.control_points_[kCP_ + 1].y = trial.y;

  bezier_curve_.CalcError(solver_.parameterization_values_,
                          solver_.data_points_, error_new);

  if (error_new.x > error_before.x) {
    bezier_curve_.control_points_[kCP_ + 1].x = tx;
  }

  if (error_new.y > error_before.y) {
    bezier_curve_.control_points_[kCP_ + 1].y = ty;
  }
}
