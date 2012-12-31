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
#include <random>
#include <cstdlib>

#include "BCDESolver.h"
#include "BezierCurve.h"
#include "Vec2.h"

BCDESolverST::BCDESolverST(const int process, BCDESolver& solver)
    : solver_(solver),
      bezier_curve_(solver.bezier_curve_.kNumberControlPoints_),
      kProcess_(process),
      dt_real_0_1_(0, 1),
      dt_int_0_2_(0, 1),
      pending_work_(false),
      finish_(false),
      work_ready_(false) {

  std::random_device rd;
  engine_.seed(rd());
  const int pop_interval = (solver_.kPopulation_ / solver_.kNProcess_);
  const int pop_start = kProcess_ * pop_interval;
  const int pop_end = pop_start + pop_interval;

  pop_start_ = pop_start;
  pop_end_ = pop_end;
  dt_int_pop_interval_from_zero_ = std::uniform_int_distribution<int>(
      0, pop_interval - 1);
  dt_int_pop_interval_ = std::uniform_int_distribution<int>(pop_start,
                                                            pop_end - 1);

  BezierCurve& bc = solver.bezier_curve_;
  for (int i = 0; i < bc.kNumberControlPoints_; i++) {
    bezier_curve_.control_points_[i].x = bc.control_points_[i].x;
    bezier_curve_.control_points_[i].y = bc.control_points_[i].y;
  }
  bezier_curve_.SetMadOptimizationCaching(solver_.parameterization_values_);
  bezier_curve_.UpdateVariableCPForMadOptimizationCaching(
      solver_.parameterization_values_, 1);
  //printf("Starting process %d\n", kProcess_);
}

BCDESolverST::~BCDESolverST() {
}

void BCDESolverST::Start() {
  //printf("Starting process %d\n", kProcess_);
  thread_ = std::thread(&BCDESolverST::Run, this);
}

void BCDESolverST::Join() {
  std::unique_lock<std::mutex> lock(mutex_);
  pending_work_ = true;
  finish_ = true;
  cond_.notify_one();
  lock.unlock();
  thread_.join();
}

void BCDESolverST::DoWork(const int control_point) {
  using namespace std;
  lock_guard<mutex> lock(mutex_);
  control_point_ = control_point;
  pending_work_ = true;
  work_ready_ = false;
  cond_.notify_one();
}

void BCDESolverST::WaitWork() {
  using namespace std;
  unique_lock<mutex> lock(work_ready_lock_);
  work_ready_cond_.wait(lock, [this]() {return this->work_ready_;});
}

void BCDESolverST::Run() {
  using namespace std;
  int work = 0;

  const int pop_interval = (solver_.kPopulation_ / solver_.kNProcess_);
  const int pop_start = kProcess_ * pop_interval;
  const int pop_end = pop_start + pop_interval;

  const std::vector<double>& parameterization_values = solver_
      .parameterization_values_;
  const std::vector<Vec2>& data_points = solver_.data_points_;
  const BezierCurve& bc = solver_.bezier_curve_;

  Vec2 error_before;
  Vec2 lowest_error_index;

  while (!finish_) {

    unique_lock<mutex> lock(mutex_);
    cond_.wait(lock, [this]() {return this->pending_work_;});
    //printf("Work! %d\n", work++);

    lock.unlock();
    if (finish_) {
      continue;
    }
    std::vector<Vec2>& parameters = solver_.parameters_[control_point_];
    std::vector<Vec2>& population_errors =
        solver_.population_errors_[control_point_];

    solver_.random_population_interval_[kProcess_] =
        dt_int_pop_interval_from_zero_(engine_);

    // update local BezierCurve
    {
      for (int i = 1; i < bc.kNumberControlPoints_ - 1; i++) {
        bezier_curve_.control_points_[i].x = bc.control_points_[i].x;
        bezier_curve_.control_points_[i].y = bc.control_points_[i].y;
      }
    }

    lowest_error_index.x = 0;
    lowest_error_index.y = 0;

    Vec2 trials;
    // just before pop processing, let's cache the BC control points
    bezier_curve_.UpdateVariableCPForMadOptimizationCaching(
        solver_.parameterization_values_, control_point_ + 1);
    bezier_curve_.CalcErrorWithMadOptimizationCaching(solver_.data_points_,
                                                      error_before);

    for (int k = pop_start; k < pop_end; k++) {
      Mutate(k, trials);
      Vec2 error_now;
      Select(trials, error_before, error_now);

      if (error_now.x <= error_before.x) {
        parameters[k].x = trials.x;
        error_before.x = error_now.x;
      }
      if (error_now.y <= error_before.y) {
        parameters[k].y = trials.y;
        error_before.y = error_now.y;
      }

      population_errors[k].x = error_before.x;
      population_errors[k].y = error_before.y;

      if (error_before.x < population_errors[lowest_error_index.x].x) {
        lowest_error_index.x = k;
      }
      if (error_before.y < population_errors[lowest_error_index.y].y) {
        lowest_error_index.y = k;
      }
    }
    solver_.lowest_error_index_[kProcess_].x = lowest_error_index.x;
    solver_.lowest_error_index_[kProcess_].y = lowest_error_index.y;
    pending_work_ = false;

    lock_guard<mutex> work_read_lock(work_ready_lock_);
    work_ready_ = true;
    work_ready_cond_.notify_one();
  }
}

void BCDESolverST::Mutate(const int actual_index, Vec2& trials) {
  const int D = 2;
  // final int NP = BCDESolver.N_POPULATION;

  // Callgrind really hated the use of "dt_int_pop_interval"

  int r1 = RN(pop_start_, pop_end_);
  while (r1 == actual_index) {
    r1 = RN(pop_start_, pop_end_);
  }
  int r2 = RN(pop_start_, pop_end_);
  while (r2 == r1 || r2 == actual_index) {
    r2 = RN(pop_start_, pop_end_);
  }
  int r3 = RN(pop_start_, pop_end_);
  while (r3 == r2 || r3 == r1 || r3 == actual_index) {
    r3 = RN(pop_start_, pop_end_);
  }

  int j = rand() % 2;
  std::vector<Vec2>& p = solver_.parameters_[control_point_];

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
  double tx = bezier_curve_.control_points_[control_point_ + 1].x;
  double ty = bezier_curve_.control_points_[control_point_ + 1].y;

  bezier_curve_.control_points_[control_point_ + 1].x = trial.x;
  bezier_curve_.control_points_[control_point_ + 1].y = trial.y;

  bezier_curve_.CalcErrorWithMadOptimizationCaching(solver_.data_points_,
                                                    error_new);

  if (error_new.x > error_before.x) {
    bezier_curve_.control_points_[control_point_ + 1].x = tx;
  }

  if (error_new.y > error_before.y) {
    bezier_curve_.control_points_[control_point_ + 1].y = ty;
  }
}

uint32_t BCDESolverST::RN(uint32_t min, uint32_t max) {
  /*double scaled = (double) rand() / RAND_MAX;
  return (max - min + 1) * scaled + min;*/
  return ( rand() % max ) + min;
}
