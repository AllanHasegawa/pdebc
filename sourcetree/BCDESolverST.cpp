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

BCDESolverST::BCDESolverST(const int process, const int population,
                           BCDESolver& solver)
    : solver_(solver),
      bezier_curve_(solver.bezier_curve_.kNumberControlPoints_),
      kProcess_(process),
      dt_real_0_1_(0, 1),
      dt_int_0_2_(0, 1),
      pending_work_(false),
      finish_(false),
      work_ready_(false),
      kPopulation_(population) {

  std::random_device rd;
  engine_.seed(rd());
  engine0_.seed(rd());
  engine1_.seed(rd());

  dt_int_pop_interval_from_zero_ = std::uniform_int_distribution<int>(
      0, kPopulation_ - 1);

  BezierCurve& bc = solver.bezier_curve_;
  for (int i = 0; i < bc.kNumberControlPoints_; i++) {
    bezier_curve_.control_points_[i].x = bc.control_points_[i].x;
    bezier_curve_.control_points_[i].y = bc.control_points_[i].y;
  }
  bezier_curve_.SetMadOptimizationCaching(solver_.data_points_,
                                          solver_.parameterization_values_);

  parameters_.resize(bezier_curve_.kNumberControlPoints_ - 2);
  for (auto i = parameters_.begin(), e = parameters_.end(); i != e; i++) {
    i->resize(kPopulation_);
  }

  population_errors_.resize(bezier_curve_.kNumberControlPoints_ - 2);
  for (auto i = population_errors_.begin(), e = population_errors_.end();
      i != e; i++) {
    i->resize(kPopulation_);
  }
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

void BCDESolverST::DoWork(const int control_point, const Vec2& error_before) {
  using namespace std;
  lock_guard<mutex> lock(mutex_);
  control_point_ = control_point;
  pending_work_ = true;
  work_ready_ = false;
  error_before_.x = error_before.x;
  error_before_.y = error_before.y;
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

  Initialize();

  /*const std::vector<double>& parameterization_values = solver_
   .parameterization_values_;
   const std::vector<Vec2>& data_points = solver_.data_points_;*/
  const BezierCurve& bc = solver_.bezier_curve_;
  r_rs_ = 0;
  GenerateRandom();
  Vec2 error_before;
  Vec2 lowest_error_index;
  //int calls = 0;
  while (!finish_) {

    unique_lock<mutex> lock(mutex_);
    cond_.wait(lock, [&]() {return this->pending_work_;});
    //printf("Work! %d\n", work++);

    lock.unlock();
    if (finish_) {
      continue;
    }
    std::vector<Vec2>& parameters = parameters_[control_point_];
    std::vector<Vec2>& population_errors = population_errors_[control_point_];

    solver_.random_population_interval_[kProcess_] =
        dt_int_pop_interval_from_zero_(engine_);

    // update local BezierCurve
    {
      //for (int i = 1; i < bc.kNumberControlPoints_ - 1; i++) {
      int cp = control_point_ - 1;
      if (cp == -1) {
        cp = bezier_curve_.kNumberControlPoints_ - 2 - 1;
      }
      bezier_curve_.control_points_[cp + 1].x = bc.control_points_[cp + 1].x;
      bezier_curve_.control_points_[cp + 1].y = bc.control_points_[cp + 1].y;
      //}
    }

    lowest_error_index.x = 0;
    lowest_error_index.y = 0;

    Vec2 trials;
    error_before.x = error_before_.x;
    error_before.y = error_before_.y;
    // just before pop processing, let's cache the BC control points
    /*bezier_curve_.UpdateVariableCPForMadOptimizationCaching(
     solver_.parameterization_values_, control_point_ + 1);*/

    //printf("Calls: %d\n", calls++);
    Vec2 error_now;
    for (int k = 0; k < kPopulation_; k++) {
      Mutate(k, trials);
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

  // Callgrind really hated the use of "dt_int_pop_interval"

  std::vector<Vec2>& p = parameters_[control_point_];
  const Vec2& tc1 = p[random_rs_[r_rs_][0]];
  const Vec2& tc2 = p[random_rs_[r_rs_][1]];
  const Vec2& tc3 = p[random_rs_[r_rs_][2]];
  int j = random_j_[r_rs_];
  for (int k = 1; k <= D; k++) {
    if (j == 0) {
      if (random_cr_[r_rs_][0] <= solver_.kDE_CR_ || k == D) {  //
        trials.x = tc1.x + solver_.kDE_F_ * (tc2.x - tc3.x);

      } else {
        trials.x = p[actual_index].x;
      }
    } else if (j == 1) {
      if (random_cr_[r_rs_][1] <= solver_.kDE_CR_ || k == D) {
        trials.y = tc1.y + solver_.kDE_F_ * (tc2.y - tc3.y);

      } else {
        trials.y = p[actual_index].y;
      }
    }
    j = (j + 1) % D;
  }
  r_rs_++;
}

void BCDESolverST::Select(const Vec2& trial, const Vec2& error_before,
                          Vec2& error_new) {
  const double tx = bezier_curve_.control_points_[control_point_ + 1].x;
  const double ty = bezier_curve_.control_points_[control_point_ + 1].y;

  bezier_curve_.control_points_[control_point_ + 1].x = trial.x;
  bezier_curve_.control_points_[control_point_ + 1].y = trial.y;

  bezier_curve_.CalcErrorWithMadOptimizationCaching(error_new);

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
  return (rand() % max) + min;
}

void BCDESolverST::Initialize() {
  const int bcn = bezier_curve_.kNumberControlPoints_;

  const int w = 128;
  const int h = 1024;

  std::uniform_real_distribution<double> distribution(-w, w);
  std::random_device rd;
  std::mt19937 engine(rd());  // Mersenne twister MT19937
  auto generator = std::bind(distribution, engine);

  // (bc-2): We ignore the first and last control points :3
  // random population for each control point...
  for (int i = 0; i < bcn - 2; i++) {
    for (int j = 0; j < kPopulation_; j++) {
      parameters_[i][j].x = generator();
      parameters_[i][j].y = generator();
    }
  }
}

void BCDESolverST::GenerateRandom() {
  const int CP = bezier_curve_.kNumberControlPoints_ - 2;
  // Hard coded generations <3 ;)
  random_rs_.resize(kPopulation_ * CP * 500);
  random_j_.resize(kPopulation_ * CP * 500);
  random_cr_.resize(kPopulation_ * CP * 500);
  for (int i = 0; i < kPopulation_ * CP * 500; i++) {
    const int actual_index = i % kPopulation_;
    int r1 = dt_int_pop_interval_from_zero_(engine_);  //pop_start_ + 1;  //  //
    int r2 = dt_int_pop_interval_from_zero_(engine0_);  //pop_start_ + 2;  //  //RN(pop_start_, pop_end_);
    int r3 = dt_int_pop_interval_from_zero_(engine1_);  //pop_start_ + 3;  //  //RN(pop_start_, pop_end_);
    while (r1 == actual_index) {
      r1 = dt_int_pop_interval_from_zero_(engine_);
    }
    while (r2 == r1 || r2 == actual_index) {
      r2 = dt_int_pop_interval_from_zero_(engine0_);
    }
    while (r3 == r2 || r3 == r1 || r3 == actual_index) {
      r3 = dt_int_pop_interval_from_zero_(engine1_);
    }
    std::array<int, 3> a = { r1, r2, r3 };
    //random_rs_.push_back(a);
    random_rs_[i] = a;
    //random_j_.push_back(dt_int_0_2_(engine_));
    random_j_[i] = dt_int_0_2_(engine_);
    std::array<float, 2> b = { (float) dt_real_0_1_(engine_),
        (float) dt_real_0_1_(engine_) };
    //random_cr_.push_back(b);
    random_cr_[i] = b;
  }
}
