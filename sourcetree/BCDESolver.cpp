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

#include "BCDESolver.h"

#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <thread>

#include "Vec2.h"
#include "BCDESolverST.h"

BCDESolver::BCDESolver(const std::vector<double>& parameterization_values,
                       const std::vector<Vec2>& data_points,
                       const int n_process, const double de_f,
                       const double de_cr, const int population,
                       BezierCurve& bezier_curve)
    : bezier_curve_(bezier_curve),
      data_points_(data_points),
      parameterization_values_(parameterization_values),
      kNProcess_(n_process),
      kDE_F_(de_f),
      kDE_CR_(de_cr),
      kPopulation_(population),
      generation_(0),
      random_population_interval_(n_process) {

  lowest_error_index_.resize(n_process);

  // place P0 in first data_point
  bezier_curve_.control_points_[0].x = data_points_[0].x;
  bezier_curve_.control_points_[0].y = data_points_[0].y;
  // place last CP in last data_point
  bezier_curve_.control_points_[bezier_curve_.control_points_.size() - 1].x =
      data_points_[data_points_.size() - 1].x;
  bezier_curve_.control_points_[bezier_curve_.control_points_.size() - 1].y =
      data_points_[data_points_.size() - 1].y;

  // this guy will hold the BCDESolverST scope
  for (int k = 0; k < kNProcess_; k++) {
    auto solver = std::shared_ptr<BCDESolverST>(
        new BCDESolverST(k, kPopulation_ / kNProcess_, *this));
    solvers_.push_back(solver);
    solver->Start();
  }

  bezier_curve_.SetMadOptimizationCaching(data_points_,
                                          parameterization_values_);
}

BCDESolver::~BCDESolver() {
  for (auto i = solvers_.begin(); i != solvers_.end(); i++) {
    (*i)->Join();
  }
}

void BCDESolver::SolveOneGeneration() {
  using namespace std;

  // -2: ignore first and last CP
  const int N_CP = bezier_curve_.kNumberControlPoints_ - 2;

  Vec2 error_before;
  for (int i = 0; i < N_CP; i++) {

    bezier_curve_.UpdateVariableCPForMadOptimizationCaching(
        parameterization_values_, i + 1);
    bezier_curve_.CalcErrorWithMadOptimizationCaching(error_before);

    for (auto k = solvers_.begin(), ke = solvers_.end(); k != ke; k++) {
      (*k)->DoWork(i, error_before);
    }

    for (auto k = solvers_.begin(), ke = solvers_.end(); k != ke; k++) {
      (*k)->WaitWork();
    }

    if (kNProcess_ > 1) {
      Migration(i);
    }

    Vec2 low_error;
    low_error.x =
        solvers_[0]->population_errors_[i][lowest_error_index_[0].x].x;
    low_error.y =
        solvers_[0]->population_errors_[i][lowest_error_index_[0].y].y;
    int low_error_x_i = 0;
    int low_error_y_i = 0;

    if (kNProcess_ > 1) {
      for (int k = 1; k < kNProcess_; k++) {
        const double dx = solvers_[k]->parameters_[i][lowest_error_index_[k].x]
            .x;
        const double dy = solvers_[k]->parameters_[i][lowest_error_index_[k].y]
            .y;
        if (dx < low_error.x) {
          low_error_x_i = k;
          low_error.x = dx;
        }
        if (dy < low_error.y) {
          low_error_y_i = k;
          low_error.y = dy;
        }
      }
    }
    bezier_curve_.control_points_[i + 1].x = solvers_[low_error_y_i]
        ->parameters_[i][lowest_error_index_[low_error_y_i].x].x;
    bezier_curve_.control_points_[i + 1].y = solvers_[low_error_y_i]
        ->parameters_[i][lowest_error_index_[low_error_y_i].y].y;
  }
  generation_++;
}

void BCDESolver::Initialize() {

}

void BCDESolver::Migration(const int control_point) {
  const double phi = 1;

  const int populationInterval = kPopulation_ / kNProcess_;

  for (int k = 0; k < kNProcess_; k++) {
//if (0 <= phi) {
    int destinyIndex = k + 1;
    if (destinyIndex >= kNProcess_) {
      destinyIndex = 0;
    }
    const int pI = random_population_interval_[k];
    solvers_[destinyIndex]->parameters_[control_point][pI].x = solvers_[k]
        ->parameters_[control_point][lowest_error_index_[k].x].x;
    solvers_[destinyIndex]->parameters_[control_point][pI].y = solvers_[k]
        ->parameters_[control_point][lowest_error_index_[k].y].y;
//}
  }
}
