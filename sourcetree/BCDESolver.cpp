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

  parameters_.resize(bezier_curve.kNumberControlPoints_ - 2);
  for (auto i = parameters_.begin(), e = parameters_.end(); i != e; i++) {
    i->resize(kPopulation_);
  }

  population_errors_.resize(bezier_curve.kNumberControlPoints_ - 2);
  for (auto i = population_errors_.begin(), e = population_errors_.end();
      i != e; i++) {
    i->resize(kPopulation_);
  }

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
    auto solver = std::shared_ptr<BCDESolverST>(new BCDESolverST(k, *this));
    solvers_.push_back(solver);
    solver->Start();
  }
  Initialize();
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

  for (int i = 0; i < N_CP; i++) {
    for (auto k = solvers_.begin(), ke = solvers_.end(); k != ke; k++) {
      (*k)->DoWork(i);
    }

    for (auto k = solvers_.begin(), ke = solvers_.end(); k != ke; k++) {
      (*k)->WaitWork();
    }

    if (kNProcess_ > 1) {
      Migration(i);
    }

    bezier_curve_.control_points_[i + 1].x =
        parameters_[i][lowest_error_index_[0].x].x;
    bezier_curve_.control_points_[i + 1].y =
        parameters_[i][lowest_error_index_[0].y].y;
  }
  generation_++;
}

void BCDESolver::Initialize() {
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
    const int pSDestiny = populationInterval * destinyIndex;
    parameters_[control_point][pSDestiny + pI].x =
        parameters_[control_point][lowest_error_index_[k].x].x;
    parameters_[control_point][pSDestiny + pI].y =
        parameters_[control_point][lowest_error_index_[k].y].y;
    //}
  }
}
