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
#include <functional>
#include <random>
#include <array>

#include "BCDESolver.h"
#include "BezierCurve.h"

class BCDESolverST {
 public:
  /*
   *  [control point][population][dimension]
   *  Here, "CP-2" is because we ignore the first and last CP
   */
  std::vector<std::vector<Vec2>> parameters_;
  std::vector<std::vector<Vec2>> population_errors_;

  BCDESolverST(const int process, const int population, BCDESolver& solver);
  virtual ~BCDESolverST();

  void Start();
  void Join();

  void Initialize();
  void DoWork(const int control_point, const Vec2& error_before);
  void WaitWork();

 private:

  const int kPopulation_;

  std::thread thread_;
  std::mutex mutex_;
  std::condition_variable cond_;
  bool pending_work_;
  bool finish_;

  bool work_ready_;
  std::mutex work_ready_lock_;
  std::condition_variable work_ready_cond_;

  BCDESolver& solver_;
  int control_point_;
  const int kProcess_;

  BezierCurve bezier_curve_;
  Vec2 error_before_;

  std::uniform_int_distribution<int> dt_int_pop_interval_from_zero_;
  std::uniform_int_distribution<int> dt_int_0_2_;
  std::uniform_real_distribution<double> dt_real_0_1_;
  std::mt19937 engine_;  // Mersenne twister MT19937
  std::mt19937 engine0_;
  std::mt19937 engine1_;

  std::vector<std::array<int,3>> random_rs_;
  int r_rs_;
  std::vector<int> random_j_;
  std::vector<std::array<float,2>> random_cr_;

  uint32_t RN(uint32_t min, uint32_t max);
  void GenerateRandom();
  void Run();
  void Mutate(const int actual_index, Vec2& trials);
  void Select(const Vec2& trial, const Vec2& error_before, Vec2& error_new);
};

#endif /* BCDESOLVERST_H_ */
