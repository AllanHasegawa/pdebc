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

#include "BCDESolver.h"
#include "BezierCurve.h"

class BCDESolverST {
 public:
  BCDESolverST(const int process, BCDESolver& solver);
  virtual ~BCDESolverST();

  void Start();
  void Join();

  void DoWork(const int control_point);
  void WaitWork();

 private:
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

  std::uniform_int_distribution<int> dt_int_pop_interval_;
  std::uniform_int_distribution<int> dt_int_0_2_;
  std::uniform_real_distribution<double> dt_real_0_1_;
  std::mt19937 engine_;  // Mersenne twister MT19937


  void Run();
  void Mutate(const int actual_index, Vec2& trials);
  void Select(const Vec2& trial, const Vec2& error_before, Vec2& error_new);
};

#endif /* BCDESOLVERST_H_ */
