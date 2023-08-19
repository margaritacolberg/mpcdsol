// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// system.h contains a list of properties needed to define a system of protein
// hard spheres and solvent particles

#ifndef MPCD_SYSTEM_H
#define MPCD_SYSTEM_H

#include "vec3.h"
#include <vector>

struct System {
  // positions of all beads at time t
  std::vector<Vec3> pos;
  // positions of all solvent particles at time t
  std::vector<Vec3> pos_sol;
  // velocities of all beads at time t
  std::vector<Vec3> vel;
  // velocities of all solvent particles at time t
  std::vector<Vec3> vel_sol;
  // bead clocks
  std::vector<double> times;
  // solvent particle clocks
  std::vector<double> times_sol;
  // collision counters of beads
  std::vector<uint64_t> counter;
  // collision counters of solvent particles
  std::vector<uint64_t> counter_sol;

  System(unsigned int nbeads, unsigned int nsol)
      : pos(nbeads), pos_sol(nsol), vel(nbeads), vel_sol(nsol), times(nbeads),
        times_sol(nsol), counter(nbeads), counter_sol(nsol) {}
};

#endif
