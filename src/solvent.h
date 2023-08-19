// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MPCD_SOLVENT_H
#define MPCD_SOLVENT_H

#include "hardspheres.h"

#if NDEBUG
#define LOG_DEBUG(x)
#else
#define LOG_DEBUG(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl
#endif

void init_pos_sol(std::vector<Vec3>& pos, std::vector<Vec3>& pos_sol, const Box& box, Random& mt, const Param& p);

void shift_grid_amount(const std::vector<Vec3>& pos_sol, const Box& box, Cells& cells, Random& mt);

void shift_grid(double dx, double dy, double dz, const std::vector<Vec3>& pos_sol, const Box& box, Cells& cells);

void collide_sol(std::vector<Vec3>& vel_sol, Random& mt, Cells& cells, unsigned int angle);

#endif
