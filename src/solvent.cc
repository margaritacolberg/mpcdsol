// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// solvent.cc simulates the multi-particle collision (MPCD) model, in which the
// solvent, which consists of point particles, is coarse-grained to eliminate
// the calculation of solvent-solvent interactions

#include "solvent.h"

// initial positions of all solvent particles
void init_pos_sol(std::vector<Vec3>& pos, std::vector<Vec3>& pos_sol, const Box& box, Random& mt, const Param& p) {
    const unsigned int nbeads = pos.size();
    const unsigned int nsol = pos_sol.size();
    std::uniform_real_distribution<> uniform_sol(0.0, box.l);

    if (pos_sol.size() < 1) return;

    unsigned int i = 0;

    while (i < nsol) {
        double x = uniform_sol(mt);
        double y = uniform_sol(mt);
        double z = uniform_sol(mt);

        // declare a variable to check if there is an overlap of
        // a solvent particle and a bead
        bool overlap = false;

        for (unsigned int j = 0; j < nbeads; j++) {
            double dx = x - pos[j].x;
            double dy = y - pos[j].y;
            double dz = z - pos[j].z;
            box.mindist(dx, dy, dz);

            // check if a solvent particle overlaps with a bead
            if ((dx*dx + dy*dy + dz*dz) < p.sigma_sb2) {
                // if there is an overlap,
                overlap = true;
                // exit without checking the other distances
                break;
            }
        }

        if (overlap) continue;

        // store only the x-, y- and z-coordinates whose solvent
        // particle does not overlap with a bead
        pos_sol[i].x = x;
        pos_sol[i].y = y;
        pos_sol[i].z = z;
        i++;
    }
}

// split the box into smaller cells, and store solvent particles in each cell
void shift_grid(double dx, double dy, double dz, const std::vector<Vec3>& pos_sol, const Box& box, Cells& cells) {
    const unsigned int nsol = pos_sol.size();

    // clear old particle indices
    for (auto &cell: cells.cells) {
        cell.clear();
    }

    // number the solvent particles
    for (unsigned int j = 0; j < nsol; j++) {
        // shift grid of cells to preserve Galilean invariance
        double x = pos_sol[j].x + dx;
        double y = pos_sol[j].y + dy;
        double z = pos_sol[j].z + dz;

        box.minpos(x, y, z);
        unsigned int icell_x = std::clamp(int(cells.inv_lcell * x), 0, int(cells.ncell - 1));
        unsigned int icell_y = std::clamp(int(cells.inv_lcell * y), 0, int(cells.ncell - 1));
        unsigned int icell_z = std::clamp(int(cells.inv_lcell * z), 0, int(cells.ncell - 1));
        // convert the 3D indices to a single index for each cell
        unsigned int icell = icell_x + icell_y * cells.ncell + icell_z * cells.ncell * cells.ncell;

        LOG_DEBUG("particle " << j << " at " << pos_sol[j].x << ", " << pos_sol[j].y << ", " << pos_sol[j].z << " in cell " << icell_x << ", " << icell_y << ", " << icell_z);
        // store the 1D index for each cell, and the solvent
        // particles contained within, in a vector of vectors
        cells.cells[icell].emplace_back(j);
    }
}

void shift_grid_amount(const std::vector<Vec3>& pos_sol, const Box& box, Cells& cells, Random& mt) {
    std::uniform_real_distribution<> uniform_shift(0.0, 1.0);

    double dx = uniform_shift(mt);
    double dy = uniform_shift(mt);
    double dz = uniform_shift(mt);

    shift_grid(dx, dy, dz, pos_sol, box, cells);
}

// update velocities of the colliding solvent particles post-collision
void collide_sol(std::vector<Vec3>& vel_sol, Random& mt, Cells& cells, unsigned int angle) {
    const double a = angle * M_PI / 180;
    unsigned int cos_a = cos(a);
    unsigned int sin_a = sin(a);

    for (unsigned int i = 0; i < cells.cells.size(); i++) {
        double u_x, u_y, u_z;
        unit_sphere(mt, u_x, u_y, u_z);

        // rotation matrix with angle "a" degrees
        const double Rxx = cos_a + u_x*u_x*(1 - cos_a);
        const double Rxy = u_x*u_y*(1 - cos_a) - u_z*sin_a;
        const double Rxz = u_x*u_z*(1 - cos_a) + u_y*sin_a;
        const double Ryx = u_y*u_x*(1 - cos_a) + u_z*sin_a;
        const double Ryy = cos_a + u_y*u_y*(1 - cos_a);
        const double Ryz = u_y*u_z*(1 - cos_a) - u_x*sin_a;
        const double Rzx = u_z*u_x*(1 - cos_a) - u_y*sin_a;
        const double Rzy = u_z*u_y*(1 - cos_a) + u_x*sin_a;
        const double Rzz = cos_a + u_z*u_z*(1 - cos_a);

        double m_vel_x = 0.0;
        double m_vel_y = 0.0;
        double m_vel_z = 0.0;
        // number of solvent particles in a cell
        const unsigned int n = cells.cells[i].size();

        for (unsigned int j = 0; j < n; j++) {
            m_vel_x += vel_sol[j].x;
            m_vel_y += vel_sol[j].y;
            m_vel_z += vel_sol[j].z;
        }

        // calculate the center of mass velocity
        const double v_cm_x = m_vel_x / n;
        const double v_cm_y = m_vel_y / n;
        const double v_cm_z = m_vel_z / n;

        for (unsigned int k = 0; k < n; k++) {
            double dv_x = vel_sol[k].x - v_cm_x;
            double dv_y = vel_sol[k].y - v_cm_y;
            double dv_z = vel_sol[k].z - v_cm_z;

            vel_sol[k].x = v_cm_x + dv_x*Rxx + dv_y*Rxy + dv_z*Rxz;
            vel_sol[k].y = v_cm_y + dv_x*Ryx + dv_y*Ryy + dv_z*Ryz;
            vel_sol[k].z = v_cm_z + dv_x*Rzx + dv_y*Rzy + dv_z*Rzz;
        }
    }
}
