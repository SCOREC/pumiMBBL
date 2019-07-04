#ifndef pumi_routines_h
#define pumi_routines_h

#include "pumiMBBL.h"

int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh);

void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight);
void pumi_locatepoint_1D(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight);
void pumi_locatepoint_uniform_1D(int *cell, double *weight, double coord, double uniform_x_left, double uniform_t0, int uniform_Nel);
void pumi_locatepoint_BL_1D(int *cell, double *weight, double coord, double x_end, double r, double t0, int local_Nel, pumi_meshflag_t pumi_flag);
void pumi_compute_covolume_1D(pumi_mesh_t *pumi_mesh, int Nel_total, double *covolume);

#endif /* pumi_routines_h */
