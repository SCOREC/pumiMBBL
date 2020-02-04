#ifndef pumi_routines_h
#define pumi_routines_h

#include "pumiMBBL.h"

int pumi_total_elements(pumi_mesh_t *pumi_mesh);
int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh);

void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight);
void pumi_locatepoint_1D(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight);
void pumi_locatepoint_uniform_1D(int *cell, double *weight, double coord, double uniform_x_left, double uniform_t0, int uniform_Nel);
void pumi_locatepoint_BL_1D(int *cell, double *weight, double coord, double x_end, double r, double t0, double log_r, double r_t0_ratio, int local_Nel, pumi_meshflag_t pumi_flag);

double pumi_compute_covolume_1D(int inode, int Nel_total, double *elemsize)__attribute__((deprecated("Use pumi_return_covolume() instead")));
void pumi_compute_elemsize_1D(pumi_mesh_t *pumi_mesh, int Nel_total, double *elemsize)__attribute__((deprecated("Use pumi_return_elemsize() instead")));
void pumi_compute_nodal_gradingratio_1D(double *elemsize, int Nel_total, double *gradingratio)__attribute__((deprecated("Use pumi_return_gradingratio() instead")));

/*!
* \brief elemsize index offset enum, possible ways to call pumi_return_elemsize()
*/
typedef enum pumi_elemsize_index_offset{
  elem_input_offset = 0, //!< no offset for direct element input
  node_input_right_elem_offset = 0, //!< no offset for node input and querying right element size
  node_input_left_elem_offset = -1, //!< -1 offset for node input and querying left element size
} pumi_elemsize_index_offset_t;

void pumi_BL_elemsize_ON(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_ON_1D(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_OFF(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_OFF_1D(pumi_mesh_t *pumi_mesh);

double pumi_return_gradingratio(pumi_mesh_t *pumi_mesh, int node);
double pumi_return_1D_gradingratio(pumi_mesh_t* pumi_mesh, int node);
double pumi_return_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
double pumi_return_1D_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
double pumi_return_covolume(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_covolume_1D(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_smallest_elemsize(pumi_mesh_t *pumi_mesh);

#endif /* pumi_routines_h */
