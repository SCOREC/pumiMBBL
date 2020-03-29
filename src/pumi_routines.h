#ifndef pumi_routines_h
#define pumi_routines_h

#include "pumiMBBL.h"

int pumi_total_elements(pumi_mesh_t *pumi_mesh);
int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh);

void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int particle_submesh, int *particle_cell, double *cell_weight);
void pumiMBBL_locatepoint_1D(pumi_mesh_t *pumi_mesh, double particle_coordinate, int particle_submesh, int *particle_cell, double *cell_weight);
void pumi_locate_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
void pumi_locate_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
void pumi_locate_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
void pumi_dummylocate(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);

double pumi_compute_covolume_1D(int inode, int Nel_total, double *elemsize)__attribute__((deprecated("Use pumi_return_covolume() instead")));
void pumi_compute_elemsize_1D(pumi_mesh_t *pumi_mesh, int Nel_total, double *elemsize)__attribute__((deprecated("Use pumi_return_elemsize() instead")));
void pumi_compute_nodal_gradingratio_1D(double *elemsize, int Nel_total, double *gradingratio)__attribute__((deprecated("Use pumi_return_gradingratio() instead")));

/*!
* \brief elemsize index offset enum, possible ways to call pumi_return_elemsize()
*/
typedef enum pumi_elemsize_index_offset{
  pumi_elem_input_offset = 0, //!< no offset for direct element input
  pumi_elem_on_right_offset = 0, //!< no offset for node input and querying right element size
  pumi_elem_on_left_offset = -1, //!< -1 offset for node input and querying left element size
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
