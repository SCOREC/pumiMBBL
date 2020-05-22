#ifndef pumi_routines_h
#define pumi_routines_h

#include "pumiMBBL.h"

int pumi_total_elements(pumi_mesh_t *pumi_mesh);
int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh);

void pumi_locate_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
void pumi_locate_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
void pumi_locate_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);

typedef void (*pumi_locate_ptr)(pumi_mesh_t*, int, double, int*, double*);
pumi_locate_ptr *pumi_locate_function;
void pumi_initialize_locate_functions(pumi_mesh_t *pumi_mesh);
void pumi_finalize_locate_functions();

int pumi_locatecell_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);

void pumi_calc_weights_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);

void pumi_calc_weights(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);

typedef int (*pumi_locatecell_ptr)(pumi_mesh_t*, int, double);
pumi_locatecell_ptr *pumi_locatecell_fnptr;
typedef void (*pumi_calc_weights_ptr)(pumi_mesh_t*, int, int, double, int*, double*);
pumi_calc_weights_ptr *pumi_calc_weights_fnptr;
void pumi_initialize_locatecell_functions(pumi_mesh_t *pumi_mesh);
void pumi_finalize_locatecell_functions();

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
double pumi_return_2D_gradingratio(pumi_mesh_t* pumi_mesh, int node);

double pumi_return_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
double pumi_return_1D_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
double pumi_return_2D_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);

double pumi_return_covolume(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_covolume_1D(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_covolume_2D(pumi_mesh_t* pumi_mesh, int inode);

double pumi_return_smallest_elemsize(pumi_mesh_t *pumi_mesh);

int pumi_locate_submesh_1D(pumi_mesh_t *pumi_mesh, double coords);
int pumi_update_submesh_1D(pumi_mesh_t *pumi_mesh, double coords, int isubmesh);

void pumi_locate_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int *submeshID, int *cellID);
void pumi_update_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int *submeshID, int *cellID);

#endif /* pumi_routines_h */
