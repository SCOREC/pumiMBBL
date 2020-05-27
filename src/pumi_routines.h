#ifndef pumi_routines_h
#define pumi_routines_h

#include "pumiMBBL.h"

#define MAX_DIM 2

int pumi_total_elements(pumi_mesh_t *pumi_mesh);
int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh);
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh);

//void pumi_locate_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
//void pumi_locate_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
//void pumi_locate_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);

//typedef void (*pumi_locate_ptr)(pumi_mesh_t*, int, double, int*, double*);
//void pumi_initialize_locate_functions(pumi_mesh_t *pumi_mesh);
//void pumi_finalize_locate_functions();

int pumi_locatecell_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);

void pumi_calc_weights_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);

void pumi_calc_weights(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);

typedef int (*pumi_locatecell_ptr)(pumi_mesh_t*, int, double);
pumi_locatecell_ptr *pumi_locatecell_fnptr;
typedef void (*pumi_calc_weights_ptr)(pumi_mesh_t*, int, int, double, int*, double*);
pumi_calc_weights_ptr *pumi_calc_weights_fnptr;
void pumi_initialize_locatecell_and_calcweights_functions(pumi_mesh_t *pumi_mesh);
void pumi_finalize_locatecell_and_calcweights_functions();

typedef int (*pumi_updatecell_ptr)(pumi_mesh_t*, int, int, double);
pumi_updatecell_ptr *pumi_updatecell_fnptr;
int pumi_updatecell_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);

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
typedef double (*pumi_gradingratio_ptr)(pumi_mesh_t*, int);
pumi_gradingratio_ptr pumi_gradingratio_fnptr[MAX_DIM];
//double (*pumi_gradingratio_fnptr[])(pumi_mesh_t*, int) = {pumi_return_1D_gradingratio, pumi_return_2D_gradingratio};

double pumi_return_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
double pumi_return_1D_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
double pumi_return_2D_elemsize(pumi_mesh_t* pumi_mesh, int index, int offset);
typedef double (*pumi_elemsize_ptr)(pumi_mesh_t*, int, int);
pumi_elemsize_ptr pumi_elemsize_fnptr[MAX_DIM];
//double (*pumi_elemsize_fnptr[])(pumi_mesh_t*, int, int) = {pumi_return_1D_elemsize, pumi_return_2D_elemsize};

double pumi_return_covolume(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_covolume_1D(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_covolume_2D(pumi_mesh_t* pumi_mesh, int inode);
typedef double (*pumi_covolume_ptr)(pumi_mesh_t*, int);
pumi_covolume_ptr pumi_covolume_fnptr[MAX_DIM];
//double (*pumi_covolume_fnptr[])(pumi_mesh_t*, int) = {pumi_return_covolume_1D, pumi_return_covolume_2D};

void pumi_initialize_multiD_functions(pumi_mesh_t *pumi_mesh);

double pumi_return_smallest_elemsize(pumi_mesh_t *pumi_mesh);

//int pumi_locate_submesh_1D(pumi_mesh_t *pumi_mesh, double coords);
//int pumi_update_submesh_1D(pumi_mesh_t *pumi_mesh, double coords, int isubmesh);

void pumi_locate_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int *submeshID, int *cellID);
void pumi_update_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int *submeshID, int *cellID);
void pumi_update_submesh_and_update_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);

#endif /* pumi_routines_h */
