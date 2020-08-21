#ifndef pumi_routines_h
#define pumi_routines_h

#include "pumiMBBL.h"

#define MAX_DIM 2

int pumi_total_elements(pumi_mesh_t *pumi_mesh);
int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh);
int pumi_total_elements_2D(pumi_mesh_t *pumi_mesh);
int pumi_total_nodes(pumi_mesh_t *pumi_mesh);
int pumi_total_nodes_1D(pumi_mesh_t *pumi_mesh);
int pumi_total_nodes_2D(pumi_mesh_t *pumi_mesh);
int pumi_submesh_total_elements_1D(pumi_mesh_t *pumi_mesh, int isubmesh);
double pumi_global_x1_min(pumi_mesh_t *pumi_mesh);
double pumi_global_x1_max(pumi_mesh_t *pumi_mesh);
double pumi_global_x2_min(pumi_mesh_t *pumi_mesh);
double pumi_global_x2_max(pumi_mesh_t *pumi_mesh);

//void pumi_locate_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
//void pumi_locate_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);
//void pumi_locate_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight);

//typedef void (*pumi_locate_ptr)(pumi_mesh_t*, int, double, int*, double*);
//void pumi_initialize_locate_functions(pumi_mesh_t *pumi_mesh);
//void pumi_finalize_locate_functions();

int pumi_locatecell_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_bottomBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);
int pumi_locatecell_in_topBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord);

void pumi_calc_weights_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_bottomBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_bottomBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_topBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights_in_topBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight);
void pumi_calc_weights(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight, int dir);

void pumi_calc_node_coords_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node);
void pumi_calc_node_coords_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node);
void pumi_calc_node_coords_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node);
void pumi_calc_node_coords_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node);
void pumi_calc_node_coords_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node);

void pumi_calc_node_coords(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node);

double pumi_calc_elem_size_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_bottomBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
double pumi_calc_elem_size_in_topBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);

double pumi_calc_elem_size(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, int dir);

int pumi_dummy_elem_node_ID(double coord_x1, double coord_x2, double dx1, double dx2, int Nel_total_x1, int *node1, int *node3);
int pumi_dummy_elem_node_ID_v2(int kcell_x1, int kcell_x2, double dx1, double dx2, int Nel_total_x1, int *node1, int *node3);


typedef int (*pumi_locatesubmesh_ptr)(pumi_mesh_t*, double);
pumi_locatesubmesh_ptr pumi_locatesubmesh_fnptr[MAX_DIM];
typedef int (*pumi_locatecell_ptr)(pumi_mesh_t*, int, double);
pumi_locatecell_ptr **pumi_locatecell_fnptr;
typedef void (*pumi_calc_weights_ptr)(pumi_mesh_t*, int, int, double, int*, double*);
pumi_calc_weights_ptr **pumi_calc_weights_fnptr;
typedef void (*pumi_calc_node_coords_ptr)(pumi_mesh_t*, int, int, double*,double*);
pumi_calc_node_coords_ptr *pumi_calc_node_coords_fnptr;
typedef double (*pumi_calc_elem_size_ptr)(pumi_mesh_t*, int, int);
pumi_calc_elem_size_ptr **pumi_calc_elem_size_fnptr;
void pumi_initialize_locatecell_and_calcweights_functions(pumi_mesh_t *pumi_mesh);
void pumi_finalize_locatecell_and_calcweights_functions();

typedef int (*pumi_updatesubmesh_ptr)(pumi_mesh_t*, double, int, int*);
pumi_updatesubmesh_ptr pumi_updatesubmesh_fnptr[MAX_DIM];
typedef int (*pumi_updatecell_ptr)(pumi_mesh_t*, int, int, double);
pumi_updatecell_ptr **pumi_updatecell_fnptr;
int pumi_updatecell_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_bottomBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_bottomBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_topBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);
int pumi_updatecell_in_topBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord);

/*!
* \brief elemsize index offset enum, possible ways to call pumi_return_elemsize()
*/
typedef enum pumi_elemsize_index_offset{
  pumi_elem_input_offset = 0, //!< no offset for direct element input
  pumi_elem_on_right_offset = 0, //!< no offset for node input and querying right element size
  pumi_elem_on_top_offset = 0,
  pumi_elem_on_left_offset = -1, //!< -1 offset for node input and querying left element size
  pumi_elem_on_bottom_offset = -1,
} pumi_elemsize_index_offset_t;

typedef enum pumi_directions{
    pumi_x1 = 0,
    pumi_x2 = 1,
} pumi_directions_t;

void pumi_BL_elemsize_ON(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_ON_1D(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_ON_2D(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_OFF(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_OFF_1D(pumi_mesh_t *pumi_mesh);
void pumi_BL_elemsize_OFF_2D(pumi_mesh_t *pumi_mesh);

double pumi_return_gradingratio(pumi_mesh_t *pumi_mesh, int node, int dir);

double pumi_return_elemsize(pumi_mesh_t *pumi_mesh, int index, int offset, int dir);

double pumi_return_covolume_1D(pumi_mesh_t* pumi_mesh, int inode);
double pumi_return_covolume_2D(pumi_mesh_t* pumi_mesh, int inp_x1, int inp_x2);
//double (*pumi_covolume_fnptr[])(pumi_mesh_t*, int) = {pumi_return_covolume_1D, pumi_return_covolume_2D};

typedef int (*pumi_nodeID_ptr)(pumi_mesh_t*, int, int, int, int, int*, int*);
pumi_nodeID_ptr **pumi_nodeID_fnptr;

bool pumi_mesh_with_no_inactive_blocks(pumi_mesh_t *pumi_mesh);
void pumi_initialize_nodeID_functions(pumi_mesh_t *pumi_mesh);
void pumi_finalize_nodeID_functions(pumi_mesh_t *pumi_mesh);

int pumi_calc_elementID_and_nodeID_typeA(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);

int pumi_calc_elementID_and_nodeID_typeB(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);
typedef void (*pumi_typeB_nodeoffset_ptr)(pumi_mesh_t*, int, int, int, int*, int*);
pumi_typeB_nodeoffset_ptr pumi_typeB_nodeoffset_fnptr[2];
void pumi_typeB_nodeoffset_expression1(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);
void pumi_typeB_nodeoffset_expression2(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);

int pumi_calc_elementID_and_nodeID_typeC(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);
typedef void (*pumi_typeC_nodeoffset_ptr)(pumi_mesh_t*, int, int, int, int*, int*);
pumi_typeC_nodeoffset_ptr pumi_typeC_nodeoffset_fnptr[2];
void pumi_typeC_nodeoffset_expression1(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);
void pumi_typeC_nodeoffset_expression2(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);

int pumi_calc_elementID_and_nodeID_typeD(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);
typedef void (*pumi_typeD_nodeoffset_ptr)(pumi_mesh_t*, int, int, int, int*, int*);
pumi_typeD_nodeoffset_ptr pumi_typeD_nodeoffset_fnptr[3];
void pumi_typeD_nodeoffset_expression1(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);
void pumi_typeD_nodeoffset_expression2(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);
void pumi_typeD_nodeoffset_expression3(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3);

int pumi_calc_elementID_and_nodeID(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);
int pumi_calc_elementID_and_nodeID_with_global_offset(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);
int pumi_calc_elementID_and_nodeID_on_fullmesh(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x1, int icell_x2, int *node1, int *node3);

bool pumi_is_node_active(pumi_mesh_t *pumi_mesh, int node_x1, int node_x2);//, int *isubmesh_x1, int *inp_x1, int *isubmesh_x2, int *inp_x2);
//int pumi_global_node_ID(pumi_mesh_t *pumi_mesh, int isubmesh_x1, int inp_x1, int isubmesh_x2, int inp_x2);
double pumi_return_smallest_elemsize(pumi_mesh_t *pumi_mesh);

void pumi_locate_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int *submeshID, int *cellID, int dir);
void pumi_update_submesh_and_update_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID, int dir);
int pumi_global_cell_ID(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell);
#endif /* pumi_routines_h */
