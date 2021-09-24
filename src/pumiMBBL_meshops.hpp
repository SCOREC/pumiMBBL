#ifndef pumiMBBLGPU_meshops_hpp
#define pumiMBBLGPU_meshops_hpp
#include "pumiMBBLGPU.hpp"

namespace pumi {

Vector3 get_rand_point_in_mesh_host(MBBL pumi_obj);

bool is_point_in_mesh_host(MBBL pumi_obj, Vector3 q);

void flatten_submeshID_and_cellID_host(MBBL pumi_obj, int isub, int icell, int jsub, int jcell, int* submeshID, int* cellID);

void locate_submesh_and_cell_x1_host(MBBL pumi_obj, double q, int* submeshID, int *cellID);

void locate_submesh_and_cell_x2_host(MBBL pumi_obj, double q, int* submeshID, int *cellID);

void update_submesh_and_cell_x1_host(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);

void update_submesh_and_cell_x2_host(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);

void calc_weights_x1_host(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2);

void calc_weights_x2_host(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2);

void calc_global_cellID_and_nodeID_fullmesh_host(MBBL pumi_obj, int kcell_x1, int kcell_x2, int *global_cell_2D, int *bottomleft_node, int *topleft_node);

void calc_global_cellID_and_nodeID_host(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node);

void get_directional_submeshID_and_cellID_host(MBBL pumi_obj, int submeshID, int cellID, int* isub, int *icell, int* jsub, int *jcell);

} // namespace pumi
#endif
