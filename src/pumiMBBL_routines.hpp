#ifndef pumiMBBLGPU_routines_hpp
#define pumiMBBLGPU_routines_hpp
#include "pumiMBBLGPU.hpp"

namespace pumi {
/*!
* \brief enum type for element index offsets
*
* To be used while querying element size using
* either element ID or associated node ID
*/
enum elemsize_index_offset{
    elem_input_offset       =  0, //!< zero offset for direct element ID input
    elem_on_max_side_offset =  0, //!< offset for node ID input and querying element to the max side
    elem_on_min_side_offset = -1, //!< offset for node ID input and querying element to the min side
};

void check_is_pumi_working();
double return_gradingratio(MBBL pumi_obj, int dir, int node);double return_elemsize(MBBL pumi_obj, int dir, int index, int offset);
double return_covolume(MBBL pumi_obj, int inode_x1);
double return_covolume_fullmesh(MBBL pumi_obj, int inode_x1, int inode_x2);
double return_covolume(MBBL pumi_obj, int inode_x1, int inode_x2);
double return_elemsize(MBBL pumi_obj, int dir, int index, int offset);
void where_is_node(MBBL pumi_obj, int knode_x1, int knode_x2, bool* on_bdry, bool* in_domain, int* bdry_tag, int* bdry_dim);
double get_global_x1_min_coord(MBBL pumi_obj);
double get_global_x1_max_coord(MBBL pumi_obj);
double get_global_x2_min_coord(MBBL pumi_obj);
double get_global_x2_max_coord(MBBL pumi_obj);
void print_mesh_skeleton(MBBL pumi_obj);
int get_total_mesh_elements(MBBL pumi_obj);
int get_total_mesh_nodes(MBBL pumi_obj);
int get_num_x1_submesh(MBBL pumi_obj);
int get_num_x1_elems_in_submesh(MBBL pumi_obj, int isubmesh);
int get_num_x1_elems_before_submesh(MBBL pumi_obj, int isubmesh);
int get_num_x2_submesh(MBBL pumi_obj);
int get_num_x2_elems_in_submesh(MBBL pumi_obj, int isubmesh);
int get_num_x2_elems_before_submesh(MBBL pumi_obj, int isubmesh);
int get_total_x1_elements(MBBL pumi_obj);
int get_total_x2_elements(MBBL pumi_obj);
int get_total_submesh_blocks(MBBL pumi_obj);
int get_total_elements_in_block(MBBL pumi_obj, int flattened_submesh_ID);
double get_mesh_volume(MBBL pumi_obj);
std::vector<double> get_bdry_normal(MBBL pumi_obj,  int iEdge);
int get_num_faces_on_bdry(MBBL pumi_obj,  int iEdge);
int get_starting_faceID_on_bdry(MBBL pumi_obj,  int iEdge);
bool is_block_active(MBBL pumi_obj, int isub, int jsub);
bool is_block_active(MBBL pumi_obj, int flattened_submesh_ID);
bool is_edge_bdry(MBBL pumi_obj,  int iEdge);
int get_total_mesh_block_edges(MBBL pumi_obj);
int get_global_nodeID(MBBL pumi_obj, int submeshID, int fullmesh_node_id);
void get_edge_info(MBBL pumi_obj,  int iEdge, int *Knp, int *next_offset, int *submeshID);
std::vector<double> get_rand_point_in_mesh(MBBL pumi_obj);
bool is_point_in_mesh(MBBL pumi_obj, std::vector<double> q);
void print_nodeIDs(MBBL pumi_obj);
void flatten_submeshID_and_cellID_host(MBBL pumi_obj, int isub, int icell, int jsub, int jcell, int* submeshID, int* cellID);

void locate_submesh_and_cell_x1_host(MBBL pumi_obj, double q, int* submeshID, int *cellID);
void locate_submesh_and_cell_x2_host(MBBL pumi_obj, double q, int* submeshID, int *cellID);
void update_submesh_and_cell_x1_host(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);
void update_submesh_and_cell_x2_host(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);
void calc_weights_x1_host(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2);
void calc_weights_x2_host(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2);
void calc_global_cellID_and_nodeID_fullmesh_host(MBBL pumi_obj, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node);
void calc_global_cellID_and_nodeID_host(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node);
void get_directional_submeshID_and_cellID_host(MBBL pumi_obj, int submeshID, int cellID, int* isub, int *icell, int* jsub, int *jcell);

} // namespace pumi
#endif
