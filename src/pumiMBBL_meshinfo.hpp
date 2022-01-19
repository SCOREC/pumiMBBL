#ifndef pumiMBBLGPU_meshinfo_hpp
#define pumiMBBLGPU_meshinfo_hpp
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

double return_gradingratio(MBBL pumi_obj, int dir, int node);double return_elemsize(MBBL pumi_obj, int dir, int index, int offset);
double return_covolume(MBBL pumi_obj, int inode_x1);
double return_covolume_fullmesh(MBBL pumi_obj, int inode_x1, int inode_x2);
double return_covolume(MBBL pumi_obj, int inode_x1, int inode_x2);
double return_elemsize(MBBL pumi_obj, int dir, int index, int offset);
void where_is_node(MBBL pumi_obj, int knode_x1, int knode_x2, bool* on_bdry, bool* in_domain, int* bdry_tag, int* bdry_dim);
bool is_fullmesh(MBBL pumi_obj);
double get_global_x1_min_coord(MBBL pumi_obj);
double get_global_x1_max_coord(MBBL pumi_obj);
double get_global_x2_min_coord(MBBL pumi_obj);
double get_global_x2_max_coord(MBBL pumi_obj);
int get_total_mesh_elements(MBBL pumi_obj);
int get_total_mesh_nodes(MBBL pumi_obj);
int get_num_x1_submesh(MBBL pumi_obj);
int get_num_x1_elems_in_submesh_host(MBBL pumi_obj, int isubmesh);
int get_num_x1_elems_before_submesh_host(MBBL pumi_obj, int isubmesh);
double get_x1_elem_size_in_submesh_host(MBBL pumi_obj, int isubmesh, int icell);
int get_num_x2_submesh(MBBL pumi_obj);
int get_num_x2_elems_in_submesh_host(MBBL pumi_obj, int isubmesh);
int get_num_x2_elems_before_submesh_host(MBBL pumi_obj, int isubmesh);
double get_x2_elem_size_in_submesh_host(MBBL pumi_obj, int isubmesh, int icell);
int get_total_x1_elements(MBBL pumi_obj);
int get_total_x2_elements(MBBL pumi_obj);
double get_x1_gradingratio_in_submesh_host(MBBL pumi_obj, int isub);
double get_x2_gradingratio_in_submesh_host(MBBL pumi_obj, int isub);
int get_x1_nodeID_at_interface_host(MBBL pumi_obj, int isub);
int get_x2_nodeID_at_interface_host(MBBL pumi_obj, int isub);
int get_x1_gradingratio_at_interface_host(MBBL pumi_obj, int isub);
int get_x2_gradingratio_at_interface_host(MBBL pumi_obj, int isub);
int get_total_submesh_blocks(MBBL pumi_obj);
int get_total_elements_in_block(MBBL pumi_obj, int flattened_submesh_ID);
bool is_horizontal_edge_host(MBBL pumi_obj, int iEdge);
int get_num_interior_nodes_on_block(MBBL pumi_obj, int isub, int jsub);
double get_mesh_volume(MBBL pumi_obj);
Vector3 get_bdry_edge_normal_host(MBBL pumi_obj,  int iEdge);
Vector3 get_bdry_vert_normal_host(MBBL pumi_obj,  int iEdge);
int get_num_interior_nodes_on_edge(MBBL pumi_obj, int iEdge);
int get_num_faces_on_edge(MBBL pumi_obj,  int iEdge);
int get_starting_faceID_on_bdry_edge(MBBL pumi_obj,  int iEdge);
int get_west_edgeID(MBBL pumi_obj, int isub, int jsub);
int get_east_edgeID(MBBL pumi_obj, int isub, int jsub);
int get_north_edgeID(MBBL pumi_obj, int isub, int jsub);
int get_south_edgeID(MBBL pumi_obj, int isub, int jsub);
bool is_block_active_host(MBBL pumi_obj, int isub, int jsub);
bool is_block_active_host(MBBL pumi_obj, int flattened_submesh_ID);
bool is_edge_bdry(MBBL pumi_obj,  int iEdge);
bool is_vert_bdry(MBBL pumi_obj,  int iVert);
int get_total_mesh_block_edges(MBBL pumi_obj);
int get_total_mesh_block_verts(MBBL pumi_obj);
int get_block_vert_submeshID_host(MBBL pumi_obj, int iVert);
int get_block_edge_submeshID_host(MBBL pumi_obj, int iEdge);
void get_block_vert_submeshIDs_host(MBBL pumi_obj, int iVert, int *isub, int *jsub);
void get_block_edge_submeshIDs_host(MBBL pumi_obj, int iEdge, int *isub, int *jsub);
int get_node_submeshID(MBBL pumi_obj, int knode_x1, int knode_x2);
int get_elem_submeshID(MBBL pumi_obj, int kcell_x1, int kcell_x2);
int get_num_block_interior_nodes(MBBL pumi_obj);
int get_num_block_edge_interior_nodes(MBBL pumi_obj);
bool check_edge_index_bounds(MBBL pumi_obj, int iEdge);
int compute_global_nodeID_2D(MBBL pumi_obj, int isubmesh, int jsubmesh, int fullmesh_node_id);
std::vector<int> get_nodes_on_bdry_edge(MBBL pumi_obj, int iEdge);
bool check_node_index_bounds(MBBL pumi_obj, int knode_x1, int knode_x2);
int get_global_nodeID_2D(MBBL pumi_obj, int knode_x1, int knode_x2);

} // namespace pumi
#endif
