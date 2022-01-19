#ifndef pumiMBBLGPU_meshutils_hpp
#define pumiMBBLGPU_meshutils_hpp
#include "pumiMBBLGPU.hpp"

namespace pumi {

void check_is_pumi_working();
void print_mesh_skeleton(MBBL pumi_obj);
void print_blockwise_nodeIDs(MBBL pumi_obj);
void print_node_submeshID(MBBL pumi_obj);
void print_fullmesh_nodeIDs(MBBL pumi_obj);
void print_2D_node_coordinates(MBBL pumi_obj);
void print_2D_node_elem_connectivity(MBBL pumi_obj);
Vector3View compute_2D_field_gradient(MBBL pumi_obj, DoubleView phi);
Vector3View compute_2D_field_gradient_v2(MBBL pumi_obj, DoubleView phi);
} // namespace pumi
#endif
