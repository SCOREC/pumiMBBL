#ifndef pumiMBBLGPU_impl_hpp
#define pumiMBBLGPU_impl_hpp

#include "pumiMBBLGPU.hpp"

namespace pumi {

/**
 * @brief Locates local cell ID of given point inside a block
 * @param[in] Submesh object (on device)
 * @param[in] point coordinate to be located
 * @return local cell ID
 */
KOKKOS_INLINE_FUNCTION
int locate_cell(DevicePointer<Submesh> submesh, double q) {
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->locate_cell(q);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->locate_cell(q);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->locate_cell(q);
        case (arbitrary) :
            return static_cast<Arbitrary_Submesh*>(submesh())->locate_cell(q);
        case (unassigned) :
            return -1;
    }
    return -1;
}

/**
 * @brief updates local cell ID of given particle coordinate
 * @param[in] Submesh object (on device)
 * @param[in] point coordinate to be located
 * @param[in] previous local cell ID of particle
 * @return upddated local cell ID
 */
KOKKOS_INLINE_FUNCTION
int update_cell(DevicePointer<Submesh> submesh, double q, int icell) {
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->update_cell(q,icell);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->update_cell(q,icell);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->update_cell(q,icell);
        case (arbitrary) :
            return static_cast<Arbitrary_Submesh*>(submesh())->update_cell(q,icell);
        case (unassigned) :
            return -1;
    }
    return -1;
}

/**
 * @brief Computes element size inside a block
 * @param[in] Submesh object (on device)
 * @param[in] local cell ID for which elemetn size is needed
 * @return Value of element size
 */
KOKKOS_INLINE_FUNCTION
double elem_size(DevicePointer<Submesh> submesh, int icell){
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->elem_size(icell);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->elem_size(icell);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->elem_size(icell);
        case (arbitrary) :
            return static_cast<Arbitrary_Submesh*>(submesh())->elem_size(icell);
        case (unassigned) :
            return -999.0;
    }
    return -999.0;
}

/**
 * @brief Computes linear 1D weights for gather/scatter operations and
 *        directional global cell ID
 * @param[in] Submesh object (on device)
 * @param[in] point coordinate inside submesh block
 * @param[in] local cell ID of point inside submesh block
 * @param[out] directional global cell ID
 * @param[out] linear 1D weight (correspoding to max-side node)
 */
KOKKOS_INLINE_FUNCTION
void calc_weights(DevicePointer<Submesh> submesh, double q, int local_cell, int *global_cell, double *Wgh2){
    switch (submesh()->meshtype){
        case (uniform) :
            static_cast<Uniform_Submesh*>(submesh())->calc_weights(q,local_cell,global_cell,Wgh2);
            return;
        case (minBL) :
            static_cast<MinBL_Submesh*>(submesh())->calc_weights(q,local_cell,global_cell,Wgh2);
            return;
        case (maxBL) :
            static_cast<MaxBL_Submesh*>(submesh())->calc_weights(q,local_cell,global_cell,Wgh2);
            return;
        case (arbitrary) :
            static_cast<Arbitrary_Submesh*>(submesh())->calc_weights(q,local_cell,global_cell,Wgh2);
            return;
        case (unassigned) :
            *global_cell = -1;
            *Wgh2 = -999.0;
            return;
    }
    *global_cell = -1;
    *Wgh2 = -999.0;
    return;
}

/**
 * @brief Fetches node coordinate inside submesh block
 * @param[in] Submesh object (on device)
 * @param[in] local node ID inside block
 * @return Value of node coordinate
 */
KOKKOS_INLINE_FUNCTION
double node_coords(DevicePointer<Submesh> submesh, int inode){
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->node_coords(inode);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->node_coords(inode);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->node_coords(inode);
        case (arbitrary) :
            return static_cast<Arbitrary_Submesh*>(submesh())->node_coords(inode);
        case (unassigned) :
            return -999.0;
    }
    return -999.0;
}

/**
 * @brief Computes grading ration inside submesh block
 * @param[in] Submesh object (on device)
 * @param[in] local node ID inside block
 * @return Value of grading ratio at given node
 */
KOKKOS_INLINE_FUNCTION
double grading_ratio(DevicePointer<Submesh> submesh, int inode){
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->grading_ratio(inode);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->grading_ratio(inode);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->grading_ratio(inode);
        case (arbitrary) :
            return static_cast<Arbitrary_Submesh*>(submesh())->grading_ratio(inode);
        case (unassigned) :
            return -999.0;
    }
    return -999.0;
}

/**
 * @brief Fetches block activity info
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @param[in] x2-submesh ID
 * @return boolean on activity status of block
 */
KOKKOS_INLINE_FUNCTION
bool is_block_active(MBBL pumi_obj, int isub, int jsub){
    return pumi_obj.mesh.isactive(isub,jsub);
}

/**
 * @brief Fetches number of x1 elements in domain
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @return number of x1 elements in domain
 */
KOKKOS_INLINE_FUNCTION
int get_num_x1_elems_in_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.submesh_x1(isubmesh)()->Nel;
}

/**
 * @brief Fetches number of x2 elements in domain
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @return number of x2 elements in domain
 */
KOKKOS_INLINE_FUNCTION
int get_num_x2_elems_in_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.submesh_x2(isubmesh)()->Nel;
}

/**
 * @brief Fetches number of x1 elements in all preceding blocks
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @return number of x1 elements in all preceding blocks
 */
KOKKOS_INLINE_FUNCTION
int get_num_x1_elems_before_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.submesh_x1(isubmesh)()->Nel_cumulative;
}

/**
 * @brief Fetches number of x2 elements in all preceding blocks
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @return number of x2 elements in all preceding blocks
 */
KOKKOS_INLINE_FUNCTION
int get_num_x2_elems_before_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.submesh_x2(isubmesh)()->Nel_cumulative;
}

/**
 * @brief Fetches number of x1 elements in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @return number of x1 elements in block
 */
KOKKOS_INLINE_FUNCTION
double get_x1_elem_size_in_submesh(MBBL pumi_obj, int isub, int icell){
    return elem_size(pumi_obj.submesh_x1(isub),icell);
}

/**
 * @brief Fetches number of x1 elements in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @return number of x1 elements in block
 */
KOKKOS_INLINE_FUNCTION
double get_x2_elem_size_in_submesh(MBBL pumi_obj, int isub, int icell){
    return elem_size(pumi_obj.submesh_x2(isub),icell);
}

/**
 * @brief Fetches x1 node-coord of node in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @param[in] local x1 node ID in block
 * @return x1 node coordinate of queried node
 */
KOKKOS_INLINE_FUNCTION
double get_x1_node_coord_in_submesh(MBBL pumi_obj, int isub, int icell){
    return node_coords(pumi_obj.submesh_x1(isub),icell);
}

/**
 * @brief Fetches x2 node-coord of node in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @param[in] local x2 node ID in block
 * @return x1 node coordinate of queried node
 */
KOKKOS_INLINE_FUNCTION
double get_x2_node_coord_in_submesh(MBBL pumi_obj, int isub, int icell){
    return node_coords(pumi_obj.submesh_x2(isub),icell);
}

KOKKOS_INLINE_FUNCTION
double get_x1_gradingratio_in_submesh(MBBL pumi_obj, int isub){
    if (pumi_obj.submesh_x1(isub)()-> meshtype & maxBL){
        return 1.0/pumi_obj.submesh_x1(isub)()->r;
    }
    else{
        return pumi_obj.submesh_x1(isub)()->r;
    }
}

/**
 * @brief Fetches x1 grading ratio in around given node
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @param[in] local x1 node ID in block
 * @return x1 grading ratio in around given node
 */
KOKKOS_INLINE_FUNCTION
double get_x1_gradingratio_in_submesh(MBBL pumi_obj, int isub, int inode){
    return grading_ratio(pumi_obj.submesh_x1(isub),inode);
}

KOKKOS_INLINE_FUNCTION
double get_x2_gradingratio_in_submesh(MBBL pumi_obj, int isub){
    if (pumi_obj.submesh_x2(isub)()-> meshtype & maxBL){
        return 1.0/pumi_obj.submesh_x2(isub)()->r;
    }
    else{
        return pumi_obj.submesh_x2(isub)()->r;
    }
}

/**
 * @brief Fetches x2 grading ratio in around given node
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @param[in] local x2 node ID in block
 * @return x1 grading ratio in around given node
 */
KOKKOS_INLINE_FUNCTION
double get_x2_gradingratio_in_submesh(MBBL pumi_obj, int isub, int inode){
    return grading_ratio(pumi_obj.submesh_x2(isub),inode);
}

/**
 * @brief Converts local cell ID to directional global ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @param[in] local x1 cell ID in block
 * @return x1 directional global ID
 */
KOKKOS_INLINE_FUNCTION
int get_x1_cellID(MBBL pumi_obj, int isub, int icell){
    return icell + pumi_obj.submesh_x1(isub)()->Nel_cumulative;
}

/**
 * @brief Converts local cell ID to directional global ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @param[in] local x2 cell ID in block
 * @return x1 directional global ID
 */
KOKKOS_INLINE_FUNCTION
int get_x2_cellID(MBBL pumi_obj, int isub, int icell){
    return icell + pumi_obj.submesh_x2(isub)()->Nel_cumulative;
}

/**
 * @brief Fetches directional global node ID at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x1 directional global node ID at block interfaces/ends
 */
KOKKOS_INLINE_FUNCTION
int get_x1_nodeID_at_interface(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.if_x1_node(if_node);
}

/**
 * @brief Fetches directional global node ID at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x2 directional global node ID at block interfaces/ends
 */
KOKKOS_INLINE_FUNCTION
int get_x2_nodeID_at_interface(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.if_x2_node(if_node);
}

/**
 * @brief Fetches directional grading ratio at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x1 grading ratio at block interfaces/ends
 */
KOKKOS_INLINE_FUNCTION
double get_x1_gradingratio_at_interface(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.if_x1_r(if_node-1);
}

/**
 * @brief Fetches directional grading ratio at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x2 grading ratio at block interfaces/ends
 */
KOKKOS_INLINE_FUNCTION
double get_x2_gradingratio_at_interface(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.if_x2_r(if_node-1);
}

/**
 * @brief Computes directional block IDs from flattened block ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] flattened block ID
 * @param[out] x1 block ID
 * @param[out] x2 block ID
 */
KOKKOS_INLINE_FUNCTION
void get_directional_submeshID(MBBL pumi_obj, int submeshID, int *isub, int *jsub){
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - pumi_obj.mesh.nsubmesh_x1*(*jsub-1) + 1;
}

/**
 * @brief Computes directional block IDs and local cell IDs from flattened block and local cell ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] flattened block ID
 * @param[in] flattened local cell ID
 * @param[out] x1 block ID
 * @param[out] x1 local cell ID
 * @param[out] x2 block ID
 * @param[out] x2 local cell ID
 */
KOKKOS_INLINE_FUNCTION
void get_directional_submeshID_and_cellID(MBBL pumi_obj, int submeshID, int cellID, int* isub, int *icell, int* jsub, int *jcell){
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - pumi_obj.mesh.nsubmesh_x1*(*jsub-1) + 1;
    *jcell = cellID/pumi_obj.submesh_x1(*isub)()->Nel;
    *icell = cellID - pumi_obj.submesh_x1(*isub)()->Nel*(*jcell);
}

/**
 * @brief Computes directional local node IDs from directional block IDs and interior node IDs
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 block ID
 * @param[in] x2 block ID
 * @param[in] interior local node ID
 * @param[out] x1 local node ID
 * @param[out] x2 local node ID
 */
KOKKOS_INLINE_FUNCTION
void get_directional_interior_nodeIDs(MBBL pumi_obj, int isub, int , int inode, int *inp, int *jnp){
    *jnp = inode/(pumi_obj.submesh_x1(isub)()->Nel-1) + 1;
    *inp = inode - (*jnp-1)*(pumi_obj.submesh_x1(isub)()->Nel-1) + 1;
}

/**
 * @brief Computes flattened submesh ID and local cell ID from directional submesh IDs and local cell IDs
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 block ID
 * @param[in] x1 local cell ID
 * @param[in] x2 block ID
 * @param[in] x2 local cell ID
 * @param[out] flattedned block ID
 * @param[out] flattedned local cell ID
 */
KOKKOS_INLINE_FUNCTION
void flatten_submeshID_and_cellID(MBBL pumi_obj, int isub, int icell, int jsub, int jcell, int* submeshID, int* cellID){
    *submeshID = (isub-1) + (jsub-1)*pumi_obj.mesh.nsubmesh_x1;
    *cellID = icell + jcell*pumi_obj.submesh_x1(isub)()->Nel;
}

/**
 * @brief Fetches edge bdry normal vector for given edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return edge bdry normal vector
 */
KOKKOS_INLINE_FUNCTION
Vector3 get_edge_normal(MBBL pumi_obj, int iEdge){
    return pumi_obj.mesh.bdry.bdry_edge_normal(iEdge);
}

/**
 * @brief Fetches bdry normal vector for given vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return vertex bdry normal vector
 */
KOKKOS_INLINE_FUNCTION
Vector3 get_vert_normal(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.bdry.bdry_vert_normal(iVert);
}

/**
 * @brief Fetches global node ID for given vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return global node ID
 */
KOKKOS_INLINE_FUNCTION
int get_block_vert_nodeID(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.vert_nodeID(iVert);
}

/**
 * @brief Fetches flattened submesh ID for given vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return flattened submesh ID
 */
KOKKOS_INLINE_FUNCTION
int get_block_vert_submeshID(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.vert_subID(iVert);
}

/**
 * @brief performs binary search
 * @param[in] array of cumulative number of nodes upto certain block/edge
 * @param[in] first index in window
 * @param[in] last index in window
 * @param[in] node ID to be located
 * @return active block/edge ID
 */
KOKKOS_INLINE_FUNCTION
int bst_search(Kokkos::View<int*> arr, int first, int last, int nodeID){
    int mid = (first+last)/2;

    if (last == first+1){
        if (arr(first) > nodeID){
            return first;
        }
        else if (arr(first) <= nodeID){
            return last;
        }
    }
    else{
        if (arr(mid) > nodeID){
            last = mid;
            return bst_search(arr, first, last, nodeID);
        }
        else if (arr(mid) < nodeID){
            first = mid;
            return bst_search(arr, first, last, nodeID);
        }
        else{
            return mid+1;
        }
    }
    return -1;
}

/**
 * @brief performs binary search to locate directional submesh ID
 * for give global directional element ID
 * @param[in] Submesh object on device
 * @param[in] first index in window
 * @param[in] last index in window
 * @param[in] entity ID to be located
 * @return directional block ID
 */
KOKKOS_INLINE_FUNCTION
int directional_bst_search(SubmeshDeviceViewPtr submesh, int first, int last, int entity_ID){
    int mid = (first+last)/2;

    if (last == first+1){
        if (submesh(first+1)()->Nel+submesh(first+1)()->Nel_cumulative > entity_ID){
            return first;
        }
        else if (submesh(first+1)()->Nel+submesh(first+1)()->Nel_cumulative <= entity_ID){
            return last;
        }
    }
    else{
        if (submesh(mid+1)()->Nel+submesh(mid+1)()->Nel_cumulative > entity_ID){
            last = mid;
            return directional_bst_search(submesh, first, last, entity_ID);
        }
        else if (submesh(mid+1)()->Nel+submesh(mid+1)()->Nel_cumulative < entity_ID){
            first = mid;
            return directional_bst_search(submesh, first, last, entity_ID);
        }
        else{
            return mid+1;
        }
    }
    return -1;
}

/**
 * @brief Locate x1 submesh and local cell ID for given x1 global cell ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 global cell ID
 * @param[out] x1 submesh ID
 * @param[out] x1 local cell ID
 */
KOKKOS_INLINE_FUNCTION
void get_x1_submeshID_and_localcellID_of_x1_elem(MBBL pumi_obj, int x1_elem_id, int *isub, int *icell){
    int nblks = pumi_obj.mesh.nsubmesh_x1;
    if (nblks>1){
        int first = 0;
        int last = nblks-1;
        *isub = directional_bst_search(pumi_obj.submesh_x1,first,last,x1_elem_id) + 1;
    }
    else{
        *isub = 1;
    }
    *icell = x1_elem_id - pumi_obj.submesh_x1(*isub)()->Nel_cumulative;
}

/**
 * @brief Locate x1 submesh and local cell ID (on the min side) for given x1 global node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 global node ID
 * @param[out] x1 submesh ID
 * @param[out] x1 local cell ID
 */
KOKKOS_INLINE_FUNCTION
void get_x1_submeshID_and_localcellID_of_x1_node(MBBL pumi_obj, int x1_node_id, int *isub, int *icell){
    int nblks = pumi_obj.mesh.nsubmesh_x1;
    if (x1_node_id == pumi_obj.mesh.Nel_tot_x1){
        *isub = nblks;
        *icell = pumi_obj.submesh_x1(*isub)()->Nel-1;
        return;
    }
    else{
        if (nblks>1){
            int first = 0;
            int last = nblks-1;
            *isub = directional_bst_search(pumi_obj.submesh_x1,first,last,x1_node_id) + 1;
        }
        else{
            *isub = 1;
        }
        *icell = x1_node_id - pumi_obj.submesh_x1(*isub)()->Nel_cumulative;
        return;
    }
}

/**
 * @brief Locate x2 submesh and local cell ID for given x2 global cell ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2 global cell ID
 * @param[out] x2 submesh ID
 * @param[out] x2 local cell ID
 */
KOKKOS_INLINE_FUNCTION
void get_x2_submeshID_and_localcellID_of_x2_elem(MBBL pumi_obj, int x2_elem_id, int *isub, int *icell){
    int nblks = pumi_obj.mesh.nsubmesh_x2;
    if (nblks>1){
        int first = 0;
        int last = nblks-1;
        *isub = directional_bst_search(pumi_obj.submesh_x2,first,last,x2_elem_id) + 1;
    }
    else{
        *isub = 1;
    }
    *icell = x2_elem_id - pumi_obj.submesh_x2(*isub)()->Nel_cumulative;
}

/**
 * @brief Locate x2 submesh and local cell ID (on the min side) for given x2 global node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2 global node ID
 * @param[out] x2 submesh ID
 * @param[out] x2 local cell ID
 */
KOKKOS_INLINE_FUNCTION
void get_x2_submeshID_and_localcellID_of_x2_node(MBBL pumi_obj, int x2_node_id, int *isub, int *icell){
    int nblks = pumi_obj.mesh.nsubmesh_x2;
    if (x2_node_id == pumi_obj.mesh.Nel_tot_x2){
        *isub = nblks;
        *icell = pumi_obj.submesh_x2(*isub)()->Nel-1;
        return;
    }
    else{
        if (nblks>1){
            int first = 0;
            int last = nblks-1;
            *isub = directional_bst_search(pumi_obj.submesh_x2,first,last,x2_node_id) + 1;
        }
        else{
            *isub = 1;
        }
        *icell = x2_node_id - pumi_obj.submesh_x2(*isub)()->Nel_cumulative;
        return;
    }
}

/**
 * @brief Locate submesh IDs and local cell IDs for given global block element ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global block element ID
 * @param[out] x1 submesh ID
 * @param[out] x2 submesh ID
 * @param[out] x2 local cell ID
 * @param[out] x1 local cell ID
 */
KOKKOS_INLINE_FUNCTION
void get_submeshIDs_and_localcellIDs_of_block_elements(MBBL pumi_obj, int ielem, int *isub, int *jsub, int *icell, int *jcell){
    int nblks = pumi_obj.mesh.bst.total_active_blocks;
    int subID;
    if (nblks>1){
        int first = 0;
        int last = nblks-1;
        subID = bst_search(pumi_obj.mesh.bst.block_elems_cumulative,first,last,ielem);
    }
    else{
        subID = 0;
    }
    int submeshID = pumi_obj.mesh.bst.active_blockID(subID);
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - (*jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;

    int ielem_loc;
    if (subID==0){
        ielem_loc = ielem;
    }
    else{
        ielem_loc = ielem - pumi_obj.mesh.bst.block_elems_cumulative(subID-1);
    }
    *jcell = ielem_loc/pumi_obj.submesh_x1(*isub)()->Nel;
    *icell = ielem_loc - (*jcell)*(pumi_obj.submesh_x1(*isub)()->Nel);
}

/**
 * @brief Locate submesh IDs for given global block-interior node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global block-interior node ID
 * @param[out] x1 submesh ID
 * @param[out] x2 submesh ID
 */
KOKKOS_INLINE_FUNCTION
void get_submeshIDs_of_block_interior_nodes(MBBL pumi_obj, int inode, int *isub, int *jsub){
    int nblks = pumi_obj.mesh.bst.total_active_blocks;
    int subID;
    if (nblks>1){
        int first = 0;
        int last = nblks-1;
        subID = bst_search(pumi_obj.mesh.bst.block_nodes_cumulative,first,last,inode);
    }
    else{
        subID = 0;
    }
    int submeshID = pumi_obj.mesh.bst.active_blockID(subID);
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - (*jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
}

/**
 * @brief Locate submesh IDs and local node IDs for given global block-interior node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global block-interior node ID
 * @param[out] x1 submesh ID
 * @param[out] x2 submesh ID
 * @param[out] x1 local node ID
 * @param[out] x2 local node ID
 */
KOKKOS_INLINE_FUNCTION
void get_submeshIDs_and_localnodeIDs_of_block_interior_nodes(MBBL pumi_obj, int inode, int *isub, int *jsub, int *inp, int *jnp){
    int nblks = pumi_obj.mesh.bst.total_active_blocks;
    int subID;
    if (nblks>1){
        int first = 0;
        int last = nblks-1;
        subID = bst_search(pumi_obj.mesh.bst.block_nodes_cumulative,first,last,inode);
    }
    else{
        subID = 0;
    }
    int submeshID = pumi_obj.mesh.bst.active_blockID(subID);
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - (*jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
    int inode_loc;
    if (subID==0){
        inode_loc = inode;
    }
    else{
        inode_loc = inode - pumi_obj.mesh.bst.block_nodes_cumulative(subID-1);
    }
    *jnp = inode_loc/(pumi_obj.submesh_x1(*isub)()->Nel-1) + 1;
    *inp = inode_loc - (*jnp-1)*(pumi_obj.submesh_x1(*isub)()->Nel-1) + 1;
}

/**
 * @brief Locate edge IDs for given global block-interior node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global block-interior node ID
 * @return edge ID
 */
KOKKOS_INLINE_FUNCTION
int get_edgeIDs_of_block_edge_interior_nodes(MBBL pumi_obj, int inode){
    int nedges = pumi_obj.mesh.bst.total_active_edges;
    int first = 0;
    int last = nedges-1;
    int edgID = bst_search(pumi_obj.mesh.bst.edge_nodes_cumulative,first,last,inode);
    return pumi_obj.mesh.bst.active_edgeID(edgID);
}

/**
 * @brief Locate edge IDs and edge-local node IDs for given global edge-interior node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global edge-interior node ID
 * @param[out] edge ID
 * @param[out] x1 submesh ID
 * @param[out] x2 submesh ID
 * @param[out] edge local node ID
 */
KOKKOS_INLINE_FUNCTION
void get_edgeIDs_submeshIDs_and_localnodeIDs_of_block_edge_interior_nodes
            (MBBL pumi_obj, int inode, int *iEdge, int *isub, int *jsub, int *inp){
    int nedges = pumi_obj.mesh.bst.total_active_edges;
    int first = 0;
    int last = nedges-1;
    int edgID = bst_search(pumi_obj.mesh.bst.edge_nodes_cumulative,first,last,inode);
    *iEdge = pumi_obj.mesh.bst.active_edgeID(edgID);
    int subID = pumi_obj.mesh.blkif.edge_subID(*iEdge);
    if (edgID == 0){
        *inp = inode;
    }
    else{
        *inp = inode - pumi_obj.mesh.bst.edge_nodes_cumulative(edgID-1);
    }
    *jsub = subID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = subID - (*jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
}

/**
 * @brief Checks if given edge ID is a edge oriented horizontally
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return true for horizontal and false for vertical edges
 */
KOKKOS_INLINE_FUNCTION
bool is_horizontal_edge(MBBL pumi_obj, int iEdge){
    int Nx = pumi_obj.mesh.nsubmesh_x1;

    int num = iEdge/(2*Nx+1);
    int rem = iEdge - num*(2*Nx+1);

    if (rem < Nx){
        return true;
    }
    else{
        return false;
    }
}

/**
 * @brief Checks if given edge ID is a boundary edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return true for bdry and false for non-bdry edges
 */
KOKKOS_INLINE_FUNCTION
bool is_edge_bdry(MBBL pumi_obj, int iEdge){
    return pumi_obj.mesh.bdry.is_bdry_edge(iEdge);
}

/**
 * @brief Checks if given vertex ID is a boundary vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return true for bdry and false for non-bdry vertices
 */
KOKKOS_INLINE_FUNCTION
bool is_vert_bdry(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.bdry.is_bdry_vert(iVert);
}

/**
* @brief Locate the submesh ID and local cell ID for a given x1-coordinate
* Uses analytical formulae to locate the input coordinate
* @param[in] Object of the wrapper mesh structure
* @param[in] x1-coordinate to be located
* @param[out] located x1-submesh ID
* @param[out] located x1-localcell ID
*/
KOKKOS_INLINE_FUNCTION
void locate_submesh_and_cell_x1(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    int nsubmesh = pumi_obj.mesh.nsubmesh_x1;
    for (isubmesh=1; isubmesh<=nsubmesh; isubmesh++){
     if (q >= (pumi_obj.submesh_x1(isubmesh)()->xmin) && q <= (pumi_obj.submesh_x1(isubmesh)()->xmax)){
         *submeshID = isubmesh;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = -1;
     *cellID = -1;
     return;
    }
    *cellID = locate_cell(pumi_obj.submesh_x1(*submeshID),q);
}

/**
* @brief Locate the submesh ID and local cell ID for a given x2-coordinate
* Uses analytical formulae to locate the input coordinate
* @param[in] Object of the wrapper mesh structure
* @param[in] x2-coordinate to be located
* @param[out] located x2-submesh ID
* @param[out] located x2-localcell ID
*/
KOKKOS_INLINE_FUNCTION
void locate_submesh_and_cell_x2(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    int nsubmesh = pumi_obj.mesh.nsubmesh_x2;
    for (isubmesh=1; isubmesh<=nsubmesh; isubmesh++){
     if (q >= (pumi_obj.submesh_x2(isubmesh)()->xmin) && q <= (pumi_obj.submesh_x2(isubmesh)()->xmax)){
         *submeshID = isubmesh;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
        *submeshID = -1;
        *cellID = -1;
        return;
    }
    *cellID  = locate_cell(pumi_obj.submesh_x2(*submeshID),q);
}

/**
 * @brief Update the submesh ID and local cell ID for a given x1-coordinate
 * based on previous submesh and cell IDs.
 * Uses adjacency search to update the IDs
 * @param[in] Object of the wrapper mesh structure
 * @param[in] new x1-coordinate
 * @param[in] old x1-submesh ID
 * @param[in] old x1-localcell ID
 * @param[out] updated x1-submesh ID
 * @param[out] updated x1-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_cell_x1(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x1(*submeshID)()->xmin)){
        *submeshID -= 1;
        prev_cellID = pumi_obj.submesh_x1(*submeshID)()->Nel - 1;
    }

    while(q > (pumi_obj.submesh_x1(*submeshID)()->xmax)){
        *submeshID += 1;
        prev_cellID = 0;
    }
    *cellID = update_cell(pumi_obj.submesh_x1(*submeshID), q, prev_cellID);
}


/**
 * @brief Update the submesh ID and local cell ID for a given x2-coordinate
 * based on previous submesh and cell IDs.
 * Uses adjacency search to update the IDs
 * @param[in] Object of the wrapper mesh structure
 * @param[in] new x2-coordinate
 * @param[in] old x2-submesh ID
 * @param[in] old x2-localcell ID
 * @param[out] updated x2-submesh ID
 * @param[out] updated x2-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_cell_x2(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x2(*submeshID)()->xmin)){
        *submeshID -= 1;
        prev_cellID = pumi_obj.submesh_x2(*submeshID)()->Nel - 1;
    }

    while(q > (pumi_obj.submesh_x2(*submeshID)()->xmax)){
        *submeshID += 1;
        prev_cellID = 0;
    }
    *cellID = update_cell(pumi_obj.submesh_x2(*submeshID), q, prev_cellID);
}

/**
 * @brief Update the submesh ID based on previous submesh and
 * locate the local cell ID
 * Uses adjacency search to update the submesh IDs and analytical
 * forumlae for cell IDs. Use this function when BL coords are
 * not stored
 * @param[in] Object of the wrapper mesh structure
 * @param[in] new x1-coordinate
 * @param[in] old x1-submesh ID
 * @param[out] updated x1-submesh ID
 * @param[out] updated x1-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_locate_cell_x1(MBBL pumi_obj, double q, int prev_submeshID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x1(*submeshID)()->xmin)){
        *submeshID -= 1;
    }

    while(q > (pumi_obj.submesh_x1(*submeshID)()->xmax)){
        *submeshID += 1;
    }

    *cellID = locate_cell(pumi_obj.submesh_x1(*submeshID),q);
}


/**
 * @brief Update the submesh ID based on previous submesh and
 * locate the local cell ID
 * Uses adjacency search to update the submesh IDs and analytical
 * forumlae for cell IDs. Use this function when BL coords are
 * not stored
 * @param[in] Object of the wrapper mesh structure
 * @param[in] new x2-coordinate
 * @param[in] old x2-submesh ID
 * @param[out] updated x2-submesh ID
 * @param[out] updated x2-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_locate_cell_x2(MBBL pumi_obj, double q, int prev_submeshID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x2(*submeshID)()->xmin)){
        *submeshID -= 1;
    }

    while(q > (pumi_obj.submesh_x2(*submeshID)()->xmax)){
        *submeshID += 1;
    }

    *cellID = locate_cell(pumi_obj.submesh_x2(*submeshID),q);
}


/**
 * @brief Computes the partial weights (correspoding to node on the max-side i.e right side)
 * for a located particle coordinate and the global directional cell ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-coordinate of the particle
 * @param[in] x1-submesh ID of the particle
 * @param[in] x1-localcell ID of the particle
 * @param[out] global cell ID in x1 direction
 * @param[out] partial weight
 */
KOKKOS_INLINE_FUNCTION
void calc_weights_x1(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2){
    calc_weights(pumi_obj.submesh_x1(isubmesh),q,icell,x1_global_cell,Wgh2);
}

/**
 * @brief Computes the partial weights (correspoding to node on the max-side i.e top side)
 * for a located particle coordinate and the global directional cell ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-coordinate of the particle
 * @param[in] x2-submesh ID of the particle
 * @param[in] x2-localcell ID of the particle
 * @param[out] global cell ID in x2 direction
 * @param[out] partial weight
 */
 KOKKOS_INLINE_FUNCTION
 void calc_weights_x2(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2){
     calc_weights(pumi_obj.submesh_x2(isubmesh),q,icell,x2_global_cell,Wgh2);
 }

 /**
  * @brief Computes the fractional 2D weights
  * @param[in] Object of the wrapper mesh structure
  * @param[in] x1-coordinate of the particle
  * @param[in] x2-coordinate of the particle
  * @param[in] x1-submesh ID of the particle
  * @param[in] x2-submesh ID of the particle
  * @param[in] x1-localcell ID of the particle
  * @param[in] x2-localcell ID of the particle
  * @param[out] global cell ID in x1 direction
  * @param[out] global cell ID in x2 direction
  * @param[out] fractional weight for bottom-left node
  * @param[out] fractional weight for bottom-right node
  * @param[out] fractional weight for top-left node
  * @param[out] fractional weight for top-left node
  */
 KOKKOS_INLINE_FUNCTION
 void calc_weights_2D(MBBL pumi_obj,
                        double q1, double q2,
                        int isub, int jsub,
                        int icell, int jcell,
                        int *x1_global_cell, int *x2_global_cell,
                        double* w1, double* w2, double* w3, double* w4){


    double x1_min_node = node_coords(pumi_obj.submesh_x1(isub),icell);
    double x2_min_node = node_coords(pumi_obj.submesh_x2(jsub),jcell);
    double x1_max_node = node_coords(pumi_obj.submesh_x1(isub),icell+1);
    double x2_max_node = node_coords(pumi_obj.submesh_x2(jsub),jcell+1);

    *x1_global_cell = pumi_obj.submesh_x1(isub)()->Nel_cumulative + icell;
    *x2_global_cell = pumi_obj.submesh_x2(jsub)()->Nel_cumulative + jcell;

    double elemsize = (x1_max_node-x1_min_node)*(x2_max_node-x2_min_node);
    *w1 = (x1_max_node-q1)*(x2_max_node-q2)/elemsize;
    *w2 = (-x1_min_node+q1)*(x2_max_node-q2)/elemsize;
    *w3 = (x1_max_node-q1)*(-x2_min_node+q2)/elemsize;
    *w4 = (-x1_min_node+q1)*(-x2_min_node+q2)/elemsize;
}

 /**
  * @brief Computes the gloabl cell ID and node ID in 2D for a full Mesh
  * with no-inactive blocks (mesh with inactive blocks will need separate implementations)
  * @param[in] Object of the wrapper mesh structure
  * @param[in] global cell ID in x1-direction
  * @param[in] global cell ID in x2-direction
  * @param[out] global cell ID in 2D
  * @param[out] global node ID of the node in left-bottom corner
  * @param[out] global node ID of the node in left-top coner
  */
  KOKKOS_INLINE_FUNCTION
  void calc_global_cellID_and_nodeID_fullmesh(MBBL pumi_obj, int kcell_x1, int kcell_x2, int *global_cell_2D, int *bottomleft_node, int *topleft_node){
      *global_cell_2D = kcell_x1 + kcell_x2*pumi_obj.mesh.Nel_tot_x1;
      *bottomleft_node = *global_cell_2D + kcell_x2;
      *topleft_node = *bottomleft_node + pumi_obj.mesh.Nel_tot_x1 + 1;
  }

/**
* @brief Computes the gloabl cell ID and node ID in 2D
* @param[in] Object of the wrapper mesh structure
* @param[in] submesh ID in x1-direction
* @param[in] submesh ID in x2-direction
* @param[in] global cell ID in x1-direction
* @param[in] global cell ID in x2-direction
* @param[out] global cell ID in 2D
* @param[out] global node ID of the node in left-bottom corner
* @param[out] global node ID of the node in left-top coner
*/
KOKKOS_INLINE_FUNCTION
void calc_global_cellID_and_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node){
    int icell_x2 = kcell_x2 - pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    int elemoffset = pumi_obj.mesh.offsets.elemoffset_start(isubmesh,jsubmesh) + icell_x2*pumi_obj.mesh.offsets.elemoffset_skip(jsubmesh);
    int fullmesh_elem = kcell_x1 + kcell_x2*pumi_obj.mesh.Nel_tot_x1;
    *global_cell_2D = fullmesh_elem - elemoffset;
    int nodeoffset_bottom = pumi_obj.mesh.offsets.nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh.offsets.nodeoffset_skip_bot(isubmesh,jsubmesh)
                    +(icell_x2-1)*pumi_obj.mesh.offsets.nodeoffset_skip_mid(isubmesh,jsubmesh);
    int nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh.offsets.nodeoffset_skip_mid(isubmesh,jsubmesh);
    if (icell_x2==0){
        nodeoffset_bottom = pumi_obj.mesh.offsets.nodeoffset_start(isubmesh,jsubmesh);
        nodeoffset_top = pumi_obj.mesh.offsets.nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh.offsets.nodeoffset_skip_bot(isubmesh,jsubmesh);
    }
    if (icell_x2==pumi_obj.submesh_x2(jsubmesh)()->Nel-1){
        nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh.offsets.nodeoffset_skip_top(isubmesh,jsubmesh);
    }
    *bottomleft_node = fullmesh_elem + kcell_x2 - nodeoffset_bottom;
    *topleft_node = fullmesh_elem + kcell_x2 + pumi_obj.mesh.Nel_tot_x1 + 1 - nodeoffset_top;
}

/**
* @brief Computes the gloabl cell ID and node ID in 2D
* @param[in] Object of the wrapper mesh structure
* @param[in] submesh ID in x1-direction
* @param[in] submesh ID in x2-direction
* @param[in] global node ID in x1-direction
* @param[in] global node ID in x2-direction
* @return global 2D node ID
*/
KOKKOS_INLINE_FUNCTION
int calc_global_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int inp, int jnp){
    int Inp = inp + pumi_obj.submesh_x1(isubmesh)()->Nel_cumulative;
    int Jnp = jnp + pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    int nodeID = Jnp*(pumi_obj.mesh.Nel_tot_x1+1) + Inp;
    int nodeoffset;
    nodeoffset = pumi_obj.mesh.offsets.nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh.offsets.nodeoffset_skip_bot(isubmesh,jsubmesh)
                    +(jnp-1)*pumi_obj.mesh.offsets.nodeoffset_skip_mid(isubmesh,jsubmesh);
    if (jnp==0){
        nodeoffset = pumi_obj.mesh.offsets.nodeoffset_start(isubmesh,jsubmesh);
    }
    if (jnp==pumi_obj.submesh_x2(jsubmesh)()->Nel){
        nodeoffset +=  (pumi_obj.mesh.offsets.nodeoffset_skip_top(isubmesh,jsubmesh)-pumi_obj.mesh.offsets.nodeoffset_skip_mid(isubmesh,jsubmesh));
    }
    return nodeID-nodeoffset;
}

/**
* @brief Computes the global node ID for nodes on a horizontal edge
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @param[in] edge-local node ID
* @return global 2D node ID
*/
KOKKOS_INLINE_FUNCTION
int calc_global_nodeID_on_horizontal_edge(MBBL pumi_obj, int iEdge, int inode){
    int nodeID = pumi_obj.mesh.blkif.edge_first_nodeID[iEdge]+inode;
    return nodeID;
}

/**
* @brief Computes the global node ID immediate north of given node (for nodes on a horizontal edge)
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @param[in] edge-local node ID
* @return global 2D node ID of immediate northern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_north_global_nodeID_to_horizontal_edge(MBBL pumi_obj, int iEdge, int inode){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jsub = iEdge/(2*Nx+1) + 1;
    int isub = iEdge - (jsub-1)*(2*Nx+1) + 1;
    int nodeID = pumi_obj.mesh.blkif.edge_first_nodeID[iEdge]+inode;
    int nodeoffset = pumi_obj.mesh.Nel_tot_x1+1-pumi_obj.mesh.offsets.nodeoffset_skip_bot(isub,jsub);
    return nodeID+nodeoffset;
}

/**
* @brief Computes the global node ID second north of given node (for nodes on a horizontal edge)
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @param[in] edge-local node ID
* @return global 2D node ID of second northern node
*/
KOKKOS_INLINE_FUNCTION
int calc_second_north_global_nodeID_to_horizontal_edge(MBBL pumi_obj, int iEdge, int inode){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jsub = iEdge/(2*Nx+1) + 1;
    int isub = iEdge - (jsub-1)*(2*Nx+1) + 1;
    int nodeID = pumi_obj.mesh.blkif.edge_first_nodeID[iEdge]+inode;
    int nodeoffset = 2*(pumi_obj.mesh.Nel_tot_x1+1)
                    -pumi_obj.mesh.offsets.nodeoffset_skip_bot(isub,jsub)
                    -pumi_obj.mesh.offsets.nodeoffset_skip_mid(isub,jsub);
    return nodeID+nodeoffset;
}

/**
* @brief Computes the global node ID immediate south of given node (for nodes on a horizontal edge)
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @param[in] edge-local node ID
* @return global 2D node ID of immediate southern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_south_global_nodeID_to_horizontal_edge(MBBL pumi_obj, int iEdge, int inode){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jsub = iEdge/(2*Nx+1) + 1;
    int isub = iEdge - (jsub-1)*(2*Nx+1) + 1;
    int nodeID = pumi_obj.mesh.blkif.edge_first_nodeID[iEdge]+inode;
    int nodeoffset = pumi_obj.mesh.Nel_tot_x1+1-pumi_obj.mesh.offsets.nodeoffset_skip_top(isub,jsub-1);
    return nodeID-nodeoffset;
}

/**
* @brief Computes the global node ID second south of given node (for nodes on a horizontal edge)
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @param[in] edge-local node ID
* @return global 2D node ID of second southern node
*/
KOKKOS_INLINE_FUNCTION
int calc_second_south_global_nodeID_to_horizontal_edge(MBBL pumi_obj, int iEdge, int inode){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jsub = iEdge/(2*Nx+1) + 1;
    int isub = iEdge - (jsub-1)*(2*Nx+1) + 1;
    int nodeID = pumi_obj.mesh.blkif.edge_first_nodeID[iEdge]+inode;
    int nodeoffset = 2*(pumi_obj.mesh.Nel_tot_x1+1)
                    -pumi_obj.mesh.offsets.nodeoffset_skip_top(isub,jsub-1)
                    -pumi_obj.mesh.offsets.nodeoffset_skip_mid(isub,jsub-1);
    return nodeID-nodeoffset;
}

/**
* @brief Computes the global node ID for nodes on a vertical edge
* @param[in] Object of the wrapper mesh structure
* @param[in] vertical edge ID
* @param[in] edge-local node ID
* @return global 2D node ID
*/
KOKKOS_INLINE_FUNCTION
int calc_global_nodeID_on_vertical_edge(MBBL pumi_obj, int iEdge, int inode){
    int nodeID = pumi_obj.mesh.blkif.edge_first_nodeID(iEdge);
    int subID = pumi_obj.mesh.blkif.edge_subID(iEdge);
    int jsub = subID/pumi_obj.mesh.nsubmesh_x1 + 1;
    int isub = subID - (jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
    nodeID += inode*(pumi_obj.mesh.Nel_tot_x1+1-pumi_obj.mesh.offsets.nodeoffset_skip_mid(isub,jsub));
    return nodeID;
}

/**
* @brief Computes the global node ID immediate north of given node (for nodes on a vertical edge)
* @param[in] Object of the wrapper mesh structure
* @param[in] vertical edge ID
* @param[in] edge-local node ID
* @return global 2D node ID of immediate northern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_north_global_nodeID_to_vertical_edge(MBBL pumi_obj, int iEdge, int inode){
    int subID = pumi_obj.mesh.blkif.edge_subID(iEdge);
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jsub = subID/Nx + 1;
    int nel_blk = pumi_obj.submesh_x2(jsub)()->Nel;
    if (inode==nel_blk-2){
        int jm1 = iEdge/(2*Nx+1);
        int im1 = iEdge - jm1*(2*Nx+1) - Nx;
        int iVert = (jm1+1)*(Nx+1)+im1;
        return pumi_obj.mesh.blkif.vert_nodeID(iVert);
    }
    else {
        return calc_global_nodeID_on_vertical_edge(pumi_obj, iEdge, inode+1);
    }
}

/**
* @brief Computes the global node ID immediate south of given node (for nodes on a vertical edge)
* @param[in] Object of the wrapper mesh structure
* @param[in] vertical edge ID
* @param[in] edge-local node ID
* @return global 2D node ID of immediate southern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_south_global_nodeID_to_vertical_edge(MBBL pumi_obj, int iEdge, int inode){
    // int subID = pumi_obj.mesh.blkif.edge_subID(iEdge);
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    if (inode==0){
        int jm1 = iEdge/(2*Nx+1);
        int im1 = iEdge - jm1*(2*Nx+1) - Nx;
        int iVert = jm1*(Nx+1)+im1;
        return pumi_obj.mesh.blkif.vert_nodeID(iVert);
    }
    else {
        return calc_global_nodeID_on_vertical_edge(pumi_obj, iEdge, inode-1);
    }
}

/**
* @brief Computes the global node ID immediate north of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of immediate northern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_north_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    int im1 = iVert - jm1*(Nx+1);
    int iEdge = jm1*(2*Nx+1)+im1+Nx;
    return calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,0);
}

/**
* @brief Computes the global node ID second north node of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of second northern node to vertex
*/
KOKKOS_INLINE_FUNCTION
int calc_second_north_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    int im1 = iVert - jm1*(Nx+1);
    int iEdge = jm1*(2*Nx+1)+im1+Nx;
    return calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,1);
}

/**
* @brief Computes the global node ID immediate south of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of immediate southern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_south_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    int im1 = iVert - jm1*(Nx+1);
    int iEdge = (jm1-1)*(2*Nx+1)+im1+Nx;
    int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jm1);
    return calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,nel_blk-2);
}

/**
* @brief Computes the global node ID second south node of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of second southern node to vertex
*/
KOKKOS_INLINE_FUNCTION
int calc_second_south_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    int im1 = iVert - jm1*(Nx+1);
    int iEdge = (jm1-1)*(2*Nx+1)+im1+Nx;
    int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jm1);
    return calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,nel_blk-3);
}

/**
* @brief Computes the global node ID immediate east of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of immediate eastern node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_east_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.vert_nodeID(iVert)+1;
}

/**
* @brief Computes the global node ID second east node of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of second eastern node to vertex
*/
KOKKOS_INLINE_FUNCTION
int calc_second_east_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.vert_nodeID(iVert)+2;
}

/**
* @brief Computes the global node ID immediate west of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of immediate western node
*/
KOKKOS_INLINE_FUNCTION
int calc_first_west_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.vert_nodeID(iVert)-1;
}

/**
* @brief Computes the global node ID second west node of given vertex node
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return global 2D node ID of second western node to vertex
*/
KOKKOS_INLINE_FUNCTION
int calc_second_west_global_nodeID_to_vertex(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.vert_nodeID(iVert)-2;
}

/**
* @brief Computes the x1 submesh ID west of vertex
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return x1 submesh ID west of vertex
*/
KOKKOS_INLINE_FUNCTION
int get_x1_submeshID_west_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    int im1 = iVert - jm1*(Nx+1);
    return im1;
}

/**
* @brief Computes the x1 submesh ID east of vertex
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return x1 submesh ID east of vertex
*/
KOKKOS_INLINE_FUNCTION
int get_x1_submeshID_east_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    int im1 = iVert - jm1*(Nx+1);
    return im1+1;
}

/**
* @brief Computes the x2 submesh ID south of vertex
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return x2 submesh ID south of vertex
*/
KOKKOS_INLINE_FUNCTION
int get_x2_submeshID_south_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    return jm1;
}

/**
* @brief Computes the x2 submesh ID north of vertex
* @param[in] Object of the wrapper mesh structure
* @param[in] vertex ID
* @return x2 submesh ID north of vertex
*/
KOKKOS_INLINE_FUNCTION
int get_x2_submeshID_north_to_vertex(MBBL pumi_obj, int iVert){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iVert/(Nx+1);
    return jm1+1;
}

/**
* @brief Computes the x2 submesh ID north of horizontal edge
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @return x2 submesh ID north of horizontal edge
*/
KOKKOS_INLINE_FUNCTION
int get_x2_submeshID_north_to_horizontal_edge(MBBL pumi_obj, int iEdge){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iEdge/(2*Nx+1);
    return jm1+1;
}

/**
* @brief Computes the x2 submesh ID south of horizontal edge
* @param[in] Object of the wrapper mesh structure
* @param[in] horizontal edge ID
* @return x2 submesh ID south of horizontal edge
*/
KOKKOS_INLINE_FUNCTION
int get_x2_submeshID_south_to_horizontal_edge(MBBL pumi_obj, int iEdge){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iEdge/(2*Nx+1);
    return jm1;
}

/**
* @brief Computes the x1 submesh ID east of vertical edge
* @param[in] Object of the wrapper mesh structure
* @param[in] vertical edge ID
* @return x1 submesh ID north of vertical edge
*/
KOKKOS_INLINE_FUNCTION
int get_x1_submeshID_east_to_vertical_edge(MBBL pumi_obj, int iEdge){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iEdge/(2*Nx+1);
    int im1 = iEdge - (jm1*(2*Nx+1)+Nx);
    return im1+1;
}

/**
* @brief Computes the x1 submesh ID east of horizontal edge
* @param[in] Object of the wrapper mesh structure
* @param[in] vertical edge ID
* @return x1 submesh ID north of vertical edge
*/
KOKKOS_INLINE_FUNCTION
int get_x1_submeshID_west_to_vertical_edge(MBBL pumi_obj, int iEdge){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int jm1 = iEdge/(2*Nx+1);
    int im1 = iEdge - (jm1*(2*Nx+1)+Nx);
    return im1;
}

/**
* @brief Perform particle push and return new coordinate
* @param[in] Object of the wrapper mesh structure
* @param[in] initial coordinate vector
* @param[in] displacement vector
* @param[in,out] x1-submesh ID
* @param[in,out] x2-submesh ID
* @param[in,out] x1-localcell ID
* @param[in,out] x2-localcell ID
* @param[out] is particle still in domain
* @param[out] boundary edge ID crossed by particle
* @param[out] fraction of push completed
* @param[out] face ID thru which particle crossed the boundary
* @return new position vector of particle
*/
KOKKOS_INLINE_FUNCTION
Vector3 push_particle(MBBL pumi_obj, Vector3 q, Vector3 dq,
                   int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain,
                   int *bdry_hit, double *fraction_done, int *faceID_on_bdry){

    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    double dq1 = dq[0];
    double dq2 = dq[1];
    double dq3 = dq[2];
    double q1_new = q1+dq1;
    double q2_new = q2+dq2;
    double q3_new = q3+dq3;
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Nxx = 2*Nx+1;
    // int Ny = pumi_obj.mesh.nsubmesh_x2;
    int case_id = (dq2>=0.0)+2*(dq1>=0.0);
    int isub = *isubmesh;
    int jsub = *jsubmesh;
    *in_domain = true;
    double eps = 1e-12;
    if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
        && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){

        *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
        *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
        *bdry_hit = -1;
        *faceID_on_bdry = -1;
        *fraction_done = 1.0;
        return Vector3(q1_new, q2_new, q3_new);
    }
    else{
        double del1, del2;
        bool located = false;
        switch (case_id) {
            case 0:
                del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                    isub--;
                    *icell = pumi_obj.submesh_x1(isub)()->Nel-1;
                }
                else{
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                    jsub--;
                    *jcell = pumi_obj.submesh_x2(jsub)()->Nel-1;
                }


                while (!located && in_domain){
                    if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                        *in_domain = false;
                        int num = *bdry_hit/Nxx;
                        int rem = *bdry_hit - num*Nxx;
                        if (rem < Nx){
                            *fraction_done = fabs(del2/dq2);
                            *isubmesh = isub;
                            *jsubmesh = jsub+1;
                            *jcell = 0;
                            q1_new = q1+(*fraction_done)*dq1;
                            q2_new = pumi_obj.submesh_x2(*jsubmesh)()->xmin + eps;
                            q3_new = q3+(*fraction_done)*dq3;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *icell;
                        }
                        else{
                            *fraction_done = fabs(del1/dq1);
                            *isubmesh = isub+1;
                            *jsubmesh = jsub;
                            *icell = 0;
                            q1_new = pumi_obj.submesh_x1(*isubmesh)()->xmin+eps;
                            q2_new = q2+(*fraction_done)*dq2;
                            q3_new = q3+(*fraction_done)*dq3;
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *jcell;
                        }
                        return Vector3(q1_new, q2_new, q3_new);
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            *faceID_on_bdry = -1;
                            *fraction_done = 1.0;
                            located = true;
                            return Vector3(q1_new, q2_new, q3_new);
                        }
                        else{
                            del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                            del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                                isub--;
                                *icell = pumi_obj.submesh_x1(isub)()->Nel-1;
                            }
                            else{
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                                jsub--;
                                *jcell = pumi_obj.submesh_x2(jsub)()->Nel-1;
                            }
                        }
                    }
                }

            case 1:
                del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                    isub--;
                    *icell = pumi_obj.submesh_x1(isub)()->Nel-1;
                }
                else{
                    *bdry_hit = jsub*(Nxx)+isub-1;
                    jsub++;
                    *jcell = 0;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                        *in_domain = false;
                        int num = *bdry_hit/Nxx;
                        int rem = *bdry_hit - num*Nxx;
                        if (rem < Nx){
                            *fraction_done = fabs(del2/dq2);
                            *isubmesh = isub;
                            *jsubmesh = jsub-1;
                            *jcell = pumi_obj.submesh_x2(*jsubmesh)()->Nel-1;
                            q1_new = q1+(*fraction_done)*dq1;
                            q2_new = pumi_obj.submesh_x2(*jsubmesh)()->xmax - eps;
                            q3_new = q3+(*fraction_done)*dq3;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *icell;
                        }
                        else{
                            *fraction_done = fabs(del1/dq1);
                            *isubmesh = isub+1;
                            *jsubmesh = jsub;
                            *icell = 0;
                            q1_new = pumi_obj.submesh_x1(*isubmesh)()->xmin+eps ;
                            q2_new = q2+(*fraction_done)*dq2;
                            q3_new = q3+(*fraction_done)*dq3;
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *jcell;
                        }
                        return Vector3(q1_new, q2_new, q3_new);
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            *faceID_on_bdry = -1;
                            located = true;
                            *fraction_done = 1.0;
                            return Vector3(q1_new, q2_new, q3_new);
                        }
                        else{
                            del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                            del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                                isub--;
                                *icell = pumi_obj.submesh_x1(isub)()->Nel-1;
                            }
                            else{
                                *bdry_hit = jsub*(Nxx)+isub-1;
                                jsub++;
                                *jcell = 0;
                            }
                        }
                    }
                }

            case 2:
                del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                    isub++;
                    *icell = 0;
                }
                else{
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                    jsub--;
                    *jcell = pumi_obj.submesh_x2(jsub)()->Nel-1;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                        *in_domain = false;
                        int num = *bdry_hit/Nxx;
                        int rem = *bdry_hit - num*Nxx;
                        if (rem < Nx){
                            *fraction_done = fabs(del2/dq2);
                            *isubmesh = isub;
                            *jsubmesh = jsub+1;
                            *jcell = 0;
                            q1_new = q1+(*fraction_done)*dq1;
                            q2_new = pumi_obj.submesh_x2(*jsubmesh)()->xmin + eps;
                            q3_new = q3+(*fraction_done)*dq3;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *icell;
                        }
                        else{
                            *fraction_done = fabs(del1/dq1);
                            *isubmesh = isub-1;
                            *jsubmesh = jsub;
                            *icell = pumi_obj.submesh_x1(*isubmesh)()->Nel-1;
                            q1_new = pumi_obj.submesh_x1(*isubmesh)()->xmax - eps;
                            q2_new = q2+(*fraction_done)*dq2;
                            q3_new = q3+(*fraction_done)*dq3;
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *jcell;
                        }
                        return Vector3(q1_new, q2_new, q3_new);
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            *faceID_on_bdry = -1;
                            located = true;
                            *fraction_done = 1.0;
                            return Vector3(q1_new, q2_new, q3_new);
                        }
                        else{
                            del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                            del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                                isub++;
                                *icell = 0;
                            }
                            else{
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                                jsub--;
                                *jcell = pumi_obj.submesh_x2(jsub)()->Nel-1;
                            }
                        }
                    }
                }

            case 3:
                del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                    isub++;
                    *icell = 0;
                }
                else{
                    *bdry_hit = jsub*(Nxx)+isub-1;
                    jsub++;
                    *jcell = 0;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                        *in_domain = false;
                        int num = *bdry_hit/Nxx;
                        int rem = *bdry_hit - num*Nxx;
                        if (rem < Nx){
                            *fraction_done = fabs(del2/dq2);
                            *isubmesh = isub;
                            *jsubmesh = jsub-1;
                            *jcell = pumi_obj.submesh_x2(*jsubmesh)()->Nel-1;
                            q1_new = q1+(*fraction_done)*dq1;
                            q2_new = pumi_obj.submesh_x2(*jsubmesh)()->xmax - eps;
                            q3_new = q3+(*fraction_done)*dq3;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *icell;
                        }
                        else{
                            *fraction_done = fabs(del1/dq1);
                            *isubmesh = isub-1;
                            *jsubmesh = jsub;
                            *icell = pumi_obj.submesh_x1(*isubmesh)()->Nel-1;
                            q1_new = pumi_obj.submesh_x1(*isubmesh)()->xmax - eps;
                            q2_new = q2+(*fraction_done)*dq2;
                            q3_new = q3+(*fraction_done)*dq3;
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *faceID_on_bdry = pumi_obj.mesh.bdry.edge_to_face(*bdry_hit) + *jcell;
                        }
                        return Vector3(q1_new, q2_new, q3_new);
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            *faceID_on_bdry = -1;
                            located = true;
                            *fraction_done = 1.0;
                            return Vector3(q1_new, q2_new, q3_new);
                        }
                        else{
                            del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                            del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                                isub++;
                                *icell = 0;
                            }
                            else{
                                *bdry_hit = jsub*(Nxx)+isub-1;
                                jsub++;
                                *jcell = 0;
                            }
                        }
                    }
                }
        }

    }

    return Vector3(-999.0,-999.0,-999.0);

}

KOKKOS_INLINE_FUNCTION
void push_particle_v2(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit, double *fraction_done){

    double q1_new = q1+dq1;
    double q2_new = q2+dq2;
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Nxx = 2*Nx+1;
    // int Ny = pumi_obj.mesh.nsubmesh_x2;
    int isub = *isubmesh;
    int jsub = *jsubmesh;

    int num_x1_crossed = 0;
    int num_x2_crossed = 0;
    int x1_sub_move = 0;
    int x2_sub_move = 0;
    while(q1_new < (pumi_obj.submesh_x1(isub)()->xmin)){
        isub--;
        num_x1_crossed++;
        x1_sub_move = 2;
        *icell = pumi_obj.submesh_x1(isub)()->Nel-1;
    }
    while(q1_new > (pumi_obj.submesh_x1(isub)()->xmax)){
        isub++;
        num_x1_crossed++;
        x1_sub_move = 1;
        *icell = 0;
    }
    while(q2_new < (pumi_obj.submesh_x2(jsub)()->xmin)){
        jsub--;
        num_x2_crossed++;
        x2_sub_move = 2;
        *jcell = pumi_obj.submesh_x2(jsub)()->Nel-1;
    }
    while(q2_new > (pumi_obj.submesh_x2(jsub)()->xmax)){
        jsub++;
        num_x2_crossed++;
        x2_sub_move = 1;
        *jcell = 0;
    }


    int case_id = x1_sub_move + 3*x2_sub_move;
    // printf("case=%d    \n",case_id );
    double del1, del2;
    int i;

    int isub_tmp = *isubmesh;
    int jsub_tmp = *jsubmesh;

    *isubmesh = isub;
    *jsubmesh = jsub;

    switch (case_id) {
        case 0:
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *isubmesh = isub;
            *jsubmesh = jsub;
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 1:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
            i=0;
            while (i<num_x1_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    *fraction_done = fabs((pumi_obj.submesh_x1(isub_tmp)()->xmax-q1)/dq1);
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d (%d,%d)\n", *bdry_hit,*isubmesh+1,*jsubmesh);
                    isub_tmp += 1;
                    *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *bdry_hit = -1;
            *in_domain = true;
            *fraction_done = 1.0;
            return;

        case 2:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
            i=0;
            while (i<num_x1_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    *fraction_done = fabs((q1-pumi_obj.submesh_x1(isub_tmp)()->xmin)/dq1);
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh-1,*jsubmesh);
                    isub_tmp -= 1;
                    *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 3:
            *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
            i=0;
            while (i<num_x2_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    *fraction_done = fabs((pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2)/dq2);
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh+1);
                    jsub_tmp += 1;
                    *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 4:
            del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
            del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                isub_tmp += 1;
            }
            else{
                *bdry_hit = jsub_tmp*(Nxx)+isub_tmp-1;
                jsub_tmp += 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    int num = *bdry_hit/Nxx;
                    int rem = *bdry_hit - num*Nxx;
                    if (rem < Nx){
                        *fraction_done = fabs(del2/dq2);
                    }
                    else{
                        *fraction_done = fabs(del1/dq1);
                    }
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
                    del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                        isub_tmp += 1;
                    }
                    else{
                        *bdry_hit = jsub_tmp*(Nxx)+isub_tmp-1;
                        jsub_tmp += 1;
                    }
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 5:
            del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
            del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                isub_tmp -= 1;
            }
            else{
                *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
                jsub_tmp += 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    int num = *bdry_hit/Nxx;
                    int rem = *bdry_hit - num*Nxx;
                    if (rem < Nx){
                        *fraction_done = fabs(del2/dq2);
                    }
                    else{
                        *fraction_done = fabs(del1/dq1);
                    }
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
                    del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                        isub_tmp -= 1;
                    }
                    else{
                        *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
                        jsub_tmp += 1;
                    }
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 6:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
            i=0;
            while (i<num_x2_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    *fraction_done = fabs((q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin)/dq2);
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh-1);
                    jsub_tmp -= 1;
                    *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 7:
            del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
            del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);
            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                isub_tmp += 1;
            }
            else{
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                jsub_tmp -= 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    int num = *bdry_hit/Nxx;
                    int rem = *bdry_hit - num*Nxx;
                    if (rem < Nx){
                        *fraction_done = fabs(del2/dq2);
                    }
                    else{
                        *fraction_done = fabs(del1/dq1);
                    }
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
                    del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);
                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                        isub_tmp += 1;
                    }
                    else{
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                        jsub_tmp -= 1;
                    }
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

        case 8:
            del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
            del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);

            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                isub_tmp -= 1;
            }
            else{
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                jsub_tmp -= 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh.bdry.is_bdry_edge(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    int num = *bdry_hit/Nxx;
                    int rem = *bdry_hit - num*Nxx;
                    if (rem < Nx){
                        *fraction_done = fabs(del2/dq2);
                    }
                    else{
                        *fraction_done = fabs(del1/dq1);
                    }
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
                    del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);

                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                        isub_tmp -= 1;
                    }
                    else{
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                        jsub_tmp -= 1;
                    }
                }
            }
            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
            *in_domain = true;
            *bdry_hit = -1;
            *fraction_done = 1.0;
            return;

    }

}

} // namespace pumi

#endif
