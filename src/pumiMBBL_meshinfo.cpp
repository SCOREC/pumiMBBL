#include "pumiMBBL_meshinfo.hpp"

namespace pumi {

/**
 * @brief Returns the grading ratio along a queried direction about a queried node
 *
 * @param[in] Object of the wrapper mesh structure
 * @param[in] direction along which grading ratio is needed
 * @param[in] global directional node ID
 * \return grading ratio about the node in the requested direction
 */
double return_gradingratio(MBBL pumi_obj, int dir, int node){

    SubmeshHostViewPtr h_submesh;
    int nsubmesh, Nel_total;
    if (dir == x1_dir){
        nsubmesh = pumi_obj.mesh.nsubmesh_x1;
        Nel_total = pumi_obj.mesh.Nel_tot_x1;
        h_submesh = pumi_obj.host_submesh_x1;
    }
    else{
        nsubmesh = pumi_obj.mesh.nsubmesh_x2;
        Nel_total = pumi_obj.mesh.Nel_tot_x2;
        h_submesh = pumi_obj.host_submesh_x2;
    }
    if (node == 0 || node == Nel_total){
        printf("Grading ratio not defined for the first and last node of the domain -- Terminating \n");
        exit(0);
    }
    else{
        for (int isubmesh=1; isubmesh<=nsubmesh; isubmesh++){
            int submesh_min_node =  h_submesh[isubmesh]->Nel_cumulative;
            int submesh_max_node = submesh_min_node + h_submesh[isubmesh]->Nel;
            if (node > submesh_min_node && node < submesh_max_node){
                if (h_submesh[isubmesh]->meshtype & maxBL){
                    return 1.0/h_submesh[isubmesh]->r;
                }
                else if (h_submesh[isubmesh]->meshtype & arbitrary){
                    int inode = node - submesh_min_node;
                    double max_elem = h_submesh[isubmesh]->host_BL_coords[inode+1]-h_submesh[isubmesh]->host_BL_coords[inode];
                    double min_elem = h_submesh[isubmesh]->host_BL_coords[inode]-h_submesh[isubmesh]->host_BL_coords[inode-1];
                    return max_elem/min_elem;
                }
                else{
                    return h_submesh[isubmesh]->r;
                }
            }
            else if (node == submesh_min_node){
                double max_elem;
                double min_elem;
                if (h_submesh[isubmesh-1]->meshtype & minBL){
                    min_elem = h_submesh[isubmesh-1]->t0 * pow(h_submesh[isubmesh-1]->r, h_submesh[isubmesh-1]->Nel-1);
                }
                else if (h_submesh[isubmesh-1]->meshtype & arbitrary){
                    int Nel = h_submesh[isubmesh-1]->Nel;
                    min_elem = h_submesh[isubmesh-1]->host_BL_coords[Nel]-h_submesh[isubmesh-1]->host_BL_coords[Nel-1];
                }
                else{
                    min_elem = h_submesh[isubmesh-1]->t0;
                }

                if (h_submesh[isubmesh]->meshtype & maxBL){
                    max_elem = h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r, h_submesh[isubmesh]->Nel-1);
                }
                else if (h_submesh[isubmesh]->meshtype & arbitrary){
                    max_elem = h_submesh[isubmesh]->host_BL_coords[1]-h_submesh[isubmesh]->host_BL_coords[0];
                }
                else{
                    max_elem = h_submesh[isubmesh]->t0;
                }
                return max_elem/min_elem;
            }
        }
    }
    return -999.0;
}

/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * @param[in] Object of the wrapper mesh structure
 * @param[in] direction along which element-size is needed
 * @param[in] global directional node ID or element ID index
 * @param[in] relevant offset enum based on provied index
 * \return element size for the queried element in the queried direction
 */
double return_elemsize(MBBL pumi_obj, int dir, int index, int offset){

    SubmeshHostViewPtr h_submesh;
    int nsubmesh, Nel_total, elem;
    if (dir == x1_dir){
        nsubmesh = pumi_obj.mesh.nsubmesh_x1;
        Nel_total = pumi_obj.mesh.Nel_tot_x1;
        h_submesh = pumi_obj.host_submesh_x1;
    }
    else{
        nsubmesh = pumi_obj.mesh.nsubmesh_x2;
        Nel_total = pumi_obj.mesh.Nel_tot_x2;
        h_submesh = pumi_obj.host_submesh_x2;
    }

    elem = index + offset;
    if (elem >= Nel_total){
        elem = Nel_total-1;
    }
    else if (elem < 0){
        elem = 0;
    }

    for (int isubmesh = 1; isubmesh <= nsubmesh; isubmesh++){
        int submesh_min_elem = h_submesh[isubmesh]->Nel_cumulative;
        int submesh_max_elem = submesh_min_elem + h_submesh[isubmesh]->Nel-1;

        if (elem >= submesh_min_elem && elem <= submesh_max_elem){
            if (h_submesh[isubmesh]->meshtype & uniform){
                return h_submesh[isubmesh]->t0;
            }
            else{
                int local_cell = elem - submesh_min_elem;
                return (h_submesh[isubmesh]->host_BL_coords[local_cell+1]-h_submesh[isubmesh]->host_BL_coords[local_cell]);
            }
            // if (h_submesh[isubmesh]->meshtype & minBL){
            //     int local_cell = elem - submesh_min_elem;
            //     return (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
            // }
            // if (h_submesh[isubmesh]->meshtype & maxBL){
            //     int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
            //     return (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
            // }
        }
    }
    return -999.0;
}

/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global node IDs along x1-direction
 * \return covolume for the requested node
 */
double return_covolume(MBBL pumi_obj, int inode_x1){

    double covolume, dx1_min, dx1_max;

    if (inode_x1 == 0){
        dx1_min = 0.0;
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
    }
    else if (inode_x1 == pumi_obj.mesh.Nel_tot_x1){
        dx1_max = 0.0;
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }
    else{
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }

    covolume = (dx1_min+dx1_max)/2.0;
    return covolume;
}


/**
 * @brief Returns the element size along a queried direction for a queried node/element
 * Only when all blocks are active
 *
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global node IDs along x1-direction
 * @param[in] global node IDs along x2-direction
 * \return covolume for the requested node
 */
double return_covolume_fullmesh(MBBL pumi_obj, int inode_x1, int inode_x2){

    double covolume, dx1_min, dx1_max, dx2_min, dx2_max;

    if (inode_x1 == 0){
        dx1_min = 0.0;
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
    }
    else if (inode_x1 == pumi_obj.mesh.Nel_tot_x1){
        dx1_max = 0.0;
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }
    else{
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }

    if (inode_x2 == 0){
        dx2_min = 0.0;
        dx2_max = return_elemsize(pumi_obj, x2_dir, elem_on_max_side_offset, inode_x2);
    }
    else if (inode_x2 == pumi_obj.mesh.Nel_tot_x2){
        dx2_max = 0.0;
        dx2_min = return_elemsize(pumi_obj, x2_dir, elem_on_min_side_offset, inode_x2);
    }
    else{
        dx2_max = return_elemsize(pumi_obj, x2_dir, elem_on_max_side_offset, inode_x2);
        dx2_min = return_elemsize(pumi_obj, x2_dir, elem_on_min_side_offset, inode_x2);
    }

    covolume = (dx1_min*dx2_min + dx1_max*dx2_min + dx1_min*dx2_max + dx1_max*dx2_max)/4.0;
    return covolume;
}

/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global node IDs along x1-direction
 * @param[in] global node IDs along x2-direction
 * \return covolume for the requested node
 */
double return_covolume(MBBL pumi_obj, int inode_x1, int inode_x2){

    int isub_min, jsub_min, isub_max, jsub_max;
    double dx1_min, dx1_max, dx2_min, dx2_max, covolume;

    dx1_min=0.0;
    dx2_min=0.0;
    dx1_max=0.0;
    dx2_max=0.0;
    covolume=0.0;
    isub_min=0;
    isub_max=0;
    jsub_min=0;
    jsub_max=0;

    SubmeshHostViewPtr h_submesh;
    int nsubmesh, Nel_total, elem;

    nsubmesh = pumi_obj.mesh.nsubmesh_x1;
    Nel_total = pumi_obj.mesh.Nel_tot_x1;
    h_submesh = pumi_obj.host_submesh_x1;

    elem = inode_x1-1;

    if (elem < 0){
        dx1_min = 0.0;
        isub_min = 0;
    }
    else{
        for (int isubmesh = 1; isubmesh <= nsubmesh; isubmesh++){
            int submesh_min_elem = h_submesh[isubmesh]->Nel_cumulative;
            int submesh_max_elem = submesh_min_elem + h_submesh[isubmesh]->Nel-1;

            if (elem >= submesh_min_elem && elem <= submesh_max_elem){
                isub_min = isubmesh;
                if (h_submesh[isubmesh]->meshtype & uniform){
                    dx1_min = h_submesh[isubmesh]->t0;
                }
                else{
                    int local_cell = elem - submesh_min_elem;
                    dx1_min = h_submesh[isubmesh]->host_BL_coords[local_cell+1]-h_submesh[isubmesh]->host_BL_coords[local_cell];
                }
                // if (h_submesh[isubmesh]->meshtype & minBL){
                //     int local_cell = elem - submesh_min_elem;
                //     dx1_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
                // if (h_submesh[isubmesh]->meshtype & maxBL){
                //     int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                //     dx1_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
            }
        }
    }

    elem = inode_x1;
    if (elem >= Nel_total){
        dx1_max = 0.0;
        isub_max = nsubmesh+1;
    }
    else{
        for (int isubmesh = 1; isubmesh <= nsubmesh; isubmesh++){
            int submesh_min_elem = h_submesh[isubmesh]->Nel_cumulative;
            int submesh_max_elem = submesh_min_elem + h_submesh[isubmesh]->Nel-1;

            if (elem >= submesh_min_elem && elem <= submesh_max_elem){
                isub_max = isubmesh;
                if (h_submesh[isubmesh]->meshtype & uniform){
                    dx1_max = h_submesh[isubmesh]->t0;
                }
                else{
                    int local_cell = elem - submesh_min_elem;
                    dx1_max = h_submesh[isubmesh]->host_BL_coords[local_cell+1]-h_submesh[isubmesh]->host_BL_coords[local_cell];
                }
                // if (h_submesh[isubmesh]->meshtype & minBL){
                //     int local_cell = elem - submesh_min_elem;
                //     dx1_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
                // if (h_submesh[isubmesh]->meshtype & maxBL){
                //     int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                //     dx1_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
            }
        }
    }

    nsubmesh = pumi_obj.mesh.nsubmesh_x2;
    Nel_total = pumi_obj.mesh.Nel_tot_x2;
    h_submesh = pumi_obj.host_submesh_x2;

    elem = inode_x2-1;

    if (elem < 0){
        dx2_min = 0.0;
        jsub_min = 0;
    }
    else{
        for (int isubmesh = 1; isubmesh <= nsubmesh; isubmesh++){
            int submesh_min_elem = h_submesh[isubmesh]->Nel_cumulative;
            int submesh_max_elem = submesh_min_elem + h_submesh[isubmesh]->Nel-1;

            if (elem >= submesh_min_elem && elem <= submesh_max_elem){
                jsub_min = isubmesh;
                if (h_submesh[isubmesh]->meshtype & uniform){
                    dx2_min = h_submesh[isubmesh]->t0;
                }
                else{
                    int local_cell = elem - submesh_min_elem;
                    dx2_min = h_submesh[isubmesh]->host_BL_coords[local_cell+1]-h_submesh[isubmesh]->host_BL_coords[local_cell];
                }
                // if (h_submesh[isubmesh]->meshtype & minBL){
                //     int local_cell = elem - submesh_min_elem;
                //     dx2_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
                // if (h_submesh[isubmesh]->meshtype & maxBL){
                //     int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                //     dx2_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
            }
        }
    }

    elem = inode_x2;
    if (elem >= Nel_total){
        dx2_max = 0.0;
        jsub_max = nsubmesh+1;
    }
    else{
        for (int isubmesh = 1; isubmesh <= nsubmesh; isubmesh++){
            int submesh_min_elem = h_submesh[isubmesh]->Nel_cumulative;
            int submesh_max_elem = submesh_min_elem + h_submesh[isubmesh]->Nel-1;

            if (elem >= submesh_min_elem && elem <= submesh_max_elem){
                jsub_max = isubmesh;
                if (h_submesh[isubmesh]->meshtype & uniform){
                    dx2_max = h_submesh[isubmesh]->t0;
                }
                else{
                    int local_cell = elem - submesh_min_elem;
                    dx2_max = h_submesh[isubmesh]->host_BL_coords[local_cell+1]-h_submesh[isubmesh]->host_BL_coords[local_cell];
                }
                // if (h_submesh[isubmesh]->meshtype & minBL){
                //     int local_cell = elem - submesh_min_elem;
                //     dx2_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
                // if (h_submesh[isubmesh]->meshtype & maxBL){
                //     int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                //     dx2_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                // }
            }
        }
    }


    covolume = (dx1_min*dx2_min*pumi_obj.mesh.host_isactive[isub_min][jsub_min] +
                dx1_max*dx2_min*pumi_obj.mesh.host_isactive[isub_max][jsub_min] +
                dx1_min*dx2_max*pumi_obj.mesh.host_isactive[isub_min][jsub_max] +
                dx1_max*dx2_max*pumi_obj.mesh.host_isactive[isub_max][jsub_max])/4.0;
    return covolume;

}

/**
 * @brief Provides info about inactivity in mesh
 *
 * @param[in] Object of the wrapper mesh structure
 * \return true for all active blocks, false otherwise
 */
bool is_fullmesh(MBBL pumi_obj){
    return pumi_obj.mesh.offsets.is_fullmesh;
}

/**
 * @brief Returns x1-coordinate of left end of the domain
 *
 * @param[in] Object of the wrapper mesh structure
 */
double get_global_x1_min_coord(MBBL pumi_obj){
    return pumi_obj.host_submesh_x1[1]->xmin;
}

/**
 * @brief Returns x1-coordinate of right end of the domain
 *
 * @param[in] Object of the wrapper mesh structure
 */
double get_global_x1_max_coord(MBBL pumi_obj){
    int nsubmesh = pumi_obj.mesh.nsubmesh_x1;
    return pumi_obj.host_submesh_x1[nsubmesh]->xmax;
}

/**
 * @brief Returns x2-coordinate of bottom end of the domain
 *
 * @param[in] Object of the wrapper mesh structure
 */
double get_global_x2_min_coord(MBBL pumi_obj){
    return pumi_obj.host_submesh_x2[1]->xmin;
}

/**
 * @brief Returns x2-coordinate of top end of the domain
 *
 * @param[in] Object of the wrapper mesh structure
 */
double get_global_x2_max_coord(MBBL pumi_obj){
    int nsubmesh = pumi_obj.mesh.nsubmesh_x2;
    return pumi_obj.host_submesh_x2[nsubmesh]->xmax;
}

/**
 * @brief Fetches number of x1 blocks in domain
 * @param[in] Object of the wrapper mesh structure
 * @return number of x1 blocks in domain
 */
int get_num_x1_submesh(MBBL pumi_obj){
    return pumi_obj.mesh.nsubmesh_x1;
}

/**
 * @brief Fetches number of x1 elements in domain
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @return number of x1 elements in domain
 */
int get_num_x1_elems_in_submesh_host(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x1[isubmesh]->Nel;
}

/**
 * @brief Fetches number of x1 elements in all preceding blocks
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @return number of x1 elements in all preceding blocks
 */
int get_num_x1_elems_before_submesh_host(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
}

/**
 * @brief Fetches number of x1 elements in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @return number of x1 elements in block
 */
double get_x1_elem_size_in_submesh_host(MBBL pumi_obj, int isubmesh, int icell){
    return pumi_obj.host_submesh_x1[isubmesh]->elem_size_host(icell);
}

/**
 * @brief Fetches number of x2 blocks in domain
 * @param[in] Object of the wrapper mesh structure
 * @return number of x2 blocks in domain
 */
int get_num_x2_submesh(MBBL pumi_obj){
    return pumi_obj.mesh.nsubmesh_x2;
}

/**
 * @brief Fetches number of x2 elements in domain
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @return number of x2 elements in domain
 */
int get_num_x2_elems_in_submesh_host(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x2[isubmesh]->Nel;
}

/**
 * @brief Fetches number of x2 elements in all preceding blocks
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @return number of x2 elements in all preceding blocks
 */
int get_num_x2_elems_before_submesh_host(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x2[isubmesh]->Nel_cumulative;
}

/**
 * @brief Fetches number of x1 elements in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @return number of x1 elements in block
 */
double get_x2_elem_size_in_submesh_host(MBBL pumi_obj, int isubmesh, int icell){
    return pumi_obj.host_submesh_x2[isubmesh]->elem_size_host(icell);
}

double get_x1_gradingratio_in_submesh_host(MBBL pumi_obj, int isub){
    if (pumi_obj.host_submesh_x1[isub]->meshtype & maxBL){
        return 1.0/pumi_obj.host_submesh_x1[isub]->r;
    }
    else{
        return pumi_obj.host_submesh_x1[isub]->r;
    }
}

/**
 * @brief Fetches x1 grading ratio in around given node
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @param[in] local x1 node ID in block
 * @return x1 grading ratio in around given node
 */
double get_x1_gradingratio_in_submesh_host(MBBL pumi_obj, int isub, int inode){
    return pumi_obj.host_submesh_x1[isub]->grading_ratio_host(inode);
}

double get_x2_gradingratio_in_submesh_host(MBBL pumi_obj, int isub){
    if (pumi_obj.host_submesh_x2[isub]->meshtype & maxBL){
        return 1.0/pumi_obj.host_submesh_x2[isub]->r;
    }
    else{
        return pumi_obj.host_submesh_x2[isub]->r;
    }
}
/**
 * @brief Fetches x2 grading ratio in around given node
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x2-submesh ID
 * @param[in] local x2 node ID in block
 * @return x1 grading ratio in around given node
 */
double get_x2_gradingratio_in_submesh_host(MBBL pumi_obj, int isub, int inode){
    return pumi_obj.host_submesh_x2[isub]->grading_ratio_host(inode);
}

/**
 * @brief Fetches directional global node ID at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x1 directional global node ID at block interfaces/ends
 */
int get_x1_nodeID_at_interface_host(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.host_if_x1_node[if_node];
}

/**
 * @brief Fetches directional global node ID at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x2 directional global node ID at block interfaces/ends
 */
int get_x2_nodeID_at_interface_host(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.host_if_x2_node[if_node];
}

/**
 * @brief Fetches directional grading ratio at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x1 grading ratio at block interfaces/ends
 */
int get_x1_gradingratio_at_interface_host(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.host_if_x1_r[if_node-1];
}

/**
 * @brief Fetches directional grading ratio at block interfaces/ends
 * @param[in] Object of the wrapper mesh structure
 * @param[in] block interface/end ID
 * @return x2 grading ratio at block interfaces/ends
 */
int get_x2_gadingratio_at_interface_host(MBBL pumi_obj, int if_node){
    return pumi_obj.mesh.blkif.host_if_x2_r[if_node-1];
}

/**
 * @brief Fetches total active mesh elements
 * @param[in] Object of the wrapper mesh structure
 * @return total active mesh elements
 */
int get_total_mesh_elements(MBBL pumi_obj){
    return pumi_obj.mesh.Nel_total;
}

/**
 * @brief Fetches total active mesh nodes
 * @param[in] Object of the wrapper mesh structure
 * @return total active mesh nodes
 */
int get_total_mesh_nodes(MBBL pumi_obj){
    return pumi_obj.mesh.Nnp_total;
}

/**
 * @brief Fetches total x1 elements
 * @param[in] Object of the wrapper mesh structure
 * @return total x1 elements
 */
int get_total_x1_elements(MBBL pumi_obj){
    return pumi_obj.mesh.Nel_tot_x1;
}

/**
 * @brief Fetches total x2 elements
 * @param[in] Object of the wrapper mesh structure
 * @return total x2 elements
 */
int get_total_x2_elements(MBBL pumi_obj){
    return pumi_obj.mesh.Nel_tot_x2;
}

/**
 * @brief Fetches total blocks (active + inactive)
 * @param[in] Object of the wrapper mesh structure
 * @return total blocks
 */
int get_total_submesh_blocks(MBBL pumi_obj){
    if (pumi_obj.mesh.ndim==1){
        return pumi_obj.mesh.nsubmesh_x1;
    }
    else if (pumi_obj.mesh.ndim==2){
        return pumi_obj.mesh.nsubmesh_x1*pumi_obj.mesh.nsubmesh_x2;
    }
    return 0;
}

/**
 * @brief Fetches total elements in a block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] flattened submesh ID
 * @return total elements in a block
 */
int get_total_elements_in_block(MBBL pumi_obj, int flattened_submesh_ID){
    if (pumi_obj.mesh.ndim==1){
        return pumi_obj.host_submesh_x1[flattened_submesh_ID]->Nel;
    }
    else if (pumi_obj.mesh.ndim==2){
        int jsub = flattened_submesh_ID/pumi_obj.mesh.nsubmesh_x1 + 1;
        int isub = flattened_submesh_ID - (jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
        int Nel = pumi_obj.host_submesh_x1[isub]->Nel * pumi_obj.host_submesh_x2[jsub]->Nel;
        return Nel;
    }
    return 0;
}

/**
 * @brief Fetches total block-interior nodes in block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 submesh ID
 * @param[in] x2 submesh ID
 * @return total block-interior nodes in block
 */
int get_num_interior_nodes_on_block(MBBL pumi_obj, int isub, int jsub){
    int Nel_x1 = pumi_obj.host_submesh_x1[isub]->Nel;
    int Nel_x2 = pumi_obj.host_submesh_x2[jsub]->Nel;
    return (Nel_x1-1)*(Nel_x2-1);
}

/**
 * @brief Fetches total volume of domain
 * @param[in] Object of the wrapper mesh structure
 * @return total volume of domain
 */
double get_mesh_volume(MBBL pumi_obj){
    if (pumi_obj.mesh.ndim==1){
        double volume = pumi_obj.host_submesh_x1[pumi_obj.mesh.nsubmesh_x1]->xmax - pumi_obj.host_submesh_x1[1]->xmin;
        return volume;
    }
    else if (pumi_obj.mesh.ndim==2){
        double volume=0.0;
        int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
        int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
        for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
            for (int jsubmesh=1; jsubmesh<=nsubmesh_x2; jsubmesh++){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    volume += pumi_obj.host_submesh_x1[isubmesh]->length * pumi_obj.host_submesh_x2[jsubmesh]->length;
                }
            }
        }
        return volume;
    }
    else {
        return 0.0;
    }
}

/**
 * @brief Fetches edge bdry normal vector for given edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return edge bdry normal vector
 */
Vector3 get_bdry_edge_normal_host(MBBL pumi_obj, int iEdge){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (iEdge>=0 && iEdge<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2){
        Vector3 nrml = pumi_obj.mesh.bdry.host_bdry_edge_normal[iEdge];
        return nrml;
    }
    else{
        std::cout << "Invalid edge ID or direction\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        std::cout << "Valid directions 0 (x1-dirxn) or 1 (x2-dirxn) or 2 (x3-dirxn)\n ";
        exit(0);
    }
}

/**
 * @brief Fetches bdry normal vector for given vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return vertex bdry normal vector
 */
Vector3 get_bdry_vert_normal_host(MBBL pumi_obj, int iVert){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (iVert>=0 && iVert<nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2+1){
        return pumi_obj.mesh.bdry.host_bdry_vert_normal[iVert];
    }
    else{
        std::cout << "Invalid vertex ID\n";
        std::cout << "Valid VertexIDs = [0,1,..," << nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2 <<"]\n";
        exit(0);
    }
}

/**
 * @brief Fetches edge bdry normal vector for given edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return edge bdry normal vector
 */
Vector3 get_edge_normal_host(MBBL pumi_obj, int iEdge){
    return pumi_obj.mesh.bdry.host_bdry_edge_normal[iEdge];
}

/**
 * @brief Fetches bdry normal vector for given vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return vertex bdry normal vector
 */
Vector3 get_vert_normal_host(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.bdry.host_bdry_vert_normal[iVert];
}

/**
 * @brief Checks if given edge ID is a edge oriented horizontally
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return true for horizontal and false for vertical edges
 */
bool is_horizontal_edge_host(MBBL pumi_obj, int iEdge){
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
 * @brief Fetches number of nodes on edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return number of nodes on edge
 */
int get_num_interior_nodes_on_edge(MBBL pumi_obj, int iEdge){
    int Nx = pumi_obj.mesh.nsubmesh_x1;

    int num = iEdge/(2*Nx+1);
    int rem = iEdge - num*(2*Nx+1);

    if (rem < Nx){
        return pumi_obj.host_submesh_x1[rem+1]->Nel-1;
    }
    else{
        return pumi_obj.host_submesh_x2[num+1]->Nel-1;
    }
}

/**
 * @brief Fetches number of faces on edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return number of faces on edge
 */
int get_num_faces_on_edge(MBBL pumi_obj, int iEdge){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    int Nx2p1 = 2*nsubmesh_x1+1;
    if (iEdge>=0 && iEdge<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2){
        int num = iEdge/Nx2p1;
        int rem = iEdge-num*Nx2p1;
        if (rem < nsubmesh_x1){
            return pumi_obj.host_submesh_x1[rem+1]->Nel;
        }
        else {
            return pumi_obj.host_submesh_x2[num+1]->Nel;
        }
    }
    else{
        std::cout << "Invalid edge ID\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        exit(0);
    }
}

/**
 * @brief Fetches starting face ID on given edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return starting face ID
 */
int get_starting_faceID_on_bdry_edge(MBBL pumi_obj, int iEdge){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (iEdge>=0 && iEdge<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2){
        return pumi_obj.mesh.bdry.host_edge_to_face[iEdge];
    }
    else{
        std::cout << "Invalid edge ID\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        exit(0);
    }
}

/**
 * @brief Fetches edge ID on west of block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 submesh ID
 * @param[in] x2 submesh ID
 * @return west edge ID
 */
int get_west_edgeID(MBBL pumi_obj, int isub, int jsub){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    return (jsub-1)*(2*Nx+1)+(isub-1)+Nx;
}

/**
 * @brief Fetches edge ID on east of block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 submesh ID
 * @param[in] x2 submesh ID
 * @return east edge ID
 */
int get_east_edgeID(MBBL pumi_obj, int isub, int jsub){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    return (jsub-1)*(2*Nx+1)+isub+Nx;
}

/**
 * @brief Fetches edge ID on north of block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 submesh ID
 * @param[in] x2 submesh ID
 * @return north edge ID
 */
int get_north_edgeID(MBBL pumi_obj, int isub, int jsub){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    return jsub*(2*Nx+1)+(isub-1);
}

/**
 * @brief Fetches edge ID on south of block
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 submesh ID
 * @param[in] x2 submesh ID
 * @return south edge ID
 */
int get_south_edgeID(MBBL pumi_obj, int isub, int jsub){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    return (jsub-1)*(2*Nx+1)+(isub-1);
}

/**
 * @brief Fetches bdry info for given edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return true for bdry edge, false otherwise
 */
bool is_edge_bdry_host(MBBL pumi_obj, int iEdge){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (iEdge>=0 && iEdge<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2){
        return pumi_obj.mesh.bdry.host_is_bdry_edge[iEdge];
    }
    else{
        std::cout << "Invalid edge ID\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        exit(0);
    }
}

/**
 * @brief Fetches bdry info for given vertex
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return true for bdry vertex, false otherwise
 */
bool is_vert_bdry_host(MBBL pumi_obj, int iVert){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (iVert>=0 && iVert<nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2+1){
        return pumi_obj.mesh.bdry.host_is_bdry_vert[iVert];
    }
    else{
        std::cout << "Invalid vertex ID\n";
        std::cout << "Valid VertexIDs = [0,1,..," << nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2 <<"]\n";
        exit(0);
    }
}

/**
 * @brief Fetches flattened active submesh ID to which given vertex belongs to
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @return flattened active submesh ID
 */
int get_block_vert_submeshID_host(MBBL pumi_obj, int iVert){
    return pumi_obj.mesh.blkif.host_vert_subID[iVert];
}

/**
 * @brief Fetches flattened active submesh ID to which given edge belongs to
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return flattened active submesh ID
 */
int get_block_edge_submeshID_host(MBBL pumi_obj, int iEdge){
    return pumi_obj.mesh.blkif.host_edge_subID[iEdge];
}

/**
 * @brief Fetches submesh IDs to which given vertex belongs to
 * @param[in] Object of the wrapper mesh structure
 * @param[in] vertex ID
 * @param[out] x1-submesh ID
 * @param[out] x2-submesh ID
 */
void get_block_vert_submeshIDs_host(MBBL pumi_obj, int iVert, int *isub, int *jsub){
    int subID = pumi_obj.mesh.blkif.host_vert_subID[iVert];
    *jsub = subID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = subID - (*jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
}

/**
 * @brief Fetches submesh IDs to which given edge belongs to
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @param[out] x1-submesh ID
 * @param[out] x2-submesh ID
 */
void get_block_edge_submeshIDs_host(MBBL pumi_obj, int iEdge, int *isub, int *jsub){
    int subID = pumi_obj.mesh.blkif.host_edge_subID[iEdge];
    *jsub = subID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = subID - (*jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
}

/**
 * @brief Fetches block activity info
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1-submesh ID
 * @param[in] x2-submesh ID
 * @return boolean on activity status of block
 */
bool is_block_active_host(MBBL pumi_obj, int isub, int jsub){
    return pumi_obj.mesh.host_isactive[isub][jsub];
}

/**
 * @brief Fetches block activity info
 * @param[in] Object of the wrapper mesh structure
 * @param[in] flattened submesh ID
 * @return boolean on activity status of block
 */
bool is_block_active_host(MBBL pumi_obj, int flattened_submesh_ID){
    int jsub = flattened_submesh_ID/pumi_obj.mesh.nsubmesh_x1 + 1;
    int isub = flattened_submesh_ID - (jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
    return pumi_obj.mesh.host_isactive[isub][jsub];
}

/**
 * @brief Fetches total edges in mesh (active + inactive)
 * @param[in] Object of the wrapper mesh structure
 * @return total edges in mesh
 */
int get_total_mesh_block_edges(MBBL pumi_obj){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    return 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2;
}

/**
 * @brief Fetches total verts in mesh (active + inactive)
 * @param[in] Object of the wrapper mesh structure
 * @return total verts in mesh
 */
int get_total_mesh_block_verts(MBBL pumi_obj){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    return nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2+1;
}

/**
 * @brief checks validity of edge ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return true for valid ID, false for Invalid id
 */
bool check_edge_index_bounds(MBBL pumi_obj, int iEdge){
    bool valid = true;
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Ny = pumi_obj.mesh.nsubmesh_x2;
    if (iEdge < 0 || iEdge >= 2*Nx*Ny+Nx+Ny)
        valid = false;

    return valid;
}

/**
 * @brief computes 2D global node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 submesh ID
 * @param[in] x2 submesh ID
 * @param[in] flattened node ID
 * @return 2D global node ID
 */
int compute_global_nodeID_2D(MBBL pumi_obj, int isubmesh, int jsubmesh, int fullmesh_node_id){
    int Jnp;
    Jnp = fullmesh_node_id/(pumi_obj.mesh.Nel_tot_x1+1);


    int nodeID = fullmesh_node_id;
    int jnp = Jnp - pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;

    int nodeoffset;
    nodeoffset = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh] +
                 pumi_obj.mesh.offsets.host_nodeoffset_skip_bot[isubmesh][jsubmesh] +
                 (jnp-1)*pumi_obj.mesh.offsets.host_nodeoffset_skip_mid[isubmesh][jsubmesh];
    if (jnp==0){
        nodeoffset = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh];
    }
    if (jnp==pumi_obj.host_submesh_x2[jsubmesh]->Nel){
        nodeoffset +=  (pumi_obj.mesh.offsets.host_nodeoffset_skip_top[isubmesh][jsubmesh]-pumi_obj.mesh.offsets.host_nodeoffset_skip_mid[isubmesh][jsubmesh]);
    }
    return nodeID-nodeoffset;
}

/**
 * @brief computes 2D list of nodes on given bdry-edge
 * @param[in] Object of the wrapper mesh structure
 * @param[in] edge ID
 * @return array of nodes on bdry-edge
 */
std::vector<int> get_nodes_on_bdry_edge(MBBL pumi_obj, int iEdge){
    bool validID = check_edge_index_bounds(pumi_obj, iEdge);
    std::vector<int> nodeIDs;
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (validID){
        int Nx2p1 = 2*nsubmesh_x1+1;

        if (is_edge_bdry(pumi_obj,iEdge)){
            int num = iEdge/Nx2p1;
            int rem = iEdge-num*Nx2p1;
            int isubmesh, jsubmesh;
            if (rem < nsubmesh_x1){
                isubmesh = rem;
                jsubmesh = num;
                if (jsubmesh >= nsubmesh_x2){
                    int Inp = pumi_obj.host_submesh_x1[isubmesh+1]->Nel_cumulative;
                    int Jnp = pumi_obj.mesh.Nel_tot_x2;
                    nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    for (int i=0; i<get_num_faces_on_edge(pumi_obj, iEdge); i++){
                        Inp++;
                        nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    }
                    return nodeIDs;
                }
                else{
                    int Inp = pumi_obj.host_submesh_x1[isubmesh+1]->Nel_cumulative;
                    int Jnp = pumi_obj.host_submesh_x2[jsubmesh+1]->Nel_cumulative;
                    nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    for (int i=0; i<get_num_faces_on_edge(pumi_obj, iEdge); i++){
                        Inp++;
                        nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    }
                    return nodeIDs;
                }
            }
            else {
                jsubmesh = num;
                isubmesh = rem-nsubmesh_x2;
                if (isubmesh >= nsubmesh_x1){
                    int Jnp = pumi_obj.host_submesh_x2[jsubmesh+1]->Nel_cumulative;
                    int Inp = pumi_obj.mesh.Nel_tot_x1;
                    nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    for (int i=0; i<get_num_faces_on_edge(pumi_obj, iEdge); i++){
                        Jnp++;
                        nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    }
                    return nodeIDs;
                }
                else{
                    int Inp = pumi_obj.host_submesh_x1[isubmesh+1]->Nel_cumulative;
                    int Jnp = pumi_obj.host_submesh_x2[jsubmesh+1]->Nel_cumulative;
                    nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    for (int i=0; i<get_num_faces_on_edge(pumi_obj, iEdge); i++){
                        Jnp++;
                        nodeIDs.push_back(get_global_nodeID_2D(pumi_obj, Inp, Jnp));
                    }
                    return nodeIDs;
                }
            }
        }
        else{
            return nodeIDs;
        }
    }
    else{
        std::cout << "Invalid edge ID\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        exit(0);
    }

    return nodeIDs;
}

/**
 * @brief checks validity of directional node IDs
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 global node ID
 * @param[in] x2 global node ID
 * @return true for valid ID, false for Invalid id
 */
bool check_node_index_bounds(MBBL pumi_obj, int knode_x1, int knode_x2){
    bool valid = true;
    if (knode_x1 < 0 || knode_x1 > pumi_obj.mesh.Nel_tot_x1){
        valid = false;
    }
    if (knode_x2 < 0 || knode_x2 > pumi_obj.mesh.Nel_tot_x2){
        valid = false;
    }
    return valid;
}

/**
 * @brief Computes 2D global node ID
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 global node ID
 * @param[in] x2 global node ID
 * @return 2D global node ID
 */
int get_global_nodeID_2D(MBBL pumi_obj, int knode_x1, int knode_x2){

    bool validID = check_node_index_bounds(pumi_obj, knode_x1, knode_x2);

    if (validID){
        int isubmesh, jsubmesh, inp, jnp;
        bool left_edge, right_edge, bottom_edge, top_edge, on_edge;
        int fullmesh_node_id = knode_x1 + knode_x2*(pumi_obj.mesh.Nel_tot_x1+1);
        on_edge = false;
        for (isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
            int submesh_min_node = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
            int submesh_max_node = pumi_obj.host_submesh_x1[isubmesh]->Nel + submesh_min_node;
            left_edge =  false;
            right_edge = false;
            if (knode_x1 >= submesh_min_node && knode_x1 <= submesh_max_node){
                inp = knode_x1 - submesh_min_node;
                if (inp == 0){
                    left_edge = true;
                    on_edge = true;
                }
                if (inp == pumi_obj.host_submesh_x1[isubmesh]->Nel){
                    right_edge = true;
                    on_edge = true;
                }
                break;
            }
        }

        for (jsubmesh=1; jsubmesh<=pumi_obj.mesh.nsubmesh_x2; jsubmesh++){
            int submesh_min_node = pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
            int submesh_max_node = pumi_obj.host_submesh_x2[jsubmesh]->Nel + submesh_min_node;
            bottom_edge =  false;
            top_edge = false;
            if (knode_x2 >= submesh_min_node && knode_x2 <= submesh_max_node){
                jnp = knode_x2 - submesh_min_node;
                if (jnp == 0){
                    bottom_edge = true;
                    on_edge = true;
                }
                if (jnp == pumi_obj.host_submesh_x2[jsubmesh]->Nel){
                    top_edge = true;
                    on_edge = true;
                }
                break;
            }
        }

        if (!on_edge){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
            }
            else{
                return -1;
            }
        }
        else{
            if (!left_edge && !right_edge && !bottom_edge && !top_edge){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                }
                else{
                    return -1;
                }
            }

            if (left_edge & !top_edge & !bottom_edge){

                if (isubmesh==1){
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }
                else{
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]) {
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }

            }

            if (left_edge & top_edge){
                if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
                    if (isubmesh==1){
                        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
                else{
                    if (isubmesh==1){
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh+1, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh+1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh+1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh+1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
            }

            if (top_edge & !left_edge & !right_edge){

                if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }
                else{
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh+1, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }
            }

            if (top_edge & right_edge){
                if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
                    if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
                else{
                    if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh+1, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh+1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh+1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh+1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
            }

            if (right_edge & !top_edge & !bottom_edge){

                if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }
                else{
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }

            }

            if (right_edge & bottom_edge){
                if (jsubmesh==1){
                    if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
                else{
                    if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh-1, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh-1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh-1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh+1, jsubmesh-1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
            }

            if (bottom_edge & !left_edge & !right_edge){

                if (jsubmesh==1){
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);            }
                    else{
                        return -1;
                    }
                }
                else{
                    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                    }
                    else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                        return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh-1, fullmesh_node_id);
                    }
                    else{
                        return -1;
                    }
                }
            }

            if (bottom_edge & left_edge){
                if (jsubmesh==1){
                    if (isubmesh==1){
                        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
                else{
                    if (isubmesh==1){
                        if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh-1, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                    else{
                        if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh-1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh-1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh-1, jsubmesh, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh-1, fullmesh_node_id);
                        }
                        else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                            return compute_global_nodeID_2D(pumi_obj, isubmesh, jsubmesh, fullmesh_node_id);
                        }
                        else{
                            return -1;
                        }
                    }
                }
            }
        }
    }
    else{
        std::cout << "Invalid node index\n";
        exit(0);
    }

    return -1;
}

/**
 * @brief Computes active block ID to which give node belongs to
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 global node ID
 * @param[in] x2 global node ID
 * @return flattened active block ID
 */
int get_node_submeshID(MBBL pumi_obj, int knode_x1, int knode_x2){

    int isubmesh, jsubmesh, inp, jnp;
    bool left_edge, right_edge, bottom_edge, top_edge;

    for (isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
        int submesh_min_node = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
        int submesh_max_node = pumi_obj.host_submesh_x1[isubmesh]->Nel + submesh_min_node;
        left_edge =  false;
        right_edge = false;
        if (knode_x1 >= submesh_min_node && knode_x1 <= submesh_max_node){
            inp = knode_x1 - submesh_min_node;
            if (inp == 0){
                left_edge = true;
            }
            if (inp == pumi_obj.host_submesh_x1[isubmesh]->Nel){
                right_edge = true;
            }
            break;
        }
    }

    for (jsubmesh=1; jsubmesh<=pumi_obj.mesh.nsubmesh_x2; jsubmesh++){
        int submesh_min_node = pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
        int submesh_max_node = pumi_obj.host_submesh_x2[jsubmesh]->Nel + submesh_min_node;
        bottom_edge =  false;
        top_edge = false;
        if (knode_x2 >= submesh_min_node && knode_x2 <= submesh_max_node){
            jnp = knode_x2 - submesh_min_node;
            if (jnp == 0){
                bottom_edge = true;
            }
            if (jnp == pumi_obj.host_submesh_x2[jsubmesh]->Nel){
                top_edge = true;
            }
            break;
        }
    }

    if (!left_edge && !right_edge && !bottom_edge && !top_edge){
        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
            return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
        }
        else{
            return -1;
        }
    }

    if (left_edge & !top_edge & !bottom_edge){

        if (isubmesh==1){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]) {
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                return (isubmesh-2)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }

    }

    if (left_edge & top_edge){
        if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
            if (isubmesh==1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                    return (isubmesh-2)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
        else{
            if (isubmesh==1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                    return (isubmesh-1)+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;;
                }
                else{
                    return -1;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                    return (isubmesh-2)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh+1]){
                    return (isubmesh-2)+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                    return (isubmesh-1)+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
    }

    if (top_edge & !left_edge & !right_edge){

        if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                return (isubmesh-1)+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }
    }

    if (top_edge & right_edge){
        if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                    return isubmesh+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
        else{
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                    return (isubmesh-1)+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                    return (isubmesh-1)+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh+1]){
                    return isubmesh+(jsubmesh)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                    return isubmesh+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
    }

    if (right_edge & !top_edge & !bottom_edge){

        if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                return isubmesh+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }

    }

    if (right_edge & bottom_edge){
        if (jsubmesh==1){
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                    return isubmesh+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
        else{
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                    return (isubmesh-1)+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                    return (isubmesh-1)+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                    return isubmesh+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh-1]){
                    return isubmesh+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
    }

    if (bottom_edge & !left_edge & !right_edge){

        if (jsubmesh==1){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;            }
            else{
                return -1;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
            }
            else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                return (isubmesh-1)+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
            }
            else{
                return -1;
            }
        }
    }

    if (bottom_edge & left_edge){
        if (jsubmesh==1){
            if (isubmesh==1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                    return (isubmesh-2)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
        else{
            if (isubmesh==1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                    return (isubmesh-1)+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh-1]){
                    return (isubmesh-2)+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                    return (isubmesh-2)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                    return (isubmesh-1)+(jsubmesh-2)*pumi_obj.mesh.nsubmesh_x1;
                }
                else if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
                }
                else{
                    return -1;
                }
            }
        }
    }

    return -1;
}

/**
 * @brief Computes active block ID to which give element belongs to
 * @param[in] Object of the wrapper mesh structure
 * @param[in] x1 global element ID
 * @param[in] x2 global element ID
 * @return flattened active block ID
 */
int get_elem_submeshID(MBBL pumi_obj, int kcell_x1, int kcell_x2){
    int isubmesh, jsubmesh;

    for (isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
        int submesh_min_elem = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
        int submesh_max_elem = pumi_obj.host_submesh_x1[isubmesh]->Nel - 1 + submesh_min_elem;
        if (kcell_x1 >= submesh_min_elem && kcell_x1 <= submesh_max_elem){
            break;
        }
    }

    for (jsubmesh=1; jsubmesh<=pumi_obj.mesh.nsubmesh_x2; jsubmesh++){
        int submesh_min_elem = pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
        int submesh_max_elem = pumi_obj.host_submesh_x2[jsubmesh]->Nel - 1 + submesh_min_elem;
        if (kcell_x2 >= submesh_min_elem && kcell_x2 <= submesh_max_elem){
            break;
        }
    }

    if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
        return (isubmesh-1)+(jsubmesh-1)*pumi_obj.mesh.nsubmesh_x1;
    }
    else{
        return -1;
    }

}

/**
 * @brief Fetches number of total block interior nodes
 * @param[in] Object of the wrapper mesh structure
 * @return number of total block interior nodes
 */
int get_num_block_interior_nodes(MBBL pumi_obj){
    return pumi_obj.mesh.bst.total_block_nodes;
}

/**
 * @brief Fetches number of total edge nodes
 * @param[in] Object of the wrapper mesh structure
 * @return number of total edge nodes
 */
int get_num_block_edge_interior_nodes(MBBL pumi_obj){
    return pumi_obj.mesh.bst.total_edge_nodes;
}

/**
 * @brief Returns node info such as if node is in active domain, if node is on a boundary
 * and boundary entity dimension (boundary vertex (dim=0) or edge (dim=1)) and entity tag
 * of the boundary
 *
 * @param[in] Object of the wrapper mesh structure
 * @param[in] global node IDs along x1-direction
 * @param[in] global node IDs along x2-direction
 * @param[out] boolean value if node is on boundary
 * @param[out] boolean value if node is on active block
 * @param[out] integer value of boundary tag
 * @param[out] integer value for boundary dimension
 */
void where_is_node(MBBL pumi_obj, int knode_x1, int knode_x2, bool* on_bdry, bool* in_domain, int* bdry_tag, int* bdry_dim){

    int isubmesh, jsubmesh, inp, jnp;
    bool left_edge, right_edge, bottom_edge, top_edge;

    for (isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
        int submesh_min_node = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
        int submesh_max_node = pumi_obj.host_submesh_x1[isubmesh]->Nel + submesh_min_node;
        left_edge =  false;
        right_edge = false;
        if (knode_x1 >= submesh_min_node && knode_x1 <= submesh_max_node){
            inp = knode_x1 - submesh_min_node;
            if (inp == 0){
                left_edge = true;
            }
            if (inp == pumi_obj.host_submesh_x1[isubmesh]->Nel){
                right_edge = true;
            }
            break;
        }
    }

    for (jsubmesh=1; jsubmesh<=pumi_obj.mesh.nsubmesh_x2; jsubmesh++){
        int submesh_min_node = pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
        int submesh_max_node = pumi_obj.host_submesh_x2[jsubmesh]->Nel + submesh_min_node;
        bottom_edge =  false;
        top_edge = false;
        if (knode_x2 >= submesh_min_node && knode_x2 <= submesh_max_node){
            jnp = knode_x2 - submesh_min_node;
            if (jnp == 0){
                bottom_edge = true;
            }
            if (jnp == pumi_obj.host_submesh_x2[jsubmesh]->Nel){
                top_edge = true;
            }
            break;
        }
    }

    if (!left_edge && !right_edge && !bottom_edge && !top_edge){
        *on_bdry = false;

        if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
            *in_domain = true;
            *bdry_tag = -1;
            *bdry_dim = -1;
            return;
        }
        else{
            *in_domain = false;
            *bdry_tag = -999;
            *bdry_dim = -1;
            return;
        }
    }

    if (left_edge & !top_edge & !bottom_edge){

        if (isubmesh==1){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = (jsubmesh-1)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + pumi_obj.mesh.nsubmesh_x1;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                *in_domain = true;
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] + pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh] == 1){
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + pumi_obj.mesh.nsubmesh_x1;
                    *bdry_dim = 1;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (left_edge & top_edge){
        if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
            if (isubmesh==1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh+1] |
                    pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh] + pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh+1] +
                        pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1] + pumi_obj.mesh.host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }

    if (top_edge & !left_edge & !right_edge){

        if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = (jsubmesh)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                *in_domain = true;
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] + pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1] == 1){
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 1;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (top_edge & right_edge){
        if (jsubmesh==pumi_obj.mesh.nsubmesh_x2){
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1] | pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh+1] |
                    pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = pumi_obj.mesh.host_isactive[isubmesh][jsubmesh+1] + pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh+1] +
                            pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh] + pumi_obj.mesh.host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = (jsubmesh)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }

    if (right_edge & !top_edge & !bottom_edge){

        if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = (jsubmesh-1)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + pumi_obj.mesh.nsubmesh_x1 + 1;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                *in_domain = true;
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] + pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh] == 1){
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + pumi_obj.mesh.nsubmesh_x1 + 1;
                    *bdry_dim = 1;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (right_edge & bottom_edge){
        if (jsubmesh==1){
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1] | pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh] |
                    pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh-1] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1] + pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh] +
                            pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh-1] + pumi_obj.mesh.host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1 + 1;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }

    if (bottom_edge & !left_edge & !right_edge){

        if (jsubmesh==1){
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = (jsubmesh-1)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                *in_domain = true;
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] + pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1] == 1){
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(2*pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (bottom_edge & left_edge){
        if (jsubmesh==1){
            if (isubmesh==1){
                if (pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==1){
                if(pumi_obj.mesh.host_isactive[isubmesh][jsubmesh] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh-1] | pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh] |
                    pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1] | pumi_obj.mesh.host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh-1] + pumi_obj.mesh.host_isactive[isubmesh-1][jsubmesh] +
                            pumi_obj.mesh.host_isactive[isubmesh][jsubmesh-1] + pumi_obj.mesh.host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = (jsubmesh-1)*(pumi_obj.mesh.nsubmesh_x1 + 1) + isubmesh-1;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }
}

} // namespace pumi
