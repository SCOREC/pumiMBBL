#include "pumiMBBL_meshinfo.hpp"

namespace pumi {

/**
 * @brief Returns the grading ratio along a queried direction about a queried node
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] direction along which grading ratio is needed
 * \param[in] global directional node ID
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
                else{
                    min_elem = h_submesh[isubmesh-1]->t0;
                }

                if (h_submesh[isubmesh]->meshtype & maxBL){
                    max_elem = h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r, h_submesh[isubmesh]->Nel-1);
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
 * \param[in] Object of the wrapper mesh structure
 * \param[in] direction along which element-size is needed
 * \param[in] global directional node ID or element ID index
 * \param[in] relevant offset enum based on provied index
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
        elem = Nel_total;
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
            if (h_submesh[isubmesh]->meshtype & minBL){
                int local_cell = elem - submesh_min_elem;
                return (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
            }
            if (h_submesh[isubmesh]->meshtype & maxBL){
                int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                return (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
            }
        }
    }
    return -999.0;
}

/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] global node IDs along x1-direction
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
 * \param[in] Object of the wrapper mesh structure
 * \param[in] global node IDs along x1-direction
 * \param[in] global node IDs along x2-direction
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
 * \param[in] Object of the wrapper mesh structure
 * \param[in] global node IDs along x1-direction
 * \param[in] global node IDs along x2-direction
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
                if (h_submesh[isubmesh]->meshtype & minBL){
                    int local_cell = elem - submesh_min_elem;
                    dx1_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
                if (h_submesh[isubmesh]->meshtype & maxBL){
                    int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                    dx1_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
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
                if (h_submesh[isubmesh]->meshtype & minBL){
                    int local_cell = elem - submesh_min_elem;
                    dx1_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
                if (h_submesh[isubmesh]->meshtype & maxBL){
                    int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                    dx1_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
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
                if (h_submesh[isubmesh]->meshtype & minBL){
                    int local_cell = elem - submesh_min_elem;
                    dx2_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
                if (h_submesh[isubmesh]->meshtype & maxBL){
                    int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                    dx2_min = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
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
                if (h_submesh[isubmesh]->meshtype & minBL){
                    int local_cell = elem - submesh_min_elem;
                    dx2_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
                if (h_submesh[isubmesh]->meshtype & maxBL){
                    int local_cell = h_submesh[isubmesh]->Nel - (elem - submesh_min_elem) - 1;
                    dx2_max = (h_submesh[isubmesh]->t0 * pow(h_submesh[isubmesh]->r , local_cell));
                }
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
 * @brief Returns x1-coordinate of left end of the domain
 *
 * \param[in] Object of the wrapper mesh structure
 */
double get_global_x1_min_coord(MBBL pumi_obj){
    return pumi_obj.host_submesh_x1[1]->xmin;
}

/**
 * @brief Returns x1-coordinate of right end of the domain
 *
 * \param[in] Object of the wrapper mesh structure
 */
double get_global_x1_max_coord(MBBL pumi_obj){
    int nsubmesh = pumi_obj.mesh.nsubmesh_x1;
    return pumi_obj.host_submesh_x1[nsubmesh]->xmax;
}

/**
 * @brief Returns x2-coordinate of bottom end of the domain
 *
 * \param[in] Object of the wrapper mesh structure
 */
double get_global_x2_min_coord(MBBL pumi_obj){
    return pumi_obj.host_submesh_x2[1]->xmin;
}

/**
 * @brief Returns x2-coordinate of top end of the domain
 *
 * \param[in] Object of the wrapper mesh structure
 */
double get_global_x2_max_coord(MBBL pumi_obj){
    int nsubmesh = pumi_obj.mesh.nsubmesh_x2;
    return pumi_obj.host_submesh_x2[nsubmesh]->xmax;
}

int get_num_x1_submesh(MBBL pumi_obj){
    return pumi_obj.mesh.nsubmesh_x1;
}

int get_num_x1_elems_in_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x1[isubmesh]->Nel;
}

int get_num_x1_elems_before_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
}

int get_num_x2_submesh(MBBL pumi_obj){
    return pumi_obj.mesh.nsubmesh_x2;
}

int get_num_x2_elems_in_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x2[isubmesh]->Nel;
}

int get_num_x2_elems_before_submesh(MBBL pumi_obj, int isubmesh){
    return pumi_obj.host_submesh_x2[isubmesh]->Nel_cumulative;
}

int get_total_mesh_elements(MBBL pumi_obj){
    return pumi_obj.mesh.Nel_total;
}

int get_total_mesh_nodes(MBBL pumi_obj){
    return pumi_obj.mesh.Nnp_total;
}

int get_total_x1_elements(MBBL pumi_obj){
    return pumi_obj.mesh.Nel_tot_x1;
}

int get_total_x2_elements(MBBL pumi_obj){
    return pumi_obj.mesh.Nel_tot_x2;
}

int get_total_submesh_blocks(MBBL pumi_obj){
    if (pumi_obj.mesh.ndim==1){
        return pumi_obj.mesh.nsubmesh_x1;
    }
    else if (pumi_obj.mesh.ndim==2){
        return pumi_obj.mesh.nsubmesh_x1*pumi_obj.mesh.nsubmesh_x2;
    }
    return 0;
}

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

std::vector<double> get_bdry_normal(MBBL pumi_obj, int iEdge){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    if (iEdge>=0 && iEdge<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2){
        // double nrml = pumi_obj.mesh.bdry.host_bdry_edge_normal[iEdge][dir];
        std::vector<double> nrml = {pumi_obj.mesh.bdry.host_bdry_edge_normal[iEdge][0],
                                    pumi_obj.mesh.bdry.host_bdry_edge_normal[iEdge][1],
                                    pumi_obj.mesh.bdry.host_bdry_edge_normal[iEdge][2]};
        return nrml;
    }
    else{
        std::cout << "Invalid edge ID or direction\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        std::cout << "Valid directions 0 (x1-dirxn) or 1 (x2-dirxn) or 2 (x3-dirxn)\n ";
        exit(0);
    }
}

int get_num_faces_on_bdry(MBBL pumi_obj, int iEdge){
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

int get_starting_faceID_on_bdry(MBBL pumi_obj, int iEdge){
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

bool is_edge_bdry(MBBL pumi_obj, int iEdge){
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

bool is_block_active(MBBL pumi_obj, int isub, int jsub){
    return pumi_obj.mesh.host_isactive[isub][jsub];
}

bool is_block_active(MBBL pumi_obj, int flattened_submesh_ID){
    int jsub = flattened_submesh_ID/pumi_obj.mesh.nsubmesh_x1 + 1;
    int isub = flattened_submesh_ID - (jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
    return pumi_obj.mesh.host_isactive[isub][jsub];
}

int get_total_mesh_block_edges(MBBL pumi_obj){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    return 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2;
}

int get_global_nodeID(MBBL pumi_obj, int submeshID, int fullmesh_node_id){
    int isubmesh, jsubmesh;
    jsubmesh = submeshID/pumi_obj.mesh.nsubmesh_x1;
    isubmesh = submeshID - jsubmesh*pumi_obj.mesh.nsubmesh_x1;

    isubmesh++;
    jsubmesh++;

    int Jnp;
    Jnp = fullmesh_node_id/(pumi_obj.mesh.Nel_tot_x1+1);


    int nodeID = fullmesh_node_id;
    int jnp = Jnp - pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
    // int nodeoffset = pumi_obj.mesh(0).nodeoffset(isubmesh,Jnp);
    int nodeoffset;
    nodeoffset = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh] + pumi_obj.mesh.offsets.host_nodeoffset_skip_bot[isubmesh][jsubmesh]
                    +(jnp-1)*pumi_obj.mesh.offsets.host_nodeoffset_skip_mid[isubmesh][jsubmesh];
    if (jnp==0){
        nodeoffset = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh];
    }
    if (jnp==pumi_obj.host_submesh_x2[jsubmesh]->Nel){
        nodeoffset +=  (pumi_obj.mesh.offsets.host_nodeoffset_skip_top[isubmesh][jsubmesh]-pumi_obj.mesh.offsets.host_nodeoffset_skip_mid[isubmesh][jsubmesh]);
    }
    return nodeID-nodeoffset;
}

void get_edge_info(MBBL pumi_obj, int iEdge, int *Knp, int *next_offset, int *submeshID){
    int nsubmesh_x1 = pumi_obj.mesh.nsubmesh_x1;
    int nsubmesh_x2 = pumi_obj.mesh.nsubmesh_x2;
    int Nx2p1 = 2*nsubmesh_x1+1;
    if (iEdge>=0 && iEdge<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2){
        if (is_edge_bdry(pumi_obj,iEdge)){
            int num = iEdge/Nx2p1;
            int rem = iEdge-num*Nx2p1;
            int isubmesh, jsubmesh;
            if (rem < nsubmesh_x1){
                *next_offset = 1;
                isubmesh = rem;
                jsubmesh = num;
                if (jsubmesh >= nsubmesh_x2){
                    jsubmesh = nsubmesh_x2-1;
                    int Inp = pumi_obj.host_submesh_x1[isubmesh+1]->Nel_cumulative;
                    int Jnp = pumi_obj.mesh.Nel_tot_x2;
                    *Knp = (pumi_obj.mesh.Nel_tot_x1+1)*Jnp + Inp;
                    *submeshID = jsubmesh*nsubmesh_x1+isubmesh;
                    return;
                }
                else{
                    int Inp = pumi_obj.host_submesh_x1[isubmesh+1]->Nel_cumulative;
                    int Jnp = pumi_obj.host_submesh_x2[jsubmesh+1]->Nel_cumulative;
                    if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh+1]){
                        *Knp = (pumi_obj.mesh.Nel_tot_x1+1)*Jnp + Inp;
                        *submeshID = jsubmesh*nsubmesh_x1+isubmesh;
                        return;
                    }
                    else{
                        *Knp = (pumi_obj.mesh.Nel_tot_x1+1)*Jnp + Inp;
                        jsubmesh--;
                        *submeshID = jsubmesh*nsubmesh_x1+isubmesh;
                        return;
                    }
                }
            }
            else {
                *next_offset = pumi_obj.mesh.Nel_tot_x1+1;
                jsubmesh = num;
                isubmesh = rem-nsubmesh_x2;
                if (isubmesh >= nsubmesh_x1){
                    isubmesh = nsubmesh_x1-1;
                    int Jnp = pumi_obj.host_submesh_x2[jsubmesh+1]->Nel_cumulative;
                    int Inp = pumi_obj.mesh.Nel_tot_x1;
                    *Knp = (pumi_obj.mesh.Nel_tot_x1+1)*Jnp + Inp;
                    *submeshID = jsubmesh*nsubmesh_x1+isubmesh;
                    return;
                }
                else{
                    int Inp = pumi_obj.host_submesh_x1[isubmesh+1]->Nel_cumulative;
                    int Jnp = pumi_obj.host_submesh_x2[jsubmesh+1]->Nel_cumulative;
                    if (pumi_obj.mesh.host_isactive[isubmesh+1][jsubmesh+1]){
                        *Knp = (pumi_obj.mesh.Nel_tot_x1+1)*Jnp + Inp;
                        *submeshID = jsubmesh*nsubmesh_x1+isubmesh;
                        return;
                    }
                    else{
                        *Knp = (pumi_obj.mesh.Nel_tot_x1+1)*Jnp + Inp;
                        isubmesh--;
                        *submeshID = jsubmesh*nsubmesh_x1+isubmesh;
                        return;
                    }
                }
            }
        }
        else{
            *Knp=-1;
            *next_offset=-1;
            *submeshID=-1;
            return;
        }

    }
    else{
        std::cout << "Invalid edge ID\n";
        std::cout << "Valid EdgeIDs = [0,1,..," << 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2-1 <<"]\n";
        exit(0);
    }
}

/**
 * @brief Returns node info such as if node is in active domain, if node is on a boundary
 * and boundary entity dimension (boundary vertex (dim=0) or edge (dim=1)) and entity tag
 * of the boundary
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] global node IDs along x1-direction
 * \param[in] global node IDs along x2-direction
 * \param[out] boolean value if node is on boundary
 * \param[out] boolean value if node is on active block
 * \param[out] integer value of boundary tag
 * \param[out] integer value for boundary dimension
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