#include "pumiMBBL_meshops.hpp"

namespace pumi {

Vector3 get_rand_point_in_mesh_host(MBBL pumi_obj){
    if (pumi_obj.mesh.ndim == 1){
        double rand_x1 = (double) rand()/RAND_MAX;
        double x1_min = get_global_x1_min_coord(pumi_obj);
        double x1_max = get_global_x1_max_coord(pumi_obj);
        double q1 = x1_min + (x1_max-x1_min)*rand_x1;
        // std::vector<double> q = {q1,0.0,0.0};
        Vector3 q = Vector3(q1,0.0,0.0);
        return q;
    }
    else if (pumi_obj.mesh.ndim == 2){
        bool part_set = false;
        double x1_min = get_global_x1_min_coord(pumi_obj);
        double x1_max = get_global_x1_max_coord(pumi_obj);
        double x2_min = get_global_x2_min_coord(pumi_obj);
        double x2_max = get_global_x2_max_coord(pumi_obj);
        double q1, q2;
        while (!part_set){
            double rand_x1 = (double) rand()/RAND_MAX;
            double rand_x2 = (double) rand()/RAND_MAX;
            q1 = x1_min + (x1_max-x1_min)*rand_x1;
            q2 = x2_min + (x2_max-x2_min)*rand_x2;

            int isub=0;
            int jsub=0;

            for (int i=1; i<=pumi_obj.mesh.nsubmesh_x1; i++){
                if (pumi_obj.host_submesh_x1[i]->xmin < q1 && pumi_obj.host_submesh_x1[i]->xmax > q1){
                    isub = i;
                    break;
                }
            }
            for (int j=1; j<=pumi_obj.mesh.nsubmesh_x2; j++){
                if (pumi_obj.host_submesh_x2[j]->xmin < q2 && pumi_obj.host_submesh_x2[j]->xmax > q2){
                    jsub = j;
                    break;
                }
            }
            if (pumi_obj.mesh.host_isactive[isub][jsub]){
                part_set = true;
            }
        }
        // std::vector<double> q = {q1,q2,0.0};
        Vector3 q = Vector3(q1,q2,0.0);
        return q;
    }
    // std::vector<double> q = {-999.0,-999.0,-999.0};
    Vector3 q = Vector3(-999.0,-999.0,-999.0);
    return q;
}

bool is_point_in_mesh_host(MBBL pumi_obj, Vector3 q){
    if (pumi_obj.mesh.ndim == 1){
        double x1_min = get_global_x1_min_coord(pumi_obj);
        double x1_max = get_global_x1_max_coord(pumi_obj);
        if (q[0]>x1_min && q[0]<x1_max){
            return true;
        }
        return false;
    }
    else if (pumi_obj.mesh.ndim == 2){
        int isub=0;
        int jsub=0;

        for (int i=1; i<=pumi_obj.mesh.nsubmesh_x1; i++){
            if (pumi_obj.host_submesh_x1[i]->xmin < q[0] && pumi_obj.host_submesh_x1[i]->xmax > q[0]){
                isub = i;
                break;
            }
        }
        for (int j=1; j<=pumi_obj.mesh.nsubmesh_x2; j++){
            if (pumi_obj.host_submesh_x2[j]->xmin < q[1] && pumi_obj.host_submesh_x2[j]->xmax > q[1]){
                jsub = j;
                break;
            }
        }
        if (pumi_obj.mesh.host_isactive[isub][jsub]){
            return true;
        }
        return false;
    }
    return false;
}

void flatten_submeshID_and_cellID_host(MBBL pumi_obj, int isub, int icell, int jsub, int jcell, int* submeshID, int* cellID){
    *submeshID = (isub-1) + (jsub-1)*pumi_obj.mesh.nsubmesh_x1;
    *cellID = icell + jcell*pumi_obj.host_submesh_x1[isub]->Nel;
}

/**
* @brief Locate the submesh ID and local cell ID for a given x1-coordinate
* Uses analytical formulae to locate the input coordinate
* \param[in] Object of the wrapper mesh structure
* \param[in] x1-coordinate to be located
* \param[out] located x1-submesh ID
* \param[out] located x1-localcell ID
*/
void locate_submesh_and_cell_x1_host(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    // int nsubmesh = pumi_obj.submesh_x1.extent(0);
    int nsubmesh = pumi_obj.mesh.nsubmesh_x1;
    for (isubmesh=2; isubmesh<=nsubmesh; isubmesh++){
     if (q < (pumi_obj.host_submesh_x1[isubmesh]->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh;
    }
    *cellID  = pumi_obj.host_submesh_x1[*submeshID]->locate_cell_host(q);
}

/**
* @brief Locate the submesh ID and local cell ID for a given x2-coordinate
* Uses analytical formulae to locate the input coordinate
* \param[in] Object of the wrapper mesh structure
* \param[in] x2-coordinate to be located
* \param[out] located x2-submesh ID
* \param[out] located x2-localcell ID
*/
void locate_submesh_and_cell_x2_host(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    // int nsubmesh = pumi_obj.submesh_x2.extent(0);
    int nsubmesh = pumi_obj.mesh.nsubmesh_x2;
    for (isubmesh=2; isubmesh<=nsubmesh; isubmesh++){
     if (q < (pumi_obj.host_submesh_x2[isubmesh]->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh;
    }
    *cellID  = pumi_obj.host_submesh_x2[*submeshID]->locate_cell_host(q);
}

/**
 * @brief Update the submesh ID and local cell ID for a given x1-coordinate
 * based on previous submesh and cell IDs.
 * Uses adjacency search to update the IDs
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x1-coordinate
 * \param[in] old x1-submesh ID
 * \param[in] old x1-localcell ID
 * \param[out] updated x1-submesh ID
 * \param[out] updated x1-localcell ID
 */
void update_submesh_and_cell_x1_host(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.host_submesh_x1[*submeshID]->xmin)){
        *submeshID -= 1;
        prev_cellID = pumi_obj.host_submesh_x1[*submeshID]->Nel - 1;
    }

    while(q > (pumi_obj.host_submesh_x1[*submeshID]->xmax)){
        *submeshID += 1;
        prev_cellID = 0;
    }
    *cellID = pumi_obj.host_submesh_x1[*submeshID]->update_cell_host(q, prev_cellID);
}


/**
 * @brief Update the submesh ID and local cell ID for a given x2-coordinate
 * based on previous submesh and cell IDs.
 * Uses adjacency search to update the IDs
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x2-coordinate
 * \param[in] old x2-submesh ID
 * \param[in] old x2-localcell ID
 * \param[out] updated x2-submesh ID
 * \param[out] updated x2-localcell ID
 */
void update_submesh_and_cell_x2_host(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.host_submesh_x2[*submeshID]->xmin)){
        *submeshID -= 1;
        prev_cellID = pumi_obj.host_submesh_x2[*submeshID]->Nel - 1;
    }

    while(q > (pumi_obj.host_submesh_x2[*submeshID]->xmax)){
        *submeshID += 1;
        prev_cellID = 0;
    }
    *cellID = pumi_obj.host_submesh_x2[*submeshID]->update_cell_host(q, prev_cellID);
}

/**
 * @brief Computes the partial weights (correspoding to node on the max-side i.e right side)
 * for a located particle coordinate and the global directional cell ID
 * \param[in] Object of the wrapper mesh structure
 * \param[in] x1-coordinate of the particle
 * \param[in] x1-submesh ID of the particle
 * \param[in] x1-localcell ID of the particle
 * \param[out] global cell ID in x1 direction
 * \param[out] partial weight
 */
void calc_weights_x1_host(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2){
    pumi_obj.host_submesh_x1[isubmesh]->calc_weights_host(q, icell, x1_global_cell, Wgh2);
}

/**
 * @brief Computes the partial weights (correspoding to node on the max-side i.e top side)
 * for a located particle coordinate and the global directional cell ID
 * \param[in] Object of the wrapper mesh structure
 * \param[in] x2-coordinate of the particle
 * \param[in] x2-submesh ID of the particle
 * \param[in] x2-localcell ID of the particle
 * \param[out] global cell ID in x2 direction
 * \param[out] partial weight
 */
 void calc_weights_x2_host(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2){
     pumi_obj.host_submesh_x2[isubmesh]->calc_weights_host(q, icell, x2_global_cell, Wgh2);
 }

/**
* @brief Computes the gloabl cell ID and node ID in 2D for a full Mesh
* with no-inactive blocks (mesh with inactive blocks will need separate implementations)
* \param[in] global cell ID in x1-direction
* \param[in] global cell ID in x2-direction
* \param[out] global cell ID in 2D
* \param[out] global node ID of the node in left-bottom corner
* \param[out] global node ID of the node in left-top coner
*/
void calc_global_cellID_and_nodeID_fullmesh_host(MBBL pumi_obj, int kcell_x1, int kcell_x2, int *global_cell_2D, int *bottomleft_node, int *topleft_node){
  *global_cell_2D = kcell_x1 + kcell_x2*pumi_obj.mesh.Nel_tot_x1;
  *bottomleft_node = *global_cell_2D + kcell_x2;
  *topleft_node = *bottomleft_node + pumi_obj.mesh.Nel_tot_x1 + 1;
}

/**
* @brief Computes the gloabl cell ID and node ID in 2D for a full Mesh
* with no-inactive blocks (mesh with inactive blocks will need separate implementations)
* \param[in] global cell ID in x1-direction
* \param[in] global cell ID in x2-direction
* \param[out] global cell ID in 2D
* \param[out] global node ID of the node in left-bottom corner
* \param[out] global node ID of the node in left-top coner
*/
void calc_global_cellID_and_nodeID_host(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node){
    int icell_x2 = kcell_x2 - pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
    int elemoffset = pumi_obj.mesh.offsets.host_elemoffset_start[isubmesh][jsubmesh] + icell_x2*pumi_obj.mesh.offsets.host_elemoffset_skip[jsubmesh];
    int fullmesh_elem = kcell_x1 + kcell_x2*pumi_obj.mesh.Nel_tot_x1;
    *global_cell_2D = fullmesh_elem - elemoffset;
    int nodeoffset_bottom = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh] + pumi_obj.mesh.offsets.host_nodeoffset_skip_bot[isubmesh][jsubmesh]
                    +(icell_x2-1)*pumi_obj.mesh.offsets.host_nodeoffset_skip_mid[isubmesh][jsubmesh];
    int nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh.offsets.host_nodeoffset_skip_mid[isubmesh][jsubmesh];
    if (icell_x2==0){
        nodeoffset_bottom = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh];
        nodeoffset_top = pumi_obj.mesh.offsets.host_nodeoffset_start[isubmesh][jsubmesh] + pumi_obj.mesh.offsets.host_nodeoffset_skip_bot[isubmesh][jsubmesh];
    }
    if (icell_x2==pumi_obj.host_submesh_x2[jsubmesh]->Nel-1){
        nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh.offsets.host_nodeoffset_skip_top[isubmesh][jsubmesh];
    }
    *bottomleft_node = fullmesh_elem + kcell_x2 - nodeoffset_bottom;
    *topleft_node = fullmesh_elem + kcell_x2 + pumi_obj.mesh.Nel_tot_x1 + 1 - nodeoffset_top;
}

void get_directional_submeshID_and_cellID_host(MBBL pumi_obj, int submeshID, int cellID, int* isub, int *icell, int* jsub, int *jcell){
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - pumi_obj.mesh.nsubmesh_x1*(*jsub-1) + 1;
    *jcell = cellID/pumi_obj.host_submesh_x1[*isub]->Nel;
    *icell = cellID - pumi_obj.host_submesh_x1[*isub]->Nel*(*jcell);
}

void get_directional_submeshID_host(MBBL pumi_obj, int submeshID, int* isub, int* jsub){
    *jsub = submeshID/pumi_obj.mesh.nsubmesh_x1 + 1;
    *isub = submeshID - pumi_obj.mesh.nsubmesh_x1*(*jsub-1) + 1;
}

} // namespace pumi
