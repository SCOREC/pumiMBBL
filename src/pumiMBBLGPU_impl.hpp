#ifndef pumiMBBLGPU_impl_hpp
#define pumiMBBLGPU_impl_hpp

#include "pumiMBBLGPU.hpp"

namespace pumi {

KOKKOS_INLINE_FUNCTION
int locate_cell(DevicePointer<Submesh> submesh, double q) {
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->locate_cell(q);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->locate_cell(q);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->locate_cell(q);
        case (unassigned) :
            return -1;
    }
    return -1;
}

KOKKOS_INLINE_FUNCTION
int update_cell(DevicePointer<Submesh> submesh, double q, int icell) {
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->update_cell(q,icell);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->update_cell(q,icell);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->update_cell(q,icell);
        case (unassigned) :
            return -1;
    }
    return -1;
}

KOKKOS_INLINE_FUNCTION
double elem_size(DevicePointer<Submesh> submesh, int icell){
    switch (submesh()->meshtype){
        case (uniform) :
            return static_cast<Uniform_Submesh*>(submesh())->elem_size(icell);
        case (minBL) :
            return static_cast<MinBL_Submesh*>(submesh())->elem_size(icell);
        case (maxBL) :
            return static_cast<MaxBL_Submesh*>(submesh())->elem_size(icell);
        case (unassigned) :
            return -999.0;
    }
    return -999.0;
}

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
        case (unassigned) :
            *global_cell = -1;
            *Wgh2 = -999.0;
            return;
    }
    *global_cell = -1;
    *Wgh2 = -999.0;
    return;
}

KOKKOS_INLINE_FUNCTION
double get_x1_elem_size_in_submesh(MBBL pumi_obj, int isub, int icell){
    return elem_size(pumi_obj.submesh_x1(isub),icell);
}

KOKKOS_INLINE_FUNCTION
double get_x2_elem_size_in_submesh(MBBL pumi_obj, int isub, int icell){
    return elem_size(pumi_obj.submesh_x2(isub),icell);
}

KOKKOS_INLINE_FUNCTION
int get_x1_cellID(MBBL pumi_obj, int isub, int icell){
    return icell + pumi_obj.submesh_x1(isub)()->Nel_cumulative;
}

KOKKOS_INLINE_FUNCTION
int get_x2_cellID(MBBL pumi_obj, int isub, int icell){
    return icell + pumi_obj.submesh_x2(isub)()->Nel_cumulative;
}

KOKKOS_INLINE_FUNCTION
void get_directional_submeshID_and_cellID(MBBL pumi_obj, int submeshID, int cellID, int* isub, int *icell, int* jsub, int *jcell){
    *jsub = submeshID/pumi_obj.mesh(0).nsubmesh_x1 + 1;
    *isub = submeshID - pumi_obj.mesh(0).nsubmesh_x1*(*jsub-1) + 1;
    *jcell = cellID/pumi_obj.submesh_x1(*isub)()->Nel;
    *icell = cellID - pumi_obj.submesh_x1(*isub)()->Nel*(*jcell);
}

KOKKOS_INLINE_FUNCTION
void flatten_submeshID_and_cellID(MBBL pumi_obj, int isub, int icell, int jsub, int jcell, int* submeshID, int* cellID){
    *submeshID = (isub-1) + (jsub-1)*pumi_obj.mesh(0).nsubmesh_x1;
    *cellID = icell + jcell*pumi_obj.submesh_x1(isub)()->Nel;
}

/**
* @brief Locate the submesh ID and local cell ID for a given x1-coordinate
* Uses analytical formulae to locate the input coordinate
* \param[in] Object of the wrapper mesh structure
* \param[in] x1-coordinate to be located
* \param[out] located x1-submesh ID
* \param[out] located x1-localcell ID
*/
KOKKOS_INLINE_FUNCTION
void locate_submesh_and_cell_x1(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    // int nsubmesh = pumi_obj.submesh_x1.extent(0);
    int nsubmesh = pumi_obj.mesh(0).nsubmesh_x1;
    for (isubmesh=2; isubmesh<=nsubmesh; isubmesh++){
     if (q < (pumi_obj.submesh_x1(isubmesh)()->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh;
    }
    // *cellID  = pumi_obj.submesh_x1(*submeshID)()->locate_cell(q);
    *cellID = locate_cell(pumi_obj.submesh_x1(*submeshID),q);
}

/**
* @brief Locate the submesh ID and local cell ID for a given x2-coordinate
* Uses analytical formulae to locate the input coordinate
* \param[in] Object of the wrapper mesh structure
* \param[in] x2-coordinate to be located
* \param[out] located x2-submesh ID
* \param[out] located x2-localcell ID
*/
KOKKOS_INLINE_FUNCTION
void locate_submesh_and_cell_x2(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    // int nsubmesh = pumi_obj.submesh_x2.extent(0);
    int nsubmesh = pumi_obj.mesh(0).nsubmesh_x2;
    for (isubmesh=2; isubmesh<=nsubmesh; isubmesh++){
     if (q < (pumi_obj.submesh_x2(isubmesh)()->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh;
    }
    // *cellID  = pumi_obj.submesh_x2(*submeshID)()->locate_cell(q);
    *cellID  = locate_cell(pumi_obj.submesh_x2(*submeshID),q);
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
    // *cellID = pumi_obj.submesh_x1(*submeshID)()->update_cell(q, prev_cellID);
    *cellID = update_cell(pumi_obj.submesh_x1(*submeshID), q, prev_cellID);
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
    // *cellID = pumi_obj.submesh_x2(*submeshID)()->update_cell(q, prev_cellID);
    *cellID = update_cell(pumi_obj.submesh_x2(*submeshID), q, prev_cellID);
}

/**
 * @brief Update the submesh ID based on previous submesh and
 * locate the local cell ID
 * Uses adjacency search to update the submesh IDs and analytical
 * forumlae for cell IDs. Use this function when BL coords are
 * not stored
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x1-coordinate
 * \param[in] old x1-submesh ID
 * \param[out] updated x1-submesh ID
 * \param[out] updated x1-localcell ID
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

    // *cellID = pumi_obj.submesh_x1(*submeshID)()->locate_cell(q);
    *cellID = locate_cell(pumi_obj.submesh_x1(*submeshID),q);
}


/**
 * @brief Update the submesh ID based on previous submesh and
 * locate the local cell ID
 * Uses adjacency search to update the submesh IDs and analytical
 * forumlae for cell IDs. Use this function when BL coords are
 * not stored
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x2-coordinate
 * \param[in] old x2-submesh ID
 * \param[out] updated x2-submesh ID
 * \param[out] updated x2-localcell ID
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

    // *cellID = pumi_obj.submesh_x2(*submeshID)()->locate_cell(q);
    *cellID = locate_cell(pumi_obj.submesh_x2(*submeshID),q);
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
KOKKOS_INLINE_FUNCTION
void calc_weights_x1(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2){
    // pumi_obj.submesh_x1(isubmesh)()->calc_weights(q, icell, x1_global_cell, Wgh2);
    calc_weights(pumi_obj.submesh_x1(isubmesh),q,icell,x1_global_cell,Wgh2);
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
 KOKKOS_INLINE_FUNCTION
 void calc_weights_x2(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2){
     // pumi_obj.submesh_x2(isubmesh)()->calc_weights(q, icell, x2_global_cell, Wgh2);
     calc_weights(pumi_obj.submesh_x2(isubmesh),q,icell,x2_global_cell,Wgh2);
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
  KOKKOS_INLINE_FUNCTION
  void calc_global_cellID_and_nodeID_fullmesh(MBBL pumi_obj, int kcell_x1, int kcell_x2, int *global_cell_2D, int *bottomleft_node, int *topleft_node){
      *global_cell_2D = kcell_x1 + kcell_x2*pumi_obj.mesh(0).Nel_tot_x1;
      *bottomleft_node = *global_cell_2D + kcell_x2;
      *topleft_node = *bottomleft_node + pumi_obj.mesh(0).Nel_tot_x1 + 1;
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
KOKKOS_INLINE_FUNCTION
void calc_global_cellID_and_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node){
    // int nodeoffset_bottom = pumi_obj.mesh(0).nodeoffset(isubmesh,kcell_x2);
    // int nodeoffset_top = pumi_obj.mesh(0).nodeoffset(isubmesh,kcell_x2+1);
    int icell_x2 = kcell_x2 - pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    int elemoffset = pumi_obj.mesh(0).elemoffset_start(isubmesh,jsubmesh) + icell_x2*pumi_obj.mesh(0).elemoffset_skip(jsubmesh);
    int fullmesh_elem = kcell_x1 + kcell_x2*pumi_obj.mesh(0).Nel_tot_x1;
    *global_cell_2D = fullmesh_elem - elemoffset;
    int nodeoffset_bottom = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh(0).nodeoffset_skip_bot(isubmesh,jsubmesh)
                    +(icell_x2-1)*pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh);
    int nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh);
    if (icell_x2==0){
        nodeoffset_bottom = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh);
        nodeoffset_top = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh(0).nodeoffset_skip_bot(isubmesh,jsubmesh);
    }
    if (icell_x2==pumi_obj.submesh_x2(jsubmesh)()->Nel-1){
        nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh(0).nodeoffset_skip_top(isubmesh,jsubmesh);
    }
    *bottomleft_node = fullmesh_elem + kcell_x2 - nodeoffset_bottom;
    *topleft_node = fullmesh_elem + kcell_x2 + pumi_obj.mesh(0).Nel_tot_x1 + 1 - nodeoffset_top;
}


KOKKOS_INLINE_FUNCTION
int calc_global_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int Inp, int Jnp){
    int nodeID = Jnp*(pumi_obj.mesh(0).Nel_tot_x1+1) + Inp;
    int jnp = Jnp - pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    // int nodeoffset = pumi_obj.mesh(0).nodeoffset(isubmesh,Jnp);
    int nodeoffset;
    nodeoffset = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh(0).nodeoffset_skip_bot(isubmesh,jsubmesh)
                    +(jnp-1)*pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh);
    if (jnp==0){
        nodeoffset = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh);
    }
    if (jnp==pumi_obj.submesh_x2(jsubmesh)()->Nel){
        nodeoffset +=  (pumi_obj.mesh(0).nodeoffset_skip_top(isubmesh,jsubmesh)-pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh));
    }
    return nodeID-nodeoffset;
}

KOKKOS_INLINE_FUNCTION
void push_particle(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit){

    double q1_new = q1+dq1;
    double q2_new = q2+dq2;
    int Nx = pumi_obj.mesh(0).nsubmesh_x1;
    int Nxx = 2*Nx+1;
    // int Ny = pumi_obj.mesh(0).nsubmesh_x2;
    int case_id = (dq2>=0.0)+2*(dq1>=0.0);
    int isub = *isubmesh;
    int jsub = *jsubmesh;
    *in_domain = true;
    if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
        && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){

        // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
        *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
        // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
        *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
        *bdry_hit = -1;
        return;
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
                }
                else{
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                    jsub--;
                }


                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                            del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                                isub--;
                            }
                            else{
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                                jsub--;
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
                }
                else{
                    *bdry_hit = jsub*(Nxx)+isub-1;
                    jsub++;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                            del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                                isub--;
                            }
                            else{
                                *bdry_hit = jsub*(Nxx)+isub-1;
                                jsub++;
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
                }
                else{
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                    jsub--;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                            del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                                isub++;
                            }
                            else{
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                                jsub--;
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
                }
                else{
                    *bdry_hit = jsub*(Nxx)+isub-1;
                    jsub++;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            // *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *icell = update_cell(pumi_obj.submesh_x1(isub),q1_new,*icell);
                            // *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *jcell = update_cell(pumi_obj.submesh_x2(jsub),q2_new,*jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                            del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                                isub++;
                            }
                            else{
                                *bdry_hit = jsub*(Nxx)+isub-1;
                                jsub++;
                            }
                        }
                    }
                }
        }

    }

}

KOKKOS_INLINE_FUNCTION
void push_particle_v2(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit){

    double q1_new = q1+dq1;
    double q2_new = q2+dq2;
    int Nx = pumi_obj.mesh(0).nsubmesh_x1;
    int Nxx = 2*Nx+1;
    // int Ny = pumi_obj.mesh(0).nsubmesh_x2;
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
            return;

        case 1:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
            i=0;
            while (i<num_x1_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
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
            return;

        case 2:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
            i=0;
            while (i<num_x1_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
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
            return;

        case 3:
            *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
            i=0;
            while (i<num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
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
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
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
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
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
            return;

        case 6:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
            i=0;
            while (i<num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
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
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
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
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
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
            return;

    }

}

} // namespace pumi

#endif
