#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumi_routines.h"

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) to compute the total number of elements in the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
int pumi_total_elements(pumi_mesh_t *pumi_mesh)
{
  int Nel_total;
  if (pumi_mesh->ndim == 1){
    Nel_total = pumi_total_elements_1D(pumi_mesh);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
  return Nel_total;
}

/*
* \brief Computes and returns the total number of elements in the mesh for 1D domain
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh)
{
  int Nel_total = 0;
  for (int isubmesh=0; isubmesh< pumi_mesh->nsubmeshes; isubmesh++){
    Nel_total += ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel;
  }
  return Nel_total;
}

/*
* \brief Computes and returns the coordinate of the left endpoint of the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh)
{
  double global_x_left = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + 0)->x_left;
  return global_x_left;
}

/*
* \brief Computes and returns the coordinate of the right endpoint of the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh)
{
  double global_x_right = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (pumi_mesh->nsubmeshes - 1))->x_right;
  return global_x_right;
}

/*
* \brief Calls the appropriate subroutine (based on the dimension of the problem) to locate the cell number of a particle
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] particle_coordinate coordinate of the particle whose cell number and local weights is to be evaluated
* \param[out] *particle_cell address of the variable where the particle cell is to be stored
* \param[out] *cell_weight address of the variable where the local weight is to be stored
*/
void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight)
{
  if (pumi_mesh->ndim == 1){
    pumi_locatepoint_1D(pumi_mesh, particle_coordinate, particle_cell, cell_weight);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
}

/*
* \brief Locates the particle in a submesh segment and calls appropriate subroutine to obtain the local cell number of the particle inside the segment
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] particle_coordinate coordinate of the particle whose cell number and local weights is to be evaluated
* \param[out] *particle_cell address of the variable where the particle cell is to be stored
* \param[out] *cell_weight address of the variable where the local weight is to be stored
*/
void pumi_locatepoint_1D(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight){
  int N_cumulative[pumi_mesh->nsubmeshes];
  N_cumulative[0] = 0;
  for (int isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    N_cumulative[isubmesh] = N_cumulative[isubmesh-1] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->submesh_total_Nel;
  }

  for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    if (particle_coordinate>=((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left && particle_coordinate<=((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right){
      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
        if (particle_coordinate <= ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_x_right && particle_coordinate>= ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left ){
          int local_lBL_cell;
          double local_lBL_weight;
          pumi_meshflag_t submeshflag = leftBL;
          pumi_locatepoint_BL_1D(&local_lBL_cell, &local_lBL_weight, particle_coordinate, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_left_r, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel, submeshflag );
          *particle_cell = N_cumulative[isubmesh]+local_lBL_cell;
          *cell_weight = local_lBL_weight;
          break;
        }
      }
      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
        if (particle_coordinate>=((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left && particle_coordinate<=((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_right){
          int local_uni_cell;
          double local_uni_weight;
          pumi_locatepoint_uniform_1D(&local_uni_cell, &local_uni_weight, particle_coordinate, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel);
          *particle_cell = N_cumulative[isubmesh] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel + local_uni_cell;
          *cell_weight = local_uni_weight;
          break;
        }
      }
      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
        if (particle_coordinate >= ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_x_left && particle_coordinate<=((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right){
          int local_rBL_cell;
          double local_rBL_weight;
          pumi_meshflag_t submeshflag = rightBL;
          pumi_locatepoint_BL_1D(&local_rBL_cell, &local_rBL_weight, particle_coordinate, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_right_r, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel, submeshflag );
          *particle_cell = N_cumulative[isubmesh]+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel +((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel - local_rBL_cell - 1;
          *cell_weight = local_rBL_weight;
          break;
        }
      }
    }
  }
}

/*
* \brief Computes the local cell number and weight of a particle located in the uniform mesh segment of a submesh block
* \param[out] *cell address of the variable where the local particle cell (w.r.t to the uniform mesh segment) is to be stored
* \param[out] *weight address of the variable where the local weight (based on linear weighting) in the located cell is to be stored
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
* \param[in] uniform_x_left left end coordinate of the uniform mesh segment in the submesh block
* \param[in] uniform_t0 size of elements in uniform mesh segment of the submesh block
* \param[in] uniform_Nel number of elements in the uniform mesh segment of the submesh block
*/
void pumi_locatepoint_uniform_1D(int *cell, double *weight, double coord, double uniform_x_left, double uniform_t0, int uniform_Nel)
{
   *cell = ((coord-uniform_x_left)/uniform_t0);
   if (*cell == uniform_Nel){ // when particle is at uniform_x_right
     *cell = *cell-1;
   }
   *weight = (coord - (uniform_x_left + uniform_t0*(*cell)))/uniform_t0; //local weight due to the charge in the cell
}

 /*
 * \brief Computes the local cell number and weight of a particle located in the left/right BL segment of a submesh block
 * \param[out] *cell address of the variable where the local particle cell (w.r.t to the left/right BL segment) is to be stored
 * \param[out] *weight address of the variable where the local weight (based on linear weighting) in the located cell is to be stored
 * \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
 * \param[in] x_end right/left end coordinate of the left/right BL segment in the submesh block respectivley
 * \param[in] r growth ratio of the elements in the left/right BL segment of the submesh block
 * \param[in] t0 size of first (leftmost/rightmost) element in the left/right BL segment inside the submesh block respectivley
 * \param[in] local_Nel number of elements in the left/right BL segment of the submesh block
 * \param[in] pumi_flag Mesh flag enum that defines the types of meshing in each segment of the submesh block
 * \details The particle cell number calculations are adjusted based on the pumi_flag that is passed to the routine. For this routine, leftBL and rightBL are the only valid inputs
 for pumi_flag
 */
void pumi_locatepoint_BL_1D(int *cell, double *weight, double coord, double x_end, double r, double t0, double log_r, double r_t0_ratio, int local_Nel, pumi_meshflag_t pumi_flag)
{
   *cell = log(1 + (fabs(x_end-coord))*r_t0_ratio)/log_r;
   if (*cell == local_Nel){ // when particle is at lBL_x_right or rBL_x_left
     *cell = *cell-1;
   }
   double r_power_cell = pow(r,*cell);

   if (pumi_flag == leftBL){
     *weight = (coord - (x_end + (r_power_cell-1.0)/r_t0_ratio))/(t0*r_power_cell);
   }
   if (pumi_flag == rightBL){
     *weight = 1 - ((x_end - (r_power_cell-1.0)/r_t0_ratio) - coord)/(t0*r_power_cell);
   }
}

 /*
 * \brief Computes the covolume at each node in the mesh
 * \param[in] *pumi_mesh pointer object to struct pumi_mesh
 * \param[in] Nel_total Total number of elements in the mesh
 * \param[out] pointer to array of nodal covolume (to be populated after this function call)

void pumi_compute_covolume_1D(pumi_mesh_t *pumi_mesh, int Nel_total, double *covolume){
   for (int inode=0; inode<Nel_total+1; inode++){
     covolume[inode] = 0.0;
   }
   int N_cumulative[pumi_mesh->nsubmeshes];
   N_cumulative[0] = 0;
   for (int isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
     N_cumulative[isubmesh] = N_cumulative[isubmesh-1] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->submesh_total_Nel;
   }

   for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
     double tmp_submesh_elem[((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel];

     tmp_submesh_elem[0] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
     for (int iCell=1; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel; iCell++){
       tmp_submesh_elem[iCell] = tmp_submesh_elem[iCell-1]*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r;
     }
     for (int iCell=0; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel; iCell++){
       tmp_submesh_elem[iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
     }
     tmp_submesh_elem[0+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-1);
     for (int iCell=1; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel; iCell++){
       tmp_submesh_elem[iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel] = tmp_submesh_elem[iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel-1]/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r;
     }

     for (int inode=0; inode<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel+1; inode++){
       if (inode==0){
         covolume[N_cumulative[isubmesh]+inode] += tmp_submesh_elem[inode]/2;
       }
       else if (inode==((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel){
         covolume[N_cumulative[isubmesh]+inode] += tmp_submesh_elem[inode-1]/2;
       }
       else{
         covolume[N_cumulative[isubmesh]+inode] = (tmp_submesh_elem[inode-1]+tmp_submesh_elem[inode])/2;
       }
     }
   }
}
*/

/*
* \brief Computes and returns the covolume for a given node in the mesh
* \param[in] node number
* \param[in] Nel_total Total number of elements in the mesh
* \param[in] pointer to array of element sizes
*/
double pumi_compute_covolume_1D(int inode, int Nel_total, double *elemsize){
  double covolume;
  if (inode == 0){
    covolume = elemsize[inode]/2.0;
  }
  else if (inode == Nel_total){
    covolume = elemsize[Nel_total-1]/2.0;
  }
  else if (inode > 0 && inode < Nel_total){
    covolume = (elemsize[inode-1]+elemsize[inode])/2.0;
  }
  else{
    printf("\tInvalid node number for covolume\n");
    exit(0);
  }
  return covolume;
}

/*
* \brief Computes the element sizes in the mesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] Nel_total Total number of elements in the mesh
* \param[out] pointer to array of element size (to be populated after this function call)
*/
void pumi_compute_elemsize_1D(pumi_mesh_t *pumi_mesh, int Nel_total, double *elemsize){
  for (int iel=0; iel<Nel_total; iel++){
    elemsize[iel] = 0.0;
  }
  int N_cumulative[pumi_mesh->nsubmeshes];
  N_cumulative[0] = 0;
  for (int isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    N_cumulative[isubmesh] = N_cumulative[isubmesh-1] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->submesh_total_Nel;
  }

  for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    elemsize[N_cumulative[isubmesh]+0] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
    for (int iCell=1; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel; iCell++){
      elemsize[N_cumulative[isubmesh]+iCell] = elemsize[N_cumulative[isubmesh]+iCell-1]*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r;
    }
    for (int iCell=0; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel; iCell++){
      elemsize[N_cumulative[isubmesh]+iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    }
    elemsize[N_cumulative[isubmesh]+0+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-1);
    for (int iCell=1; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel; iCell++){
      elemsize[N_cumulative[isubmesh]+iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel] = elemsize[N_cumulative[isubmesh]+iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel-1]/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r;
    }
  }
}


/*
* \brief Computes the grading ratio array
* \param[in] pointer to array of element size
* \param[in] Nel_total Total number of elements in the mesh
* \param[out] pointer to array of grading ratios (to be populated after this function call)
*/
void pumi_compute_nodal_gradingratio_1D(double *elemsize, int Nel_total, double *gradingratio){
  for (int i=0; i<Nel_total-1; i++){
    gradingratio[i] = elemsize[i+1]/elemsize[i];
  }
}
