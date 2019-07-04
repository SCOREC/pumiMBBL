#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumi_routines.h"

int pumi_total_elements_1D(pumi_mesh_t *pumi_mesh)
{
  int Nel_total = 0;
  for (int isubmesh=0; isubmesh< pumi_mesh->nsubmeshes; isubmesh++){
    Nel_total += ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel;
  }
  return Nel_total;
}

double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh)
{
  double global_x_left = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + 0)->x_left;
  return global_x_left;
}

double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh)
{
  double global_x_right = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (pumi_mesh->nsubmeshes - 1))->x_right;
  return global_x_right;
}

void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int *particle_cell, double *cell_weight)
{
  if (pumi_mesh->ndim == 1){
    pumi_locatepoint_1D(pumi_mesh, particle_coordinate, particle_cell, cell_weight);
  }
}

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
          pumi_locatepoint_BL_1D(&local_lBL_cell, &local_lBL_weight, particle_coordinate, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel, submeshflag );
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
          pumi_locatepoint_BL_1D(&local_rBL_cell, &local_rBL_weight, particle_coordinate, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0, ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel, submeshflag );
          *particle_cell = N_cumulative[isubmesh]+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel +((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel - local_rBL_cell - 1;
          *cell_weight = local_rBL_weight;
          break;
        }
      }
    }
  }
}

void pumi_locatepoint_uniform_1D(int *cell, double *weight, double coord, double uniform_x_left, double uniform_t0, int uniform_Nel){
   *cell = ((coord-uniform_x_left)/uniform_t0);
   if (*cell == uniform_Nel){ // when particle is at uniform_x_right
     *cell = *cell-1;
   }
   *weight = (coord - (uniform_x_left + uniform_t0*(*cell)))/uniform_t0; //local weight due to the charge in the cell
 }

void pumi_locatepoint_BL_1D(int *cell, double *weight, double coord, double x_end, double r, double t0, int local_Nel, pumi_meshflag_t pumi_flag){
   *cell = log(1 + (fabs(x_end-coord))*(r-1)/t0)/log(r);
   if (*cell == local_Nel){ // when particle is at lBL_x_right or rBL_x_left
     *cell = *cell-1;
   }

   if (pumi_flag == leftBL){
     *weight = (coord - (x_end + t0*(pow(r,*cell)-1)/(r-1)))/(t0*pow(r,*cell)); //local weight due to the charge in the cell
   }
   if (pumi_flag == rightBL){
     *weight = 1 - ((x_end - t0*(pow(r,*cell)-1)/(r-1)) - coord)/(t0*pow(r,*cell)); //local weight due to the charge in the cell
   }
 }

void pumi_compute_covolume_1D(pumi_mesh_t *pumi_mesh, int N_total, double *covolume){
   for (int inode=0; inode<N_total+1; inode++){
     covolume[inode] = 0.0;
   }
   int N_cumulative[pumi_mesh->nsubmeshes];
   N_cumulative[0] = 0;
   for (int isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
     N_cumulative[isubmesh] = N_cumulative[isubmesh-1] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->submesh_total_Nel;
   }

   for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
     double tmp_submesh_elem[((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel];

     for (int iCell=0; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel; iCell++){
       tmp_submesh_elem[iCell] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r,iCell);
     }
     for (int iCell=0; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel; iCell++){
       tmp_submesh_elem[iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
     }
     for (int iCell=0; iCell<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel; iCell++){
       tmp_submesh_elem[iCell+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel] = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-iCell-1);
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
