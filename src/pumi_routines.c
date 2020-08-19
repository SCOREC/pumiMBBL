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
  int Nel_total = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (pumi_mesh->nsubmeshes_x1-1))->submesh_Nel + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (pumi_mesh->nsubmeshes_x1-1))->Nel_cumulative;
  return Nel_total;
}

/*
* \brief Computes and returns the total number of elements in a submesh block for 1D domain
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param ID of the submesh
*/
int pumi_submesh_total_elements_1D(pumi_mesh_t *pumi_mesh, int isubmesh)
{
  return ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel;
}

/*
* \brief Computes and returns the coordinate of the left endpoint of the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
double pumi_global_x_left_1D(pumi_mesh_t *pumi_mesh)
{
  double global_x1_left = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + 0)->coord_min;
  return global_x1_left;
}

/*
* \brief Computes and returns the coordinate of the right endpoint of the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
double pumi_global_x_right_1D(pumi_mesh_t *pumi_mesh)
{
  double global_x1_right = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (pumi_mesh->nsubmeshes_x1 - 1))->coord_max;
  return global_x1_right;
}


/*
* \brief Computes and returns the covolume for a given node in the mesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] node number
*/
double pumi_return_covolume_1D(pumi_mesh_t* pumi_mesh, int inode){
  double covolume;
  if (inode == 0){
    covolume = pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_right_offset, pumi_x1)/2.0;
  }
  else if (inode == pumi_mesh->pumi_Nel_total_x1){
    covolume = pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_left_offset, pumi_x1)/2.0;
  }
  else if (inode > 0 && inode < pumi_mesh->pumi_Nel_total_x1){
    covolume = pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_left_offset, pumi_x1)/2.0 + pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_right_offset, pumi_x1)/2.0;
  }
  else{
    printf("\tInvalid node number for covolume\n");
    exit(0);
  }
  return covolume;
}

double pumi_return_covolume_2D(pumi_mesh_t* pumi_mesh, int inp_x1, int inp_x2){
    double dx1_left, dx1_right, dx2_bottom, dx2_top;
    if (inp_x1 == 0){
      dx1_left = 0.0;
      dx1_right = pumi_return_elemsize(pumi_mesh, inp_x1, pumi_elem_on_right_offset, pumi_x1);
    }
    else if (inp_x1 == pumi_mesh->pumi_Nel_total_x1){
      dx1_left = pumi_return_elemsize(pumi_mesh, inp_x1, pumi_elem_on_left_offset, pumi_x1);
      dx1_right = 0.0;
    }
    else if (inp_x1 > 0 && inp_x1 < pumi_mesh->pumi_Nel_total_x1){
      dx1_left = pumi_return_elemsize(pumi_mesh, inp_x1, pumi_elem_on_left_offset, pumi_x1);
      dx1_right = pumi_return_elemsize(pumi_mesh, inp_x1, pumi_elem_on_right_offset, pumi_x1);
    }
    else{
      printf("\tInvalid x1 - node number for covolume\n");
      exit(0);
    }

    if (inp_x2 == 0){
      dx2_bottom = 0.0;
      dx2_top = pumi_return_elemsize(pumi_mesh, inp_x2, pumi_elem_on_top_offset, pumi_x2);
    }
    else if (inp_x2 == pumi_mesh->pumi_Nel_total_x2){
      dx2_bottom = pumi_return_elemsize(pumi_mesh, inp_x2, pumi_elem_on_bottom_offset, pumi_x2);
      dx2_top = 0.0;
    }
    else if (inp_x2 > 0 && inp_x2 < pumi_mesh->pumi_Nel_total_x2){
      dx2_bottom = pumi_return_elemsize(pumi_mesh, inp_x2, pumi_elem_on_bottom_offset, pumi_x2);
      dx2_top = pumi_return_elemsize(pumi_mesh, inp_x2, pumi_elem_on_top_offset, pumi_x2);
    }
    else{
      printf("\tInvalid x2 - node number for covolume\n");
      exit(0);
    }

    double covolume = (dx1_left*dx2_bottom + dx1_right*dx2_bottom +  dx1_left*dx2_top + dx1_right*dx2_top)/4.0;
    return covolume;

}

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) to compute BL element sizes
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_ON(pumi_mesh_t *pumi_mesh){
  if (pumi_mesh->ndim == 1){
    pumi_BL_elemsize_ON_1D(pumi_mesh);
  }
  else {
    pumi_BL_elemsize_ON_2D(pumi_mesh);
    //printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    //exit(0);
  }
}

/*
* \brief Computes and stores BL element size in submesh stuct member
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_ON_1D(pumi_mesh_t *pumi_mesh){
    int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
      int left_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize = (double*) malloc(left_Nel*sizeof(double));
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords = (double*) malloc((left_Nel+1)*sizeof(double));
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 1) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min+((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
      int iCell;
      for (iCell=1; iCell<left_Nel-1; iCell++){
        *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
        *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell)) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell);
      }
      iCell = left_Nel-1;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max;
    }
    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
      int right_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize = (double*) malloc(right_Nel*sizeof(double));
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords = (double*) malloc((right_Nel+1)*sizeof(double));
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0*pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel-1);
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 1) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 0) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + 0);
      int iCell;
      for (iCell=1; iCell<right_Nel-1; iCell++){
        *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
        *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell)) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell);
      }
      iCell = right_Nel-1;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
      *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max;
    }
  }
}

/*
* \brief Computes and stores BL element size in submesh stuct member
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_ON_2D(pumi_mesh_t *pumi_mesh){
    int isubmesh, jsubmesh;
    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
            int left_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel;
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize = (double*) malloc(left_Nel*sizeof(double));
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords = (double*) malloc((left_Nel+1)*sizeof(double));
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 1) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min+((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
            int iCell;
            for (iCell=1; iCell<left_Nel-1; iCell++){
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell)) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell);
            }
            iCell = left_Nel-1;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max;
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
            int right_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel;
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize = (double*) malloc(right_Nel*sizeof(double));
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords = (double*) malloc((right_Nel+1)*sizeof(double));
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0*pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel-1);
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 1) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + 0) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + 0);
            int iCell;
            for (iCell=1; iCell<right_Nel-1; iCell++){
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell)) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell);
            }
            iCell = right_Nel-1;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + (iCell-1))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (iCell+1)) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max;
        }
    }

    for (jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & bottomBL){
            int bottom_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->submesh_Nel;
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize = (double*) malloc(bottom_Nel*sizeof(double));
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords = (double*) malloc((bottom_Nel+1)*sizeof(double));
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->t0;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->coord_min;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + 1) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->coord_min+((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->t0;
            int iCell;
            for (iCell=1; iCell<bottom_Nel-1; iCell++){
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + (iCell-1))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->r;
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + (iCell+1)) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + (iCell)) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + iCell);
            }
            iCell = bottom_Nel-1;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + (iCell-1))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->r;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + (iCell+1)) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->coord_max;
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & topBL){
            int top_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->submesh_Nel;
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize = (double*) malloc(top_Nel*sizeof(double));
            ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords = (double*) malloc((top_Nel+1)*sizeof(double));
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->t0*pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->r,((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->submesh_Nel-1);
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + 0) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->coord_min;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + 1) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + 0) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + 0);
            int iCell;
            for (iCell=1; iCell<top_Nel-1; iCell++){
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + (iCell-1))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->r;
                *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + (iCell+1)) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + (iCell)) + *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + iCell);
            }
            iCell = top_Nel-1;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + iCell) = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize + (iCell-1))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->r;
            *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords + (iCell+1)) = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->coord_max;
        }
    }
}

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) to free allocated memory to compute BL elem sizes
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_OFF(pumi_mesh_t *pumi_mesh){
  if (pumi_mesh->ndim == 1){
    pumi_BL_elemsize_OFF_1D(pumi_mesh);
  }
  else {
    pumi_BL_elemsize_OFF_2D(pumi_mesh);
    //printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    //exit(0);
  }
}

/*
* \brief Frees allocated memory to compute 1D BL elem sizes for all submeshes
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_OFF_1D(pumi_mesh_t *pumi_mesh){
    pumi_mesh->BL_elem_coords_cache_flag = 0;
    int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
      free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize);
      free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords);
    }
    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
      free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize);
      free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords);
    }
  }
}

/*
* \brief Frees allocated memory to compute 1D BL elem sizes for all submeshes
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_OFF_2D(pumi_mesh_t *pumi_mesh){
    pumi_mesh->BL_elem_coords_cache_flag = 0;
    int isubmesh, jsubmesh;
    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize);
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords);
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize);
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords);
        }
    }

    for (jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & bottomBL){
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize);
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords);
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & topBL){
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_elemsize);
            free(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->BL_coords);
        }
    }
}


double pumi_return_gradingratio(pumi_mesh_t *pumi_mesh, int node, int dir){
    pumi_submesh_t *submeshes;
    int Nel_total, nsubmeshes;
    if (dir == pumi_x1){
        submeshes = (pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1;
        Nel_total = pumi_mesh->pumi_Nel_total_x1;
        nsubmeshes = pumi_mesh->nsubmeshes_x1;
    }
    else{
        submeshes = (pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2;
        Nel_total = pumi_mesh->pumi_Nel_total_x2;
        nsubmeshes = pumi_mesh->nsubmeshes_x2;
    }

    if (node == 0 || node == Nel_total){
        printf("Grading ratio not defined for the first and last node of the domain -- Terminating \n");
        exit(0);
    }
    else{
        int isubmesh;
        for (isubmesh=0; isubmesh<nsubmeshes; isubmesh++){
            int submesh_min_node = (submeshes + isubmesh)->Nel_cumulative;
            int submesh_max_node = submesh_min_node + (submeshes + isubmesh)->submesh_Nel;
            if (node > submesh_min_node && node < submesh_max_node){
                return (submeshes + isubmesh)->r;
            }
            else if (node == (submeshes + isubmesh)->Nel_cumulative){
                double min_elem_size;
                double max_elem_size;
                // On min of the node
                if ((submeshes + isubmesh-1)->pumi_flag & rightBL){
                    min_elem_size = (submeshes + isubmesh-1)->t0;
                }
                if ((submeshes + isubmesh-1)->pumi_flag & uniform){
                    min_elem_size = (submeshes + isubmesh-1)->t0;
                }
                if ((submeshes + isubmesh-1)->pumi_flag & leftBL){
                    double t0 = (submeshes + isubmesh-1)->t0;
                    double r = (submeshes + isubmesh-1)->r;
                    double Nel = (submeshes + isubmesh-1)->submesh_Nel;
                    min_elem_size = t0*pow(r,Nel-1);
                }

                // On max of the node
                if ((submeshes + isubmesh)->pumi_flag & leftBL){
                    max_elem_size = (submeshes + isubmesh)->t0;
                }
                if ((submeshes + isubmesh)->pumi_flag & uniform){
                    max_elem_size = (submeshes + isubmesh)->t0;
                }
                if ((submeshes + isubmesh)->pumi_flag & rightBL){
                    double t0 = (submeshes + isubmesh)->t0;
                    double r = (submeshes + isubmesh)->r;
                    double Nel = (submeshes + isubmesh)->submesh_Nel;
                    max_elem_size = t0*pow(r,Nel-1);
                }

                return max_elem_size/min_elem_size;
            }
        }
    }
}



double pumi_return_elemsize(pumi_mesh_t *pumi_mesh, int index, int offset, int dir){
    pumi_submesh_t *submeshes;
    int Nel_total, nsubmeshes, elem, isubmesh;
    if (dir == pumi_x1){
        submeshes = (pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1;
        Nel_total = pumi_mesh->pumi_Nel_total_x1;
        nsubmeshes = pumi_mesh->nsubmeshes_x1;
    }
    else{
        submeshes = (pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2;
        Nel_total = pumi_mesh->pumi_Nel_total_x2;
        nsubmeshes = pumi_mesh->nsubmeshes_x2;
    }

    elem = index+offset;

    if (elem < 0){
        elem = 0;
    }
    if (elem >= Nel_total){
        elem = Nel_total-1;
    }

    for (isubmesh=0; isubmesh<nsubmeshes; isubmesh++){
        int submesh_min_elem = (submeshes + isubmesh)->Nel_cumulative;
        int submesh_max_elem = (submeshes + isubmesh)->Nel_cumulative + (submeshes + isubmesh)->submesh_Nel-1;

        if (elem >= submesh_min_elem && elem <= submesh_max_elem){
            int local_cell = elem - submesh_min_elem;
            //return (pumi_calc_elem_size_x2(pumi_mesh, isubmesh, local_cell));
            return (pumi_calc_elem_size_fnptr[dir][isubmesh](pumi_mesh, isubmesh, local_cell));
        }
    }
}


/*
* \brief Returns smallest element size in the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
double pumi_return_smallest_elemsize(pumi_mesh_t *pumi_mesh){
  double smallest_elemsize;
  int isubmesh = 0;

  if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
    smallest_elemsize = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (isubmesh))->t0;
  }
  else{
    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
      smallest_elemsize = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (isubmesh))->t0;
    }
    else{
      smallest_elemsize = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
    }
  }
  for (isubmesh=1; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
    double new_smallest_elemsize;
    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
      new_smallest_elemsize = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (isubmesh))->t0;
    }
    else{
      if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
        new_smallest_elemsize = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + (isubmesh))->t0;
      }
      else{
        new_smallest_elemsize = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
      }
    }

    if (new_smallest_elemsize < smallest_elemsize){
      smallest_elemsize = new_smallest_elemsize;
    }
  }

  return smallest_elemsize;
}


/*
* \brief Locates submesh ID and local cell ID of a initialized particle with global search and analytically
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param coords - new coordinate of the particle
* \param prev_submeshID - old submesh ID of the pushed particle (i.e. before the push)
* \param[out] pointers to new submesh ID and new local cell which will be populated inside the routine
*/
void pumi_locate_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int *submeshID, int *cellID, int dir){
    //if (pumi_mesh->nsubmeshes_x1 == 1){
    //    *submeshID = 0;
    //}
    //else{
        int isubmesh;
        int submesh_located = 0;
        for (isubmesh=1; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
            if (coords < ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min){
                *submeshID = isubmesh-1;
                submesh_located++;
                break;
            }
        }
        if (!(submesh_located)){
            *submeshID = pumi_mesh->nsubmeshes_x1-1;
        }
    //}

    *cellID = pumi_locatecell_fnptr[dir][*submeshID](pumi_mesh, *submeshID, coords);

}

/*
* \brief Calculates new submesh ID and new local cell ID of a pushed particle using adjacency searches
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param coords - new coordinate of the particle
* \param prev_submeshID - old submesh ID of the pushed particle (i.e. before the push)
* \param prev_cellID - old local cell ID of the pushed particle (i.e. before the push)
* \param[out] pointers to new submesh ID and new local cell which will be populated inside the routine
*/
void pumi_update_submesh_and_update_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID, int dir){
    //if (pumi_mesh->nsubmeshes_x1 == 1){
    //    *submeshID = 0;
    //}
    //else{
        *submeshID = prev_submeshID;
        while(coords<((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + *submeshID)->coord_min){
            *submeshID -= 1;
            prev_cellID = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + *submeshID)->submesh_Nel - 1;
        }

        while(coords>((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + *submeshID)->coord_max){
            *submeshID += 1;
            prev_cellID = 0;
        }
    //}
    *cellID = pumi_updatecell_fnptr[dir][*submeshID](pumi_mesh, *submeshID, prev_cellID, coords);
}

/*
* \brief subroutine to call relevant weight-calculate routine using function pointers
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the block
* \param[in] local cell ID in the block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
* \param[out] pointers to global cell and weight2 which will be populated inside the routine
*/
void pumi_calc_weights(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight, int dir){
    pumi_calc_weights_fnptr[dir][isubmesh](pumi_mesh, isubmesh, local_cell, coord, global_cell, weight);
}

/*
* \brief subroutine to locate the cell number of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_locatecell_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh_x1, double coord_x1){
    int icell_x1 = (coord_x1 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh_x1)->coord_min)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh_x1)->t0;
    return icell_x1;
}

int pumi_locatecell_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh_x2, double coord_x2){
    int icell_x2 = (coord_x2 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->coord_min)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->t0;
    return icell_x2;
}

/*
* \brief subroutine to update the cell number of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] local cell ID of particle in the uniform block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_updatecell_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = (coord - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
    return icell;
}

int pumi_updatecell_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = (coord - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_min)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0;
    return icell;
}

/*
* \brief subroutine to calculate golbal cell ID and weight contribution of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] local cell ID of particle in the uniform block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
void pumi_calc_weights_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    *weight = (coord - (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0*local_cell))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
    *global_cell = local_cell + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative;
}

void pumi_calc_weights_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    *weight = (coord - (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_min + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0*local_cell))/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0;
    *global_cell = local_cell + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->Nel_cumulative;
}

/*
* \brief subroutine to locate the cell number of a particle inside a leftBL block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_locatecell_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->log_r;
    return icell;
}

int pumi_locatecell_in_bottomBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_min-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->log_r;
    return icell;
}

/*
* \brief subroutine to update the local cell number of a particle inside a leftBL block with adjacency search
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] local cell ID in the leftBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
*/
int pumi_updatecell_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    while(coord < *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + icell)){
        icell -= 1;
    }

    while(coord > *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (icell+1))){
        icell += 1;
    }
    return icell;
}

int pumi_updatecell_in_bottomBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    while(coord < *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_coords + icell)){
        icell -= 1;
    }

    while(coord > *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_coords + (icell+1))){
        icell += 1;
    }
    return icell;
}
/*
* \brief subroutine to update the local cell number of a particle inside a leftBL block analytically
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] local cell ID in the leftBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
*/
int pumi_updatecell_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->log_r;
    return icell;
}

int pumi_updatecell_in_bottomBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_min-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->log_r;
    return icell;
}

/*
* \brief subroutine to calculate golbal cell ID and weight contribution of a particle inside a leftBL block analytically
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] local cell ID in the leftBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
* \param[out] pointers to global cell and weight2 which will be populated inside the routine
*/
void pumi_calc_weights_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,local_cell);
    *global_cell = local_cell + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative;
    *weight = (coord - (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min + (r_power_cell-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio))/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0*r_power_cell);
}

void pumi_calc_weights_in_bottomBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r,local_cell);
    *global_cell = local_cell + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->Nel_cumulative;
    *weight = (coord - (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_min + (r_power_cell-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r_t0_ratio))/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0*r_power_cell);
}

/*
* \brief subroutine to calculate golbal cell ID and weight contribution of a particle inside a leftBL block using cached elementsize and nodal coords
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] local cell ID in the leftBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
* \param[out] pointers to global cell and weight2 which will be populated inside the routine
*/
void pumi_calc_weights_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    *global_cell = local_cell + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative;
//    *weight = (coord - (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->x_left + (*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->leftBL_elemsize + local_cell) - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->lBL_t0)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->left_r - 1.0)))/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->leftBL_elemsize + local_cell));
    *weight = (coord - *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + local_cell))/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + local_cell));
}

void pumi_calc_weights_in_bottomBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    *global_cell = local_cell + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->Nel_cumulative;
//    *weight = (coord - (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->x_left + (*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->leftBL_elemsize + local_cell) - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->lBL_t0)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->left_r - 1.0)))/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->leftBL_elemsize + local_cell));
    *weight = (coord - *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_coords + local_cell))/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_elemsize + local_cell));
}

/*
* \brief subroutine to locate the cell number of a particle inside a rightBL block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_locatecell_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->log_r;
    return ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel - icell - 1;
}

int pumi_locatecell_in_topBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_max-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->log_r;
    return ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->submesh_Nel - icell - 1;
}

/*
* \brief subroutine to update the local cell number of a particle inside a rightBL block with adjacency search
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] local cell ID in the rightBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
*/
int pumi_updatecell_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    while(coord < *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + icell)){
        icell -= 1;
    }

    while(coord > *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (icell+1))){
        icell += 1;
    }
    return icell;
}

int pumi_updatecell_in_topBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    while(coord < *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_coords + icell)){
        icell -= 1;
    }

    while(coord > *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_coords + (icell+1))){
        icell += 1;
    }
    return icell;
}

/*
* \brief subroutine to update the local cell number of a particle inside a rightBL block analytically
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
*/
int pumi_updatecell_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->log_r;
    return ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel - icell - 1;
}

int pumi_updatecell_in_topBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = log(1 + (fabs(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_max-coord))*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r_t0_ratio)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->log_r;
    return ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->submesh_Nel - icell - 1;
}

/*
* \brief subroutine to calculate golbal cell ID and weight contribution of a particle inside a rightBL block analytically
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] local cell ID in the rightBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
* \param[out] pointers to global cell and weight2 which will be populated inside the routine
*/
void pumi_calc_weights_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,local_cell);
    *global_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel - local_cell - 1 ;
    *weight = 1 - ((((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max - (r_power_cell-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio) - coord)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0*r_power_cell);
}

void pumi_calc_weights_in_topBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->submesh_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r,local_cell);
    *global_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->Nel_cumulative + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->submesh_Nel - local_cell - 1 ;
    *weight = 1 - ((((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->coord_max - (r_power_cell-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r_t0_ratio) - coord)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0*r_power_cell);
}

/*
* \brief subroutine to calculate golbal cell ID and weight contribution of a particle inside a rightBL block using cached elementsize and nodal coords
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] local cell ID in the rightBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
* \param[out] pointers to global cell and weight2 which will be populated inside the routine
*/
void pumi_calc_weights_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
//    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->right_Nel - local_cell - 1;
    *global_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative + local_cell ;
//    *weight = 1 - ((((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->x_right - (*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->rightBL_elemsize + local_cell) - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->rBL_t0)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->right_r-1.0)) - coord)/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->rightBL_elemsize + local_cell));
    *weight = (coord - *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + local_cell))/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + local_cell));
}

void pumi_calc_weights_in_topBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
//    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->right_Nel - local_cell - 1;
    *global_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->Nel_cumulative + local_cell ;
//    *weight = 1 - ((((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->x_right - (*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->rightBL_elemsize + local_cell) - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->rBL_t0)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->right_r-1.0)) - coord)/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->rightBL_elemsize + local_cell));
    *weight = (coord - *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_coords + local_cell))/(*(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_elemsize + local_cell));
}

/*
* \brief subroutine that returns global cell number from submesh ID and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
int pumi_global_cell_ID(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative + local_cell);
}

/*
* \brief subroutine that calculates the element node coordinates from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[out] pointer to variable where left node coord is to be stored
* \param[out] pointer to variable where right node coord is to be stored
*/
void pumi_calc_node_coords(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node){
    pumi_calc_node_coords_fnptr[isubmesh](pumi_mesh, isubmesh, local_cell, left_node, right_node);
}

/*
* \brief subroutine that calculates the element node coordinates in uniform block from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[out] pointer to variable where left node coord is to be stored
* \param[out] pointer to variable where right node coord is to be stored
*/
void pumi_calc_node_coords_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node){
    *left_node = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0*local_cell;
    *right_node = *left_node + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
}

/*
* \brief subroutine that calculates the element node coordinates in leftBL block using cached node coords, submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[out] pointer to variable where left node coord is to be stored
* \param[out] pointer to variable where right node coord is to be stored
*/
void pumi_calc_node_coords_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node){
    *left_node = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + local_cell);
    *right_node = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (local_cell+1));
}

/*
* \brief subroutine that calculates the element node coordinates in leftBL analytically from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[out] pointer to variable where left node coord is to be stored
* \param[out] pointer to variable where right node coord is to be stored
*/
void pumi_calc_node_coords_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node){
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,local_cell);
    *left_node = (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min + (r_power_cell-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio);
    *right_node = *left_node + r_power_cell*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
}

/*
* \brief subroutine that calculates the element node coordinates in rightBL block using cached node coords, submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[out] pointer to variable where left node coord is to be stored
* \param[out] pointer to variable where right node coord is to be stored
*/
void pumi_calc_node_coords_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node){
    *left_node = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + local_cell);
    *right_node = *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_coords + (local_cell+1));
}

/*
* \brief subroutine that calculates the element node coordinates in rightBL analytically from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
* \param[out] pointer to variable where left node coord is to be stored
* \param[out] pointer to variable where right node coord is to be stored
*/
void pumi_calc_node_coords_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double *left_node, double *right_node){
    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,local_cell);
    *right_node = (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max - (r_power_cell-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r_t0_ratio);
    *left_node = *right_node - r_power_cell*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0;
}

/*
* \brief subroutine that calculates the element size from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
double pumi_calc_elem_size(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, int dir){
    return (pumi_calc_elem_size_fnptr[dir][isubmesh](pumi_mesh, isubmesh, local_cell));
}


/*
* \brief subroutine that calculates the element size in uniform block from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
double pumi_calc_elem_size_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0);
}

double pumi_calc_elem_size_in_uni_x1(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0);
}

double pumi_calc_elem_size_in_uni_x2(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0);
}

/*
* \brief subroutine that calculates the element size in leftBL block using cached element sizes, submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
double pumi_calc_elem_size_in_leftBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + local_cell);
}

double pumi_calc_elem_size_in_bottomBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_elemsize + local_cell);
}

/*
* \brief subroutine that calculates the element size in leftBL analytically from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
double pumi_calc_elem_size_in_leftBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,local_cell);
    return (r_power_cell*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0);
}

double pumi_calc_elem_size_in_bottomBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r,local_cell);
    return (r_power_cell*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0);
}

/*
* \brief subroutine that calculates the element size in rightBL block using cached element sizes, submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
double pumi_calc_elem_size_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->BL_elemsize + local_cell);
}

double pumi_calc_elem_size_in_topBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return *(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->BL_elemsize + local_cell);
}

/*
* \brief subroutine that calculates the element size in rightBL analytically from submesh and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
double pumi_calc_elem_size_in_rightBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->r,local_cell);
    return (r_power_cell*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->t0);
}

double pumi_calc_elem_size_in_topBL_analytic(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    local_cell = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->submesh_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->r,local_cell);
    return (r_power_cell*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh)->t0);
}


/*
* \brief Assigns the appropriate subroutine for each submesh block to locate/update the local cell number of a particle in that submesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \details Allocates memory to function pointer *pumi_locatecell_fnptr, *pumi_updatecell_fnptr, *pumi_calc_weights_fnptr
and assign relevant locate/update/weight-calculate pumi routines to each submesh during pumi mesh initialization
*/
/*
void pumi_initialize_locatecell_and_calcweights_functions(pumi_mesh_t *pumi_mesh){
    pumi_locatecell_fnptr = malloc(pumi_mesh->nsubmeshes_x1*sizeof(pumi_locatecell_ptr));
    pumi_updatecell_fnptr = malloc(pumi_mesh->nsubmeshes_x1*sizeof(pumi_updatecell_ptr));
    pumi_calc_weights_fnptr = malloc(pumi_mesh->nsubmeshes_x1*sizeof(pumi_calc_weights_ptr));
    pumi_calc_node_coords_fnptr = malloc(pumi_mesh->nsubmeshes_x1*sizeof(pumi_calc_node_coords_ptr));
    pumi_calc_elem_size_fnptr = (pumi_calc_elem_size_ptr**) malloc(MAX_DIM * sizeof(pumi_calc_elem_size_ptr *));
    pumi_calc_elem_size_fnptr[0] = (pumi_calc_elem_size_ptr*) malloc(pumi_mesh->nsubmeshes_x1 * sizeof(pumi_calc_elem_size_ptr));
    pumi_calc_elem_size_fnptr[1] = NULL;
    int isubmesh;
    for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
            pumi_locatecell_fnptr[isubmesh] = &pumi_locatecell_in_leftBL;
            printf("submesh=%d -- leftBL locate cell routine initialized\n",isubmesh );
            if (pumi_mesh->BL_elem_coords_cache_flag){
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_leftBL_cached;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_leftBL_cached;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_leftBL_cached;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_leftBL_cached;
                printf("submesh=%d -- leftBL calc weight and cell update (with cache) routines initialized\n",isubmesh );
            }
            else{
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_leftBL_analytic;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_leftBL_analytic;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_leftBL_analytic;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_leftBL_analytic;
                printf("submesh=%d -- leftBL calc weight and cell update (without cache) routines initialized\n",isubmesh );
            }
        }
        else if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
            pumi_locatecell_fnptr[isubmesh] = &pumi_locatecell_in_rightBL;
            printf("submesh=%d -- rightBL locate cell routine initialized\n",isubmesh );
            if (pumi_mesh->BL_elem_coords_cache_flag){
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_rightBL_cached;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_rightBL_cached;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_rightBL_cached;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_rightBL_cached;
                printf("submesh=%d -- rightBL calc weights and cell update (with cache) routines initialized\n",isubmesh );
            }
            else{
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_rightBL_analytic;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_rightBL_analytic;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_rightBL_analytic;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_rightBL_analytic;
                printf("submesh=%d -- rightBL calc weight and cell update (without cache) routines initialized\n",isubmesh );
            }
        }
        else if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & uniform){
            pumi_locatecell_fnptr[isubmesh] = &pumi_locatecell_in_uni;
            pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_uni;
            pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_uni;
            pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_uni;
            pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_uni;
            printf("submesh=%d -- uniform routines initialized\n",isubmesh );
        }
        else{
            printf("Error in meshtype for submesh %d \n", isubmesh);
            exit(0);
        }
    }
}
*/
void pumi_initialize_locatecell_and_calcweights_functions(pumi_mesh_t *pumi_mesh){
    pumi_locatecell_fnptr = (pumi_locatecell_ptr**) malloc(MAX_DIM*sizeof(pumi_locatecell_ptr *));
    pumi_locatecell_fnptr[0] = (pumi_locatecell_ptr*) malloc(pumi_mesh->nsubmeshes_x1*sizeof(pumi_locatecell_ptr));

    pumi_calc_elem_size_fnptr = (pumi_calc_elem_size_ptr**) malloc(MAX_DIM * sizeof(pumi_calc_elem_size_ptr *));
    pumi_calc_elem_size_fnptr[0] = (pumi_calc_elem_size_ptr*) malloc(pumi_mesh->nsubmeshes_x1 * sizeof(pumi_calc_elem_size_ptr));

    pumi_updatecell_fnptr = (pumi_updatecell_ptr**) malloc(MAX_DIM * sizeof(pumi_updatecell_ptr *));
    pumi_updatecell_fnptr[0] = (pumi_updatecell_ptr*) malloc(pumi_mesh->nsubmeshes_x1 * sizeof(pumi_updatecell_ptr));

    pumi_calc_weights_fnptr = (pumi_calc_weights_ptr**) malloc(MAX_DIM * sizeof(pumi_calc_weights_ptr *));
    pumi_calc_weights_fnptr[0] = (pumi_calc_weights_ptr*) malloc(pumi_mesh->nsubmeshes_x1 * sizeof(pumi_calc_weights_ptr));

    if (pumi_mesh->ndim == 1){
        pumi_locatecell_fnptr[1] = NULL;
        pumi_calc_elem_size_fnptr[1] = NULL;
        pumi_updatecell_fnptr[1] = NULL;
        pumi_calc_weights_fnptr[1] = NULL;
    }
    else{
        pumi_locatecell_fnptr[1] = (pumi_locatecell_ptr*) malloc(pumi_mesh->nsubmeshes_x2 * sizeof(pumi_locatecell_ptr));
        pumi_calc_elem_size_fnptr[1] = (pumi_calc_elem_size_ptr*) malloc(pumi_mesh->nsubmeshes_x2 * sizeof(pumi_calc_elem_size_ptr));
        pumi_updatecell_fnptr[1] = (pumi_updatecell_ptr*) malloc(pumi_mesh->nsubmeshes_x2 * sizeof(pumi_updatecell_ptr));
        pumi_calc_weights_fnptr[1] = (pumi_calc_weights_ptr*) malloc(pumi_mesh->nsubmeshes_x2 * sizeof(pumi_calc_weights_ptr));
    }

    pumi_calc_node_coords_fnptr = malloc(pumi_mesh->nsubmeshes_x1 * sizeof(pumi_calc_node_coords_ptr));

    int isubmesh, jsubmesh;
    for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
            pumi_locatecell_fnptr[0][isubmesh] = &pumi_locatecell_in_leftBL;
            printf("submesh=%d -- leftBL locate cell routine initialized\n",isubmesh );
            if (pumi_mesh->BL_elem_coords_cache_flag){
                pumi_calc_weights_fnptr[0][isubmesh] = &pumi_calc_weights_in_leftBL_cached;
                pumi_updatecell_fnptr[0][isubmesh] = &pumi_updatecell_in_leftBL_cached;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_leftBL_cached;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_leftBL_cached;
                printf("submesh=%d -- leftBL calc weight and cell update (with cache) routines initialized\n",isubmesh );
            }
            else{
                pumi_calc_weights_fnptr[0][isubmesh] = &pumi_calc_weights_in_leftBL_analytic;
                pumi_updatecell_fnptr[0][isubmesh] = &pumi_updatecell_in_leftBL_analytic;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_leftBL_analytic;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_leftBL_analytic;
                printf("submesh=%d -- leftBL calc weight and cell update (without cache) routines initialized\n",isubmesh );
            }
        }
        else if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
            pumi_locatecell_fnptr[0][isubmesh] = &pumi_locatecell_in_rightBL;
            printf("submesh=%d -- rightBL locate cell routine initialized\n",isubmesh );
            if (pumi_mesh->BL_elem_coords_cache_flag){
                pumi_calc_weights_fnptr[0][isubmesh] = &pumi_calc_weights_in_rightBL_cached;
                pumi_updatecell_fnptr[0][isubmesh] = &pumi_updatecell_in_rightBL_cached;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_rightBL_cached;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_rightBL_cached;
                printf("submesh=%d -- rightBL calc weights and cell update (with cache) routines initialized\n",isubmesh );
            }
            else{
                pumi_calc_weights_fnptr[0][isubmesh] = &pumi_calc_weights_in_rightBL_analytic;
                pumi_updatecell_fnptr[0][isubmesh] = &pumi_updatecell_in_rightBL_analytic;
                pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_rightBL_analytic;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_rightBL_analytic;
                printf("submesh=%d -- rightBL calc weight and cell update (without cache) routines initialized\n",isubmesh );
            }
        }
        else if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & uniform){
            pumi_locatecell_fnptr[0][isubmesh] = &pumi_locatecell_in_uni_x1;
            pumi_updatecell_fnptr[0][isubmesh] = &pumi_updatecell_in_uni_x1;
            pumi_calc_weights_fnptr[0][isubmesh] = &pumi_calc_weights_in_uni_x1;
            pumi_calc_elem_size_fnptr[0][isubmesh] = &pumi_calc_elem_size_in_uni_x1;
            pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_uni;
            printf("submesh=%d -- uniform routines initialized\n",isubmesh );
        }
        else{
            printf("Error in meshtype for submesh %d \n", isubmesh);
            exit(0);
        }
    }
    if (pumi_mesh->ndim > 1){
        for(jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){
            if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & bottomBL){
                pumi_locatecell_fnptr[1][jsubmesh] = &pumi_locatecell_in_bottomBL;
                printf("submesh=%d -- bottomBL locate cell routine initialized\n",jsubmesh );
                if (pumi_mesh->BL_elem_coords_cache_flag){
                    pumi_calc_weights_fnptr[1][jsubmesh] = &pumi_calc_weights_in_bottomBL_cached;
                    pumi_updatecell_fnptr[1][jsubmesh] = &pumi_updatecell_in_bottomBL_cached;
                    pumi_calc_elem_size_fnptr[1][jsubmesh] = &pumi_calc_elem_size_in_bottomBL_cached;
                    printf("submesh=%d -- bottomBL calc weight and cell update (with cache) routines initialized\n",jsubmesh );
                }
                else{
                    pumi_calc_weights_fnptr[1][jsubmesh] = &pumi_calc_weights_in_bottomBL_analytic;
                    pumi_updatecell_fnptr[1][jsubmesh] = &pumi_updatecell_in_bottomBL_analytic;
                    pumi_calc_elem_size_fnptr[1][jsubmesh] = &pumi_calc_elem_size_in_bottomBL_analytic;
                    printf("submesh=%d -- bottomBL calc weight and cell update (without cache) routines initialized\n",jsubmesh );
                }
            }
            else if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & topBL){
                pumi_locatecell_fnptr[1][jsubmesh] = &pumi_locatecell_in_topBL;
                printf("submesh=%d -- topBL locate cell routine initialized\n",jsubmesh );
                if (pumi_mesh->BL_elem_coords_cache_flag){
                    pumi_calc_weights_fnptr[1][jsubmesh] = &pumi_calc_weights_in_topBL_cached;
                    pumi_updatecell_fnptr[1][jsubmesh] = &pumi_updatecell_in_topBL_cached;
                    pumi_calc_elem_size_fnptr[1][jsubmesh] = &pumi_calc_elem_size_in_topBL_cached;
                    printf("submesh=%d -- topBL calc weights and cell update (with cache) routines initialized\n",jsubmesh );
                }
                else{
                    pumi_calc_weights_fnptr[1][jsubmesh] = &pumi_calc_weights_in_topBL_analytic;
                    pumi_updatecell_fnptr[1][jsubmesh] = &pumi_updatecell_in_topBL_analytic;
                    pumi_calc_elem_size_fnptr[1][jsubmesh] = &pumi_calc_elem_size_in_topBL_analytic;
                    printf("submesh=%d -- topBL calc weight and cell update (without cache) routines initialized\n",jsubmesh );
                }
            }
            else if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->pumi_flag & uniform){
                pumi_locatecell_fnptr[1][jsubmesh] = &pumi_locatecell_in_uni_x2;
                pumi_updatecell_fnptr[1][jsubmesh] = &pumi_updatecell_in_uni_x2;
                pumi_calc_weights_fnptr[1][jsubmesh] = &pumi_calc_weights_in_uni_x2;
                pumi_calc_elem_size_fnptr[1][jsubmesh] = &pumi_calc_elem_size_in_uni_x2;
                printf("submesh=%d -- uniform routines initialized\n",jsubmesh );
            }
            else{
                printf("Error in meshtype for submesh %d \n", jsubmesh);
                exit(0);
            }
        }
    }
}

/*!
* \brief Deallocates/Frees the memory allocated to the funciton pointers pumi_*_fnptr
*/
/*
void pumi_finalize_locatecell_and_calcweights_functions(){
    free(pumi_locatecell_fnptr);
    free(pumi_updatecell_fnptr);
    free(pumi_calc_weights_fnptr);
    free(pumi_calc_node_coords_fnptr);
    free(pumi_calc_elem_size_fnptr[0]);
    free(pumi_calc_elem_size_fnptr[1]);
    free(pumi_calc_elem_size_fnptr);
}
*/
void pumi_finalize_locatecell_and_calcweights_functions(){
    free(pumi_locatecell_fnptr[0]);
    free(pumi_locatecell_fnptr[1]);
    free(pumi_locatecell_fnptr);
    free(pumi_updatecell_fnptr[0]);
    free(pumi_updatecell_fnptr[1]);
    free(pumi_updatecell_fnptr);
    free(pumi_calc_weights_fnptr[0]);
    free(pumi_calc_weights_fnptr[1]);
    free(pumi_calc_weights_fnptr);
    free(pumi_calc_elem_size_fnptr[0]);
    free(pumi_calc_elem_size_fnptr[1]);
    free(pumi_calc_elem_size_fnptr);
    free(pumi_calc_node_coords_fnptr);
}

int pumi_calc_elementID_and_nodeID_typeA(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    int nodeoffset1, nodeoffset3, elemoffset, elemID, icell_x2;
    icell_x2 = kcell_x2 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->Nel_cumulative;
    nodeoffset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    nodeoffset3 = nodeoffset1 + pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    elemoffset = pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2];
    elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2;
    *node1 = elemID + kcell_x2 - nodeoffset1;
    *node3 = elemID + kcell_x2 + pumi_mesh->pumi_Nnp_total_x1 - nodeoffset3;
    return (elemID-elemoffset);
}

int pumi_calc_elementID_and_nodeID_typeB(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    int nodeoffset1, nodeoffset3, index, elemoffset, elemID, icell_x2;
    icell_x2 = kcell_x2 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->Nel_cumulative;
    index = (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel_minus_1-icell_x2)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel_minus_1);
    pumi_typeB_nodeoffset_fnptr[index](pumi_mesh, isubmesh_x1, isubmesh_x2, icell_x2, &nodeoffset1, &nodeoffset3);
    elemoffset = pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2];
    elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2;
    *node1 = elemID + kcell_x2 - nodeoffset1;
    *node3 = elemID + kcell_x2 + pumi_mesh->pumi_Nnp_total_x1 - nodeoffset3;
    //if (icell_x2 <= 1){
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2];
    //}
    //else{
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2] + (icell_x2-1)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    //}
    //elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2 - (pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2]);
    return (elemID-elemoffset);
}

void pumi_typeB_nodeoffset_expression2(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2];
}

void pumi_typeB_nodeoffset_expression1(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2] + (icell_x2-1)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
}

int pumi_calc_elementID_and_nodeID_typeC(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    int nodeoffset1, nodeoffset3, index, elemoffset, elemID, icell_x2;
    icell_x2 = kcell_x2 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->Nel_cumulative;
    index = icell_x2/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel_minus_1);
    pumi_typeC_nodeoffset_fnptr[index](pumi_mesh, isubmesh_x1, isubmesh_x2, icell_x2, &nodeoffset1, &nodeoffset3);
    elemoffset = pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2];
    elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2;
    *node1 = elemID + kcell_x2 - nodeoffset1;
    *node3 = elemID + kcell_x2 + pumi_mesh->pumi_Nnp_total_x1 - nodeoffset3;
    //if (icell_x2 == ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel){
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + (icell_x2-1)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_top[isubmesh_x1][isubmesh_x2];
    //}
    //else{
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    //}
    //elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2 - (pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2]);
    return (elemID-elemoffset);
}

void pumi_typeC_nodeoffset_expression1(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
}

void pumi_typeC_nodeoffset_expression2(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_top[isubmesh_x1][isubmesh_x2];
}

int pumi_calc_elementID_and_nodeID_typeD(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    int nodeoffset1, nodeoffset3, index, elemoffset, elemID, icell_x2;
    icell_x2 = kcell_x2 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->Nel_cumulative;
    index = icell_x2/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel_minus_1) + 1 -
    (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel-1-icell_x2)/(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel-1);
    pumi_typeD_nodeoffset_fnptr[index](pumi_mesh, isubmesh_x1, isubmesh_x2, icell_x2, &nodeoffset1, &nodeoffset3);
    elemoffset = pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2];
    elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2;
    *node1 = elemID + kcell_x2 - nodeoffset1;
    *node3 = elemID + kcell_x2 + pumi_mesh->pumi_Nnp_total_x1 - nodeoffset3;
    //if (icell_x2 <= 1){
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2];
    //}
    //else if (icell_x2 == ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->submesh_Nel){
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2] + (icell_x2-2)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_top[isubmesh_x1][isubmesh_x2];
    //}
    //else{
    //    nodeoffset = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2] + (icell_x2-1)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    //}
    //elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2 - (pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2]);
    return (elemID-elemoffset);
}

void pumi_typeD_nodeoffset_expression1(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2];
}

void pumi_typeD_nodeoffset_expression2(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2] + (icell_x2-1)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
}

void pumi_typeD_nodeoffset_expression3(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int icell_x2, int *offset1, int *offset3){
    *offset1 = pumi_mesh->nodeoffset_start[isubmesh_x1][isubmesh_x2] + pumi_mesh->nodeoffset_skip_bottom[isubmesh_x1][isubmesh_x2] + (icell_x2-1)*pumi_mesh->nodeoffset_skip_mid[isubmesh_x1][isubmesh_x2];
    *offset3 = *offset1 + pumi_mesh->nodeoffset_skip_top[isubmesh_x1][isubmesh_x2];
}

int pumi_calc_elementID_and_nodeID(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    return pumi_nodeID_fnptr[isubmesh_x1][isubmesh_x2](pumi_mesh, isubmesh_x1, isubmesh_x2, kcell_x1, kcell_x2, node1, node3);
}

int pumi_calc_elementID_and_nodeID_with_global_offset(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    int nodeoffset1, nodeoffset3, elemoffset, elemID, icell_x2;
    nodeoffset1 = pumi_mesh->global_nodeoffset[isubmesh_x1][kcell_x2];
    nodeoffset3 = pumi_mesh->global_nodeoffset[isubmesh_x1][kcell_x2+1];
    icell_x2 = kcell_x2 - ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + isubmesh_x2)->Nel_cumulative;
    elemoffset = pumi_mesh->elemoffset_start[isubmesh_x1][isubmesh_x2] + icell_x2*pumi_mesh->elemoffset_skip[isubmesh_x2];
    elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2;
    *node1 = elemID + kcell_x2 - nodeoffset1;
    *node3 = elemID + kcell_x2 + pumi_mesh->pumi_Nnp_total_x1 - nodeoffset3;
    return (elemID-elemoffset);
}

int pumi_calc_elementID_and_nodeID_on_fullmesh(pumi_mesh_t* pumi_mesh, int isubmesh_x1, int isubmesh_x2, int kcell_x1, int kcell_x2, int *node1, int *node3){
    int elemID = kcell_x1 + pumi_mesh->pumi_Nel_total_x1*kcell_x2;
    *node1 = elemID + kcell_x2;
    *node3 = *node1 + pumi_mesh->pumi_Nnp_total_x1;
    return elemID;
}

int pumi_dummy_elem_node_ID(double coord_x1, double coord_x2, double dx1, double dx2, int Nel_total_x1, int *node1, int *node3){
    int kcell_x1 = floor(coord_x1/dx1);
    int kcell_x2 = floor(coord_x2/dx2);
    int kcell = kcell_x2*Nel_total_x1+kcell_x1;
    *node1 = kcell + kcell_x2;
    *node3 = *node1 + Nel_total_x1 + 1;
    return kcell;
}

int pumi_dummy_elem_node_ID_v2(int kcell_x1, int kcell_x2, double dx1, double dx2, int Nel_total_x1, int *node1, int *node3){
    int kcell = kcell_x2*Nel_total_x1+kcell_x1;
    *node1 = kcell + kcell_x2;
    *node3 = *node1 + Nel_total_x1 + 1;
    return kcell;
}

bool pumi_mesh_with_no_inactive_blocks(pumi_mesh_t *pumi_mesh){
    bool is_fullmesh = true;
    int isubmesh, jsubmesh;
    for (jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){
        for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
            if(!(pumi_mesh->isactive[isubmesh][jsubmesh])){
                is_fullmesh = false;
            }
        }
    }
    return is_fullmesh;
}

bool pumi_is_node_active(pumi_mesh_t *pumi_mesh, int inp_x1, int inp_x2){
    int isubmesh, jsubmesh, local_x1_node, local_x2_node;
    bool left_edge, right_edge, bottom_edge, top_edge;

    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){

        int submesh_left_node = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative;
        int submesh_right_node = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->Nel_cumulative + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel;
        left_edge =  false;
        right_edge = false;
        if (inp_x1 >= submesh_left_node && inp_x1 <= submesh_right_node){
            local_x1_node = inp_x1 - submesh_left_node;
            if (local_x1_node == 0){
                left_edge = true;
            }
            if (local_x1_node == ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_Nel){
                right_edge = true;
            }
            break;
        }
    }

    for (jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){

        int submesh_bottom_node = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->Nel_cumulative;
        int submesh_top_node = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->Nel_cumulative + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->submesh_Nel;
        bottom_edge =  false;
        top_edge = false;
        if (inp_x2 >= submesh_bottom_node && inp_x2 <= submesh_top_node){
            local_x2_node = inp_x2 - submesh_bottom_node;
            if (local_x2_node == 0){
                bottom_edge = true;
            }
            if (local_x2_node == ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x2 + jsubmesh)->submesh_Nel){
                top_edge = true;
            }
            break;
        }
    }

    if (pumi_mesh->isactive[isubmesh][jsubmesh]){
        return true;
    }
    else{
        if (!left_edge && !right_edge && !bottom_edge && !top_edge){
            return false;
        }
        else{
            if (left_edge){
                if (isubmesh==0){
                    return false;
                }
                else{
                    if (pumi_mesh->isactive[isubmesh-1][jsubmesh]){
                        return true;
                    }
                    else{
                        return false;
                    }
                }
            }

            if (right_edge){
                if (isubmesh==pumi_mesh->nsubmeshes_x1-1){
                    return false;
                }
                else{
                    if (pumi_mesh->isactive[isubmesh+1][jsubmesh]){
                        return true;
                    }
                    else{
                        return false;
                    }
                }
            }

            if (bottom_edge){
                if (jsubmesh==0){
                    return false;
                }
                else{
                    if (pumi_mesh->isactive[isubmesh][jsubmesh-1]){
                        return true;
                    }
                    else{
                        return false;
                    }
                }
            }

            if (top_edge){
                if (jsubmesh==pumi_mesh->nsubmeshes_x2-1){
                    return false;
                }
                else{
                    if (pumi_mesh->isactive[isubmesh][jsubmesh+1]){
                        return true;
                    }
                    else{
                        return false;
                    }
                }
            }

        }
    }

}

void pumi_initialize_nodeID_functions(pumi_mesh_t *pumi_mesh){
    pumi_nodeID_fnptr = (pumi_nodeID_ptr**) malloc(pumi_mesh->nsubmeshes_x1 * sizeof(pumi_nodeID_ptr *));
    int isubmesh, jsubmesh;
    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        pumi_nodeID_fnptr[isubmesh] = (pumi_nodeID_ptr*) malloc(pumi_mesh->nsubmeshes_x2 * sizeof(pumi_nodeID_ptr));
    }

    if (pumi_mesh_with_no_inactive_blocks(pumi_mesh)){
        printf("No INACTIVE blocks in mesh detected -- Intializing element/node ID routines without offsets\n\n");
        for (jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){
            for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
                pumi_nodeID_fnptr[isubmesh][jsubmesh] = pumi_calc_elementID_and_nodeID_on_fullmesh;
            }
        }
    }
    else{
        printf("INACTIVE blocks in mesh detected -- Intializing element/node ID routines with offsets\n\n");
        for (jsubmesh=0; jsubmesh<pumi_mesh->nsubmeshes_x2; jsubmesh++){
            for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
                if (pumi_mesh->isactive[isubmesh][jsubmesh]){
                    if (pumi_mesh->nodeoffset_cache_flag){
                        pumi_nodeID_fnptr[isubmesh][jsubmesh] = pumi_calc_elementID_and_nodeID_with_global_offset;
                    }
                    else{
                        if (pumi_mesh->blocktype[isubmesh][jsubmesh] == type_A){
                            pumi_nodeID_fnptr[isubmesh][jsubmesh] = &pumi_calc_elementID_and_nodeID_typeA;
                        }
                        else if (pumi_mesh->blocktype[isubmesh][jsubmesh] == type_B){
                            pumi_nodeID_fnptr[isubmesh][jsubmesh] = &pumi_calc_elementID_and_nodeID_typeB;
                        }
                        else if (pumi_mesh->blocktype[isubmesh][jsubmesh] == type_C){
                            pumi_nodeID_fnptr[isubmesh][jsubmesh] = &pumi_calc_elementID_and_nodeID_typeC;
                        }
                        else if (pumi_mesh->blocktype[isubmesh][jsubmesh] == type_D){
                            pumi_nodeID_fnptr[isubmesh][jsubmesh] = &pumi_calc_elementID_and_nodeID_typeD;
                        }
                    }
                }
            }
        }

        pumi_typeB_nodeoffset_fnptr[0] = &pumi_typeB_nodeoffset_expression1;
        pumi_typeB_nodeoffset_fnptr[1] = &pumi_typeB_nodeoffset_expression2;

        pumi_typeC_nodeoffset_fnptr[0] = &pumi_typeC_nodeoffset_expression1;
        pumi_typeC_nodeoffset_fnptr[1] = &pumi_typeC_nodeoffset_expression2;

        pumi_typeD_nodeoffset_fnptr[0] = &pumi_typeD_nodeoffset_expression1;
        pumi_typeD_nodeoffset_fnptr[1] = &pumi_typeD_nodeoffset_expression2;
        pumi_typeD_nodeoffset_fnptr[2] = &pumi_typeD_nodeoffset_expression3;
    }


}

void pumi_finalize_nodeID_functions(pumi_mesh_t *pumi_mesh){
    int isubmesh;
    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        free(pumi_nodeID_fnptr[isubmesh]);
    }
    free(pumi_nodeID_fnptr);
}
