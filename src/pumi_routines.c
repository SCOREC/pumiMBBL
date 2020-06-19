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
  int Nel_total = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (pumi_mesh->nsubmeshes-1))->submesh_total_Nel + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (pumi_mesh->nsubmeshes-1))->Nel_cumulative;
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
/*
void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int particle_submesh, int *particle_cell, double *cell_weight)
{
  if (pumi_mesh->ndim == 1){
    pumiMBBL_locatepoint_1D(pumi_mesh, particle_coordinate, particle_submesh, particle_cell, cell_weight);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
}*/

/*
* \brief subroutine to locate the cell number of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] particle_coordinate coordinate of the particle whose cell number and local weights is to be evaluated
* \param[out] *particle_cell address of the variable where the particle cell is to be stored
* \param[out] *cell_weight address of the variable where the local weight is to be stored
*/
/*
void pumi_locate_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){
    *cell = (coord - ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0*(*cell)))/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    *cell += ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
}
*/
/*
* \brief subroutine to locate the cell number of a particle inside a leftBL block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] particle_coordinate coordinate of the particle whose cell number and local weights is to be evaluated
* \param[out] *particle_cell address of the variable where the particle cell is to be stored
* \param[out] *cell_weight address of the variable where the local weight is to be stored
*/
/*
void pumi_locate_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){
    *cell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_left_r;
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r,*cell);
    *cell += ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left + (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio))/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0*r_power_cell);
}
*/
/*
* \brief subroutine to locate the cell number of a particle inside a rightBL block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] particle_coordinate coordinate of the particle whose cell number and local weights is to be evaluated
* \param[out] *particle_cell address of the variable where the particle cell is to be stored
* \param[out] *cell_weight address of the variable where the local weight is to be stored
*/
/*
void pumi_locate_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){
    *cell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_right_r;
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,*cell);
    *cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel - (*cell) - 1 ;
    *weight = 1 - ((((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right - (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio) - coord)/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*r_power_cell);
}
*/
/*
* \brief Assigns the appropriate subroutine for each submesh block to locate the local cell number of a particle in that submesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \details Allocates memory to function pointer *pumi_locate_function and assign relevant locate pumi routines
during pumi mesh initialization
*/
/*
void pumi_initialize_locate_functions(pumi_mesh_t *pumi_mesh){
    pumi_locate_function = malloc(pumi_mesh->nsubmeshes*sizeof(pumi_locate_ptr));
    int isubmesh;
    for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
            pumi_locate_function[isubmesh] = &pumi_locate_in_leftBL;
            printf("submesh=%d -- leftBL routine initialized\n",isubmesh );
        }
        else if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
            pumi_locate_function[isubmesh] = &pumi_locate_in_rightBL;
            printf("submesh=%d -- rightBL routine initialized\n",isubmesh );
        }
        else if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
            pumi_locate_function[isubmesh] = &pumi_locate_in_uni;
            printf("submesh=%d -- uniform routine initialized\n",isubmesh );
        }
        else{
            printf("Error in meshtype for submesh %d \n", isubmesh);
            exit(0);
        }
    }
}*/

/*!
* \brief Deallocates/Frees the memory allocated to the funciton pointer *pumi_locate_function
*/
/*
void pumi_finalize_locate_functions(){
    free(pumi_locate_function);
}
*/
/*
* \brief Computes and returns the covolume for a given node in the mesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] node number
*/
double pumi_return_covolume_1D(pumi_mesh_t* pumi_mesh, int inode){
  double covolume;
  if (inode == 0){
    covolume = pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_right_offset)/2.0;
  }
  else if (inode == pumi_mesh->pumi_Nel_total){
    covolume = pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_left_offset)/2.0;
  }
  else if (inode > 0 && inode < pumi_mesh->pumi_Nel_total){
    covolume = pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_left_offset)/2.0 + pumi_return_elemsize(pumi_mesh, inode, pumi_elem_on_right_offset)/2.0;
  }
  else{
    printf("\tInvalid node number for covolume\n");
    exit(0);
  }
  return covolume;
}

/*
* \brief Computes and returns the covolume for a given node in the mesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] node number
*/
double pumi_return_covolume_2D(pumi_mesh_t* pumi_mesh, int inode){
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
}

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) and returns covolume for given node
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param node number
*/
double pumi_return_covolume(pumi_mesh_t* pumi_mesh, int inode){
    return (pumi_covolume_fnptr[pumi_mesh->ndim-1](pumi_mesh, inode));
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
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
}

/*
* \brief Computes and stores BL element size in submesh stuct member
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_ON_1D(pumi_mesh_t *pumi_mesh){
    int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize_calc_flag) = 1;
      int left_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel;
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize = (double*) malloc(left_Nel*sizeof(double));
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords = (double*) malloc((left_Nel+1)*sizeof(double));
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + 0) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + 0) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left;
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + 1) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left+((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
      int iCell;
      for (iCell=1; iCell<left_Nel; iCell++){
        *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + iCell) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + (iCell-1))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r;
        *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + (iCell+1)) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + (iCell)) + *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + iCell);
      }
    }
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize_calc_flag) = 1;
      int right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize = (double*) malloc(right_Nel*sizeof(double));
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords = (double*) malloc((right_Nel+1)*sizeof(double));
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + 0) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-1);
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + 0) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_x_left;
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + 1) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + 0) + *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + 0);
      int iCell;
      for (iCell=1; iCell<right_Nel; iCell++){
        *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + iCell) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + (iCell-1))/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r;
        *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + (iCell+1)) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + (iCell)) + *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + iCell);
      }
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
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
}

/*
* \brief Frees allocated memory to compute 1D BL elem sizes for all submeshes
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
void pumi_BL_elemsize_OFF_1D(pumi_mesh_t *pumi_mesh){
    pumi_mesh->BL_elem_coords_cache_flag = 0;
    int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize_calc_flag = 0;
      free(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize);
      free(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords);
    }
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize_calc_flag = 0;
      free(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize);
      free(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords);
    }
  }
}

/*
* \brief Returns grading ratio for given node
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param node number
*/
double pumi_return_1D_gradingratio(pumi_mesh_t *pumi_mesh, int node){

  if (node == 0 || node == pumi_mesh->pumi_Nel_total){
    printf("Grading ratio not defined for the first and last node of the domain -- Terminating \n");
    exit(0);
  }
  else{
      int isubmesh;
    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){

      int submesh_left_node = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
      int submesh_right_node = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel;

      if (node >= (submesh_left_node + 1) && node <= (submesh_right_node - 1)){

        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
          int submesh_left_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel;
          int submesh_right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
          if (node >= submesh_left_Nel+submesh_left_node && node <= submesh_right_node-submesh_right_Nel){
            return 1.0;
            break;
          }
        }
        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
          int submesh_left_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel;
          if (node >= submesh_left_node+1 && node <= submesh_left_node+submesh_left_Nel-1){
            return ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r;
            break;
          }
        }
        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
          int submesh_right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
          if (node <= submesh_right_node-1 && node >= submesh_right_node-submesh_right_Nel+1){
            return ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r;
            break;
          }
        }


      }

      else if (node == ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative){
        double left_elem_size;
        double right_elem_size;
        // On LHS of the node
        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->pumi_flag & rightBL){
          left_elem_size = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->rBL_t0;
        }
        else{
          if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->pumi_flag & uniform){
            left_elem_size = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->uniform_t0;
          }
          else{
            double t0 = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->lBL_t0;
            double r = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->left_r;
            double Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->left_Nel;
            left_elem_size = t0*pow(r,Nel-1);
          }
        }
        // On RHS of the node
        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->pumi_flag & leftBL){
          right_elem_size = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->lBL_t0;
        }
        else{
          if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->pumi_flag & uniform){
            right_elem_size = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->uniform_t0;
          }
          else{
            double t0 = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->rBL_t0;
            double r = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->right_r;
            double Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->right_Nel;
            right_elem_size = t0*pow(r,Nel-1);
          }
        }
        return right_elem_size/left_elem_size;
        break;
      }
    }
  }

}

/*
* \brief Returns grading ratio for given node
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param node number
*/
double pumi_return_2D_gradingratio(pumi_mesh_t *pumi_mesh, int node){
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
    return 0.0;
}

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) and returns grading ratio for given node
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param node number
*/
double pumi_return_gradingratio(pumi_mesh_t *pumi_mesh, int node){
  return (pumi_gradingratio_fnptr[pumi_mesh->ndim-1](pumi_mesh, node));
}

/*
* \brief Returns element size for given index (and offset)
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param index - node or element index
* \param offset - If node index is inpute, index will be 0 or -1 depending on right/left elementsize to the node to be returned
*/
double pumi_return_1D_elemsize(pumi_mesh_t *pumi_mesh, int index, int offset){
  int elem = index+offset;

  if (elem < 0){
    elem = 0;
  }
  if (elem > pumi_mesh->pumi_Nel_total){
    elem = pumi_mesh->pumi_Nel_total;
  }
  int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    int submesh_left_elem = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
    int submesh_right_elem = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel-1;

    if (elem >= submesh_left_elem && elem <= submesh_right_elem){

      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
        int submesh_left_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel;
        int submesh_right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
        if (elem >= submesh_left_Nel+submesh_left_elem && elem <= submesh_right_elem-submesh_right_Nel){
          return ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
          break;
        }
      }
      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
        int submesh_left_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel;
        if (elem >= submesh_left_elem && elem <= submesh_left_elem+submesh_left_Nel-1){
          int local_elem = elem - submesh_left_elem;
          if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize_calc_flag){
            return *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + local_elem);
            break;
          }
          else {
            double t0 = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->lBL_t0;
            double r = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->left_r;
            return t0*pow(r,local_elem);
            break;
          }
        }
      }
      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
        int submesh_right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
        if (elem <= submesh_right_elem && elem >= submesh_right_elem-submesh_right_Nel+1){
          int local_elem = submesh_right_Nel-(submesh_right_elem-elem+1);
          if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize_calc_flag){
            return *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + local_elem);
            break;
          }
          else {
            double t0 = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->rBL_t0;
            double r = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->right_r;
            return t0*pow(r,submesh_right_Nel-1-local_elem);
            break;
          }
        }
      }

    }
  }
}


/*
* \brief Returns element size for given index (and offset)
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param index - node or element index
* \param offset - If node index is inpute, index will be 0 or -1 depending on right/left elementsize to the node to be returned
*/
double pumi_return_2D_elemsize(pumi_mesh_t *pumi_mesh, int index, int offset){
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
}

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) and returns element size for given index (and offset)
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param index - node or element index
* \param offset - If node index is inpute, index will be 0 or -1 depending on right/left elementsize to the node to be returned
*/
double pumi_return_elemsize(pumi_mesh_t *pumi_mesh, int index, int offset){
    return (pumi_elemsize_fnptr[pumi_mesh->ndim-1](pumi_mesh, index, offset));
}

void pumi_initialize_multiD_functions(pumi_mesh_t *pumi_mesh){

    pumi_gradingratio_fnptr[0] = &pumi_return_1D_gradingratio;
    pumi_gradingratio_fnptr[1] = &pumi_return_2D_gradingratio;

    pumi_elemsize_fnptr[0] = &pumi_return_1D_elemsize;
    pumi_elemsize_fnptr[1] = &pumi_return_2D_elemsize;

    pumi_covolume_fnptr[0] = &pumi_return_covolume_1D;
    pumi_covolume_fnptr[1] = &pumi_return_covolume_2D;
}

/*
* \brief Returns smallest element size in the mesh
* \param *pumi_mesh pointer object to struct pumi_mesh
*/
double pumi_return_smallest_elemsize(pumi_mesh_t *pumi_mesh){
  double smallest_elemsize;
  int isubmesh = 0;

  if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
    smallest_elemsize = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->lBL_t0;
  }
  else{
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      smallest_elemsize = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->rBL_t0;
    }
    else{
      smallest_elemsize = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    }
  }
  for (isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    double new_smallest_elemsize;
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      new_smallest_elemsize = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->lBL_t0;
    }
    else{
      if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
        new_smallest_elemsize = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh))->rBL_t0;
      }
      else{
        new_smallest_elemsize = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
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
void pumi_locate_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int *submeshID, int *cellID){
    if (pumi_mesh->nsubmeshes == 1){
        *submeshID = 0;
    }
    else{
        int isubmesh;
        int submesh_located = 0;
        for (isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
            if (coords < ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Length_cumulative){
                *submeshID = isubmesh-1;
                submesh_located++;
                break;
            }
        }
        if (!(submesh_located)){
            *submeshID = pumi_mesh->nsubmeshes-1;
        }
    }

    *cellID = pumi_locatecell_fnptr[*submeshID](pumi_mesh, *submeshID, coords);

}

/*
* \brief Calculates new submesh ID with adjacency search and new local cell ID of a pushed particle analytically
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param coords - new coordinate of the particle
* \param prev_submeshID - old submesh ID of the pushed particle (i.e. before the push)
* \param[out] pointers to new submesh ID and new local cell which will be populated inside the routine
*/
void pumi_update_submesh_and_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int *submeshID, int *cellID){
    if (pumi_mesh->nsubmeshes == 1){
        *submeshID = 0;
    }
    else{
        *submeshID = prev_submeshID;
        while(coords<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + *submeshID)->x_left){
            *submeshID -= 1;
        }

        while(coords>((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + *submeshID)->x_right){
            *submeshID += 1;
        }
    }
    *cellID = pumi_locatecell_fnptr[*submeshID](pumi_mesh, *submeshID, coords);
}

/*
* \brief Calculates new submesh ID and new local cell ID of a pushed particle using adjacency searches
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param coords - new coordinate of the particle
* \param prev_submeshID - old submesh ID of the pushed particle (i.e. before the push)
* \param prev_cellID - old local cell ID of the pushed particle (i.e. before the push)
* \param[out] pointers to new submesh ID and new local cell which will be populated inside the routine
*/
void pumi_update_submesh_and_update_cell(pumi_mesh_t *pumi_mesh, double coords, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    if (pumi_mesh->nsubmeshes == 1){
        *submeshID = 0;
    }
    else{
        *submeshID = prev_submeshID;
        while(coords<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + *submeshID)->x_left){
            *submeshID -= 1;
            prev_cellID = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + *submeshID)->submesh_total_Nel - 1;
        }

        while(coords>((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + *submeshID)->x_right){
            *submeshID += 1;
            prev_cellID = 0;
        }
    }
    *cellID = pumi_updatecell_fnptr[*submeshID](pumi_mesh, *submeshID, prev_cellID, coords);
}

/*
* \brief subroutine to call relevant weight-calculate routine using function pointers
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the block
* \param[in] local cell ID in the block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
* \param[out] pointers to global cell and weight2 which will be populated inside the routine
*/
void pumi_calc_weights(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    pumi_calc_weights_fnptr[isubmesh](pumi_mesh, isubmesh, local_cell, coord, global_cell, weight);
}

/*
* \brief subroutine to locate the cell number of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_locatecell_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = (coord - ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    return icell;
}

/*
* \brief subroutine to update the cell number of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] local cell ID of particle in the uniform block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_updatecell_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    icell = (coord - ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    return icell;
}

/*
* \brief subroutine to calculate golbal cell ID and weight contribution of a particle inside a uniform block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the uniform block
* \param[in] local cell ID of particle in the uniform block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
void pumi_calc_weights_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell, double coord, int *global_cell, double* weight){
    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0*local_cell))/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    *global_cell = local_cell + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
}

/*
* \brief subroutine to locate the cell number of a particle inside a leftBL block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the leftBL block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_locatecell_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_left_r;
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
    while(coord < *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + icell)){
        icell -= 1;
    }

    while(coord > *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + (icell+1))){
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
    icell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_left_r;
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
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r,local_cell);
    *global_cell = local_cell + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left + (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio))/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0*r_power_cell);
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
    *global_cell = local_cell + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative;
//    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left + (*(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + local_cell) - ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0)/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r - 1.0)))/(*(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + local_cell));
    *weight = (coord - *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + local_cell))/(*(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + local_cell));
}

/*
* \brief subroutine to locate the cell number of a particle inside a rightBL block
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] coord coordinate of the particle whose cell number and local weights is to be evaluated
*/
int pumi_locatecell_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord){
    int icell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_right_r;
    return ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel - icell - 1;
}

/*
* \brief subroutine to update the local cell number of a particle inside a rightBL block with adjacency search
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] local cell ID in the rightBL block
* \param[in] coord coordinate of the particle nodal weight is to be evaluated
*/
int pumi_updatecell_in_rightBL_cached(pumi_mesh_t *pumi_mesh, int isubmesh, int icell, double coord){
    while(coord < *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + icell)){
        icell -= 1;
    }

    while(coord > *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + (icell+1))){
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
    icell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_right_r;
    return icell;
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
    local_cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,local_cell);
    *global_cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel - local_cell - 1 ;
    *weight = 1 - ((((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right - (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio) - coord)/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*r_power_cell);
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
//    local_cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel - local_cell - 1;
    *global_cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative + local_cell ;
//    *weight = 1 - ((((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right - (*(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + local_cell) - ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0)/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r-1.0)) - coord)/(*(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + local_cell));
    *weight = (coord - *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + local_cell))/(*(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + local_cell));
}

/*
* \brief subroutine that returns global cell number from submesh ID and local cell ID
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] submesh ID of the rightBL block
* \param[in] icell local cell ID in the rightBL block
*/
int pumi_global_cell_ID(pumi_mesh_t *pumi_mesh, int isubmesh, int local_cell){
    return (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->Nel_cumulative + local_cell);
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
    *left_node = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0*local_cell;
    *right_node = *left_node + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
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
    *left_node = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + local_cell);
    *right_node = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_coords + (local_cell+1));
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
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r,local_cell);
    *left_node = (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left + (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio);
    *right_node = *left_node + r_power_cell*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
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
    *left_node = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + local_cell);
    *right_node = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_coords + (local_cell+1));
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
    local_cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel - local_cell - 1;
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,local_cell);
    *right_node = (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right - (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio);
    *left_node = *right_node - r_power_cell*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0;
}

/*
* \brief Assigns the appropriate subroutine for each submesh block to locate/update the local cell number of a particle in that submesh
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \details Allocates memory to function pointer *pumi_locatecell_fnptr, *pumi_updatecell_fnptr, *pumi_calc_weights_fnptr
and assign relevant locate/update/weight-calculate pumi routines to each submesh during pumi mesh initialization
*/
void pumi_initialize_locatecell_and_calcweights_functions(pumi_mesh_t *pumi_mesh){
    pumi_locatecell_fnptr = malloc(pumi_mesh->nsubmeshes*sizeof(pumi_locatecell_ptr));
    pumi_updatecell_fnptr = malloc(pumi_mesh->nsubmeshes*sizeof(pumi_updatecell_ptr));
    pumi_calc_weights_fnptr = malloc(pumi_mesh->nsubmeshes*sizeof(pumi_calc_weights_ptr));
    pumi_calc_node_coords_fnptr = malloc(pumi_mesh->nsubmeshes*sizeof(pumi_calc_node_coords_ptr));
    int isubmesh;
    for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
        if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
            pumi_locatecell_fnptr[isubmesh] = &pumi_locatecell_in_leftBL;
            printf("submesh=%d -- leftBL locate cell routine initialized\n",isubmesh );
            if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize_calc_flag){
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_leftBL_cached;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_leftBL_cached;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_leftBL_cached;
                printf("submesh=%d -- leftBL calc weight and cell update (with cache) routines initialized\n",isubmesh );
            }
            else{
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_leftBL_analytic;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_leftBL_analytic;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_leftBL_analytic;
                printf("submesh=%d -- leftBL calc weight and cell update (without cache) routines initialized\n",isubmesh );
            }
        }
        else if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
            pumi_locatecell_fnptr[isubmesh] = &pumi_locatecell_in_rightBL;
            printf("submesh=%d -- rightBL locate cell routine initialized\n",isubmesh );
            if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize_calc_flag){
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_rightBL_cached;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_rightBL_cached;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_rightBL_cached;
                printf("submesh=%d -- rightBL calc weights and cell update (with cache) routines initialized\n",isubmesh );
            }
            else{
                pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_rightBL_analytic;
                pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_rightBL_analytic;
                pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_rightBL_analytic;
                printf("submesh=%d -- rightBL calc weight and cell update (without cache) routines initialized\n",isubmesh );
            }
        }
        else if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
            pumi_locatecell_fnptr[isubmesh] = &pumi_locatecell_in_uni;
            pumi_updatecell_fnptr[isubmesh] = &pumi_updatecell_in_uni;
            pumi_calc_weights_fnptr[isubmesh] = &pumi_calc_weights_in_uni;
            pumi_calc_node_coords_fnptr[isubmesh] = &pumi_calc_node_coords_in_uni;
            printf("submesh=%d -- uniform routines initialized\n",isubmesh );
        }
        else{
            printf("Error in meshtype for submesh %d \n", isubmesh);
            exit(0);
        }
    }
}

/*!
* \brief Deallocates/Frees the memory allocated to the funciton pointers pumi_*_fnptr
*/
void pumi_finalize_locatecell_and_calcweights_functions(){
    free(pumi_locatecell_fnptr);
    free(pumi_updatecell_fnptr);
    free(pumi_calc_weights_fnptr);
}
