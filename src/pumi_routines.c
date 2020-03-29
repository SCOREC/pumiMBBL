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
void pumi_locatepoint(pumi_mesh_t *pumi_mesh, double particle_coordinate, int particle_submesh, int *particle_cell, double *cell_weight)
{
  if (pumi_mesh->ndim == 1){
    pumiMBBL_locatepoint_1D(pumi_mesh, particle_coordinate, particle_submesh, particle_cell, cell_weight);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
}

void pumi_locate_in_uni(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){
    *cell = (coord - ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0*(*cell)))/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
}

void pumi_locate_in_leftBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){
    *cell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_left_r;
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r,*cell);
    *weight = (coord - (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left + (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r_lBL_t0_ratio))/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0*r_power_cell);
}

void pumi_locate_in_rightBL(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){
    *cell = log(1 + (fabs(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right-coord))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->log_right_r;
    double r_power_cell = pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,*cell);
    *weight = 1 - ((((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right - (r_power_cell-1.0)/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r_rBL_t0_ratio) - coord)/(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*r_power_cell);
}

void pumi_dummylocate(pumi_mesh_t *pumi_mesh, int isubmesh, double coord, int *cell, double *weight){

}
//                                                                          0x00            0x01                0x02                    0x03                0x04
void (*pumi_locate_arr[])(pumi_mesh_t*, int, double, int*, double*) = {pumi_dummylocate, pumi_locate_in_uni, pumi_locate_in_leftBL, pumi_dummylocate, pumi_locate_in_rightBL};

/*
* \brief Calls appropriate subroutine to obtain the local cell number of the particle inside the segment
* \param[in] *pumi_mesh pointer object to struct pumi_mesh
* \param[in] particle_coordinate coordinate of the particle whose cell number and local weights is to be evaluated
* \param[in] submesh of the particle whose cell number and local weights is to be evaluated
* \param[out] *particle_cell address of the variable where the particle cell is to be stored
* \param[out] *cell_weight address of the variable where the local weight is to be stored
*/
void pumiMBBL_locatepoint_1D(pumi_mesh_t *pumi_mesh, double particle_coordinate, int particle_submesh, int *particle_cell, double *cell_weight){
    int local_cell;
    double local_weight;
    (*pumi_locate_arr[((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + particle_submesh)->pumi_flag])(pumi_mesh, particle_submesh, particle_coordinate, &local_cell, &local_weight);
    *particle_cell = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + particle_submesh)->Nel_cumulative + local_cell;
    *cell_weight = local_weight;
}


/*
* \brief Call appropriate subroutine (based on the dimension of the problem) and returns covolume for given node
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param node number
*/
double pumi_return_covolume(pumi_mesh_t* pumi_mesh, int inode){
  if (pumi_mesh->ndim == 1){
    return pumi_return_covolume_1D(pumi_mesh, inode);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
}


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
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + 0) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
      int iCell;
      for (iCell=1; iCell<left_Nel; iCell++){
        *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + iCell) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + (iCell-1))*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r;
      }
    }
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize_calc_flag) = 1;
      int right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize = (double*) malloc(right_Nel*sizeof(double));
      *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + 0) = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-1);
      int iCell;
      for (iCell=1; iCell<right_Nel; iCell++){
        *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + iCell) = *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + (iCell-1))/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r;
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
    int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      /*
      int left_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel;
      for (int i=0; i<left_Nel; i++){
        printf("Submesh %d: leftBL element %d size is %2.4e \n", isubmesh, i, *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize + i) );
      }*/
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize_calc_flag = 0;
      free(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->leftBL_elemsize);
      // printf("\n");
    }
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      /*
      int right_Nel = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel;
      for (int i=0; i<right_Nel; i++){
        printf("Submesh %d: rightBL element %d size is %2.4e \n", isubmesh, i, *(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize + i) );
      }*/
      ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize_calc_flag = 0;
      free(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rightBL_elemsize);
      // printf("\n");
    }
    // printf("\n\n");
  }
}

/*
* \brief Call appropriate subroutine (based on the dimension of the problem) and returns grading ratio for given node
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param node number
*/
double pumi_return_gradingratio(pumi_mesh_t *pumi_mesh, int node){
  if (pumi_mesh->ndim == 1){
    return pumi_return_1D_gradingratio(pumi_mesh, node);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
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
* \brief Call appropriate subroutine (based on the dimension of the problem) and returns element size for given index (and offset)
* \param *pumi_mesh pointer object to struct pumi_mesh
* \param index - node or element index
* \param offset - If node index is inpute, index will be 0 or -1 depending on right/left elementsize to the node to be returned
*/
double pumi_return_elemsize(pumi_mesh_t *pumi_mesh, int index, int offset){
  if (pumi_mesh->ndim == 1){
    return pumi_return_1D_elemsize(pumi_mesh, index, offset);
  }
  else {
    printf("Multi dimension pumi mesh not implemented -- Terminating\n");
    exit(0);
  }
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
