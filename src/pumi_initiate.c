#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumi_initiate.h"

pumi_mesh_t* pumi_initiate(int PUMI_INITIATE_FLAG, pumi_initiate_input_t pumi_inputs){

  pumi_mesh_t* pumi_mesh = (pumi_mesh_t*) malloc(sizeof(pumi_mesh_t));

  if (PUMI_INITIATE_FLAG == 0){
    int dimension;
    int submesh_num;
    double **submesh_params;
    unsigned int *submesh_flag;
    pumi_getmeshparameters_from_terminal(&dimension,&submesh_num,&submesh_params,&submesh_flag);
    pumi_mesh->ndim = dimension;
    pumi_mesh->nsubmeshes = submesh_num;
    if (pumi_mesh->ndim == 1){
      pumi_mesh->pumi_submeshes = (void*) malloc(pumi_mesh->nsubmeshes * sizeof(pumi_submesh1D_t));
    }
    for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
      pumi_setsubmesh(pumi_mesh, isubmesh, submesh_params[isubmesh][0], submesh_params[isubmesh][1], submesh_flag[isubmesh], (int) submesh_params[isubmesh][2], submesh_params[isubmesh][3], submesh_params[isubmesh][4], (int) submesh_params[isubmesh][5], submesh_params[isubmesh][6], submesh_params[isubmesh][7], (int) submesh_params[isubmesh][8]);
    }
  }
  else if (PUMI_INITIATE_FLAG == 1){
    //only use this if pumi_initiate_input struct members are populated accordingly
    pumi_mesh->ndim = pumi_inputs.ndim;
    pumi_mesh->nsubmeshes = pumi_inputs.nsubmeshes;
    if (pumi_mesh->ndim == 1){
      pumi_mesh->pumi_submeshes = (void*) malloc(pumi_mesh->nsubmeshes * sizeof(pumi_submesh1D_t));
    }
    for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
      char flagstring[30];
      strcpy(flagstring, pumi_inputs.type_flag[isubmesh]);
      int n = strlen(flagstring);
      flagstring[n]='L'; // add this dummy letter to the string because the string is one charecter smaller than usual i.e. pumi_getmeshparameters_from_terminal  
      unsigned int submesh_flag = pumi_getsubmeshflag(flagstring);
      pumi_setsubmesh(pumi_mesh, isubmesh, *(pumi_inputs.x_left + isubmesh), *(pumi_inputs.x_right + isubmesh), submesh_flag, *(pumi_inputs.uniform_Nel + isubmesh), *(pumi_inputs.left_T + isubmesh), *(pumi_inputs.left_r + isubmesh), *(pumi_inputs.left_Nel + isubmesh), *(pumi_inputs.right_T + isubmesh), *(pumi_inputs.right_r + isubmesh), *(pumi_inputs.right_Nel + isubmesh));
    }
  }

  return (pumi_mesh);
}

void pumi_setsubmesh(pumi_mesh_t *pumi_mesh, int isubmesh, double xleft, double xright, unsigned int submeshflag, int N_uniform, double T_left, double r_left, int N_left, double T_right, double r_right, int N_right){
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag = submeshflag;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left = xleft;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right = xright;

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel = N_uniform;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left = xleft + T_left; // (dependent variable)
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_right = xright - T_right; // (dependent variable)
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0 = ((xright - T_right) - (xleft + T_left))/N_uniform; // (dependent variable)

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_T = T_left;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r = r_left;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel = N_left;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_x_right = xleft + T_left; // (dependent variable)
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0 = T_left*(r_left-1)/(pow(r_left,N_left)-1); // (dependent variable)

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_T = T_right;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r = r_right;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel = N_right;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_x_left = xright - T_right;// (dependent variable)
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0 = T_right*(r_right-1)/(pow(r_right,N_right)-1);// (dependent variable)

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel = N_uniform + N_left + N_right; // (dependent variable)
}

void pumi_getmeshparameters_from_terminal (int *dimension, int *submesh_num, double ***submesh_params, unsigned int **submesh_flag){
  printf("\nEnter the dimension of the domain: ");
  scanf("%d", dimension);
  printf("\nEnter the number of submeshes available: ");
  scanf("%d", submesh_num); //user supplied
  int m = (int) *submesh_num;
  double **tmp_param = (double **) malloc( m * sizeof(**tmp_param));
  unsigned int *tmp_flag = malloc (m * sizeof(*tmp_flag));

  for (int isubmesh=0; isubmesh<m; isubmesh++){
      tmp_param[isubmesh] = (double*) malloc(9*sizeof(double));
      char flagstring[30];
      for (int j=2; j<9; j++){
        tmp_param[isubmesh][j] = 0.0; // default value for all struct members
      }
      printf("\n\nSUBMESH %d :",isubmesh+1 );
      printf("\nCoordinates of left end point of the submesh (%d) : ",isubmesh+1);
      scanf("%lf", &tmp_param[isubmesh][0]); //user supplied
      printf("\nCoordinates of right end point of the submesh (%d) : ",isubmesh+1);
      scanf("%lf", &tmp_param[isubmesh][1]); //user supplied
      printf("\nEnter the active flags (uniform, leftBL, rightBL -- separated by '&') for submesh (%d) : ",isubmesh+1);
      fgetc(stdin);
      fgets(flagstring, sizeof flagstring, stdin);
      tmp_flag[isubmesh] = pumi_getsubmeshflag(flagstring); //user supplied

      if (tmp_flag[isubmesh] & uniform){
        printf("\nNumber of elements in uniform mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][2]); //user supplied
      }
      if (tmp_flag[isubmesh] & leftBL){
        printf("\nThickness of (left) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][3]); //user supplied
        printf("\nGrading ratio in the (left) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][4]); //user supplied
        printf("\nNumber of elements in the (left) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][5]); //user supplied
      }
      if (tmp_flag[isubmesh] & rightBL){
        printf("\nThickness of (right) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][6]); //user supplied
        printf("\nGrading ratio in the (right) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][7]); //user supplied
        printf("\nNumber of elements in the (right) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][8]); //user supplied
      }
  }
  *submesh_flag = tmp_flag;
  *submesh_params = tmp_param;

  // deallocate_tmp_flag(tmp_flag);
  // deallocate_tmp_param(m, tmp_param);

}

unsigned int pumi_getsubmeshflag(char flagstring[30]){

  int n = strlen(flagstring);
  flagstring[n-1]='&';

  char newflagstring[3][10];
  char flagtypes[3][10] = {"uniform","leftBL","rightBL"};
  int j=0;
  int numstring=0;
  for (int i=0;i<(strlen(flagstring));i++){
    if (flagstring[i]=='&'){
      newflagstring[numstring][j] = '\0';
      numstring++;
      j=0;
    }
    else{
      newflagstring[numstring][j] = flagstring[i];
      j++;
    }
  }

  unsigned int intflag = 0x00;
  for (int i=0; i<numstring; i++){
    int l1 = strcmp(newflagstring[i],flagtypes[0]);
    int l2 = strcmp(newflagstring[i],flagtypes[1]);
    int l3 = strcmp(newflagstring[i],flagtypes[2]);
    if (l1*l2*l3 != 0){
      printf("Invalid flag input -- Terminating\n");
      printf("Valid input(s): leftBL&uniform or rightBL&uniform or leftBL&rightBL or uniform&leftBL&rightBL\n");
      exit(0);
    }
    if (strcmp(newflagstring[i],flagtypes[0])==0){
      intflag =  intflag + 0x01;
    }
    if (strcmp(newflagstring[i],flagtypes[1])==0){
      intflag = intflag + 0x02;
    }
    if (strcmp(newflagstring[i],flagtypes[2])==0){
      intflag = intflag + 0x04;
    }
  }
    return intflag;
}
/*
void deallocate_tmp_flag(unsigned int *tmp_flag){
  free(tmp_flag);
}

void deallocate_tmp_param(int dim, double **tmp_param){
  int i;
  for (i=0; i<dim; i++){
      free(tmp_param[i]);
  }
  free(tmp_param);
}
*/
