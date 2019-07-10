#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumi_initiate.h"

pumi_mesh_t* pumi_initiate(pumi_initiate_flag_t pumi_input_initiate_flag, pumi_initiate_input_t *pumi_inputs){

  pumi_mesh_t* pumi_mesh = (pumi_mesh_t*) malloc(sizeof(pumi_mesh_t));

  if (pumi_input_initiate_flag == initiate_from_terminal){
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
    else{
      printf("Multi dimension pumi mesh not implemented -- Terminating\n");
      exit(0);
    }
    for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
      pumi_setsubmesh(pumi_mesh, isubmesh, submesh_params[isubmesh][0], submesh_params[isubmesh][1], submesh_flag[isubmesh], (int) submesh_params[isubmesh][2], submesh_params[isubmesh][3], submesh_params[isubmesh][4], (int) submesh_params[isubmesh][5], submesh_params[isubmesh][6], submesh_params[isubmesh][7], (int) submesh_params[isubmesh][8]);
    }
    pumi_freemeshparameters_from_terminal(pumi_mesh->nsubmeshes, submesh_params, submesh_flag);

  }
  else if (pumi_input_initiate_flag == initiate_from_commandline_inputs){
    //only use this if pumi_initiate_input struct members are populated accordingly
    pumi_mesh->ndim = pumi_inputs->ndim;
    pumi_mesh->nsubmeshes = pumi_inputs->nsubmeshes;
    if (pumi_mesh->ndim == 1){
      pumi_mesh->pumi_submeshes = (void*) malloc(pumi_mesh->nsubmeshes * sizeof(pumi_submesh1D_t));
    }
    else{
      printf("Multi dimension pumi mesh not implemented -- Terminating\n");
      exit(0);
    }
    for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
      char flagstring[SUBMESH_FLAGSTRING_LENGTH];
      strcpy(flagstring, pumi_inputs->type_flag[isubmesh]);
      unsigned int submesh_flag = pumi_getsubmeshflag(flagstring);
      pumi_setsubmesh(pumi_mesh, isubmesh, *(pumi_inputs->x_left + isubmesh), *(pumi_inputs->x_right + isubmesh), submesh_flag, *(pumi_inputs->uniform_Nel + isubmesh), *(pumi_inputs->left_T + isubmesh), *(pumi_inputs->left_r + isubmesh), *(pumi_inputs->left_Nel + isubmesh), *(pumi_inputs->right_T + isubmesh), *(pumi_inputs->right_r + isubmesh), *(pumi_inputs->right_Nel + isubmesh));
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
      char flagstring[SUBMESH_FLAGSTRING_LENGTH];

      printf("\n\nSUBMESH %d :",isubmesh+1 );
      printf("\nCoordinates of left end point of the submesh (%d) : ",isubmesh+1);
      scanf("%lf", &tmp_param[isubmesh][0]); //user supplied
      printf("\nCoordinates of right end point of the submesh (%d) : ",isubmesh+1);
      scanf("%lf", &tmp_param[isubmesh][1]); //user supplied
      printf("\nEnter the active flags (uniform, leftBL, rightBL -- separated by '&') for submesh (%d) : ",isubmesh+1);
      scanf("%s", flagstring);
      tmp_flag[isubmesh] = pumi_getsubmeshflag(flagstring); //user supplied

      if (tmp_flag[isubmesh] & uniform){
        printf("\nNumber of elements in uniform mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][2]); //user supplied
      }
      else {
        tmp_param[isubmesh][2] = DEFAULT_PARAM_VAL;
      }
      if (tmp_flag[isubmesh] & leftBL){
        printf("\nThickness of (left) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][3]); //user supplied
        printf("\nGrading ratio in the (left) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][4]); //user supplied
        printf("\nNumber of elements in the (left) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][5]); //user supplied
      }
      else {
        tmp_param[isubmesh][3] = DEFAULT_PARAM_VAL;
        tmp_param[isubmesh][4] = DEFAULT_PARAM_VAL;
        tmp_param[isubmesh][5] = DEFAULT_PARAM_VAL;
      }
      if (tmp_flag[isubmesh] & rightBL){
        printf("\nThickness of (right) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][6]); //user supplied
        printf("\nGrading ratio in the (right) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][7]); //user supplied
        printf("\nNumber of elements in the (right) graded mesh region: ");
        scanf("%lf", &tmp_param[isubmesh][8]); //user supplied
      }
      else {
        tmp_param[isubmesh][6] = DEFAULT_PARAM_VAL;
        tmp_param[isubmesh][7] = DEFAULT_PARAM_VAL;
        tmp_param[isubmesh][8] = DEFAULT_PARAM_VAL;
      }
  }
  *submesh_flag = tmp_flag;
  *submesh_params = tmp_param;
}

unsigned int pumi_getsubmeshflag(char flagstring[SUBMESH_FLAGSTRING_LENGTH]){
  char newflagstring[SUBMESH_MAX_SEGMENTS][SEGMENT_STRING_LENGTH];
  char flagtypes[SUBMESH_MAX_SEGMENTS][SEGMENT_STRING_LENGTH] = {"uniform","leftBL","rightBL"};

  char *tok = strtok(flagstring, "&");
  int numstring=0;
  while (tok != NULL){
    strcpy (newflagstring[numstring], tok);
    tok = strtok(NULL, "&");
    numstring++;
  }

  unsigned int intflag = unassigned;
  for (int i=0; i<numstring; i++){
    int l1 = strcmp(newflagstring[i],flagtypes[0]);
    int l2 = strcmp(newflagstring[i],flagtypes[1]);
    int l3 = strcmp(newflagstring[i],flagtypes[2]);
    int l123 = l1*l2*l3;
    if (l123 != 0){
      printf("Invalid flag input -- Terminating\n");
      printf("Valid input(s):\nuniform\nleftBL\nrightBL\nleftBL&uniform\nrightBL&uniform\nleftBL&rightBL\nuniform&leftBL&rightBL\n");
      exit(0);
    }
    if (strcmp(newflagstring[i],flagtypes[0])==0){
      intflag =  intflag + uniform;
    }
    if (strcmp(newflagstring[i],flagtypes[1])==0){
      intflag = intflag + leftBL;
    }
    if (strcmp(newflagstring[i],flagtypes[2])==0){
      intflag = intflag + rightBL;
    }
  }
    return intflag;
}

void pumi_initiate_allocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes){
  pumi_inputs->x_left = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->x_right = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->type_flag = (char**) malloc(nsubmeshes*sizeof(char*));
  for (int i=0; i<nsubmeshes; i++){
    pumi_inputs->type_flag[i] = (char*) malloc(SUBMESH_FLAGSTRING_LENGTH);
  }
  pumi_inputs->left_T = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->left_r = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->left_Nel = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->right_T = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->right_r = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->right_Nel = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->uniform_Nel = malloc(nsubmeshes*sizeof(int));
}

void pumi_initiate_deallocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes){
  free(pumi_inputs->x_left);
  free(pumi_inputs->x_right);
  free(pumi_inputs->left_T);
  free(pumi_inputs->left_r);
  free(pumi_inputs->left_Nel);
  free(pumi_inputs->right_T);
  free(pumi_inputs->right_r);
  free(pumi_inputs->right_Nel);
  free(pumi_inputs->uniform_Nel);
  for (int i=0; i<nsubmeshes; i++){
    free(pumi_inputs->type_flag[i]);
  }
  free(pumi_inputs->type_flag);
  free(pumi_inputs);
}

void pumi_finalize(pumi_mesh_t* pumi_mesh){
  free(pumi_mesh->pumi_submeshes);
  free(pumi_mesh);
}

void pumi_freemeshparameters_from_terminal(int nsubmeshes, double **submesh_params, unsigned int *submesh_flag){
  free(submesh_flag);
  for (int i=0; i<nsubmeshes; i++){
      free(submesh_params[i]);
  }
  free(submesh_params);
}
