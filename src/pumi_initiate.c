#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumi_initiate.h"

/*!
* \brief Initiates the pumi mesh and populates the members of the structures pumi_mesh and pumi_submesh1D and returns the object to pumi_mesh
* \param pumi_input_initiate_flag initiate flag enum that specifies how to read the mesh parameters from user
* \param *pumi_inputs pointer object to pumi_initiate_input structure
* \details If initiate flag is set to initiate_from_terminal the members of pumi_initiate_input need not be assigned any value. The user will be prompted to supply the inputs in the terminal.
If the flag is set to initiate_from_commandline_inputs, the members of pumi_initiate_input has to populated by the user (based on command line arguments of hpic main code) before this function call.

*/
pumi_mesh_t* pumi_initiate(pumi_initiate_flag_t pumi_input_initiate_flag, pumi_initiate_input_t *pumi_inputs, int BL_caching_flag){

  pumi_mesh_t* pumi_mesh = (pumi_mesh_t*) malloc(sizeof(pumi_mesh_t));

  if (pumi_input_initiate_flag == initiate_from_terminal){
    int dimension;
    int submesh_num;
    double **submesh_params;
    unsigned int *submesh_flag;
    pumi_getmeshparameters_from_terminal(&dimension,&submesh_num,&submesh_params,&submesh_flag);
    pumi_mesh->ndim = dimension;
    pumi_mesh->nsubmeshes_x = submesh_num;
    if (pumi_mesh->ndim == 1){
      pumi_mesh->pumi_submeshes_x = (void*) malloc(pumi_mesh->nsubmeshes_x * sizeof(pumi_submesh_t));
    }
    else{
      printf("Multi dimension pumi mesh not implemented -- Terminating\n");
      exit(0);
    }
    int isubmesh;
    for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x; isubmesh++){
      pumi_setsubmesh_x(pumi_mesh, isubmesh, submesh_params[isubmesh][0], submesh_params[isubmesh][1], submesh_flag[isubmesh], (int) submesh_params[isubmesh][2], submesh_params[isubmesh][3], submesh_params[isubmesh][4], (int) submesh_params[isubmesh][5], submesh_params[isubmesh][6], submesh_params[isubmesh][7], (int) submesh_params[isubmesh][8]);
    }
    pumi_freemeshparameters_from_terminal(pumi_mesh->nsubmeshes_x, submesh_params, submesh_flag);

  }
  else if (pumi_input_initiate_flag == initiate_from_commandline_inputs){
    //only use this if pumi_initiate_input struct members are populated accordingly
    pumi_mesh->ndim = pumi_inputs->ndim;

    if (pumi_mesh->ndim == 1){
      pumi_mesh->nsubmeshes_x = pumi_inputs->nsubmeshes;
      pumi_mesh->pumi_Nel_total_x = 0;
      pumi_mesh->pumi_submeshes_x = (void*) malloc(pumi_mesh->nsubmeshes_x * sizeof(pumi_submesh_t));
      int isubmesh;
      for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x; isubmesh++){
        char flagstring[SUBMESH_FLAGSTRING_LENGTH];
        strcpy(flagstring, pumi_inputs->type_flag[isubmesh]);
        unsigned int submesh_flag = pumi_getsubmeshflag(flagstring);
        pumi_setsubmesh_x(pumi_mesh, isubmesh, *(pumi_inputs->x_left + isubmesh), *(pumi_inputs->x_right + isubmesh), submesh_flag, *(pumi_inputs->uniform_Nel + isubmesh), *(pumi_inputs->left_T + isubmesh), *(pumi_inputs->left_r + isubmesh), *(pumi_inputs->left_Nel + isubmesh), *(pumi_inputs->right_T + isubmesh), *(pumi_inputs->right_r + isubmesh), *(pumi_inputs->right_Nel + isubmesh));
      }
      printf("PUMI mesh parameter info :\n\n");
      printf("\tTotal elements in mesh = %d\n\n", pumi_mesh->pumi_Nel_total_x);

      for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x; isubmesh++){
        printf("\tSUBMESH %d parameters:\n", isubmesh+1 );

        printf("\n\t submeshflag = ");
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & leftBL){
          printf("leftBL\n");
          printf("\t left_t0     = %2.4e \t [m] Cell size of first/leftmost cell in left BL segment\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0);
          printf("\t left_T      = %2.4e \t [m] Left boundary layer (left BL) thickness\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T);
          printf("\t left_r      = %2.4e \t Grading ratio in left BL mesh\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r);
          printf("\t left_Nel    = %d    \t\t Number of Cells in left BL mesh region\n\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel);
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & uniform){
          printf("uniform\n");
          printf("\t uniform_Nel = %d    \t\t Number of Cells in uniform mesh region\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel);
          printf("\t uniform_dx  = %2.4e \t [m] Cell size in uniform mesh segment\n\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0);
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & rightBL){
          printf("rightBL\n");
          printf("\t right_t0    = %2.4e \t [m] Cell size of last/rightmost cell in right BL segment\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0);
          printf("\t right_T     = %2.4e \t [m] Right boundary layer (right BL) thickness\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T);
          printf("\t right_r     = %2.4e \t Grading ratio in right BL mesh\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r);
          printf("\t right_Nel   = %d    \t\t Number of Cells in right BL mesh region\n\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel);
        }
      }
      pumi_initialize_multiD_functions(pumi_mesh);
      //pumi_initialize_locate_functions(pumi_mesh);
      pumi_initialize_locatecell_and_calcweights_functions(pumi_mesh);
    }
    else{
      printf("Multi dimension pumi mesh not implemented -- Terminating\n");
      exit(0);
    }

  }

  pumi_verify_params(pumi_mesh);
  pumi_print_node_coordinates(pumi_mesh);
  if (BL_caching_flag){
      pumi_mesh->BL_elem_coords_cache_flag = 1;
      pumi_BL_elemsize_ON(pumi_mesh);
  }
  else{
      pumi_mesh->BL_elem_coords_cache_flag = 0;
  }

  return (pumi_mesh);
}

/*!
* \brief Assigns the values to x-members of struct pumi_submesh based on user inputs
* \param[out] *pumi_mesh pointer object of struct pumi_mesh whose members are assigned their corresponding values in this routine
* \param[in] isubmesh The index number of the submesh whose members (i.e submesh parameters) are being assigned their values
* \param[in] xmin Left end coordinate of the submesh block
* \param[in] xmax Right end coordinate of the submesh block
* \param[in] submeshflag Mesh type flag of submesh block converted into a unsigned int using pumi_getsubmeshflag
* \param[in] N_uniform Number of elements in the uniform mesh segment of the submesh block
* \param[in] T_minBL Thickness of the left BL segment of the submesh block
* \param[in] r_minBL Growth ratio of elements in the left BL segment of the submesh block
* \param[in] N_minBL Number of elements in the left BL segment of the submesh block
* \param[in] T_maxBL Thickness of the right BL segment of the submesh block
* \param[in] r_maxBL Growth ratio of elements in the right BL segment of the submesh block
* \param[in] N_maxBL Number of elements in the right BL segment of the submesh block
* \details The submesh parameters (used to define the submesh along with the dependent submesh paramaters) defined in the struct pumi_submesh1D are assigned their corresponding
values based on inputs from command line of hpic code or based on values supplied by the user directly from terminal prompt
*/
void pumi_setsubmesh_x(pumi_mesh_t *pumi_mesh, int isubmesh, double xmin, double xmax, unsigned int submeshflag, int N_uniform, double T_minBL, double r_minBL, int N_minBL, double T_maxBL, double r_maxBL, int N_maxBL){
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag = submeshflag;
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->x_min = xmin;
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->x_max = xmax;
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T = xmax-xmin;
  if (submeshflag & uniform){
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel = N_uniform;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0 = (xmax-xmin)/N_uniform;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r = 1.0;
  }

  if (submeshflag & leftBL){
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel = N_minBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0 = T_minBL*(r_minBL-1.0)/(pow(r_minBL,N_minBL)-1.0);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r = r_minBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->log_r = log(r_minBL);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r_t0_ratio = (r_minBL-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0;
  }

  if (submeshflag & rightBL){
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel = N_maxBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0 = T_maxBL*(r_maxBL-1)/(pow(r_maxBL,N_maxBL)-1.0);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r = r_maxBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->log_r = log(r_maxBL);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r_t0_ratio = (r_maxBL-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0;
  }

  if (isubmesh==0){
    ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->Nel_cumulative = 0;
  }
  else{
    ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->Nel_cumulative = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + (isubmesh-1))->Nel_cumulative + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + (isubmesh-1))->submesh_Nel;
  }
  pumi_mesh->pumi_Nel_total_x += ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel;
}

/*!
* \brief Assigns the values to y-members of struct pumi_submesh based on user inputs
* \param[out] *pumi_mesh pointer object of struct pumi_mesh whose members are assigned their corresponding values in this routine
* \param[in] isubmesh The index number of the submesh whose members (i.e submesh parameters) are being assigned their values
* \param[in] xmin Left end coordinate of the submesh block
* \param[in] xmax Right end coordinate of the submesh block
* \param[in] submeshflag Mesh type flag of submesh block converted into a unsigned int using pumi_getsubmeshflag
* \param[in] N_uniform Number of elements in the uniform mesh segment of the submesh block
* \param[in] T_minBL Thickness of the left BL segment of the submesh block
* \param[in] r_minBL Growth ratio of elements in the left BL segment of the submesh block
* \param[in] N_minBL Number of elements in the left BL segment of the submesh block
* \param[in] T_maxBL Thickness of the right BL segment of the submesh block
* \param[in] r_maxBL Growth ratio of elements in the right BL segment of the submesh block
* \param[in] N_maxBL Number of elements in the right BL segment of the submesh block
* \details The submesh parameters (used to define the submesh along with the dependent submesh paramaters) defined in the struct pumi_submesh1D are assigned their corresponding
values based on inputs from command line of hpic code or based on values supplied by the user directly from terminal prompt
*/
void pumi_setsubmesh_y(pumi_mesh_t *pumi_mesh, int isubmesh, double xmin, double xmax, unsigned int submeshflag, int N_uniform, double T_minBL, double r_minBL, int N_minBL, double T_maxBL, double r_maxBL, int N_maxBL){
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->pumi_flag = submeshflag;
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->x_min = xmin;
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->x_max = xmax;
  ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->submesh_T = xmax-xmin;
  if (submeshflag & uniform){
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->submesh_Nel = N_uniform;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->t0 = (xmax-xmin)/N_uniform;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->r = 1.0;
  }

  if (submeshflag & bottomBL){
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->submesh_Nel = N_minBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->t0 = T_minBL*(r_minBL-1.0)/(pow(r_minBL,N_minBL)-1.0);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->r = r_minBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->log_r = log(r_minBL);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->r_t0_ratio = (r_minBL-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->t0;
  }

  if (submeshflag & topBL){
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->submesh_Nel = N_maxBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->t0 = T_maxBL*(r_maxBL-1)/(pow(r_maxBL,N_maxBL)-1.0);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->r = r_maxBL;
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->log_r = log(r_maxBL);
      ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->r_t0_ratio = (r_maxBL-1.0)/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->t0;
  }

  if (isubmesh==0){
    ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->Nel_cumulative = 0;
  }
  else{
    ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + isubmesh)->Nel_cumulative = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + (isubmesh-1))->Nel_cumulative + ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_y + (isubmesh-1))->submesh_Nel;
  }
  pumi_mesh->pumi_Nel_total_y += ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel;
}

/*!
* \brief Prompts the user to supply values to the parameters that defines the mesh (and its submeshes)
* \param[out] *dimension address of the variable where number of dimensions of the problem space is stored
* \param[out] *submesh_num address of the variable where number of submeshes to be used for the problem is stored
* \param[out] ***submesh_params address of a matrix which stores the values to the submesh parameters for all submesh blocks
* \param[out] **submesh_flag address to the array which stores the mesh type flag (converetd to unsigned int) of each submesh block
* \details Based on the mesh type flag assigned for each submesh the terminal will only prompt the user to supply values relevant to the submesh.
Default values will be assigned to rest of members of pumi_submesh1D. For example, the terminal will not prompt the user to supply values to
the three members defining the right BL segment if the mesh type for the submesh is "uniform&leftBL"
*/
void pumi_getmeshparameters_from_terminal (int *dimension, int *submesh_num, double ***submesh_params, unsigned int **submesh_flag){
  printf("\nEnter the dimension of the domain: ");
  scanf("%d", dimension);
  printf("\nEnter the number of submeshes available: ");
  scanf("%d", submesh_num); //user supplied
  int m = (int) *submesh_num;
  double **tmp_param = (double **) malloc( m * sizeof(**tmp_param));
  unsigned int *tmp_flag = malloc (m * sizeof(*tmp_flag));
  int isubmesh;
  for (isubmesh=0; isubmesh<m; isubmesh++){
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
        tmp_param[isubmesh][2] = (int) DEFAULT_PARAM_VAL;
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
        tmp_param[isubmesh][5] = (int) DEFAULT_PARAM_VAL;
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
        tmp_param[isubmesh][8] = (int) DEFAULT_PARAM_VAL;
      }
  }
  *submesh_flag = tmp_flag;
  *submesh_params = tmp_param;
}

/*!
* \brief Reads the mesh type flag string supplied by the user (in command line inputs or from terminal prompt) and returns a unique
unsigned int number (i.e a bit flag based on enum pumi_meshflag) corresponding to that string
* \param flagstring mesh type flag string supplied by the user (in command line inputs or from terminal prompt)
*/
unsigned int pumi_getsubmeshflag(char flagstring[SUBMESH_FLAGSTRING_LENGTH]){
  char newflagstring[SUBMESH_MAX_SEGMENTS][SEGMENT_STRING_LENGTH];
  char flagtypes[SUBMESH_MAXTYPES][SEGMENT_STRING_LENGTH] = {"uniform","leftBL","rightBL","bottomBL","topBL"};

  char *tok = strtok(flagstring, "&");
  int numstring=0;
  while (tok != NULL){
    strcpy (newflagstring[numstring], tok);
    tok = strtok(NULL, "&");
    numstring++;
  }

  unsigned int intflag = unassigned;
  int i;
  for (i=0; i<numstring; i++){
/*    int l1 = strcmp(newflagstring[i],flagtypes[0]);
    int l2 = strcmp(newflagstring[i],flagtypes[1]);
    int l3 = strcmp(newflagstring[i],flagtypes[2]);
    int l123 = l1*l2*l3;
    if (l123 != 0){
      printf("\n\tInvalid flag input\n");
      printf("\tList of valid input(s):\n\tuniform\n\tleftBL\n\trightBL\n\tleftBL&uniform\n\trightBL&uniform\n\tleftBL&rightBL\n\tuniform&leftBL&rightBL\n\n");
      printf("\tAssigning \"leftBL&uniform&rightBL\" (default input) to submesh flag...\n\n");
      intflag = uniform+leftBL+rightBL;
  }*/
    //else{
      if (strcmp(newflagstring[i],flagtypes[0])==0){
        intflag =  intflag + uniform;
      }
      if (strcmp(newflagstring[i],flagtypes[1])==0 || strcmp(newflagstring[i],flagtypes[3])==0){
        intflag = intflag + leftBL;
      }
      if (strcmp(newflagstring[i],flagtypes[2])==0 || strcmp(newflagstring[i],flagtypes[4])==0){
        intflag = intflag + rightBL;
      }
    //}
  }
    return intflag;
}

/*!
* \brief Allocates the memory to members of struct pumi_initiate_input based on the number of submesh blocks
* \param *pumi_inputs pointer object to struct pumi_initiate_input
* \param nsubmeshes number of submesh blocks
*/
//void pumi_inputs_allocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes){
pumi_initiate_input_t* pumi_inputs_allocate(int nsubmeshes){
  pumi_initiate_input_t* pumi_inputs = (pumi_initiate_input_t*) malloc(sizeof(pumi_initiate_input_t));
  // 1D allocs
  pumi_inputs->p1_i = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->Nel_i = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->p2max_i = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->p2min_i = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->x_left = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->x_right = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->type_flag = (char**) malloc(nsubmeshes*sizeof(char*));
  int i;
  for (i=0; i<nsubmeshes; i++){
    pumi_inputs->type_flag[i] = (char*) malloc(SUBMESH_FLAGSTRING_LENGTH);
  }
  pumi_inputs->left_T = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->left_r = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->left_Nel = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->right_T = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->right_r = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->right_Nel = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->uniform_Nel = malloc(nsubmeshes*sizeof(int));
  // 2D allocs
  pumi_inputs->p1_i_x = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->p1_i_y = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->Nel_i_x = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->Nel_i_y = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->p2max_i_x = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->p2max_i_y = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->p2min_i_x = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->p2min_i_y = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->y_bottom = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->y_top = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->uniform_Nel_x = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->uniform_Nel_y = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->bottom_T = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->bottom_r = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->bottom_Nel = malloc(nsubmeshes*sizeof(int));
  pumi_inputs->top_T = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->top_r = malloc(nsubmeshes*sizeof(double));
  pumi_inputs->top_Nel = malloc(nsubmeshes*sizeof(int));
  return pumi_inputs;
}

/*!
* \brief Deallocates/Frees the memory allocated to members and object of struct pumi_initiate_input in pumi_initiate_allocate()
* \param *pumi_inputs pointer object to struct pumi_initiate_input
*/
void pumi_inputs_deallocate(pumi_initiate_input_t *pumi_inputs){
  //1D params
  free(pumi_inputs->p1_i);
  free(pumi_inputs->Nel_i);
  free(pumi_inputs->p2max_i);
  free(pumi_inputs->p2min_i);
  free(pumi_inputs->x_left);
  free(pumi_inputs->x_right);
  free(pumi_inputs->left_T);
  free(pumi_inputs->left_r);
  free(pumi_inputs->left_Nel);
  free(pumi_inputs->right_T);
  free(pumi_inputs->right_r);
  free(pumi_inputs->right_Nel);
  free(pumi_inputs->uniform_Nel);
  //2D params
  free(pumi_inputs->p1_i_x);
  free(pumi_inputs->p1_i_y);
  free(pumi_inputs->Nel_i_x);
  free(pumi_inputs->Nel_i_y);
  free(pumi_inputs->p2max_i_x);
  free(pumi_inputs->p2max_i_y);
  free(pumi_inputs->p2min_i_x);
  free(pumi_inputs->p2min_i_y);
  free(pumi_inputs->y_bottom);
  free(pumi_inputs->y_top);
  free(pumi_inputs->bottom_T);
  free(pumi_inputs->bottom_r);
  free(pumi_inputs->bottom_Nel);
  free(pumi_inputs->top_T);
  free(pumi_inputs->top_r);
  free(pumi_inputs->top_Nel);
  free(pumi_inputs->uniform_Nel_x);
  free(pumi_inputs->uniform_Nel_y);
  // common params
  int i;
  for (i=0; i<pumi_inputs->nsubmeshes; i++){
    free(pumi_inputs->type_flag[i]);
  }
  free(pumi_inputs->type_flag);
  free(pumi_inputs);
}

/*!
* \brief Deallocates/Frees the memory allocated to members and object of struct pumi_mesh in pumi_initiate()
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_finalize(pumi_mesh_t* pumi_mesh){
  //pumi_finalize_locate_functions();
  pumi_finalize_locatecell_and_calcweights_functions();
  if (pumi_mesh->BL_elem_coords_cache_flag){
      pumi_BL_elemsize_OFF(pumi_mesh);
  }
  free(pumi_mesh->pumi_submeshes_x);
  free(pumi_mesh);
}
/*!
* \brief Deallocates/Frees the memory allocated to variables in pumi_getmeshparameters_from_terminal()
* \param nsubmeshes number of submesh blocks
* \param submesh_params Matrix which contains the input submesh parameters for all submesh blocks
* \param submesh_flag Array which contains the bit flag representation of the mesh type flags for each submesh block

*/
void pumi_freemeshparameters_from_terminal(int nsubmeshes, double **submesh_params, unsigned int *submesh_flag){
  free(submesh_flag);
  int i;
  for (i=0; i<nsubmeshes; i++){
      free(submesh_params[i]);
  }
  free(submesh_params);
}

/*!
* \brief Computes a new grading ratio for left/right BL blocks
* \param BL_T thickness of left/right BL blocks
* \param BL_t0 first layer thickness in BL
* \param BL_Nel number of elements in the left/right BL block
* \details This routine will only be called inside the hpic initialize code
*/
double pumi_compute_grading_ratio_new(double BL_T, double BL_t0, int BL_Nel){
    double tol = 1e-5;
    double r = 1.5;
    double del_r = 1.0;
    int iter = 0;
    int max_iter = 10000;
    double f, f_r;
    while (fabs(del_r) > tol){
        f = BL_t0*pow(r,BL_Nel) - BL_T*r + (BL_T-BL_t0);
        f_r = BL_Nel*BL_t0*pow(r,BL_Nel-1) - BL_T;
        del_r = f/f_r;
        r -= del_r;
        iter++;
        if (iter > max_iter){
          printf("Cannot compute grading ratio. Solver not converging\n" );
          exit(0);
        }
    }
    return r;
}


/*!
* \brief Calls relevant routine based on dimension to perform checks on pumi params and verfies their validity (for diagnostics purposes)
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_verify_params(pumi_mesh_t *pumi_mesh){
  if (pumi_mesh->ndim == 1){
      pumi_verify_params_1D(pumi_mesh);
  }
  else{
      exit(0);
  }
}

/*!
* \brief Peforms checks on 1D pumi params and verfies their validity (for diagnostics purposes)
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_verify_params_1D(pumi_mesh_t *pumi_mesh){
  printf("\tVerifying valdity of 1D pumi parameters for\n");
  int flag = 0;
  int isubmesh;
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x; isubmesh++){
    printf("\tSUBMESH %d:\n", isubmesh+1 );

    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & leftBL){
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel > 0)){
        printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel);
        flag++;
      }
      else{
        printf("\t\t left_Nel    -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r > 1.0)){
        printf("\t\t left_r = %2.4e is not a valid input. It has to be greater than 1.0 for a graded BL mesh.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r);
        flag++;
      }
      else{
        printf("\t\t left_r      -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T > 0.0)){
        printf("\t\t left_T = %2.4e is not a valid input. It has to be positive real number.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T);
        flag++;
      }
      else{
        printf("\t\t left_T      -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0 > 0.0)){
        printf("\t\t left_t0 = %2.4e is not a valid calculated parameter. It has to be positive real number.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0);
        flag++;
      }
      else{
        printf("\t\t left_t0     -- verified...\n");
      }
    }

    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & uniform){
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel > 0)){
        printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel);
        flag++;
      }
      else{
        printf("\t\t uniform_Nel -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0 > 0.0)){
        printf("\t\t uniform_dx = %2.4e is not a valid calculated parameter. It has to be positive real number.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0);
        flag++;
      }
      else{
        printf("\t\t uniform_dx  -- verified...\n");
      }
    }

    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & rightBL){
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel > 0)){
        printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel);
        flag++;
      }
      else{
        printf("\t\t right_Nel   -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r > 1.0)){
        printf("\t\t right_r = %2.4e is not a valid input. It has to be greater than 1.0 for a graded BL mesh.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r);
        flag++;
      }
      else{
        printf("\t\t right_r     -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T > 0.0)){
        printf("\t\t right_T = %2.4e is not a valid input. It has to be positive real number.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_T);
        flag++;
      }
      else{
        printf("\t\t right_T     -- verified...\n");
      }
      if (!(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0 > 0.0)){
        printf("\t\t right_t0 = %2.4e is not a valid calculated parameter. It has to be positive real number.\n", ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0);
        flag++;
      }
      else{
        printf("\t\t right_t0    -- verified...\n");
      }
    }
  }

  if (flag == 0){
    printf("\n\tThe input mesh parameters and the calculated mesh parameters are all valid and verified\n\n");
  }
  else{
    printf("\t\nERROR: One or more input/calculated mesh paramater is not valid. Abort\n");
    pumi_finalize(pumi_mesh);
    exit(0);
  }
}
/*!
* \brief Prints the coordinates of the nodes (for diagnostics purposes)
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_print_node_coordinates(pumi_mesh_t *pumi_mesh){
    if (pumi_mesh->ndim == 1){
        pumi_print_node_coordinates_1D(pumi_mesh);
    }
    else{
        exit(0);
    }
}

/*!
* \brief Prints the coordinates of the nodes (for diagnostics purposes)
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_print_node_coordinates_1D(pumi_mesh_t *pumi_mesh){
  int isubmesh;
  int inode = 0;
  double coord, elem_size;
  printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
  for (isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x; isubmesh++){
    printf("SUBMESH %d:\n", isubmesh+1 );

    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & leftBL){
      FILE *lBL_fptr;
      char lBL_coord_file[30];
      sprintf(lBL_coord_file,"submesh%d_coord_leftBL.txt",isubmesh+1);
      lBL_fptr = fopen(lBL_coord_file,"w");
      printf("\tLeft BL segment:\n");
      inode = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->Nel_cumulative+1;
      coord = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->x_min;
      elem_size = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0;
      printf("\t\tNode %6d: %2.8e\n", inode, coord );
      fprintf(lBL_fptr, "%.16e\n", coord );
      int i;
      for (i=0; i<((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel; i++ ){
        inode++;
        coord += elem_size;
        printf("\t\tNode %6d: %2.8e\n", inode, coord );
        fprintf(lBL_fptr, "%.16e\n", coord );
        elem_size = elem_size*((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r;
      }
      printf("\tLeftBL coordinates written to the file \"%s\"\n\n",lBL_coord_file );
      fclose(lBL_fptr);
    }

    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & uniform){
      FILE *uni_fptr;
      char uni_coord_file[30];
      sprintf(uni_coord_file,"submesh%d_coord_uniform.txt",isubmesh+1);
      uni_fptr = fopen(uni_coord_file,"w");
      printf("\tUniform segment:\n");
      inode = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->Nel_cumulative+1;
      coord = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->x_min;
      double dx_uniform = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0;
      printf("\t\tNode %6d: %2.8e\n", inode, coord );
      fprintf(uni_fptr, "%.16e\n", coord );
      int i;
      for (i=0; i<((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel; i++ ){
        inode++;
        coord += dx_uniform;
        printf("\t\tNode %6d: %2.8e\n", inode, coord );
        fprintf(uni_fptr, "%.16e\n", coord );
      }
      printf("\tUniform segemnt coordinates written to the file \"%s\"\n\n",uni_coord_file );
      fclose(uni_fptr);
    }

    if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->pumi_flag & rightBL){
      FILE *rBL_fptr;
      char rBL_coord_file[30];
      sprintf(rBL_coord_file,"submesh%d_coord_rightBL.txt",isubmesh+1);
      rBL_fptr = fopen(rBL_coord_file,"w");
      printf("\tRight BL segment:\n");
      inode = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->Nel_cumulative+1;
      coord = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->x_min;
      elem_size = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->t0*pow(((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r,((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel-1);
      printf("\t\tNode %6d: %2.8e\n", inode, coord );
      fprintf(rBL_fptr, "%.16e\n", coord );
      int i;
      for (i=((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->submesh_Nel-1; i>=0; i-- ){
        inode++;
        coord += elem_size;
        printf("\t\tNode %6d: %2.8e\n", inode, coord );
        fprintf(rBL_fptr, "%.16e\n", coord );
        elem_size = elem_size/((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x + isubmesh)->r;
      }
      printf("\tRightBL coordinates written to the file \"%s\"\n\n",rBL_coord_file );
      fclose(rBL_fptr);
    }
    printf("\n");
  }
}
