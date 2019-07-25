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
  printf("PUMI mesh parameter info :\n\n");
  for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    printf("\tSUBMESH %d parameters:\n", isubmesh+1 );

    printf("\n\t submeshflag = ");
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      printf("leftBL&");
    }
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
      printf("uniform");
    }
    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      printf("&rightBL");
    }

    printf("\n\n\t uniform_Nel = %d    \t\t Number of Cells in uniform mesh region\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel);
    printf("\t uniform_dx  = %2.4e \t [m] Cell size in uniform mesh segment\n\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0);

    printf("\t left_t0     = %2.4e \t [m] Cell size of first/leftmost cell in left BL segment\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0);
    printf("\t left_T      = %2.4e \t [m] Left boundary layer (left BL) thickness\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_T);
    printf("\t left_r      = %2.4e \t Grading ratio in left BL mesh\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r);
    printf("\t left_Nel    = %d    \t\t Number of Cells in left BL mesh region\n\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel);

    printf("\t right_t0    = %2.4e \t [m] Cell size of last/rightmost cell in right BL segment\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0);
    printf("\t right_T     = %2.4e \t [m] Right boundary layer (right BL) thickness\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_T);
    printf("\t right_r     = %2.4e \t Grading ratio in right BL mesh\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r);
    printf("\t right_Nel   = %d    \t\t Number of Cells in right BL mesh region\n\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel);
  }
  pumi_verify_params(pumi_mesh);
  pumi_print_node_coordinates(pumi_mesh);
  return (pumi_mesh);
}

/*!
* \brief Assigns the values to members of struct pumi_submesh1D based on user inputs
* \param[out] *pumi_mesh pointer object of struct pumi_mesh whose members are assigned their corresponding values in this routine
* \param[in] isubmesh The index number of the submesh whose members (i.e submesh parameters) are being assigned their values
* \param[in] xleft Left end coordinate of the submesh block
* \param[in] xright Right end coordinate of the submesh block
* \param[in] submeshflag Mesh type flag of submesh block converted into a unsigned int using pumi_getsubmeshflag
* \param[in] N_uniform Number of elements in the uniform mesh segment of the submesh block
* \param[in] T_left Thickness of the left BL segment of the submesh block
* \param[in] r_left Growth ratio of elements in the left BL segment of the submesh block
* \param[in] N_left Number of elements in the left BL segment of the submesh block
* \param[in] T_right Thickness of the right BL segment of the submesh block
* \param[in] r_right Growth ratio of elements in the right BL segment of the submesh block
* \param[in] N_right Number of elements in the right BL segment of the submesh block
* \details The submesh parameters (used to define the submesh along with the dependent submesh paramaters) defined in the struct pumi_submesh1D are assigned their corresponding
values based on inputs from command line of hpic code or based on values supplied by the user directly from terminal prompt
*/
void pumi_setsubmesh(pumi_mesh_t *pumi_mesh, int isubmesh, double xleft, double xright, unsigned int submeshflag, int N_uniform, double T_left, double r_left, int N_left, double T_right, double r_right, int N_right){
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag = submeshflag;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left = xleft;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right = xright;

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel = N_uniform;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left = xleft + T_left; // (dependent variable)
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_right = xright - T_right; // (dependent variable)
  if (submeshflag & uniform){
    ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0 = ((xright - T_right) - (xleft + T_left))/N_uniform; // (dependent variable)
  }
  else{
    ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0 = 0.0; // (dependent variable)
  }

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_T = T_left;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r = r_left;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel = N_left;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_x_right = xleft + T_left; // (dependent variable)
  if (submeshflag & leftBL){
    ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0 = T_left*(r_left-1)/(pow(r_left,N_left)-1); // (dependent variable)
  }
  else{
    ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0 = 0.0; // (dependent variable)
  }

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_T = T_right;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r = r_right;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel = N_right;
  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_x_left = xright - T_right;// (dependent variable)
  if (submeshflag & rightBL){
    ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0 = T_right*(r_right-1)/(pow(r_right,N_right)-1);// (dependent variable)
  }
  else{
    ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0 = 0.0;// (dependent variable)
  }

  ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->submesh_total_Nel = N_uniform + N_left + N_right; // (dependent variable)
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
      printf("\n\tInvalid flag input\n");
      printf("\tList of valid input(s):\n\tuniform\n\tleftBL\n\trightBL\n\tleftBL&uniform\n\trightBL&uniform\n\tleftBL&rightBL\n\tuniform&leftBL&rightBL\n\n");
      printf("\tAssigning \"leftBL&uniform&rightBL\" (default input) to submesh flag...\n\n");
      intflag = uniform+leftBL+rightBL;
    }
    else{
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
  }
    return intflag;
}

/*!
* \brief Allocates the memory to members of struct pumi_initiate_input based on the number of submesh blocks
* \param *pumi_inputs pointer object to struct pumi_initiate_input
* \param nsubmeshes number of submesh blocks
*/
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

/*!
* \brief Deallocates/Frees the memory allocated to members and object of struct pumi_initiate_input in pumi_initiate_allocate()
* \param *pumi_inputs pointer object to struct pumi_initiate_input
* \param nsubmeshes number of submesh blocks
*/
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

/*!
* \brief Deallocates/Frees the memory allocated to members and object of struct pumi_mesh in pumi_initiate()
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_finalize(pumi_mesh_t* pumi_mesh){
  free(pumi_mesh->pumi_submeshes);
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
  for (int i=0; i<nsubmeshes; i++){
      free(submesh_params[i]);
  }
  free(submesh_params);
}

/*!
* \brief Computes a new grading ratio for left/right BL segment if a valid BL mesh with grading ratio 1.2 cannot be constructed with hpic inputs
* \param p1_lr number of debye lenghts in the left/right BL segment
* \param p2 number of grid cells per debye length in the uniform region
* \param BL_nel number of elements in the left/right BL segment
* \details This routine will only be called inside the hpic code when the inputs to hpic does
not correspond to a valid BL mesh with grading ratio 1.2. The routine returns a new grading ratio that will satisfy
the hpic inputs.
*/
double pumi_compute_grading_ratio(int p1_lr, int p2, int BL_Nel){
  double tol = 1e-5;
  double r = 1.2;
  if (BL_Nel == 0){
    r = 0.0;
    return r;
  }
  else{
    double del_r = 1.0;
    int iter = 0;
    int max_iter = 10000;
    double f, f_r;
    while (fabs(del_r) > tol){
      f = (double) (p1_lr*p2-1.0)*pow(r,BL_Nel) - (p1_lr*p2)*pow(r,BL_Nel-1) + 1.0;
      f_r = (double) BL_Nel*(p1_lr*p2-1.0)*pow(r,BL_Nel-1) - (BL_Nel-1)*(p1_lr*p2)*pow(r,BL_Nel-2);
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
}
/*!
* \brief Peforms checks on pumi params and verfies their validity (for diagnostics purposes)
* \param *pumi_mesh pointer object to struct pumi_initiate
*/
void pumi_verify_params(pumi_mesh_t *pumi_mesh){
  printf("\tVerifying valdity of pumi parameters for\n");
  int flag = 0;
  for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    printf("\tSUBMESH %d:\n", isubmesh+1 );

    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel > 0)){
        printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel);
        flag++;
      }
      else{
        printf("\t\t left_Nel    -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r > 1.0)){
        printf("\t\t left_r = %2.4e is not a valid input. It has to be greater than 1.0 for a graded BL mesh.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r);
        flag++;
      }
      else{
        printf("\t\t left_r      -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_T > 0.0)){
        printf("\t\t left_T = %2.4e is not a valid input. It has to be positive real number.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_T);
        flag++;
      }
      else{
        printf("\t\t left_T      -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0 > 0.0)){
        printf("\t\t left_t0 = %2.4e is not a valid calculated parameter. It has to be positive real number.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0);
        flag++;
      }
      else{
        printf("\t\t left_t0     -- verified...\n");
      }
    }

    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel > 0)){
        printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel);
        flag++;
      }
      else{
        printf("\t\t uniform_Nel -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0 > 0.0)){
        printf("\t\t uniform_dx = %2.4e is not a valid calculated parameter. It has to be positive real number.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0);
        flag++;
      }
      else{
        printf("\t\t uniform_dx  -- verified...\n");
      }
    }

    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel > 0)){
        printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel);
        flag++;
      }
      else{
        printf("\t\t right_Nel   -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r > 1.0)){
        printf("\t\t right_r = %2.4e is not a valid input. It has to be greater than 1.0 for a graded BL mesh.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r);
        flag++;
      }
      else{
        printf("\t\t right_r     -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_T > 0.0)){
        printf("\t\t right_T = %2.4e is not a valid input. It has to be positive real number.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_T);
        flag++;
      }
      else{
        printf("\t\t right_T     -- verified...\n");
      }
      if (!(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0 > 0.0)){
        printf("\t\t right_t0 = %2.4e is not a valid calculated parameter. It has to be positive real number.\n", ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0);
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
  int N_cumulative[pumi_mesh->nsubmeshes];
  N_cumulative[0] = 0;
  for (int isubmesh=1; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    N_cumulative[isubmesh] = N_cumulative[isubmesh-1] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + (isubmesh-1))->submesh_total_Nel;
  }

  int inode = 0;
  double coord, elem_size;
  printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
  for (int isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
    printf("SUBMESH %d:\n", isubmesh+1 );

    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & leftBL){
      FILE *lBL_fptr;
      char lBL_coord_file[] = "coord_leftBL.txt";
      lBL_fptr = fopen(lBL_coord_file,"w");
      printf("\tLeft BL segment:\n");
      inode = N_cumulative[isubmesh]+1;
      coord = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left;
      elem_size = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->lBL_t0;
      printf("\t\tNode %6d: %2.4e\n", inode, coord );
      fprintf(lBL_fptr, "%2.4e\n", coord );
      for (int i=0; i<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel; i++ ){
        inode++;
        coord += elem_size;
        printf("\t\tNode %6d: %2.4e\n", inode, coord );
        fprintf(lBL_fptr, "%2.4e\n", coord );
        elem_size = elem_size*((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_r;
      }
      printf("\tLeftBL coordinates written to the file \"%s\"\n\n",lBL_coord_file );
      fclose(lBL_fptr);
    }

    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & uniform){
      FILE *uni_fptr;
      char uni_coord_file[] = "coord_uniform.txt";
      uni_fptr = fopen(uni_coord_file,"w");
      printf("\tUniform segment:\n");
      inode = N_cumulative[isubmesh] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel + 1;
      coord = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_x_left;
      double dx_uniform = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_t0;
      printf("\t\tNode %6d: %2.4e\n", inode, coord );
      fprintf(uni_fptr, "%2.4e\n", coord );
      for (int i=0; i<((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel; i++ ){
        inode++;
        coord += dx_uniform;
        printf("\t\tNode %6d: %2.4e\n", inode, coord );
        fprintf(uni_fptr, "%2.4e\n", coord );
      }
      printf("\tUniform segemnt coordinates written to the file \"%s\"\n\n",uni_coord_file );
      fclose(uni_fptr);
    }

    if (((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->pumi_flag & rightBL){
      FILE *rBL_fptr;
      char rBL_coord_file[] = "coord_rightBL.txt";
      rBL_fptr = fopen(rBL_coord_file,"w");
      printf("\tRight BL segment:\n");
      inode = N_cumulative[isubmesh] + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->left_Nel + ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->uniform_Nel + 1;
      coord = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_x_left;
      elem_size = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->rBL_t0*pow(((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r,((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-1);
      printf("\t\tNode %6d: %2.4e\n", inode, coord );
      fprintf(rBL_fptr, "%2.4e\n", coord );
      for (int i=((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_Nel-1; i>=0; i-- ){
        inode++;
        coord += elem_size;
        printf("\t\tNode %6d: %2.4e\n", inode, coord );
        fprintf(rBL_fptr, "%2.4e\n", coord );
        elem_size = elem_size/((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->right_r;
      }
      printf("\tRightBL coordinates written to the file \"%s\"\n\n",rBL_coord_file );
      fclose(rBL_fptr);
    }
    printf("\n");
  }
}
