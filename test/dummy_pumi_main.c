#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumiMBBL.h"

int main(int argc, char *argv[])
{
  pumi_initiate_input_t    *pumi_inputs;
  //pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_terminal, pumi_inputs);

  ///*
  if (argc < 8){
    printf("Execute the code with the following command line arguments -- \n\n" );
    printf("\t ./install/bin/pumiMBBL_Demo p1 p2 N \"typeflag_1,..,typeflag_N\" \"p1_1,..,p1_N\" \"Nel_max_FLAG_1,..,Nel_max_FLAG_N\" \"alpha_1/Nel_max_1,..,alpha_N/Nel_max_N\"\n\n");
    printf("\t p1    \t\t Number of Debye Lengths In Domain \n");
    printf("\t p2    \t\t Number of Grid Points Per Debye Length \n");
    printf("\t N     \t\t Number of submeshes in the domain \n");
    printf("\t \"typeflag_i\" \t Active mesh type segments in i-th submesh\n" );
    printf("\t p1_i  \t\t Number of Debye Lengths In i-th submesh \n");
    printf("\t Nel_max_FLAG_i \t Flag to specify type of input for maximum number of elements allowed in i-th submesh\n\n");
    printf("\t if (Nel_max_FLAG_i = 0)\n");
    printf("\t \talpha_i \t Multiplicative factor to determine maximum number of elements allowed in i-th submesh (has to be greater then 1.0)\n ");
    printf("\t \t        \t This sets Nel_max to be floor(alpha_i*p1_i*p2)\n");
    printf("\t if (Nel_max_FLAG_i = 1)\n");
    printf("\t \tNel_max_i  Maximum number of elements allowed in i-th submesh (has to be greater than p1_i*p2)\n\n");
    printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
    printf("  E.g.\n\n");
    printf("    ./install/bin/pumiMBBL_Demo 50 1 2 \"uniform&leftBL&rightBL,uniform&rightBL\" \"30,20\" \"0,0\" \"1.3,1.5\"\n");
    printf("    ./install/bin/pumiMBBL_Demo 50 1 2 \"uniform&leftBL,uniform&rightBL\" \"30,20\" \"0,1\" \"1.3,35\"\n");
    exit(0);
  }

  // necessary hpic parameters needed for pumi mesh
  int NumberDebyeLengthsInDomain = atoi ( argv[1] );
  int NumberPointsPerDebyeLength = atoi ( argv[2] );
  int Ncells = NumberDebyeLengthsInDomain*NumberPointsPerDebyeLength;

  pumi_inputs = malloc(sizeof(pumi_initiate_input_t));
  pumi_inputs->ndim = 1; // Fixed pumi input
  pumi_inputs->nsubmeshes = atoi( argv[3] );

  char all_submesh_flag[1000];
  char each_submesh_flag[10][100];
  strcpy(all_submesh_flag, argv[4]);

  char *tok = strtok(all_submesh_flag, ",");
  int i_submesh=0;
  while (tok != NULL){
    strcpy (each_submesh_flag[i_submesh], tok);
    tok = strtok(NULL, ",");
    i_submesh++;
  }

  if (i_submesh != pumi_inputs->nsubmeshes){
      printf("ERROR: Number of typeflag arguments not equal to number of submeshes...\n");
      exit(0);
  }

  char all_p1_submesh[50];
  char each_p1_submesh[10][5];
  strcpy(all_p1_submesh, argv[5]);

  tok = strtok(all_p1_submesh, ",");
  i_submesh=0;
  while (tok != NULL){
    strcpy (each_p1_submesh[i_submesh], tok);
    tok = strtok(NULL, ",");
    i_submesh++;
  }

  if (i_submesh != pumi_inputs->nsubmeshes){
      printf("ERROR: Number of p1_i arguments not equal to number of submeshes...\n");
      exit(0);
  }

  int p1_check=0;
  for (int j=0; j<i_submesh; j++){
      p1_check += atoi(each_p1_submesh[j]);
  }
  if (p1_check != NumberDebyeLengthsInDomain){
      printf("ERROR: Submesh Debye lengths not adding up to total Debye length of domain...\n");
      exit(0);
  }

  char all_Nel_max_flag[30];
  char each_Nel_max_flag[10][3];
  strcpy(all_Nel_max_flag, argv[6]);

  tok = strtok(all_Nel_max_flag, ",");
  i_submesh=0;
  while (tok != NULL){
    strcpy (each_Nel_max_flag[i_submesh], tok);
    tok = strtok(NULL, ",");
    i_submesh++;
  }

  if (i_submesh != pumi_inputs->nsubmeshes){
      printf("ERROR: Number of Nel_max_flag arguments not equal to number of submeshes...\n");
      exit(0);
  }

  char all_Nel_max_alpha[100];
  char each_Nel_max[10][10];
  char each_alpha[10][10];
  strcpy(all_Nel_max_alpha, argv[7]);

  tok = strtok(all_Nel_max_alpha, ",");
  i_submesh=0;
  while (tok != NULL){
    if (atoi(each_Nel_max_flag[i_submesh])){
        strcpy(each_Nel_max[i_submesh], tok);
    }
    else{
        strcpy(each_alpha[i_submesh], tok);
    }
    tok = strtok(NULL, ",");
    i_submesh++;
  }

  if (i_submesh != pumi_inputs->nsubmeshes){
      printf("ERROR: Number of Nel_max/alpha arguments not equal to number of submeshes...\n");
      exit(0);
  }

  pumi_inputs_allocate(pumi_inputs, pumi_inputs->nsubmeshes);

  double lambda_D = 3.3246e-04;
  double x1_min = 0.0;
  double x1_max = lambda_D*NumberDebyeLengthsInDomain;

  int isubmesh;

  for(isubmesh=0; isubmesh<pumi_inputs->nsubmeshes; isubmesh++){
      strcpy(pumi_inputs->type_flag[isubmesh], each_submesh_flag[isubmesh]);
      *(pumi_inputs->p1_i + i_submesh) = atoi( each_p1_submesh[isubmesh]);
      *(pumi_inputs->Nel_max_FLAG + isubmesh) = atoi( each_Nel_max_flag[isubmesh] );
      if (*(pumi_inputs->Nel_max_FLAG + isubmesh) == 0){ // read alpha and compute Nel_max
        *(pumi_inputs->alpha + isubmesh) = atof( each_alpha[isubmesh] );
        *(pumi_inputs->Nel_max + isubmesh) = floor( *(pumi_inputs->alpha + isubmesh) * *(pumi_inputs->p1_i + i_submesh) * NumberPointsPerDebyeLength );
        //printf("For SUBMESH %d: Nel_max = %d\n",isubmesh,  *(pumi_inputs->Nel_max + isubmesh));
      }
      else if (*(pumi_inputs->Nel_max_FLAG + isubmesh) == 1){// read Nel_max directly
        *(pumi_inputs->Nel_max + isubmesh) = atoi( each_Nel_max[isubmesh] );
        *(pumi_inputs->alpha + isubmesh) = (double) *(pumi_inputs->Nel_max + isubmesh)/(*(pumi_inputs->p1_i + i_submesh) * NumberPointsPerDebyeLength );
        //printf("For SUBMESH %d: Nel_max = %d\n",isubmesh,  *(pumi_inputs->Nel_max + isubmesh));
      }
      else{
        *(pumi_inputs->alpha + isubmesh) = 1.1;
        *(pumi_inputs->Nel_max + isubmesh) = floor( *(pumi_inputs->alpha + isubmesh) * *(pumi_inputs->p1_i + i_submesh) * NumberPointsPerDebyeLength );
        //printf("For SUBMESH %d: Nel_max = %d\n",isubmesh,  *(pumi_inputs->Nel_max + isubmesh));
      }


      if (*(pumi_inputs->alpha + isubmesh) <= 1.0){ //take default values
        if (*(pumi_inputs->Nel_max + isubmesh) == 0){
          printf("WARNING: Boundary layer mesh not possible for given alpha. Assigning default value for alpha as 1.1 and computing Nel_max...\n\n");
        }
        if (*(pumi_inputs->Nel_max + isubmesh) == 1){
          printf("WARNING: Boundary layer mesh not possible for given Nel_max. Assigning default value for alpha as 1.1 and computing Nel_max...\n\n");
        }
        *(pumi_inputs->alpha + isubmesh) = 1.1;
        *(pumi_inputs->Nel_max + isubmesh) = floor( *(pumi_inputs->alpha + isubmesh) * *(pumi_inputs->p1_i + i_submesh) * NumberPointsPerDebyeLength );
      }



      // calculating all pumi_inputs to initiate the mesh
      if (isubmesh == 0){
          *(pumi_inputs->x_left + isubmesh) = x1_min;
          *(pumi_inputs->x_right + isubmesh) = *(pumi_inputs->x_left + isubmesh) + lambda_D* *(pumi_inputs->p1_i + i_submesh);
      }
      else if (isubmesh>0 && isubmesh<(pumi_inputs->nsubmeshes-1)){
          *(pumi_inputs->x_left + isubmesh) = *(pumi_inputs->x_right + (isubmesh-1));
          *(pumi_inputs->x_right + isubmesh) = *(pumi_inputs->x_left + isubmesh) + lambda_D* *(pumi_inputs->p1_i + i_submesh);
      }
      else{
          *(pumi_inputs->x_left + isubmesh) = *(pumi_inputs->x_right + (isubmesh-1));
          *(pumi_inputs->x_right + isubmesh) = x1_max;
      }

      char tmp_typeflag_string[SUBMESH_FLAGSTRING_LENGTH];
      strcpy(tmp_typeflag_string, pumi_inputs->type_flag[isubmesh]);
      unsigned int typeflag = pumi_getsubmeshflag(tmp_typeflag_string);
      if (!(typeflag & uniform)){
        printf("\t\"uniform\" mesh flag needs to be specified. The mesh will not span the full domain otherwise\n" );
        printf("\tAdding \"uniform\" flag to the typeflag input\n\n");
        strcat(pumi_inputs->type_flag[isubmesh], "&uniform");
        strcpy(tmp_typeflag_string, pumi_inputs->type_flag[isubmesh]);
        typeflag = pumi_getsubmeshflag(tmp_typeflag_string);
      }

      double subdomain_L = *(pumi_inputs->x_right + isubmesh) - *(pumi_inputs->x_left + isubmesh);
      double dx1 = lambda_D/NumberPointsPerDebyeLength;

      // Initialize submesh parameters to be 0
      pumi_inputs->p1_l = 0;
      double left_r = 0.0;
      double left_T = 0.0;
      double left_t0 = 0.0;
      int left_Nel = 0;

      pumi_inputs->p1_r = 0;
      double right_r = 0.0;
      double right_T = 0.0;
      double right_t0 = 0.0;
      int right_Nel = 0;

      double uniform_dx1 = dx1;
      double uniform_L = 0.0;
      int uniform_Nel = 0;

      double beta = 10.0; // factor that determines the lower limits of left_t0 and right_t0

      if (typeflag & leftBL){// set values if leftBL flag is active
        pumi_inputs->p1_l = 10;
        left_r = 1.25;
        left_T = pumi_inputs->p1_l*lambda_D;
      }

      if (typeflag & rightBL){// set values if rightBL flag is active
        pumi_inputs->p1_r = 10;
        right_r = 1.25;
        right_T = pumi_inputs->p1_r*lambda_D;
      }

      if (typeflag & uniform){
        uniform_L = subdomain_L-left_T-right_T;
        uniform_Nel = floor (uniform_L/uniform_dx1);
        left_t0 = left_r*uniform_dx1 - (left_r-1.0)*left_T;
        right_t0 = right_r*uniform_dx1 - (right_r-1.0)*right_T;

        if (typeflag & leftBL){ //leftBL&uniform
          left_Nel = 1 + floor( log(uniform_dx1/left_t0)/log(left_r) );
          left_r = pumi_compute_grading_ratio(pumi_inputs->p1_l, NumberPointsPerDebyeLength, left_Nel); //recompute grading ratio for rounded off left_Nel
        }
        if (typeflag & rightBL){ //rightBL&uniform
          right_Nel = 1 + floor( log(uniform_dx1/right_t0)/log(right_r) );
          right_r = pumi_compute_grading_ratio(pumi_inputs->p1_r, NumberPointsPerDebyeLength, right_Nel); //recompute grading ratio for rounded off right_Nel
        }
        int Nel_total = uniform_Nel + right_Nel + left_Nel;
        if (Nel_total >  *(pumi_inputs->Nel_max + isubmesh)){ // case (ii) in pumi mesh parameters document
          //printf("For SUBMESH %d: Computed Nel_total = %d is greater than input Nel_max %d\n", isubmesh, Nel_total, *(pumi_inputs->Nel_max + isubmesh));
          int Nel_BL_total = *(pumi_inputs->Nel_max + isubmesh) - uniform_Nel;
          left_Nel = floor( (double) (pumi_inputs->p1_l*Nel_BL_total/(pumi_inputs->p1_l+pumi_inputs->p1_r)) );
          right_Nel = floor( (double) (pumi_inputs->p1_r*Nel_BL_total/(pumi_inputs->p1_l+pumi_inputs->p1_r)) );
          left_r = pumi_compute_grading_ratio(pumi_inputs->p1_l, NumberPointsPerDebyeLength, left_Nel);
          right_r = pumi_compute_grading_ratio(pumi_inputs->p1_r, NumberPointsPerDebyeLength, right_Nel);
          left_t0 = left_T*(left_r-1)/(pow(left_r,left_Nel)-1.0);
          right_t0 = right_T*(right_r-1)/(pow(right_r,right_Nel)-1.0);
        }
      }

      if (left_t0 < 0.0 || right_t0 < 0.0){ //case (iii) in pumi mesh parameters document
        printf("For SUBMESH %d: mesh will not span the submesh domain. Recalculating parameters...\n", isubmesh);
        int Nel_BL_total = *(pumi_inputs->Nel_max + isubmesh) - uniform_Nel;
        left_Nel = floor( (double) (pumi_inputs->p1_l*Nel_BL_total/(pumi_inputs->p1_l+pumi_inputs->p1_r)) );
        right_Nel = floor( (double) (pumi_inputs->p1_r*Nel_BL_total/(pumi_inputs->p1_l+pumi_inputs->p1_r)) );
        left_r = pumi_compute_grading_ratio(pumi_inputs->p1_l, NumberPointsPerDebyeLength, left_Nel);
        right_r = pumi_compute_grading_ratio(pumi_inputs->p1_r, NumberPointsPerDebyeLength, right_Nel);
        left_t0 = left_T*(left_r-1)/(pow(left_r,left_Nel)-1.0);
        right_t0 = right_T*(right_r-1)/(pow(right_r,right_Nel)-1.0);
      }

      if (typeflag & leftBL){ //case (iv) in pumi mesh parameters document
        if (left_t0 < lambda_D/beta ){
          printf("\tWARNING: For the given inputs left_t0 is less than lambda_D/10...\n" );
          printf("\t\tRecomputing left_r and left_Nel with the new constraint -- left_t0 = lambda_D/10...\n" );
          left_r = NumberPointsPerDebyeLength*(beta*pumi_inputs->p1_l - 1.0)/(beta*pumi_inputs->p1_l*NumberPointsPerDebyeLength -beta);
          left_Nel = 1 + floor( log(beta/NumberPointsPerDebyeLength)/log(left_r) );
          left_r = pumi_compute_grading_ratio(pumi_inputs->p1_l, NumberPointsPerDebyeLength, left_Nel);
        }
      }
      if (typeflag & rightBL){ //case (iv) in pumi mesh parameters document
        if (right_t0 < lambda_D/beta ){
          printf("\tWARNING: For the given inputs right_t0 is less than lambda_D/10...\n" );
          printf("\t\tRecomputing right_r and right_Nel with the new constraint -- right_t0 = lambda_D/10...\n" );
          right_r = NumberPointsPerDebyeLength*(beta*pumi_inputs->p1_r - 1.0)/(beta*pumi_inputs->p1_r*NumberPointsPerDebyeLength -beta);
          right_Nel = 1 + floor( log(beta/NumberPointsPerDebyeLength)/log(right_r) );
          right_r = pumi_compute_grading_ratio(pumi_inputs->p1_r, NumberPointsPerDebyeLength, right_Nel);
        }
      }

      // all pumi_inputs are calculated
      *(pumi_inputs->uniform_Nel + isubmesh) = uniform_Nel;
      *(pumi_inputs->left_T + isubmesh)      = left_T;
      *(pumi_inputs->left_r + isubmesh)      = left_r;
      *(pumi_inputs->left_Nel + isubmesh)    = left_Nel;
      *(pumi_inputs->right_T + isubmesh)     = right_T;
      *(pumi_inputs->right_r + isubmesh)     = right_r;
      *(pumi_inputs->right_Nel + isubmesh)   = right_Nel;
  }



  // the pumi_input object NEEDS TO BE POPULATED before initializing pumi_mesh
  pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs);
  // deallocate memory allocated to pumi_inputs -- Always do this IMMEDIATELY AFTER pumi_initiate()
  pumi_inputs_deallocate(pumi_inputs, pumi_inputs->nsubmeshes);
  
  // Call this function if BL element sizes are to be precomputed
  // HIGHLY RECOMMENDED to call this function to ensure log and power functions are not used to locate particle in BL
  pumi_BL_elemsize_ON(pumi_mesh);

  // calculating total number of elements in the mesh
  int Nel_total = pumi_total_elements(pumi_mesh);
  double X_LEFT = pumi_global_x_left_1D(pumi_mesh); //Global x_left
  double X_RIGHT = pumi_global_x_right_1D(pumi_mesh); //Global x_right

  printf("Smallest element size is %2.4e\n", pumi_return_smallest_elemsize(pumi_mesh));

  //double *elemsize;
  //double *gradingratio;
  //pumi_compute_elemsize_1D(pumi_mesh, Nel_total, elemsize);
  //pumi_compute_covolume_1D(0, Nel_total, elemsize);
  //pumi_compute_nodal_gradingratio_1D(elemsize, Nel_total, gradingratio);

  for (int i=1; i<Nel_total; i++){
    // pumi_return_gradingratio() returns the grading ratio for given node 'i'
    double r = pumi_return_gradingratio(pumi_mesh, i);
    printf("Node %d = %2.4e\n", i+1, r);
  }

  int num_particles;
  printf("\nEnter number of particles : ");
  scanf("%d", &num_particles); //user supplied

  double *coordinates = (double*) malloc(num_particles * sizeof(double));

  srand48(time(NULL));
  for (int i=0; i<num_particles; i++){
    coordinates[i] = X_LEFT + (X_RIGHT-X_LEFT)*drand48(); //generates random coordinates within the global domain for each particle
  }

  double *grid_weights = (double*) malloc((Nel_total+1)*sizeof(double));
  for (int i=0; i<= Nel_total; i++){
    grid_weights[i] = 0.0; //intialize charge to be zero at all nodes
  }
  for (int iparticle=0; iparticle<num_particles; iparticle++){ //loop over all particles
    int kcell;
    double Wgh1, Wgh2;
    pumi_locatepoint(pumi_mesh, coordinates[iparticle], &kcell, &Wgh2); // computes paricle cell and weight (based on linear weighting)
    Wgh1 = 1.0 - Wgh2;
    grid_weights[kcell]   += Wgh1;
    grid_weights[kcell+1] += Wgh2; //accumulate the weights for each particle
  }

  for (int i=0; i<(Nel_total+1); i++){
    printf("node=%d   gridweight=%lf\n", i, grid_weights[i]);
  }

  for (int i=0; i<Nel_total; i++){
    // pumi_return_elemsize()  returns elemtent size for given element idex 'i' [offset=0].
    // nodal index can also be specified but appropriate offset to be supplied with it
    printf("Element size of element %d is %2.4e\n", i, pumi_return_elemsize(pumi_mesh, i, 0) );
  }

  double charge_density[Nel_total+1];
  for (int i=0; i<(Nel_total+1); i++){
    // pumi_return_covolume_1D() retruns covolume of a given node
    charge_density[i] = grid_weights[i]/pumi_return_covolume_1D(pumi_mesh, i);
    printf("Charge density at node %d is %2.4e\n", i+1, charge_density[i]);
  }

  // always call this function if pumi_BL_elemsize_ON() is previously calculated
  // This function is to be called towards the end of the code (after the simulation is completed and before
  // hpic_finalize() is called)
  pumi_BL_elemsize_OFF(pumi_mesh);

  // free allocated memory
  free(coordinates);
  free(grid_weights);

  pumi_finalize(pumi_mesh); //deallocates pumi_mesh object


  return 0;
}
