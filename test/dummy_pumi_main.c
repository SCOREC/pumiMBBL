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

  if (argc < 6){
    printf("Execute the code with the following command line arguments -- \n\n" );
    printf("\t ./install/bin/pumiMBBL_Demo p1 p2 \"typeflag\" Nel_max_FLAG alpha-or-Nel_max\n\n");
    printf("\t p1    \t\t Number of Debye Lengths In Domain \n");
    printf("\t p2    \t\t Number of Grid Points Per Debye Length \n");
    printf("\t \"typeflag\" \t Active mesh type segments in the submesh\n" );
    printf("\t Nel_max_FLAG \t Flag to specify type of input for maximum number of elements allowed in the domain\n\n");
    printf("\t if (Nel_max_FLAG = 0)\n");
    printf("\t \talpha \t Multiplicative factor to determine maximum number of elements allowed in the domain (has to be greater then 1.0)\n ");
    printf("\t \t      \t This sets Nel_max to be floor(alpha*p1*p2)\n");
    printf("\t if (Nel_max_FLAG = 1)\n");
    printf("\t \tNel_max  Maximum number of elements allowed in the domain (has to be greater than p1*p2)\n");
    printf("  E.g.\n\n");
    printf("    ./install/bin/pumiMBBL_Demo  50 1 \"uniform&leftBL&rightBL\" 0 1.3\n");
    printf("    ./install/bin/pumiMBBL_Demo  50 2 \"uniform&leftBL&rightBL\" 1 120\n");
    exit(0);
  }

  int NumberDebyeLengthsInDomain = atoi ( argv[1] );
  int NumberPointsPerDebyeLength = atoi ( argv[2] );
  int Ncells = NumberDebyeLengthsInDomain*NumberPointsPerDebyeLength;

  pumi_inputs = malloc(sizeof(pumi_initiate_input_t));
  pumi_inputs->ndim = 1; // Fixed pumi input
  pumi_inputs->nsubmeshes = 1; // Fixed pumi input
  pumi_initiate_allocate(pumi_inputs, pumi_inputs->nsubmeshes);
  int isubmesh = 0;
  strcpy(pumi_inputs->type_flag[isubmesh], argv[3]);
  pumi_inputs->Nel_max_FLAG = atoi( argv[4] );
  if (pumi_inputs->Nel_max_FLAG == 0){ // read alpha and compute Nel_max
    pumi_inputs->alpha = atof( argv[5] );
    pumi_inputs->Nel_max = floor( pumi_inputs->alpha*NumberDebyeLengthsInDomain*NumberPointsPerDebyeLength );
  }
  else if (pumi_inputs->Nel_max_FLAG == 1){// read Nel_max directly
    pumi_inputs->Nel_max = atoi( argv[5] );
    pumi_inputs->alpha = (double) pumi_inputs->Nel_max/(NumberDebyeLengthsInDomain*NumberPointsPerDebyeLength );
  }
  else{
    pumi_inputs->alpha = 1.1;
    pumi_inputs->Nel_max = floor( pumi_inputs->alpha*NumberDebyeLengthsInDomain*NumberPointsPerDebyeLength );
  }


  if (pumi_inputs->alpha <= 1.0 || pumi_inputs->Nel_max <= Ncells){ //take default values
    if (pumi_inputs->Nel_max_FLAG == 0){
      printf("WARNING: Boundary layer mesh not possible for given alpha. Assigning default value for alpha as 1.1 and computing Nel_max...\n\n");
    }
    if (pumi_inputs->Nel_max_FLAG == 1){
      printf("WARNING: Boundary layer mesh not possible for given Nel_max. Assigning default value for alpha as 1.1 and computing Nel_max...\n\n");
    }
    pumi_inputs->alpha = 1.1;
    pumi_inputs->Nel_max = floor( pumi_inputs->alpha*NumberDebyeLengthsInDomain*NumberPointsPerDebyeLength );
  }

  double lambda_D = 3.3246e-04;
  double x1_min = 0.0;
  double x1_max = lambda_D*NumberDebyeLengthsInDomain;

  *(pumi_inputs->x_left + isubmesh) = x1_min;
  *(pumi_inputs->x_right + isubmesh) = x1_max;
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

  double domain_L = x1_max - x1_min;
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
    uniform_L = domain_L-left_T-right_T;
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
    if (Nel_total >  pumi_inputs->Nel_max){ // case (ii) in pumi mesh parameters document
      int Nel_BL_total = pumi_inputs->Nel_max - uniform_Nel;
      left_Nel = floor( (double) (pumi_inputs->p1_l*Nel_BL_total/(pumi_inputs->p1_l+pumi_inputs->p1_r)) );
      right_Nel = floor( (double) (pumi_inputs->p1_r*Nel_BL_total/(pumi_inputs->p1_l+pumi_inputs->p1_r)) );
      left_r = pumi_compute_grading_ratio(pumi_inputs->p1_l, NumberPointsPerDebyeLength, left_Nel);
      right_r = pumi_compute_grading_ratio(pumi_inputs->p1_r, NumberPointsPerDebyeLength, right_Nel);
    }
  }

  if (left_t0 < 0.0 || right_t0 < 0.0){ //case (iii) in pumi mesh parameters document
    int Nel_BL_total = pumi_inputs->Nel_max - uniform_Nel;
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

  *(pumi_inputs->uniform_Nel + isubmesh) = uniform_Nel;
  *(pumi_inputs->left_T + isubmesh)      = left_T;
  *(pumi_inputs->left_r + isubmesh)      = left_r;
  *(pumi_inputs->left_Nel + isubmesh)    = left_Nel;
  *(pumi_inputs->right_T + isubmesh)     = right_T;
  *(pumi_inputs->right_r + isubmesh)     = right_r;
  *(pumi_inputs->right_Nel + isubmesh)   = right_Nel;

  // the pumi_initiate_input struct NEEDS TO BE POPULATED before calling this function
  pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs);
  pumi_initiate_deallocate(pumi_inputs, pumi_inputs->nsubmeshes);


  // calculating total number of elements in the mesh
  int Nel_total = pumi_total_elements(pumi_mesh);
  double X_LEFT = pumi_global_x_left_1D(pumi_mesh); //Global x_left
  double X_RIGHT = pumi_global_x_right_1D(pumi_mesh); //Global x_right

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

  double *covolume = (double*) malloc((Nel_total+1)*sizeof(double));
  pumi_compute_covolume_1D(pumi_mesh, Nel_total, covolume); //populates the variable covolume

  double charge_density[Nel_total+1];
  for (int i=0; i<(Nel_total+1); i++){
    charge_density[i] = grid_weights[i]/covolume[i];
    printf("Charge density at node %d is %lf\n", i, charge_density[i]);
  }

  free(coordinates);
  free(grid_weights);
  free(covolume);
  pumi_finalize(pumi_mesh); //deallocates pumi_mesh object


  return 0;
}
