#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumiMBBL.h"

int main()
{
  pumi_initiate_input_t *pumi_inputs;
  //pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_terminal, pumi_inputs);

  // assigning vals to struct pumi_initiate_input members
  pumi_inputs = malloc(sizeof(pumi_initiate_input_t));
  pumi_inputs->ndim = 1;
  pumi_inputs->nsubmeshes = 1;
  pumi_initiate_allocate(pumi_inputs, pumi_inputs->nsubmeshes);

  *(pumi_inputs->x_left + 0) = 0.0;
  *(pumi_inputs->x_right + 0) = 1.0;
  strcpy(pumi_inputs->type_flag[0], "uniform&rightBL");
  *(pumi_inputs->left_T + 0) = 0.0; //value zero needs to be assigned for inactive flag parameters
  *(pumi_inputs->left_r + 0) = 0.0;
  *(pumi_inputs->left_Nel + 0) = 0;
  *(pumi_inputs->right_T + 0) = 0.5;
  *(pumi_inputs->right_r + 0) = 1.2;
  *(pumi_inputs->right_Nel + 0) = 10;
  *(pumi_inputs->uniform_Nel + 0) = 10;

  // the pumi_initiate_input struct NEEDS TO BE POPULATED before calling this function
  pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs);
  pumi_initiate_deallocate(pumi_inputs, pumi_inputs->nsubmeshes);


  // calculating total number of elements in the mesh
  int Nel_total = pumi_total_elements_1D(pumi_mesh);
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
