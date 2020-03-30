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
    ///*
    if (argc < 6){
      printf("Execute the code with the following command line arguments -- \n\n" );
      printf("\t ./install/bin/pumiMBBL_Demo N \"typeflag_i\" \"p1_i\" \"p2_max_i\" \"p2_min_i\"\n\n\n");
      printf("\t N     \t\t\t Total Number of submeshes in the domain \n\n");
      printf("\t \"typeflag_i\" \t\t Active mesh type segment in i-th submesh\n\n" );
      printf("\t \"p1_i\"  \t\t Number of Debye Lengths In i-th submesh \n\n");
      printf("\t \"p2_max_i\"  \t\t For leftBL/rightBL, Number of maximum size cells in a Debye Length for i-th submesh \n");
      printf("\t \t  \t\t For uniform, Number of cells in a Debye Length for i-th submesh \n\n");
      printf("\t \"p2_min_i\"  \t\t For leftBL/rightBL, Number of minimum size cells in a Debye Length for i-th submesh \n");
      printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
      printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
      printf("  E.g.\n\n");
      printf("    ./install/bin/pumiMBBL_Demo 3 \"leftBL,uniform,rightBL\" \"10,25,15\" \"1,1,1\" \"10,0,10\"\n");
      exit(0);
    }


    pumi_inputs = malloc(sizeof(pumi_initiate_input_t));
    pumi_inputs->ndim = 1; // Fixed pumi input
    pumi_inputs->nsubmeshes = atoi( argv[1] );

    // reading submesh meshtypes
    char all_submesh_flag[MAX_SUBMESHES*SEGMENT_STRING_LENGTH];
    char each_submesh_flag[MAX_SUBMESHES][SEGMENT_STRING_LENGTH];
    strcpy(all_submesh_flag, argv[2]);

    char *tok = strtok(all_submesh_flag, ",");
    int isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_flag[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes){
        printf("ERROR: Number of typeflag arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh lengths
    char all_p1_submesh[MAX_SUBMESHES*5];
    char each_p1_submesh[MAX_SUBMESHES][5];
    strcpy(all_p1_submesh, argv[3]);

    tok = strtok(all_p1_submesh, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes){
        printf("ERROR: Number of p1_i arguments not equal to number of submeshes...\n");
        exit(0);
    }


    //reading submesh max elemsize
    char all_p2max_submesh[MAX_SUBMESHES*2];
    char each_p2max_submesh[MAX_SUBMESHES][2];
    strcpy(all_p2max_submesh, argv[4]);

    tok = strtok(all_p2max_submesh, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2max_submesh[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes){
        printf("ERROR: Number of p2max_i arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh min elemsize
    char all_p2min_submesh[MAX_SUBMESHES*2];
    char each_p2min_submesh[MAX_SUBMESHES][2];
    strcpy(all_p2min_submesh, argv[5]);

    tok = strtok(all_p2min_submesh, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes){
        printf("ERROR: Number of p2max_i arguments not equal to number of submeshes...\n");
        exit(0);
    }

    pumi_inputs_allocate(pumi_inputs, pumi_inputs->nsubmeshes);

    double lambda_D = 3.3246e-04;
    double x1_min = 0.0;
    int NumberDebyeLengthsInDomain=0;
    for(isubmesh=0; isubmesh<pumi_inputs->nsubmeshes; isubmesh++){
        strcpy(pumi_inputs->type_flag[isubmesh], each_submesh_flag[isubmesh]);
        *(pumi_inputs->p1_i + isubmesh) = atoi( each_p1_submesh[isubmesh]);
        NumberDebyeLengthsInDomain += *(pumi_inputs->p1_i + isubmesh);
        *(pumi_inputs->p2max_i + isubmesh) = atoi( each_p2max_submesh[isubmesh]);
        *(pumi_inputs->p2min_i + isubmesh) = atoi( each_p2min_submesh[isubmesh]);

        // calculating all pumi_inputs to initiate the mesh
        if (isubmesh == 0){
            *(pumi_inputs->x_left + isubmesh) = x1_min;
            *(pumi_inputs->x_right + isubmesh) = *(pumi_inputs->x_left + isubmesh) + lambda_D* *(pumi_inputs->p1_i + isubmesh);
        }
        else {
            *(pumi_inputs->x_left + isubmesh) = *(pumi_inputs->x_right + (isubmesh-1));
            *(pumi_inputs->x_right + isubmesh) = *(pumi_inputs->x_left + isubmesh) + lambda_D* *(pumi_inputs->p1_i + isubmesh);
        }

        char tmp_typeflag_string[SUBMESH_FLAGSTRING_LENGTH];
        strcpy(tmp_typeflag_string, pumi_inputs->type_flag[isubmesh]);
        unsigned int typeflag = pumi_getsubmeshflag(tmp_typeflag_string);

        double subdomain_L = *(pumi_inputs->x_right + isubmesh) - *(pumi_inputs->x_left + isubmesh);

        // Initialize submesh parameters to be 0

        double left_r = 0.0;
        double left_T = 0.0;
        double left_t0 = 0.0;
        int left_Nel = 0;

        double right_r = 0.0;
        double right_T = 0.0;
        double right_t0 = 0.0;
        int right_Nel = 0;

        double uniform_L = 0.0;
        double uniform_dx1 = 0.0;
        int uniform_Nel = 0;

        if (typeflag & leftBL){// set values if leftBL flag is active
          left_T = *(pumi_inputs->p1_i + isubmesh)*lambda_D;
          double left_tmax = lambda_D/(*(pumi_inputs->p2max_i + isubmesh));
          left_t0 = lambda_D/(*(pumi_inputs->p2min_i + isubmesh));
          left_r = (left_T-left_t0)/(left_T-left_tmax);
          left_Nel = 1 + floor( log(left_tmax/left_t0)/log(left_r) );
          left_r = pumi_compute_grading_ratio(*(pumi_inputs->p1_i + isubmesh), *(pumi_inputs->p2max_i + isubmesh), left_Nel);
        }

        if (typeflag & rightBL){// set values if leftBL flag is active
          right_T = *(pumi_inputs->p1_i + isubmesh)*lambda_D;
          double right_tmax = lambda_D/(*(pumi_inputs->p2max_i + isubmesh));
          right_t0 = lambda_D/(*(pumi_inputs->p2min_i + isubmesh));
          right_r = (right_T-right_t0)/(right_T-right_tmax);
          right_Nel = 1 + floor( log(right_tmax/right_t0)/log(right_r) );
          right_r = pumi_compute_grading_ratio(*(pumi_inputs->p1_i + isubmesh), *(pumi_inputs->p2max_i + isubmesh), right_Nel);
        }

        if (typeflag & uniform){
          uniform_L = *(pumi_inputs->p1_i + isubmesh)*lambda_D;
          uniform_dx1 = lambda_D/(*(pumi_inputs->p2max_i + isubmesh));
          uniform_Nel = floor (uniform_L/uniform_dx1);
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

    int Nel_total = pumi_total_elements(pumi_mesh);

    int num_particles_per_debyelength;
    printf("\nEnter number of particles per Debye Length : ");
    scanf("%d", &num_particles_per_debyelength); //user supplied
    int num_particles = num_particles_per_debyelength*NumberDebyeLengthsInDomain;
    double *coordinates = (double*) malloc(num_particles * sizeof(double));
    int *particle_isactive = (int*) malloc(num_particles * sizeof(int));

    srand48(time(NULL));
    int iparticle=0;
    for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes; isubmesh++){
        double submesh_xleft = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_left;
        double submesh_xright = ((pumi_submesh1D_t*) pumi_mesh->pumi_submeshes + isubmesh)->x_right;
        double submesh_L = submesh_xright-submesh_xleft;
        double deblenth = submesh_L/lambda_D;
        int submesh_debyelengths = (int) deblenth;
        int numparticle_submesh = num_particles_per_debyelength*submesh_debyelengths;
        int j;
        for(j=0; j<numparticle_submesh;j++){
            coordinates[iparticle] = submesh_xleft + submesh_L*drand48();
            //particle_isactive[iparticle] = isubmesh+1;
            particle_isactive[iparticle] = pumi_locate_submesh_1D(pumi_mesh, coordinates[iparticle]);
            iparticle++;
        }
    }

    double *grid_weights = (double*) malloc((Nel_total+1)*sizeof(double));
    int i;
    for (i=0; i<= Nel_total; i++){
      grid_weights[i] = 0.0; //intialize charge to be zero at all nodes
    }

    for (iparticle=0; iparticle<num_particles; iparticle++){ //loop over all particles
      int kcell;
      double Wgh1, Wgh2;
      pumiMBBL_locatepoint_1D(pumi_mesh, coordinates[iparticle], particle_isactive[iparticle], &kcell, &Wgh2); // computes paricle cell and weight (based on linear weighting)
      Wgh1 = 1.0 - Wgh2;
      grid_weights[kcell]   += Wgh1;
      grid_weights[kcell+1] += Wgh2; //accumulate the weights for each particle
    }
    double sumweights=0;
    for (i=0; i<(Nel_total+1); i++){
      printf("node=%d   gridweight=%lf\n", i, grid_weights[i]);
      sumweights += grid_weights[i];
    }

    printf("Total weight = %lf\n\n", sumweights );
    free(coordinates);
    free(particle_isactive);
    free(grid_weights);


    pumi_finalize(pumi_mesh); //deallocates pumi_mesh object


    return 0;
}