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
      printf("\t \"p1_i\"  \t\t Number of Debye Lengths in i-th submesh \n\n");
      printf("\t \"Nel_i\" \t\t Number of elements in i-th submesh \n");
      printf("\t \"p2_min_i\"  \t\t For leftBL/rightBL, Number of minimum size cells in a Debye Length for i-th submesh \n");
      printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
      printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
      printf("  E.g.\n\n");
      printf("    ./install/bin/pumiMBBL1D_Demo 3 \"leftBL,uniform,rightBL\" \"20,60,20\" \"12,23,12\" \"5,0,5\"\n");
      exit(0);
    }


    //pumi_inputs = malloc(sizeof(pumi_initiate_input_t));
    int nsubmesh = atoi( argv[1] );
    pumi_inputs = pumi_inputs_allocate(nsubmesh);
    pumi_inputs->ndim = 1; // Fixed pumi input
    pumi_inputs->nsubmeshes = nsubmesh;

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
    char all_Nel_submesh[MAX_SUBMESHES*4];
    char each_Nel_submesh[MAX_SUBMESHES][4];
    strcpy(all_Nel_submesh, argv[4]);

    tok = strtok(all_Nel_submesh, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_Nel_submesh[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes){
        printf("ERROR: Number of Nel_i arguments not equal to number of submeshes...\n");
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

    //pumi_inputs_allocate(pumi_inputs, pumi_inputs->nsubmeshes);

    double lambda_D = 4.7016e-05;
    double x1_min = 0.0;
    int NumberDebyeLengthsInDomain=0;
    for(isubmesh=0; isubmesh<pumi_inputs->nsubmeshes; isubmesh++){
        strcpy(pumi_inputs->type_flag[isubmesh], each_submesh_flag[isubmesh]);
        *(pumi_inputs->p1_i + isubmesh) = atoi( each_p1_submesh[isubmesh]);
        NumberDebyeLengthsInDomain += *(pumi_inputs->p1_i + isubmesh);
        *(pumi_inputs->Nel_i + isubmesh) = atoi( each_Nel_submesh[isubmesh]);
        *(pumi_inputs->p2min_i + isubmesh) = atof( each_p2min_submesh[isubmesh]);

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
        double uniform_dx11 = 0.0;
        int uniform_Nel = 0;

        if (typeflag & leftBL){// set values if leftBL flag is active
          left_T = *(pumi_inputs->p1_i + isubmesh)*lambda_D;
          left_t0 = lambda_D/(*(pumi_inputs->p2min_i + isubmesh));
          left_Nel = *(pumi_inputs->Nel_i + isubmesh);
          left_r = pumi_compute_grading_ratio_new(left_T, left_t0, left_Nel);
        }

        if (typeflag & rightBL){// set values if leftBL flag is active
            right_T = *(pumi_inputs->p1_i + isubmesh)*lambda_D;
            right_t0 = lambda_D/(*(pumi_inputs->p2min_i + isubmesh));
            right_Nel = *(pumi_inputs->Nel_i + isubmesh);
            right_r = pumi_compute_grading_ratio_new(right_T, right_t0, right_Nel);
        }

        if (typeflag & uniform){
          uniform_L = *(pumi_inputs->p1_i + isubmesh)*lambda_D;
          uniform_Nel = *(pumi_inputs->Nel_i + isubmesh);
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


    pumi_initiate_mesh_options_t pumi_initiate_options;
    pumi_initiate_options.BL_cache_flag = pumi_cache_BL_elemsize_OFF;
    //pumi_initiate_options.nodeoffset_cache_flag = nodeoffset_caching_flag;

    // the pumi_input object NEEDS TO BE POPULATED before initializing pumi_mesh
    pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs, pumi_initiate_options);

    // deallocate memory allocated to pumi_inputs -- Always do this IMMEDIATELY AFTER pumi_initiate()
    pumi_inputs_deallocate(pumi_inputs);

    int Nel_total = pumi_total_elements(pumi_mesh);

    int num_particles_per_debyelength;
    printf("\nEnter number of particles per Debye Length : ");
    scanf("%d", &num_particles_per_debyelength); //user supplied
    int num_particles = num_particles_per_debyelength*NumberDebyeLengthsInDomain;
    double *coordinates = (double*) malloc(num_particles * sizeof(double));
    //int *particle_isactive = (int*) malloc(num_particles * sizeof(int));

    srand48(time(NULL));
    int iparticle=0;
    for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        double submesh_x1left = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_min;
        double submesh_x1right = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->coord_max;
        double submesh_L = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->submesh_T;
        double deblenth = submesh_L/lambda_D;
        int submesh_debyelengths = (int) deblenth;
        int numparticle_submesh = num_particles_per_debyelength*submesh_debyelengths;
        int j;
        for(j=0; j<numparticle_submesh;j++){
            coordinates[iparticle] = submesh_x1left + submesh_L*drand48();
            //particle_isactive[iparticle] = isubmesh+1;
            //particle_isactive[iparticle] = pumi_locate_submesh_1D(pumi_mesh, coordinates[iparticle]);
            iparticle++;
        }
    }

    /*for(isubmesh=0; isubmesh<pumi_mesh->nsubmeshes_x1; isubmesh++){
        int submesh_total_Nel, icell;
        double left_node, right_node, elem_size, err;
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & leftBL){
            submesh_total_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->left_Nel;
            for(icell=0; icell<submesh_total_Nel; icell++){
                pumi_calc_node_coords(pumi_mesh, isubmesh, icell, &left_node, &right_node);
                elem_size = pumi_calc_elem_size(pumi_mesh, isubmesh, icell);
                err = right_node-left_node-elem_size;
                printf("icell=%d left_node=%2.8e right_node=%2.8e err=%2.8e\n", icell, left_node, right_node, err);
            }
            printf("\n\n");
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & uniform){
            submesh_total_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->uniform_Nel;
            for(icell=0; icell<submesh_total_Nel; icell++){
                pumi_calc_node_coords(pumi_mesh, isubmesh, icell, &left_node, &right_node);
                elem_size = pumi_calc_elem_size(pumi_mesh, isubmesh, icell);
                err = right_node-left_node-elem_size;
                printf("icell=%d left_node=%2.8e right_node=%2.8e err=%2.8e\n", icell, left_node, right_node, err);
            }
            printf("\n\n");
        }
        if (((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->pumi_flag & rightBL){
            submesh_total_Nel = ((pumi_submesh_t*) pumi_mesh->pumi_submeshes_x1 + isubmesh)->right_Nel;
            for(icell=0; icell<submesh_total_Nel; icell++){
                pumi_calc_node_coords(pumi_mesh, isubmesh, icell, &left_node, &right_node);
                elem_size = pumi_calc_elem_size(pumi_mesh, isubmesh, icell);
                err = right_node-left_node-elem_size;
                printf("icell=%d left_node=%2.8e right_node=%2.8e err=%2.8e\n", icell, left_node, right_node, err);
            }
            printf("\n\n");
        }
    }*/

    double *grid_weights = (double*) malloc((Nel_total+1)*sizeof(double));
    int i;
    for (i=0; i<= Nel_total; i++){
      grid_weights[i] = 0.0; //intialize charge to be zero at all nodes
    }
    int kcell;
    double Wgh1, Wgh2;
    int submesh_icell;
    for (iparticle=0; iparticle<num_particles; iparticle++){ //loop over all particles
      //isubmesh = particle_isactive[iparticle];
      //pumiMBBL_locatepoint_1D(pumi_mesh, coordinates[iparticle], isubmesh, &kcell, &Wgh2); // computes paricle cell and weight (based on linear weighting)
      //pumi_locate_function[isubmesh](pumi_mesh, isubmesh, coordinates[iparticle], &kcell, &Wgh2);
      pumi_locate_submesh_and_cell(pumi_mesh, coordinates[iparticle], &isubmesh, &submesh_icell, pumi_x1);
      //printf("particle %d located at submesh %d and local cell %d\n",iparticle, isubmesh, submesh_icell );
      pumi_calc_weights(pumi_mesh, isubmesh, submesh_icell, coordinates[iparticle], &kcell, &Wgh2, pumi_x1);
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
    //free(particle_isactive);
    free(grid_weights);

    pumi_finalize(pumi_mesh); //deallocates pumi_mesh object

    return 0;
}
