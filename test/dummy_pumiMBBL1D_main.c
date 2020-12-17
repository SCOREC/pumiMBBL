#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumiMBBL.h"

#define NUMSTEPS 10

void write2file(double *field, int size, int ID){
    FILE *field_fptr;
    char field_file[30];
    sprintf(field_file,"Field_%d.txt",ID);
    field_fptr = fopen(field_file,"w");
    int i;
    for (i=0; i<size; i++ ){
        fprintf(field_fptr, "%.16e\n", field[i] );
    }
}

int main(int argc, char *argv[])
{
    clock_t total_run_time;
    total_run_time = clock();
    pumi_initiate_input_t    *pumi_inputs;
    ///*
    if (argc < 7){
      printf("Execute the code with the following command line arguments -- \n\n" );
      printf("\t ./install/bin/pumiMBBL_Demo N \"typeflag_i\" \"p1_i\" \"p2_max_i\" \"p2_min_i\"\n\n\n");
      printf("\t N     \t\t\t Total Number of submeshes in the domain \n\n");
      printf("\t \"typeflag_i\" \t\t Active mesh type segment in i-th submesh\n\n" );
      printf("\t \"p1_i\"  \t\t Number of Debye Lengths in i-th submesh \n\n");
      printf("\t \"Nel_i\" \t\t Number of elements in i-th submesh \n");
      printf("\t \"p2_min_i\"  \t\t For leftBL/rightBL, Number of minimum size cells in a Debye Length for i-th submesh \n");
      printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
      printf("\t P     \t\t\t B-spline polynomial order for charge distirbution\n");
      printf("\t \t  \t\t Input will be ignored if bspline flag is turned off in main function\n\n");
      printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
      printf("  E.g.\n\n");
      printf("    ./install/bin/pumiMBBL1D_Demo 3 \"leftBL,uniform,rightBL\" \"20,60,20\" \"12,23,12\" \"5,0,5\" 2\n");
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

    int P_spline = atoi( argv[6] );
    pumi_inputs->P_spline = P_spline;
    //pumi_inputs_allocate(pumi_inputs, pumi_inputs->nsubmeshes);

    double lambda_D = 1.0;
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
    pumi_initiate_options.BL_cache_flag = pumi_cache_BL_elemsize_ON;
    //pumi_initiate_options.nodeoffset_cache_flag = nodeoffset_caching_flag;
    pumi_initiate_options.bspline_flag = pumi_bspline_ON;
    pumi_initiate_options.periodic_mesh_flag = pumi_periodic_mesh_ON;

    // the pumi_input object NEEDS TO BE POPULATED before initializing pumi_mesh
    pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs, pumi_initiate_options);

    // deallocate memory allocated to pumi_inputs -- Always do this IMMEDIATELY AFTER pumi_initiate()
    pumi_inputs_deallocate(pumi_inputs);
    int p = pumi_mesh->P_spline;
    int i,j,k;

    int num_particles;
    printf("\nEnter number of particles in domain : ");
    scanf("%d", &num_particles); //user supplied

    double *coords;
    coords = (double*) malloc(num_particles * sizeof(double));

    int *part_isubmesh, *part_icell;
    part_isubmesh = (int*) malloc( num_particles * sizeof(int));
    part_icell = (int*) malloc( num_particles * sizeof(int));

    double *field1, *field2;
    int Nnp = pumi_total_nodes(pumi_mesh);
    int Nel = pumi_total_elements(pumi_mesh);
    field1 = (double *) malloc(Nnp * sizeof(double));
    field2 = (double *) malloc(Nnp * sizeof(double));
    int inp, jnp;
    for (inp=0; inp<Nnp; inp++){
        field1[inp] = 0.0;
        field2[inp] = 0.0;
    }

    x1_min = pumi_global_x1_min(pumi_mesh);
    double x1_max = pumi_global_x1_max(pumi_mesh);
    if (pumi_mesh->bspline_flag){
        int iparticle, icell, kcell, node_left, node_right;
        double Q_macro_particle = 1.0;
        double Wgh1, Wgh2;
        int ispline;
        pumi_reset_Qspl_coeffs(pumi_mesh);
        // double delx_particle = (x1_max-0.0001)/(num_particles-1);
        for(iparticle=0; iparticle<num_particles; iparticle++){
            coords[iparticle] = (1.0*x1_min + 0.0*x1_max) + 1.0*(x1_max-x1_min)*drand48();
            // coords[iparticle] = iparticle*delx_particle;
            double q0 = coords[iparticle];

            pumi_locate_submesh_and_cell(pumi_mesh, q0, &isubmesh, &icell, pumi_x1);
            part_isubmesh[iparticle] = isubmesh;
            part_icell[iparticle] = icell;

            pumi_calc_weights(pumi_mesh, isubmesh, icell, q0, &kcell, &Wgh2, pumi_x1);
            Wgh1 = 1.0 - Wgh2;

            pumi_compute_Qspl_coeffs_periodic(pumi_mesh, Wgh2, kcell, Q_macro_particle);

            node_left = kcell;
            node_right = (kcell+1)%Nel;
            field2[node_left] += Q_macro_particle*Wgh1;
            field2[node_right] += Q_macro_particle*Wgh2;
        }
        pumi_compute_bspline_nodal_density_periodic(pumi_mesh, pumi_x1, field1);
        int inode[1];
        double cov;
        for (i=0; i<Nnp; i++){
            inode[0] = i;
            cov = pumi_return_covolume_periodic_1D(pumi_mesh, inode);
            field2[i] /= cov;
        }
        write2file(field1,pumi_mesh->pumi_Nnp_total_x1,1);
        write2file(field2,pumi_mesh->pumi_Nnp_total_x1,2);
    }
    else{
        int iparticle, icell, kcell, node_left, node_right;
        double Q_macro_particle = 1.0;
        double Wgh1, Wgh2;
        // double *qhat_spline = (double*) malloc(pumi_mesh->N_spline * sizeof(double));
        int ispline;

        for(iparticle=0; iparticle<num_particles; iparticle++){
            coords[iparticle] = (1.0*x1_min + 0.0*x1_max) + 1.0*(x1_max-x1_min)*drand48();
            // coords[iparticle] = iparticle*delx_particle;
            double q0 = coords[iparticle];

            pumi_locate_submesh_and_cell(pumi_mesh, q0, &isubmesh, &icell, pumi_x1);
            part_isubmesh[iparticle] = isubmesh;
            part_icell[iparticle] = icell;

            pumi_calc_weights(pumi_mesh, isubmesh, icell, q0, &kcell, &Wgh2, pumi_x1);
            Wgh1 = 1.0 - Wgh2;
            node_left = kcell;
            node_right = (kcell+1)%Nel;
            field2[node_left] += Q_macro_particle*Wgh1;
            field2[node_right] += Q_macro_particle*Wgh2;
        }

        int inode[1];
        double cov, gr;
        for (i=0; i<Nnp; i++){
            inode[0] = i;
            cov = pumi_return_covolume_periodic_1D(pumi_mesh, inode);
            gr = pumi_return_gradingratio_periodic(pumi_mesh, inode[0], pumi_x1);
            field2[i] /= cov;
            printf("Density[%3d] = %3.4f\tgr[%3d] = %3.4f\n",i,field2[i],i,gr );
        }
    }

    free(field1);
    free(field2);
    free(coords);
    free(part_isubmesh);
    free(part_icell);

    pumi_finalize(pumi_mesh);

    return 0;
}
