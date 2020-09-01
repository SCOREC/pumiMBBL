#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "pumiMBBL.h"

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

void write2file_uni(double *field, int size, int ID){
    FILE *field_fptr;
    char field_file[30];
    sprintf(field_file,"Field_uni_%d.txt",ID);
    field_fptr = fopen(field_file,"w");
    int i;
    for (i=0; i<size; i++ ){
        fprintf(field_fptr, "%.16e\n", field[i] );
    }
}

void std_hpic_calc(double q0, double q1, double dx1, double dx2, int Nel_x1, double *Wgh2_x1, double *Wgh2_x2, int *node1, int *node3){
    int icell, jcell, kcell;
    icell = floor(q0/dx1);
    jcell = floor(q1/dx2);
    kcell = jcell*Nel_x1 + icell;
    *node1 = kcell+jcell;
    *node3 = kcell+jcell+Nel_x1+1;

    *Wgh2_x1 = (q0-icell*dx1)/dx1;
    *Wgh2_x2 = (q1-jcell*dx2)/dx2;
}

#define NUMSTEPS 10
int main(int argc, char *argv[])
{
    clock_t total_run_time;
    total_run_time = clock();
    pumi_initiate_input_t    *pumi_inputs;

    ///*
    if (argc != 12){
      printf("Execute the code with the following command line arguments -- \n\n" );
      printf("\t ./install/bin/pumiMBBL2D_Demo N_x1 \"typeflag_i_x1\" \"p1_i_x1\" \"Nel_i_x1\" \"p2min_i_x1\" N_x2 \"typeflag_i_x2\" \"p1_i_x2\" \"Nel_i_x2\" \"p2min_i_x2\" \"block_isactive\"\n\n\n");
      printf("\t N_x1     \t\t Total Number of submeshes along the x1-direction \n");
      printf("\t \"typeflag_i_x1\" \t Active mesh type segment in i-th submesh along the x1-direction \n" );
      printf("\t \"p1_i_x1\"  \t\t Number of Debye Lengths in i-th submesh along the x1-direction \n");
      printf("\t \"Nel_i_x1\" \t\t Number of elements in i-th submesh along the x1-direction \n");
      printf("\t \"p2min_i_x1\"  \t\t For leftBL/rightBL, Number of minimum size cells in a Debye Length for i-th submesh along the x1-direction \n");
      printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
      printf("\t N_x2     \t\t Total Number of submeshes along the x2-direction \n");
      printf("\t \"typeflag_i_x2\" \t Active mesh type segment in i-th submesh along the x2-direction \n" );
      printf("\t \"p1_i_x2\"  \t\t Number of Debye Lengths in i-th submesh along the x2-direction \n");
      printf("\t \"Nel_i_x2\" \t\t Number of elements in i-th submesh along the x2-direction \n");
      printf("\t \"p2min_i_x2\"  \t\t For bottomBL/topBL, Number of minimum size cells in a Debye Length for i-th submesh along the x2-direction \n");
      printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
      printf("\t block_isactive \t Activity info of each submesh-block (N_x1*N_x2 inputs required)\n" );
      printf("\t \t  \t\t 0 is inactive \n\n");
      printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
      printf("  E.g.#1\n\n");
      printf("    ./install/bin/pumiMBBL2D_Demo 4 \"leftBL,uniform,uniform,rightBL\" \"20,30,30,20\" \"10,10,10,10\" \"1.0,0,0,1.0\" 3 \"bottomBL,uniform,topBL\" \"10,30,10\" \"5,10,5\" \"1.0,0,1.0\" \"1,0,0,1,1,1,1,1,1,0,0,1\"\n\n");
      printf("  E.g.#2\n\n");
      printf("    ./install/bin/pumiMBBL2D_Demo 3 \"leftBL,uniform,rightBL\" \"30,90,30\" \"15,25,15\" \"1.0,0.0,1.0\" 2 \"topBL,bottomBL\" \"50,50\" \"25,25\" \"1.0,1.0\" \"1,1,1,1,1,1\"\n\n");
      exit(0);
    }

    int nsubmesh_x1 = atoi( argv[1] );
    int nsubmesh_x2 = atoi( argv[6] );
    int nsubmesh = nsubmesh_x1+nsubmesh_x2;
    pumi_inputs = pumi_inputs_allocate(nsubmesh);
    pumi_inputs->ndim = 2; // Fixed pumi input
    pumi_inputs->nsubmeshes_x1 = nsubmesh_x1;
    pumi_inputs->nsubmeshes_x2 = nsubmesh_x2;
    // reading submesh meshtypes
    char all_submesh_flag_x1[MAX_SUBMESHES*SEGMENT_STRING_LENGTH];
    char each_submesh_flag_x1[MAX_SUBMESHES][SEGMENT_STRING_LENGTH];
    strcpy(all_submesh_flag_x1, argv[2]);

    char *tok = strtok(all_submesh_flag_x1, ",");
    int isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_flag_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x1){
        printf("ERROR: Number of typeflag_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh lengths
    char all_p1_submesh_x1[MAX_SUBMESHES*5];
    char each_p1_submesh_x1[MAX_SUBMESHES][5];
    strcpy(all_p1_submesh_x1, argv[3]);

    tok = strtok(all_p1_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x1){
        printf("ERROR: Number of p1_i_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh max elemsize
    char all_Nel_submesh_x1[MAX_SUBMESHES*4];
    char each_Nel_submesh_x1[MAX_SUBMESHES][4];
    strcpy(all_Nel_submesh_x1, argv[4]);

    tok = strtok(all_Nel_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_Nel_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x1){
        printf("ERROR: Number of Nel_i_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh min elemsize
    char all_p2min_submesh_x1[MAX_SUBMESHES*4];
    char each_p2min_submesh_x1[MAX_SUBMESHES][4];
    strcpy(all_p2min_submesh_x1, argv[5]);

    tok = strtok(all_p2min_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x1){
        printf("ERROR: Number of p2max_i_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    // reading submesh meshtypes
    char all_submesh_flag_x2[MAX_SUBMESHES*SEGMENT_STRING_LENGTH];
    char each_submesh_flag_x2[MAX_SUBMESHES][SEGMENT_STRING_LENGTH];
    strcpy(all_submesh_flag_x2, argv[7]);

    tok = strtok(all_submesh_flag_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_flag_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x2){
        printf("ERROR: Number of typeflag_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh lengths
    char all_p1_submesh_x2[MAX_SUBMESHES*5];
    char each_p1_submesh_x2[MAX_SUBMESHES][5];
    strcpy(all_p1_submesh_x2, argv[8]);

    tok = strtok(all_p1_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x2){
        printf("ERROR: Number of p1_i_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh max elemsize
    char all_Nel_submesh_x2[MAX_SUBMESHES*4];
    char each_Nel_submesh_x2[MAX_SUBMESHES][4];
    strcpy(all_Nel_submesh_x2, argv[9]);

    tok = strtok(all_Nel_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_Nel_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x2){
        printf("ERROR: Number of Nel_i_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh min elemsize
    char all_p2min_submesh_x2[MAX_SUBMESHES*4];
    char each_p2min_submesh_x2[MAX_SUBMESHES][4];
    strcpy(all_p2min_submesh_x2, argv[10]);

    tok = strtok(all_p2min_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x2){
        printf("ERROR: Number of p2max_i_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh activity info
    char all_submesh_isactive[MAX_SUBMESHES*2];
    char each_submesh_isactive[MAX_SUBMESHES][2];
    strcpy(all_submesh_isactive, argv[11]);

    tok = strtok(all_submesh_isactive, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_isactive[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmeshes_x2*pumi_inputs->nsubmeshes_x1){
        printf("ERROR: Number of block_isactive arguments not equal to number of submeshes...\n");
        exit(0);
    }


    double lambda_D = 1.0;
    double x1_min = 0.0;
    int NumberDebyeLengthsInDomain_x1=0;
    for(isubmesh=0; isubmesh<pumi_inputs->nsubmeshes_x1; isubmesh++){
        strcpy(pumi_inputs->type_flag[isubmesh], each_submesh_flag_x1[isubmesh]);
        *(pumi_inputs->p1_i_x1 + isubmesh) = atoi( each_p1_submesh_x1[isubmesh]);
        NumberDebyeLengthsInDomain_x1 += *(pumi_inputs->p1_i_x1 + isubmesh);
        *(pumi_inputs->Nel_i_x1 + isubmesh) = atoi( each_Nel_submesh_x1[isubmesh]);
        *(pumi_inputs->p2min_i_x1 + isubmesh) = atof( each_p2min_submesh_x1[isubmesh]);

        // calculating all pumi_inputs to initiate the mesh
        if (isubmesh == 0){
            *(pumi_inputs->x_left + isubmesh) = x1_min;
            *(pumi_inputs->x_right + isubmesh) = *(pumi_inputs->x_left + isubmesh) + lambda_D* *(pumi_inputs->p1_i_x1 + isubmesh);
        }
        else {
            *(pumi_inputs->x_left + isubmesh) = *(pumi_inputs->x_right + (isubmesh-1));
            *(pumi_inputs->x_right + isubmesh) = *(pumi_inputs->x_left + isubmesh) + lambda_D* *(pumi_inputs->p1_i_x1 + isubmesh);
        }

        char tmp_typeflag_string[SUBMESH_FLAGSTRING_LENGTH];
        strcpy(tmp_typeflag_string, pumi_inputs->type_flag[isubmesh]);
        unsigned int typeflag = pumi_getsubmeshflag(tmp_typeflag_string);

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
        int uniform_Nel = 0;

        if (typeflag & leftBL){// set values if leftBL flag is active
          left_T = *(pumi_inputs->p1_i_x1 + isubmesh)*lambda_D;
          left_t0 = lambda_D/(*(pumi_inputs->p2min_i_x1 + isubmesh));
          left_Nel = *(pumi_inputs->Nel_i_x1 + isubmesh);
          left_r = pumi_compute_grading_ratio_new(left_T, left_t0, left_Nel);
        }

        if (typeflag & rightBL){// set values if leftBL flag is active
            right_T = *(pumi_inputs->p1_i_x1 + isubmesh)*lambda_D;
            right_t0 = lambda_D/(*(pumi_inputs->p2min_i_x1 + isubmesh));
            right_Nel = *(pumi_inputs->Nel_i_x1 + isubmesh);
            right_r = pumi_compute_grading_ratio_new(right_T, right_t0, right_Nel);
        }

        if (typeflag & uniform){
          uniform_L = *(pumi_inputs->p1_i_x1 + isubmesh)*lambda_D;
          uniform_Nel = *(pumi_inputs->Nel_i_x1 + isubmesh);
        }

        // all pumi_inputs are calculated
        *(pumi_inputs->uniform_Nel_x1 + isubmesh) = uniform_Nel;
        *(pumi_inputs->left_T + isubmesh)      = left_T;
        *(pumi_inputs->left_r + isubmesh)      = left_r;
        *(pumi_inputs->left_Nel + isubmesh)    = left_Nel;
        *(pumi_inputs->right_T + isubmesh)     = right_T;
        *(pumi_inputs->right_r + isubmesh)     = right_r;
        *(pumi_inputs->right_Nel + isubmesh)   = right_Nel;
    }

    double x2_min = 0.0;
    int NumberDebyeLengthsInDomain_x2=0;
    int jsubmesh;
    for(jsubmesh=0; jsubmesh<pumi_inputs->nsubmeshes_x2; jsubmesh++){
        isubmesh = jsubmesh + pumi_inputs->nsubmeshes_x1;
        strcpy(pumi_inputs->type_flag[isubmesh], each_submesh_flag_x2[jsubmesh]);
        *(pumi_inputs->p1_i_x2 + isubmesh) = atoi( each_p1_submesh_x2[jsubmesh]);
        NumberDebyeLengthsInDomain_x2 += *(pumi_inputs->p1_i_x2 + isubmesh);
        *(pumi_inputs->Nel_i_x2 + isubmesh) = atoi( each_Nel_submesh_x2[jsubmesh]);
        *(pumi_inputs->p2min_i_x2 + isubmesh) = atof( each_p2min_submesh_x2[jsubmesh]);

        // calculating all pumi_inputs to initiate the mesh
        if (jsubmesh == 0){
            *(pumi_inputs->y_bottom + isubmesh) = x2_min;
            *(pumi_inputs->y_top + isubmesh) = *(pumi_inputs->y_bottom + isubmesh) + lambda_D* *(pumi_inputs->p1_i_x2 + isubmesh);
        }
        else {
            *(pumi_inputs->y_bottom + isubmesh) = *(pumi_inputs->y_top + (isubmesh-1));
            *(pumi_inputs->y_top + isubmesh) = *(pumi_inputs->y_bottom + isubmesh) + lambda_D* *(pumi_inputs->p1_i_x2 + isubmesh);
        }

        char tmp_typeflag_string[SUBMESH_FLAGSTRING_LENGTH];
        strcpy(tmp_typeflag_string, pumi_inputs->type_flag[isubmesh]);
        unsigned int typeflag = pumi_getsubmeshflag(tmp_typeflag_string);

        // Initialize submesh parameters to be 0

        double bottom_r = 0.0;
        double bottom_T = 0.0;
        double bottom_t0 = 0.0;
        int bottom_Nel = 0;

        double top_r = 0.0;
        double top_T = 0.0;
        double top_t0 = 0.0;
        int top_Nel = 0;

        double uniform_L = 0.0;
        int uniform_Nel = 0;

        if (typeflag & bottomBL){// set values if leftBL flag is active
          bottom_T = *(pumi_inputs->p1_i_x2 + isubmesh)*lambda_D;
          bottom_t0 = lambda_D/(*(pumi_inputs->p2min_i_x2 + isubmesh));
          bottom_Nel = *(pumi_inputs->Nel_i_x2 + isubmesh);
          bottom_r = pumi_compute_grading_ratio_new(bottom_T, bottom_t0, bottom_Nel);
        }

        if (typeflag & topBL){// set values if leftBL flag is active
            top_T = *(pumi_inputs->p1_i_x2 + isubmesh)*lambda_D;
            top_t0 = lambda_D/(*(pumi_inputs->p2min_i_x2 + isubmesh));
            top_Nel = *(pumi_inputs->Nel_i_x2 + isubmesh);
            top_r = pumi_compute_grading_ratio_new(top_T, top_t0, top_Nel);
        }

        if (typeflag & uniform){
          uniform_L = *(pumi_inputs->p1_i_x2 + isubmesh)*lambda_D;
          uniform_Nel = *(pumi_inputs->Nel_i_x2 + isubmesh);
        }

        // all pumi_inputs are calculated
        *(pumi_inputs->uniform_Nel_x2 + isubmesh) = uniform_Nel;
        *(pumi_inputs->bottom_T + isubmesh)      = bottom_T;
        *(pumi_inputs->bottom_r + isubmesh)      = bottom_r;
        *(pumi_inputs->bottom_Nel + isubmesh)    = bottom_Nel;
        *(pumi_inputs->top_T + isubmesh)     = top_T;
        *(pumi_inputs->top_r + isubmesh)     = top_r;
        *(pumi_inputs->top_Nel + isubmesh)   = top_Nel;
    }

    int ksubmesh=0;
    for(jsubmesh=0; jsubmesh<pumi_inputs->nsubmeshes_x2; jsubmesh++){
        for(isubmesh=0; isubmesh<pumi_inputs->nsubmeshes_x1; isubmesh++){
            int val = atoi(each_submesh_isactive[ksubmesh]);
            ksubmesh++;
            if (val){
                pumi_inputs->isactive[isubmesh][jsubmesh] = true;
            }
            else{
                pumi_inputs->isactive[isubmesh][jsubmesh] = false;
            }
        }
    }

    // the pumi_input object NEEDS TO BE POPULATED before initializing pumi_mesh
    pumi_initiate_mesh_options_t pumi_initiate_options;
    //pumi_initiate_options.BL_cache_flag = pumi_cache_BL_elemsize_OFF;
    //pumi_initiate_options.nodeoffset_cache_flag = pumi_cache_nodeoffset_ON;
    pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs, pumi_initiate_options);
    // deallocate memory allocated to pumi_inputs -- Always do this IMMEDIATELY AFTER pumi_initiate()
    pumi_inputs_deallocate(pumi_inputs);

    if (!(pumi_mesh_with_no_inactive_blocks(pumi_mesh))){
        printf("Particle locate/update not implemented for mesh with inactive blocks -- Terminating...\n");
        exit(0);
    }

    int num_particles;
    printf("\nEnter number of particles in domain : ");
    scanf("%d", &num_particles); //user supplied

    double **coords;
    coords = (double**) malloc(2 * sizeof(double *));
    coords[0] = (double*) malloc(num_particles * sizeof(double));
    coords[1] = (double*) malloc(num_particles * sizeof(double));

    int **part_isubmesh, **part_icell;
    part_isubmesh = (int**) malloc (2 * sizeof(int *));
    part_isubmesh[0] = (int*) malloc( num_particles * sizeof(int));
    part_isubmesh[1] = (int*) malloc( num_particles * sizeof(int));
    part_icell = (int**) malloc (2 * sizeof(int *));
    part_icell[0] = (int*) malloc( num_particles * sizeof(int));
    part_icell[1] = (int*) malloc( num_particles * sizeof(int));

    double *field;
    int Nnp_2D = pumi_total_nodes(pumi_mesh);
    field = (double *) malloc(Nnp_2D * sizeof(double));
    int inp, jnp;
    for (inp=0; inp<Nnp_2D; inp++){
        field[inp] = 0.0;
    }

    x1_min = pumi_global_x1_min(pumi_mesh);
    double x1_max = pumi_global_x1_max(pumi_mesh);
    x2_min = pumi_global_x2_min(pumi_mesh);
    double x2_max = pumi_global_x2_max(pumi_mesh);

    //printf("x1_min = %2.8f  x1_max = %2.8f\n", x1_min, x1_max);
    //printf("x2_min = %2.8f  x2_max = %2.8f\n", x2_min, x2_max);
    // particle initiate
    srand48(time(NULL));
    int iparticle, icell, jcell, kcell_x1, kcell_x2, kcell, node1, node3;
    double Wgh1_x1, Wgh2_x1, Wgh1_x2, Wgh2_x2;
    clock_t time_pumi, time_hpic;
    clock_t time_pumi_loop, time_pumi_loop_var;
    time_pumi_loop = 0.0;
    int istep, nodeID;

    time_pumi = clock();
    for(iparticle=0; iparticle<num_particles; iparticle++){
        coords[0][iparticle] = (0.9*x1_min + 0.1*x1_max) + 0.25*(x1_max-x1_min)*drand48();
        coords[1][iparticle] = (0.9*x2_min + 0.1*x2_max) + 0.25*(x2_max-x2_min)*drand48();

        double q0 = coords[0][iparticle];
        double q1 = coords[1][iparticle];

        pumi_locate_submesh_and_cell(pumi_mesh, q0, &isubmesh, &icell, pumi_x1);
        part_isubmesh[0][iparticle] = isubmesh;
        part_icell[0][iparticle] = icell;

        pumi_locate_submesh_and_cell(pumi_mesh, q1, &jsubmesh, &jcell, pumi_x2);
        part_isubmesh[1][iparticle] = jsubmesh;
        part_icell[1][iparticle] = jcell;

        pumi_calc_weights(pumi_mesh, isubmesh, icell, q0, &kcell_x1, &Wgh2_x1, pumi_x1);
        Wgh1_x1 = 1.0 - Wgh2_x1;

        pumi_calc_weights(pumi_mesh, jsubmesh, jcell, q1, &kcell_x2, &Wgh2_x2, pumi_x2);
        Wgh1_x2 = 1.0 - Wgh2_x2;

        kcell = pumi_calc_elementID_and_nodeID(pumi_mesh, isubmesh, jsubmesh, kcell_x1, kcell_x2, &node1, &node3);

        field[node1]   += Wgh1_x1*Wgh1_x2;
        field[node1+1] += Wgh2_x1*Wgh1_x2;
        field[node3]   += Wgh1_x1*Wgh2_x2;
        field[node3+1] += Wgh2_x1*Wgh2_x2;
    }

    nodeID = 0;
    for (jnp=0; jnp<pumi_mesh->pumi_Nnp_total_x2; jnp++){
        for (inp=0; inp<pumi_mesh->pumi_Nnp_total_x1; inp++){
            if (pumi_is_node_active(pumi_mesh, inp, jnp)){
                field[nodeID] /= pumi_return_covolume_2D(pumi_mesh, inp, jnp);
                nodeID++;
            }
        }
    }

    write2file(field,Nnp_2D,0);

    //particle push
    for (istep=1; istep<=NUMSTEPS; istep++){
        printf("TIME STEP %d\n", istep);
        for (inp=0; inp<Nnp_2D; inp++){
            field[inp] = 0.0;
        }
        time_pumi_loop_var = clock();
        for (iparticle=0; iparticle<num_particles; iparticle++){
            double q0 = coords[0][iparticle];
            double q1 = coords[1][iparticle];

            isubmesh = part_isubmesh[0][iparticle];
            jsubmesh = part_isubmesh[1][iparticle];
            icell = part_icell[0][iparticle];
            jcell = part_icell[1][iparticle];

            pumi_calc_weights(pumi_mesh, isubmesh, icell, q0, &kcell_x1, &Wgh2_x1, pumi_x1);
            Wgh1_x1 = 1.0 - Wgh2_x1;

            pumi_calc_weights(pumi_mesh, jsubmesh, jcell, q1, &kcell_x2, &Wgh2_x2, pumi_x2);
            Wgh1_x2 = 1.0 - Wgh2_x2;

            kcell = pumi_calc_elementID_and_nodeID(pumi_mesh, isubmesh, jsubmesh, kcell_x1, kcell_x2, &node1, &node3);

            //Ex = E1[node1]*Wgh1_x1*Wgh1_x2 + E1[node1+1]*Wgh2_x1*Wgh1_x2 + E1[node3]*Wgh1_x1*Wgh2_x2 + E1[node3+1]*Wgh2_x1*Wgh2_x2;
            //Ey = E2[node1]*Wgh1_x1*Wgh1_x2 + E2[node1+1]*Wgh2_x1*Wgh1_x2 + E2[node3]*Wgh1_x1*Wgh2_x2 + E2[node3+1]*Wgh2_x1*Wgh2_x2;

            //coords[0][iparticle] += 3.0;
            //coords[1][iparticle] += 1.0;
            coords[0][iparticle] += 3.0*drand48();
            coords[1][iparticle] += 1.0*drand48();

            q0 = coords[0][iparticle];
            q1 = coords[1][iparticle];

            pumi_update_submesh_and_update_cell(pumi_mesh, q0, isubmesh, icell, &isubmesh, &icell, pumi_x1);
            part_isubmesh[0][iparticle] = isubmesh;
            part_icell[0][iparticle] = icell;

            pumi_update_submesh_and_update_cell(pumi_mesh, q1, jsubmesh, jcell, &jsubmesh, &jcell, pumi_x2);
            part_isubmesh[1][iparticle] = jsubmesh;
            part_icell[1][iparticle] = jcell;

            pumi_calc_weights(pumi_mesh, isubmesh, icell, q0, &kcell_x1, &Wgh2_x1, pumi_x1);
            Wgh1_x1 = 1.0 - Wgh2_x1;

            pumi_calc_weights(pumi_mesh, jsubmesh, jcell, q1, &kcell_x2, &Wgh2_x2, pumi_x2);
            Wgh1_x2 = 1.0 - Wgh2_x2;

            kcell = pumi_calc_elementID_and_nodeID(pumi_mesh, isubmesh, jsubmesh, kcell_x1, kcell_x2, &node1, &node3);

            field[node1]   += Wgh1_x1*Wgh1_x2;
            field[node1+1] += Wgh2_x1*Wgh1_x2;
            field[node3]   += Wgh1_x1*Wgh2_x2;
            field[node3+1] += Wgh2_x1*Wgh2_x2;

        }
        time_pumi_loop_var = clock() - time_pumi_loop_var;
        time_pumi_loop += time_pumi_loop_var;

        nodeID = 0;
        for (jnp=0; jnp<pumi_mesh->pumi_Nnp_total_x2; jnp++){
            for (inp=0; inp<pumi_mesh->pumi_Nnp_total_x1; inp++){
                if (pumi_is_node_active(pumi_mesh, inp, jnp)){
                    field[nodeID] /= pumi_return_covolume_2D(pumi_mesh, inp, jnp);
                    nodeID++;
                }
            }
        }
        write2file(field,Nnp_2D,istep);
    }
    time_pumi = clock() - time_pumi;

    int Nel_x1 = pumi_mesh->pumi_Nel_total_x1;
    int Nel_x2 = pumi_mesh->pumi_Nel_total_x2;
    double dx1 = (x1_max - x1_min)/(double) Nel_x1;
    double dx2 = (x2_max - x2_min)/(double) Nel_x2;
    for (inp=0; inp<Nnp_2D; inp++){
        field[inp] = 0.0;
    }

    time_hpic = clock();
    for(iparticle=0; iparticle<num_particles; iparticle++){
        coords[0][iparticle] = (0.9*x1_min + 0.1*x1_max) + 0.25*(x1_max-x1_min)*drand48();
        coords[1][iparticle] = (0.9*x2_min + 0.1*x2_max) + 0.25*(x2_max-x2_min)*drand48();

        double q0 = coords[0][iparticle];
        double q1 = coords[1][iparticle];

        icell = floor(q0/dx1);
        jcell = floor(q1/dx2);
        kcell = jcell*Nel_x1 + icell;

        Wgh2_x1 = (q0-icell*dx1)/dx1;
        Wgh1_x1 = 1.0 - Wgh2_x1;

        Wgh2_x2 = (q1-jcell*dx2)/dx2;
        Wgh1_x2 = 1.0 - Wgh2_x2;


        field[kcell+jcell]   += Wgh1_x1*Wgh1_x2;
        field[kcell+jcell+1] += Wgh2_x1*Wgh1_x2;
        field[kcell+jcell+Nel_x1+1] += Wgh1_x1*Wgh2_x2;
        field[kcell+jcell+Nel_x1+2] += Wgh2_x1*Wgh2_x2;
    }

    nodeID=0;
    double x1_cov, x2_cov;
    for (jnp=0; jnp<=Nel_x2; jnp++){
        for (inp=0; inp<=Nel_x1; inp++){
            if (inp==0 || inp==Nel_x1){
                x1_cov = dx1/2.0;
            }
            else{
                x1_cov = dx1;
            }
            if (jnp==0 || jnp==Nel_x2){
                x2_cov = dx2/2.0;
            }
            else{
                x2_cov = dx2;
            }
            field[nodeID] /= x1_cov*x2_cov;
            //field[nodeID] /= pumi_return_covolume_2D(pumi_mesh, inp, jnp);
            nodeID++;
        }
    }

    write2file_uni(field, Nnp_2D, 0);
    clock_t time_hpic_loop, time_hpic_loop_var;
    time_hpic_loop = 0.0;
    for (istep=1; istep<=NUMSTEPS; istep++){
        printf("TIME STEP %d\n", istep);
        for (inp=0; inp<Nnp_2D; inp++){
            field[inp] = 0.0;
        }
        time_hpic_loop_var = clock();
        for (iparticle=0; iparticle<num_particles; iparticle++){
            double q0 = coords[0][iparticle];
            double q1 = coords[1][iparticle];

            std_hpic_calc(q0, q1, dx1, dx2, Nel_x1, &Wgh2_x1, &Wgh2_x2, &node1, &node3);
            Wgh1_x1 = 1.0 - Wgh2_x1;
            Wgh1_x2 = 1.0 - Wgh2_x2;

            //Ex = E1[node1]*Wgh1_x1*Wgh1_x2 + E1[node1+1]*Wgh2_x1*Wgh1_x2 + E1[node3]*Wgh1_x1*Wgh2_x2 + E1[node3+1]*Wgh2_x1*Wgh2_x2;
            //Ey = E2[node1]*Wgh1_x1*Wgh1_x2 + E2[node1+1]*Wgh2_x1*Wgh1_x2 + E2[node3]*Wgh1_x1*Wgh2_x2 + E2[node3+1]*Wgh2_x1*Wgh2_x2;

            //coords[0][iparticle] += 3.0;
            //coords[1][iparticle] += 1.0;
            coords[0][iparticle] += 3.0*drand48();
            coords[1][iparticle] += 1.0*drand48();

            q0 = coords[0][iparticle];
            q1 = coords[1][iparticle];


            std_hpic_calc(q0, q1, dx1, dx2, Nel_x1, &Wgh2_x1, &Wgh2_x2, &node1, &node3);
            Wgh1_x1 = 1.0 - Wgh2_x1;
            Wgh1_x2 = 1.0 - Wgh2_x2;

            field[node1]   += Wgh1_x1*Wgh1_x2;
            field[node1+1] += Wgh2_x1*Wgh1_x2;
            field[node3]   += Wgh1_x1*Wgh2_x2;
            field[node3+1] += Wgh2_x1*Wgh2_x2;

        }
        time_hpic_loop_var = clock() - time_hpic_loop_var;
        time_hpic_loop += time_hpic_loop_var;
        nodeID=0;
        double x1_cov, x2_cov;
        for (jnp=0; jnp<=Nel_x2; jnp++){
            for (inp=0; inp<=Nel_x1; inp++){
                if (inp==0 || inp==Nel_x1){
                    x1_cov = dx1/2.0;
                }
                else{
                    x1_cov = dx1;
                }
                if (jnp==0 || jnp==Nel_x2){
                    x2_cov = dx2/2.0;
                }
                else{
                    x2_cov = dx2;
                }
                field[nodeID] /= x1_cov*x2_cov;
                //field[nodeID] /= pumi_return_covolume_2D(pumi_mesh, inp, jnp);
                nodeID++;
            }
        }
        write2file_uni(field, Nnp_2D, istep);
    }
    time_hpic = clock() - time_hpic;

    double t_p, t_h, t_t;
    t_p = ((double) time_pumi_loop)/CLOCKS_PER_SEC;
    t_h = ((double) time_hpic_loop)/CLOCKS_PER_SEC;

    printf("pumi = %2.8f s    hpic = %2.8f s       slowdown = %2.2fx\n",t_p, t_h, t_p/t_h );

    free(coords[0]);
    free(coords[1]);
    free(coords);

    free(part_isubmesh[0]);
    free(part_isubmesh[1]);
    free(part_isubmesh);

    free(part_icell[0]);
    free(part_icell[1]);
    free(part_icell);

    pumi_finalize(pumi_mesh);
    total_run_time = clock() - total_run_time;
    t_t = ((double) total_run_time)/CLOCKS_PER_SEC;
    printf("Total Run Time = %2.8f s\n",t_t );
    return 0;
}
