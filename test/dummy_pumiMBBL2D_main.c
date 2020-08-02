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
      printf("    ./install/bin/pumiMBBL2D_Demo 4 \"leftBL,uniform,uniform,rightBL\" \"20,30,30,20\" \"10,10,10,10\" \"1.0,0,0,1.0\" 3 \"bottomBL,uniform,topBL\" \"10,30,10\" \"5,10,5\" \"1.0,0,1.0\" \"1,1,1,1,1,1,1,1,1,1,1,1,\"\n\n");
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
    pumi_mesh_t *pumi_mesh = pumi_initiate(initiate_from_commandline_inputs, pumi_inputs, pumi_cache_BL_elemsize_ON);

    // deallocate memory allocated to pumi_inputs -- Always do this IMMEDIATELY AFTER pumi_initiate()
    pumi_inputs_deallocate(pumi_inputs);

    //pumi_finalize(pumi_mesh);
    return 0;
}
