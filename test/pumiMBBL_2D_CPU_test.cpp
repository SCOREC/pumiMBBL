#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_parsed_inputs(pumi::Mesh_Inputs *pumi_inputs);
void print_usage();
void write2file(pumi::ParticleDataCPU* hp, int N_part, int nstep);

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );
  {
      if (argc != 14)
      {
          print_usage();
      }
    pumi::check_is_pumi_working();
    int nsubmesh_x1 = atoi( argv[1] );
    double domain_x1_min = atof( argv[2] );
    int nsubmesh_x2 = atoi( argv[7] );
    double domain_x2_min = atof( argv[8] );
    // int nsubmesh = nsubmesh_x1+nsubmesh_x2;
    pumi::Mesh_Inputs *pumi_inputs = pumi::inputs_allocate();

    pumi_inputs->ndim = 2; // Fixed pumi input
    pumi_inputs->nsubmesh_x1 = nsubmesh_x1;
    pumi_inputs->nsubmesh_x2 = nsubmesh_x2;
    pumi_inputs->domain_x1_min = domain_x1_min;
    pumi_inputs->domain_x2_min = domain_x2_min;

    parse_inputs(argc, argv, pumi_inputs);


    pumi::Mesh_Options pumi_options;
    pumi_options.BL_storage_option = pumi::store_BL_coords_ON;
    pumi_options.print_node_option = pumi::print_node_coords_OFF;

    pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);

    printf("Mesh volume = %2.2f\n",pumi::get_mesh_volume(pumi_obj) );
    pumi::inputs_deallocate(pumi_inputs);

    pumi::print_mesh_skeleton(pumi_obj);

    int N_part = 1000;
    int N_step = 10;

    pumi::ParticleDataCPU* Partdata = new pumi::ParticleDataCPU[N_part];

    bool test0 = true;

    std::srand((unsigned)(std::time(nullptr)));

    double x1_min = pumi::get_global_x1_min_coord(pumi_obj);
    double x1_max = pumi::get_global_x1_max_coord(pumi_obj);
    double L_x1 = x1_max-x1_min;
    double x2_min = pumi::get_global_x2_min_coord(pumi_obj);
    double x2_max = pumi::get_global_x2_max_coord(pumi_obj);
    double L_x2 = x2_max-x2_min;
    double dist_factor = 20.0;
    int num_push = 0;

    if (test0){
        printf("push-test-0\n");
        for (int ipart=0; ipart<N_part; ipart++){
            pumi::Vector3 q = pumi::get_rand_point_in_mesh_host(pumi_obj);
            Partdata[ipart] = pumi::ParticleDataCPU(q[0],q[1]);
        }

        for (int ipart=0; ipart<N_part; ipart++){
            int isub, jsub, icell, jcell, submeshID, cellID;
            double q1 = Partdata[ipart].x1;
            double q2 = Partdata[ipart].x2;
            pumi::locate_submesh_and_cell_x1_host(pumi_obj, q1, &isub, &icell);
            pumi::locate_submesh_and_cell_x2_host(pumi_obj, q2, &jsub, &jcell);
            pumi::flatten_submeshID_and_cellID_host(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
            Partdata[ipart] = pumi::ParticleDataCPU(q1,q2,submeshID,cellID,true,-1);
        }

        write2file(Partdata, N_part, 0);

        for (int istep=0; istep<N_step; istep++){
            for (int p=0; p<N_part; p++){
                bool part_active = Partdata[p].part_active;
                num_push += part_active;
            }
            for (int ipart=0; ipart<N_part; ipart++){
                bool part_active = Partdata[ipart].part_active;
                if (part_active){
                    int isub, jsub, icell, jcell, kcell_x1, kcell_x2, submeshID, cellID;
                    int global_cell, topleft_node, bottomleft_node;
                    double q1 = Partdata[ipart].x1;
                    double q2 = Partdata[ipart].x2;
                    submeshID = Partdata[ipart].submeshID;
                    cellID = Partdata[ipart].cellID;

                    pumi::get_directional_submeshID_and_cellID_host(pumi_obj,submeshID,cellID,&isub,&icell,&jsub,&jcell);

                    double Wgh2_x1, Wgh2_x2, Wgh1_x1, Wgh1_x2;
                    pumi::calc_weights_x1_host(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                    pumi::calc_weights_x2_host(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                    Wgh1_x1 = 1.0-Wgh2_x1;
                    Wgh1_x2 = 1.0-Wgh2_x2;
                    pumi::calc_global_cellID_and_nodeID_fullmesh_host(pumi_obj, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                    double dq1 = L_x1/dist_factor;
                    double dq2 = -L_x2/dist_factor;

                    q1 += dq1;
                    q2 += dq2;

                    if (q1 < x1_min || q1 > x1_max || q2 < x2_min || q2 > x2_max){
                        submeshID=-1;
                        cellID=-1;
                        Partdata[ipart] = pumi::ParticleDataCPU(q1,q2,submeshID,cellID,false,-1);
                    }
                    else {
                        pumi::update_submesh_and_cell_x1_host(pumi_obj, q1, isub, icell, &isub, &icell);
                        pumi::update_submesh_and_cell_x2_host(pumi_obj, q2, jsub, jcell, &jsub, &jcell);
                        pumi::calc_weights_x1_host(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                        pumi::calc_weights_x2_host(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                        Wgh1_x1 = 1.0-Wgh2_x1;
                        Wgh1_x2 = 1.0-Wgh2_x2;
                        pumi::calc_global_cellID_and_nodeID_fullmesh_host(pumi_obj, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        pumi::flatten_submeshID_and_cellID_host(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
                        Partdata[ipart] = pumi::ParticleDataCPU(q1,q2,submeshID,cellID,true,-1,Wgh1_x1,Wgh2_x1,Wgh1_x2,Wgh2_x2);
                    }

                }
            }
            write2file(Partdata, N_part, istep+1);
        }
        printf("Total number of particle pushes executed in Test-0 = %d\n",num_push );
    }

  }
  Kokkos::finalize();

  return 0;
}

void parse_inputs(int , char* argv[], pumi::Mesh_Inputs *pumi_inputs)
{
    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int nsubmesh_x2 = pumi_inputs->nsubmesh_x2;
    // reading submesh meshtypes
    char all_submesh_flag_x1[MAX_SUBMESHES*10];
    char each_submesh_flag_x1[MAX_SUBMESHES][10];
    strcpy(all_submesh_flag_x1, argv[3]);

    char *tok = strtok(all_submesh_flag_x1, ",");
    int isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_flag_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of typeflag_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh lengths
    char all_p1_submesh_x1[MAX_SUBMESHES*10];
    char each_p1_submesh_x1[MAX_SUBMESHES][10];
    strcpy(all_p1_submesh_x1, argv[4]);

    tok = strtok(all_p1_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of block_length_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh max elemsize
    char all_p2max_submesh_x1[MAX_SUBMESHES*10];
    char each_p2max_submesh_x1[MAX_SUBMESHES][10];
    strcpy(all_p2max_submesh_x1, argv[5]);

    tok = strtok(all_p2max_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2max_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of Nel_i_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh min elemsize
    char all_p2min_submesh_x1[MAX_SUBMESHES*10];
    char each_p2min_submesh_x1[MAX_SUBMESHES][10];
    strcpy(all_p2min_submesh_x1, argv[6]);

    tok = strtok(all_p2min_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of max_elem_size_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading arbitrary block elemsize files
    char all_arb_submesh_x1[MAX_SUBMESHES*100];
    char each_arb_submesh_x1[MAX_SUBMESHES][100];
    strcpy(all_arb_submesh_x1, argv[7]);

    tok = strtok(all_arb_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_arb_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of elem-size-file arguments not equal to number of submeshes...\n");
        exit(0);
    }

    // reading submesh meshtypes
    char all_submesh_flag_x2[MAX_SUBMESHES*10];
    char each_submesh_flag_x2[MAX_SUBMESHES][10];
    strcpy(all_submesh_flag_x2, argv[10]);

    tok = strtok(all_submesh_flag_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_flag_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
        printf("ERROR: Number of typeflag_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh lengths
    char all_p1_submesh_x2[MAX_SUBMESHES*10];
    char each_p1_submesh_x2[MAX_SUBMESHES][10];
    strcpy(all_p1_submesh_x2, argv[11]);

    tok = strtok(all_p1_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
        printf("ERROR: Number of block_length_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh max elemsize
    char all_p2max_submesh_x2[MAX_SUBMESHES*10];
    char each_p2max_submesh_x2[MAX_SUBMESHES][10];
    strcpy(all_p2max_submesh_x2, argv[12]);

    tok = strtok(all_p2max_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2max_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
        printf("ERROR: Number of Nel_i_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh min elemsize
    char all_p2min_submesh_x2[MAX_SUBMESHES*10];
    char each_p2min_submesh_x2[MAX_SUBMESHES][10];
    strcpy(all_p2min_submesh_x2, argv[13]);

    tok = strtok(all_p2min_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
        printf("ERROR: Number of max_elem_size_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading arbitrary block elemsize files
    char all_arb_submesh_x2[MAX_SUBMESHES*100];
    char each_arb_submesh_x2[MAX_SUBMESHES][100];
    strcpy(all_arb_submesh_x2, argv[14]);

    tok = strtok(all_arb_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_arb_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
        printf("ERROR: Number of elem-size-file arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh activity info
    char all_submesh_isactive[MAX_SUBMESHES*2];
    char each_submesh_isactive[MAX_SUBMESHES][2];
    strcpy(all_submesh_isactive, argv[15]);

    tok = strtok(all_submesh_isactive, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_submesh_isactive[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2*pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of block_isactive arguments not equal to number of submeshes...\n");
        exit(0);
    }

    for (isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
        (pumi_inputs->meshtype_x1).push_back(each_submesh_flag_x1[isubmesh]);
        (pumi_inputs->block_length_x1).push_back(atof(each_p1_submesh_x1[isubmesh]));
        (pumi_inputs->max_elem_size_x1).push_back(atof(each_p2max_submesh_x1[isubmesh]));
        (pumi_inputs->min_elem_size_x1).push_back(atof(each_p2min_submesh_x1[isubmesh]));
        std::string x1_elemsize_file(each_arb_submesh_x1[isubmesh]);
        (pumi_inputs->arbitrary_x1_elemsize_file).push_back(x1_elemsize_file);
    }

    for (isubmesh=0; isubmesh<nsubmesh_x2; isubmesh++){
        (pumi_inputs->meshtype_x2).push_back(each_submesh_flag_x2[isubmesh]);
        (pumi_inputs->block_length_x2).push_back(atof(each_p1_submesh_x2[isubmesh]));
        (pumi_inputs->max_elem_size_x2).push_back(atof(each_p2max_submesh_x2[isubmesh]));
        (pumi_inputs->min_elem_size_x2).push_back(atof(each_p2min_submesh_x2[isubmesh]));
        std::string x2_elemsize_file(each_arb_submesh_x2[isubmesh]);
        (pumi_inputs->arbitrary_x2_elemsize_file).push_back(x2_elemsize_file);
    }

    int ksubmesh=0;
    for(int jsubmesh=0; jsubmesh<pumi_inputs->nsubmesh_x2; jsubmesh++){
        for(isubmesh=0; isubmesh<pumi_inputs->nsubmesh_x1; isubmesh++){
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

    print_parsed_inputs(pumi_inputs);
}

void print_parsed_inputs(pumi::Mesh_Inputs *pumi_inputs)
{
    printf("Printing the parsed commandline inputs\n\n");
    printf("Number of blocks in x1-direction = %d\n",pumi_inputs->nsubmesh_x1 );
    printf("Domain min-side x1-coordinate    = %2.4f\n",pumi_inputs->domain_x1_min );
    for(int i=0; i<pumi_inputs->nsubmesh_x1; i++){
        printf("\n[X1-BLOCK #%d]\n",i+1 );
        std::cout << "\t type          = " << pumi_inputs->meshtype_x1[i] << "\n";
        printf("\t length        = %2.4e\n",pumi_inputs->block_length_x1[i] );
        printf("\t max-elem-size = %2.4e\n",pumi_inputs->max_elem_size_x1[i] );
        printf("\t min-elem-size = %2.4e\n",pumi_inputs->min_elem_size_x1[i] );
    }
    printf("\n");
    printf("Number of blocks in x2-direction = %d\n",pumi_inputs->nsubmesh_x2 );
    printf("Domain min-side x2-coordinate    = %2.4f\n",pumi_inputs->domain_x2_min );

    for(int i=0; i<pumi_inputs->nsubmesh_x2; i++){
        printf("\n[X2-BLOCK #%d]\n",i+1 );
        std::cout << "\t type          = " << pumi_inputs->meshtype_x2[i] << "\n";
        printf("\t length        = %2.4e\n",pumi_inputs->block_length_x2[i] );
        printf("\t max-elem-size = %2.4e\n",pumi_inputs->max_elem_size_x2[i] );
        printf("\t min-elem-size = %2.4e\n",pumi_inputs->min_elem_size_x2[i] );
    }
}

void print_usage()
{
    printf("Execute the code with the following command line arguments -- \n\n" );
    printf("\t ./install/bin/pumiMBBL2D_Demo N_x1 domain_x1_min \"typeflag_i_x1\" \"block_length_x1\" \"Nel_i_x1\" \"min_elem_size_x1\" N_x2 domain_x2_min \"typeflag_i_x2\" \"block_length_x2\" \"Nel_i_x2\" \"min_elem_size_x2\" \"block_isactive\"\n\n\n");
    printf("\t N_x1     \t\t Total Number of submeshes along the x1-direction \n");
    printf("\t domain_x1_min \t\t Starting x1 coordinate of the domain");
    printf("\t \"typeflag_i_x1\" \t Active mesh type segment in i-th submesh along the x1-direction \n" );
    printf("\t \"block_length_x1\"  \t\t Number of Debye Lengths in i-th submesh along the x1-direction \n");
    printf("\t \"max_elem_size_x1\" \t\t Maximum cell size in Debye lengths for i-th submesh along the x1-direction \n");
    printf("\t \"min_elem_size_x1\"  \t\t For leftBL/rightBL, Minimum cell size in Debye lengths for i-th submesh for i-th submesh along the x1-direction \n");
    printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
    printf("\t \"elem_size_file_x1\" \t\t For arbitrary, list of cell sizes in Debye lengths for i-th submesh along the x1-direction \n");
    printf("\t \t  \t\t For uniform/maxBL/minBL, the inputs will be ignored \n\n");
    printf("\t N_x2     \t\t Total Number of submeshes along the x2-direction \n");
    printf("\t domain_x2_min \t\t Starting x2 coordinate of the domain");
    printf("\t \"typeflag_i_x2\" \t Active mesh type segment in i-th submesh along the x2-direction \n" );
    printf("\t \"block_length_x2\"  \t\t Number of Debye Lengths in i-th submesh along the x2-direction \n");
    printf("\t \"max_elem_size_x2\" \t\t Maximum cell size in Debye lengths for i-th submesh along the x2-direction \n");
    printf("\t \"min_elem_size_x2\"  \t\t For bottomBL/topBL, Minimum cell size in Debye lengths for i-th submesh for i-th submesh along the x2-direction \n");
    printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
    printf("\t \"elem_size_file_x2\" \t\t For arbitrary, list of cell sizes in Debye lengths for i-th submesh along the x2-direction \n");
    printf("\t \t  \t\t For uniform/maxBL/minBL, the inputs will be ignored \n\n");
    printf("\t block_isactive \t Activity info of each submesh-block (N_x1*N_x2 inputs required)\n" );
    printf("\t \t  \t\t 0 is inactive \n\n");
    printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
    // printf("  E.g.#1 [On-HOST]\n\n");
    // printf("    ./pumi-test.host 3 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" 3 \"maxBL,uniform,minBL\" \"50.0,20.0,50.0\" \"4.0,1.0,4.0\" \"1.0,2.0,1.0\" \"1,1,1,1,1,1,1,1,1\" \n\n");
    printf("  E.g.#1 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 3 0.0 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" 3 1.0 \"maxBL,uniform,minBL\" \"50.0,20.0,50.0\" \"4.0,1.0,4.0\" \"1.0,2.0,1.0\" \"1,1,1,1,1,1,1,1,1\" \n\n");
    printf("  E.g.#2 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 4 0.0 \"minBL,uniform,uniform,maxBL\" \"10.0,5.0,5.0,10.0\" \"3.0,1.0,1.0,3.0\" \"1.0,1.0,1.0,1.0\" 4 1.0 \"maxBL,uniform,uniform,minBL\" \"20.0,20.0,20.0,20.0\" \"4.0,1.0,1.0,4.0\" \"1.0,2.0,2.0,1.0\" \"1,0,1,1,1,0,0,1,1,1,1,1,0,1,0,1\" \n\n");
    printf("  E.g.#3 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 1 1.0 \"uniform\" \"20.0\" \"1.0\" \"1.0\" 1 1.0 \"uniform\" \"20.0\" \"1.0\" \"1.0\" \"1\"\n\n");
    printf("  E.g.#4 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 4 0.0 \"minBL,uniform,arbitrary,maxBL\" \"10.0,5.0,5.0,10.0\" \"3.0,1.0,1.0,3.0\" \"1.0,1.0,1.0,1.0\" \"NA,NA,x1-size.dat,NA\" 4 1.0 \"maxBL,arbitrary,uniform,minBL\" \"20.0,20.0,20.0,20.0\" \"4.0,1.0,1.0,4.0\" \"1.0,2.0,2.0,1.0\" \"NA,x2-size.dat,NA,NA\" \"1,0,1,1,1,0,0,1,1,1,1,1,0,1,0,1\" \n\n");
    Kokkos::finalize();
    exit(0);
}

void write2file(pumi::ParticleDataCPU* hp, int N_part, int nstep){
    FILE *part_file;
    char part_filename[30];
    sprintf(part_filename,"part_coords_t%d.dat",nstep);
    part_file = fopen(part_filename,"w");
    int skip = 1;
    for (int i=0; i<N_part; i=i+skip){
        int part_active = hp[i].part_active;
        int subID = hp[i].submeshID;
        int exit_faceID = hp[i].exit_faceID;
        fprintf(part_file, "%d %.5e %.5e %d %d\n", part_active, hp[i].x1, hp[i].x2, subID, exit_faceID);
    }

    fclose(part_file);
}
