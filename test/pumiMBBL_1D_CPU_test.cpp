#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_parsed_inputs(pumi::Mesh_Inputs *pumi_inputs);
void print_usage();

// void print_partdata(pumi::ParticleDataCPU* pd, int Npart){
//     for (int ipart=0; ipart<Npart; ipart++){
//         printf("Particle-%3d - q=%2.4f - sub=%2d - cell=%3d \n",ipart, pd[ipart].x1, pd[ipart].submeshID, pd[ipart].cellID );
//     }
// }

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );
  {
      if (argc != 8)
      {
          print_usage();
      }

    int nsubmesh_x1 = atoi( argv[1] );
    double domain_x1_min = atof( argv[2] );

    pumi::Mesh_Inputs *pumi_inputs = pumi::inputs_allocate();

    pumi_inputs->ndim = 1; // Fixed pumi input
    pumi_inputs->nsubmesh_x1 = nsubmesh_x1;
    pumi_inputs->domain_x1_min = domain_x1_min;

    parse_inputs(argc, argv, pumi_inputs);

    pumi::Mesh_Options pumi_options;
    pumi_options.BL_storage_option = pumi::store_BL_coords_ON;
    pumi_options.print_node_option = pumi::print_node_coords_OFF;

    pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);

    printf("Mesh volume = %2.2f\n",pumi::get_mesh_volume(pumi_obj) );
    pumi::inputs_deallocate(pumi_inputs);

    int N_part = 100;
    int N_step = 20;

    pumi::ParticleDataCPU* Partdata = new pumi::ParticleDataCPU[N_part];

    bool test0 = true;

    std::srand((unsigned)(std::time(nullptr)));

    double x1_min = pumi::get_global_x1_min_coord(pumi_obj);
    double x1_max = pumi::get_global_x1_max_coord(pumi_obj);
    double L_x1 = x1_max-x1_min;
    double dist_factor = 100.0;
    int num_push = 0;

    if (test0){
        printf("push-test-0\n");
        for (int ipart=0; ipart<N_part; ipart++){
            bool part_set = false;
            while (!part_set){
                double rand_x1 = (double) rand()/RAND_MAX;

                double q1 = x1_min + L_x1*rand_x1;

                Partdata[ipart] = pumi::ParticleDataCPU(q1,0.0);
                part_set = true;
            }
        }

        for (int ipart=0; ipart<N_part; ipart++){
            int isub, icell;
            double q1 = Partdata[ipart].x1;
            pumi::locate_submesh_and_cell_x1_host(pumi_obj, q1, &isub, &icell);
            Partdata[ipart] = pumi::ParticleDataCPU(q1,0.0,isub,icell,true,-1);
        }

        for (int istep=0; istep<N_step; istep++){
            for (int p=0; p<N_part; p++){
                bool part_active = Partdata[p].part_active;
                num_push += part_active;
            }
            for (int ipart=0; ipart<N_part; ipart++){
                bool part_active = Partdata[ipart].part_active;
                if (part_active){
                    int isub, icell, kcell_x1;
                    double q1 = Partdata[ipart].x1;
                    isub = Partdata[ipart].submeshID;
                    icell = Partdata[ipart].cellID;

                    double Wgh2_x1, Wgh1_x1;
                    pumi::calc_weights_x1_host(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                    Wgh1_x1 = 1.0-Wgh2_x1;

                    double dq1 = L_x1/dist_factor;

                    q1 += dq1;

                    if (q1 < x1_min || q1 > x1_max){
                        int exit_faceID = (q1 < x1_min) + 2*(q1>x1_max) - 1;
                        isub = -1;
                        icell = -1;
                        Partdata[ipart] = pumi::ParticleDataCPU(q1,0.0,isub,icell,false,exit_faceID);
                    }
                    else {
                        pumi::update_submesh_and_cell_x1_host(pumi_obj, q1, isub, icell, &isub, &icell);
                        pumi::calc_weights_x1_host(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                        Wgh1_x1 = 1.0-Wgh2_x1;
                        Partdata[ipart] = pumi::ParticleDataCPU(q1,0.0,isub,icell,true,-1,Wgh1_x1,Wgh2_x1,0.0,0.0);
                    }

                }
            }
        }
        printf("Total number of particle pushes executed in Test-0 = %d\n",num_push );
    }

    // print_partdata(Partdata, N_part);

  }
  Kokkos::finalize();

  return 0;
}

void parse_inputs(int , char* argv[], pumi::Mesh_Inputs *pumi_inputs)
{
    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;

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

    for (isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
        (pumi_inputs->meshtype_x1).push_back(each_submesh_flag_x1[isubmesh]);
        (pumi_inputs->block_length_x1).push_back(atof(each_p1_submesh_x1[isubmesh]));
        (pumi_inputs->max_elem_size_x1).push_back(atof(each_p2max_submesh_x1[isubmesh]));
        (pumi_inputs->min_elem_size_x1).push_back(atof(each_p2min_submesh_x1[isubmesh]));
        std::string x1_elemsize_file(each_arb_submesh_x1[isubmesh]);
        (pumi_inputs->arbitrary_x1_elemsize_file).push_back(x1_elemsize_file);
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
}

void print_usage()
{
    printf("Execute the code with the following command line arguments -- \n\n" );
    printf("\t ./install/bin/pumiMBBL1D_Demo N_x1 domain_x1_min \"typeflag_i_x1\" \"block_length_x1\" \"Nel_i_x1\" \"min_elem_size_x1\" \n\n\n");
    printf("\t N_x1     \t\t Total Number of submeshes along the x1-direction \n");
    printf("\t domain_x1_min \t\t Starting x1 coordinate of the domain");
    printf("\t \"typeflag_i_x1\" \t Active mesh type segment in i-th submesh along the x1-direction \n" );
    printf("\t \"block_length_x1\"  \t\t Number of Debye Lengths in i-th submesh along the x1-direction \n");
    printf("\t \"max_elem_size_x1\" \t\t Maximum cell size in Debye lengths for i-th submesh along the x1-direction \n");
    printf("\t \"min_elem_size_x1\"  \t\t For leftBL/rightBL, Minimum cell size in Debye lengths for i-th submesh for i-th submesh along the x1-direction \n");
    printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
    printf("\t \"elem_size_file_x1\" \t\t For arbitrary, list of cell sizes in Debye lengths for i-th submesh along the x1-direction \n");
    printf("\t \t  \t\t For uniform/maxBL/minBL, the inputs will be ignored \n\n");
    printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
    printf("  E.g.#1 \n\n");
    printf("    ./install/bin/pumiMBBL1D_Demo 3 1.0 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" \"NA,NA,NA\" \n\n");
    printf("  E.g.#2 \n\n");
    printf("    ./install/bin/pumiMBBL1D_Demo 3 1.0 \"minBL,uniform,arbitrary\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" \"NA,NA,x1-size.dat\" \n\n");
    Kokkos::finalize();
    exit(0);
}
