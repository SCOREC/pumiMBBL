#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_usage();

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );
  {
      if (argc != 6)
      {
          print_usage();
      }

    int nsubmesh_x1 = atoi( argv[1] );
    int nsubmesh = nsubmesh_x1;
    pumi::Mesh_Inputs *pumi_inputs;
    pumi_inputs = pumi::inputs_allocate(nsubmesh);

    pumi_inputs->ndim = 1; // Fixed pumi input
    pumi_inputs->nsubmesh_x1 = nsubmesh_x1;

    parse_inputs(argc, argv, pumi_inputs);

    pumi::MeshDeviceViewPtr mesh;
    pumi::SubmeshDeviceViewPtr submesh_x1;
    pumi::SubmeshHostViewPtr host_submesh_x1;
    pumi::Mesh_Options pumi_options;
    pumi_options.BL_storage_option = pumi::store_BL_coords_ON;

    submesh_x1 = pumi::submesh_initialize(pumi_inputs, pumi_options, pumi::x1_dir, &host_submesh_x1);
    mesh = pumi::mesh_initialize(pumi_inputs, submesh_x1, host_submesh_x1);

    pumi::MBBL pumi_obj(mesh, submesh_x1, host_submesh_x1);

    pumi::inputs_deallocate(pumi_inputs);

    pumi::MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(mesh);
    Kokkos::deep_copy(h_pumi_mesh, mesh);

    // for (int j=0; j<=h_pumi_mesh(0).Nel_tot_x1; j++){
    //     double cv = pumi::return_covolume(pumi_obj, j);
    //     printf("cv[%3d] = %2.4f \n", j, cv);
    // }


    int N_part = 100000;
    int N_step = 20;
    Kokkos::View<double**> part_coords("particle-coordinates",N_part,3);

    bool test0 = true;

    std::srand((unsigned)(std::time(nullptr)));

    double x1_min = pumi::get_global_left_coord(pumi_obj);
    double x1_max = pumi::get_global_right_coord(pumi_obj);
    double L_x1 = x1_max-x1_min;
    double dist_factor = 100.0;
    int num_push = 0;
    Kokkos::View<double**>::HostMirror h_part_coords = Kokkos::create_mirror_view(part_coords);

    if (test0){
        printf("push-test-0\n");
        for (int ipart=0; ipart<N_part; ipart++){
            bool part_set = false;
            while (!part_set){
                double rand_x1 = (double) rand()/RAND_MAX;

                double q1 = x1_min + L_x1*rand_x1;

                h_part_coords(ipart,0) = q1;
                part_set = true;
            }

        }

        Kokkos::deep_copy(part_coords, h_part_coords);


        Kokkos::parallel_for("particle-locate-0", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, icell;
            double q1 = part_coords(ipart,0);
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            part_coords(ipart,1) = isub;
            part_coords(ipart,2) = icell;
        });
        Kokkos::deep_copy(h_part_coords, part_coords);

        Kokkos::Profiling::pushRegion("push_test_0");
        for (int istep=0; istep<N_step; istep++){
            for (int p=0; p<N_part; p++){
                int part_active = h_part_coords(p,1);
                num_push += (part_active+1 != 0);
            }
            Kokkos::parallel_for("particle-push-test-0", 1, KOKKOS_LAMBDA (const int) {
                for (int ipart=0; ipart<N_part; ipart++){
                    int part_active = part_coords(ipart,1);
                    if (part_active+1){
                        int isub, icell, kcell_x1;
                        double q1 = part_coords(ipart,0);
                        isub = part_coords(ipart,1);
                        icell = part_coords(ipart,2);

                        double Wgh2_x1, Wgh1_x1;
                        pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                        Wgh1_x1 = 1.0-Wgh2_x1;

                        double dq1 = L_x1/dist_factor;

                        q1 += dq1;

                        if (q1 < x1_min || q1 > x1_max){
                            isub = -1;
                            icell = -1;
                        }
                        else {
                            pumi::update_submesh_and_cell_x1(pumi_obj, q1, isub, icell, &isub, &icell);
                            pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                            Wgh1_x1 = 1.0-Wgh2_x1;
                        }

                        part_coords(ipart,1) = isub;
                        part_coords(ipart,2) = icell;
                        part_coords(ipart,0) = q1;
                    }
                }
            });
            Kokkos::deep_copy(h_part_coords, part_coords);
        }
        Kokkos::Profiling::popRegion();
        printf("Total number of particle pushes executed in Test-0 = %d\n",num_push );
    }

  }
  Kokkos::finalize();

  return 0;
}

void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs)
{
    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int nsubmesh_x2 = pumi_inputs->nsubmesh_x2;
    // reading submesh meshtypes
    char all_submesh_flag_x1[MAX_SUBMESHES*10];
    char each_submesh_flag_x1[MAX_SUBMESHES][10];
    strcpy(all_submesh_flag_x1, argv[2]);

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
    strcpy(all_p1_submesh_x1, argv[3]);

    tok = strtok(all_p1_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of p1_i_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh max elemsize
    char all_p2max_submesh_x1[MAX_SUBMESHES*10];
    char each_p2max_submesh_x1[MAX_SUBMESHES][10];
    strcpy(all_p2max_submesh_x1, argv[4]);

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
    strcpy(all_p2min_submesh_x1, argv[5]);

    tok = strtok(all_p2min_submesh_x1, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh_x1[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of p2max_i_x1 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    for (isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
        (pumi_inputs->meshtype).push_back(each_submesh_flag_x1[isubmesh]);
        *(pumi_inputs->p1_i_x1 + isubmesh) = atof(each_p1_submesh_x1[isubmesh]);
        *(pumi_inputs->p2max_i_x1 + isubmesh) = atof( each_p2max_submesh_x1[isubmesh]);
        *(pumi_inputs->p2min_i_x1 + isubmesh) = atof( each_p2min_submesh_x1[isubmesh]);
    }

}

void print_usage()
{
    printf("Execute the code with the following command line arguments -- \n\n" );
    printf("\t ./install/bin/pumiMBBL1D_Demo N_x1 \"typeflag_i_x1\" \"p1_i_x1\" \"Nel_i_x1\" \"p2min_i_x1\" \n\n\n");
    printf("\t N_x1     \t\t Total Number of submeshes along the x1-direction \n");
    printf("\t \"typeflag_i_x1\" \t Active mesh type segment in i-th submesh along the x1-direction \n" );
    printf("\t \"p1_i_x1\"  \t\t Number of Debye Lengths in i-th submesh along the x1-direction \n");
    printf("\t \"p2max_i_x1\" \t\t Maximum cell size in Debye lengths for i-th submesh along the x1-direction \n");
    printf("\t \"p2min_i_x1\"  \t\t For leftBL/rightBL, Minimum cell size in Debye lengths for i-th submesh for i-th submesh along the x1-direction \n");
    printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
    printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
    printf("  E.g.#1 [On-DEVICE]\n\n");
    printf("    ./install/bin/pumiMBBL1D_Demo 3 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" \n\n");
    Kokkos::finalize();
    exit(0);
}
