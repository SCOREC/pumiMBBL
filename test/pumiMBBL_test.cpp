#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_usage();
void write2file(Kokkos::View<double**>::HostMirror hp_coords, int N_part, int nstep);

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );
  {
      if (argc != 12)
      {
          print_usage();
      }

    int nsubmesh_x1 = atoi( argv[1] );
    int nsubmesh_x2 = atoi( argv[6] );
    int nsubmesh = nsubmesh_x1+nsubmesh_x2;
    pumi::Mesh_Inputs *pumi_inputs;
    pumi_inputs = pumi::inputs_allocate(nsubmesh);

    pumi_inputs->ndim = 2; // Fixed pumi input
    pumi_inputs->nsubmesh_x1 = nsubmesh_x1;
    pumi_inputs->nsubmesh_x2 = nsubmesh_x2;

    parse_inputs(argc, argv, pumi_inputs);


    pumi::MeshDeviceViewPtr mesh;
    pumi::SubmeshDeviceViewPtr submesh_x1;
    pumi::SubmeshHostViewPtr host_submesh_x1;
    pumi::SubmeshDeviceViewPtr submesh_x2;
    pumi::SubmeshHostViewPtr host_submesh_x2;
    pumi::Mesh_Options pumi_options;
    pumi_options.BL_storage_option = pumi::store_BL_coords_ON;

    submesh_x1 = pumi::submesh_initialize(pumi_inputs, pumi_options, pumi::x1_dir, &host_submesh_x1);
    submesh_x2 = pumi::submesh_initialize(pumi_inputs, pumi_options, pumi::x2_dir, &host_submesh_x2);
    mesh = pumi::mesh_initialize(pumi_inputs, submesh_x1, host_submesh_x1, submesh_x2, host_submesh_x2);

    pumi::MBBL pumi_obj(mesh, submesh_x1, host_submesh_x1, submesh_x2, host_submesh_x2);

    pumi::inputs_deallocate(pumi_inputs);

    pumi::MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(mesh);
    Kokkos::deep_copy(h_pumi_mesh, mesh);

    pumi::print_mesh_skeleton(pumi_obj);


    Kokkos::parallel_for("bdry-test-1", 1, KOKKOS_LAMBDA (int j) {
        int Nx = pumi_obj.mesh(0).nsubmesh_x1;
        int Ny = pumi_obj.mesh(0).nsubmesh_x2;
        for (int i=0; i<2*Nx*Ny+Nx+Ny; i++){
            printf("edge-%2d -- isbdry-%d\n",i,pumi_obj.mesh(0).is_bdry(i) );
        }
    });

    int N_part = 100000;
    int N_step = 20;
    Kokkos::View<double**> part_coords("particle-coordinates",N_part,6);

    bool test0 = true;
    bool test1 = true;
    bool test2 = false;

    std::srand((unsigned)(std::time(nullptr)));

    double x1_min = pumi::get_global_left_coord(pumi_obj);
    double x1_max = pumi::get_global_right_coord(pumi_obj);
    double L_x1 = x1_max-x1_min;
    double x2_min = pumi::get_global_bottom_coord(pumi_obj);
    double x2_max = pumi::get_global_top_coord(pumi_obj);
    double L_x2 = x2_max-x2_min;
    double dist_factor = 100.0;
    int num_push = 0;
    Kokkos::View<double**>::HostMirror h_part_coords = Kokkos::create_mirror_view(part_coords);

    if (test0){
        printf("push-test-0\n");
        for (int ipart=0; ipart<N_part; ipart++){
            bool part_set = false;
            while (!part_set){
                double rand_x1 = (double) rand()/RAND_MAX;
                double rand_x2 = (double) rand()/RAND_MAX;

                double q1 = x1_min + L_x1*rand_x1;
                double q2 = x2_min + L_x2*rand_x2;

                int isub, jsub;

                for (int i=1; i<=h_pumi_mesh(0).nsubmesh_x1; i++){
                    if (pumi_obj.host_submesh_x1[i].xmin < q1 && pumi_obj.host_submesh_x1[i].xmax > q1){
                        isub = i;
                        break;
                    }
                }
                for (int j=1; j<=h_pumi_mesh(0).nsubmesh_x2; j++){
                    if (pumi_obj.host_submesh_x2[j].xmin < q2 && pumi_obj.host_submesh_x2[j].xmax > q2){
                        jsub = j;
                        break;
                    }
                }

                if (h_pumi_mesh(0).host_isactive[isub][jsub]){
                    h_part_coords(ipart,0) = q1;
                    h_part_coords(ipart,1) = q2;
                    // h_part_activity(ipart) = true;
                    part_set = true;
                }
            }

        }

        Kokkos::deep_copy(part_coords, h_part_coords);


        Kokkos::parallel_for("particle-locate-0", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, jsub, icell, jcell;
            double q1 = part_coords(ipart,0);
            double q2 = part_coords(ipart,1);
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
            part_coords(ipart,2) = isub;
            part_coords(ipart,3) = jsub;
            part_coords(ipart,4) = icell;
            part_coords(ipart,5) = jcell;
        });
        Kokkos::deep_copy(h_part_coords, part_coords);
        // write2file(h_part_coords, N_part, 0);

        Kokkos::Profiling::pushRegion("push_test_0");
        for (int istep=0; istep<N_step; istep++){
            for (int p=0; p<N_part; p++){
                int part_active = h_part_coords(p,2);
                num_push += (part_active+1 != 0);
            }
            Kokkos::parallel_for("particle-push-test-0", 1, KOKKOS_LAMBDA (int j) {
                for (int ipart=0; ipart<N_part; ipart++){
                    int part_active = part_coords(ipart,2);
                    if (part_active+1){
                        int isub, jsub, icell, jcell, kcell_x1, kcell_x2;
                        int global_cell, topleft_node, bottomleft_node;
                        double q1 = part_coords(ipart,0);
                        double q2 = part_coords(ipart,1);
                        isub = part_coords(ipart,2);
                        jsub = part_coords(ipart,3);
                        icell = part_coords(ipart,4);
                        jcell = part_coords(ipart,5);

                        double Wgh2_x1, Wgh2_x2, Wgh1_x1, Wgh1_x2;
                        pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                        pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                        Wgh1_x1 = 1.0-Wgh2_x1;
                        Wgh1_x2 = 1.0-Wgh2_x2;
                        pumi::calc_global_cellID_and_nodeID_fullmesh(pumi_obj, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        // double dq1 = -10.0 + 20.0*part_x1_disp(istep,ipart);
                        // double dq2 = -5.0 + 10.0*part_x2_disp(istep,ipart);
                        double dq1 = L_x1/dist_factor;
                        double dq2 = -L_x2/dist_factor;

                        q1 += dq1;
                        q2 += dq2;

                        if (q1 < x1_min || q1 > x1_max || q2 < x2_min || q2 > x2_max){
                            isub = -1;
                            icell = -1;
                            jsub = -1;
                            jcell = -1;
                        }
                        else {
                            pumi::update_submesh_and_cell_x1(pumi_obj, q1, isub, icell, &isub, &icell);
                            pumi::update_submesh_and_cell_x2(pumi_obj, q2, jsub, jcell, &jsub, &jcell);
                            pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                            pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                            Wgh1_x1 = 1.0-Wgh2_x1;
                            Wgh1_x2 = 1.0-Wgh2_x2;
                            pumi::calc_global_cellID_and_nodeID_fullmesh(pumi_obj, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        }

                        part_coords(ipart,2) = isub;
                        part_coords(ipart,3) = jsub;
                        part_coords(ipart,4) = icell;
                        part_coords(ipart,5) = jcell;
                        part_coords(ipart,0) = q1;
                        part_coords(ipart,1) = q2;
                    }
                }
            });
            Kokkos::deep_copy(h_part_coords, part_coords);
            // write2file(h_part_coords, N_part, istep+1);
        }
        Kokkos::Profiling::popRegion();
        printf("Total number of particle pushes executed in Test-0 = %d\n",num_push );
    }
    if (test1){
        printf("push-test-1\n");
        for (int ipart=0; ipart<N_part; ipart++){
            bool part_set = false;
            while (!part_set){
                double rand_x1 = (double) rand()/RAND_MAX;
                double rand_x2 = (double) rand()/RAND_MAX;

                double q1 = x1_min + L_x1*rand_x1;
                double q2 = x2_min + L_x2*rand_x2;

                int isub, jsub;

                for (int i=1; i<=h_pumi_mesh(0).nsubmesh_x1; i++){
                    if (pumi_obj.host_submesh_x1[i].xmin < q1 && pumi_obj.host_submesh_x1[i].xmax > q1){
                        isub = i;
                        break;
                    }
                }
                for (int j=1; j<=h_pumi_mesh(0).nsubmesh_x2; j++){
                    if (pumi_obj.host_submesh_x2[j].xmin < q2 && pumi_obj.host_submesh_x2[j].xmax > q2){
                        jsub = j;
                        break;
                    }
                }

                if (h_pumi_mesh(0).host_isactive[isub][jsub]){
                    h_part_coords(ipart,0) = q1;
                    h_part_coords(ipart,1) = q2;
                    part_set = true;
                }
            }

        }

        Kokkos::deep_copy(part_coords, h_part_coords);

        Kokkos::parallel_for("particle-locate-1", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, jsub, icell, jcell;
            double q1 = part_coords(ipart,0);
            double q2 = part_coords(ipart,1);
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
            part_coords(ipart,2) = isub;
            part_coords(ipart,3) = jsub;
            part_coords(ipart,4) = icell;
            part_coords(ipart,5) = jcell;
        });
        Kokkos::deep_copy(h_part_coords, part_coords);
        // write2file(h_part_coords, N_part, N_step+1);

        num_push = 0;
        Kokkos::Profiling::pushRegion("push_test_1");
        for (int istep=0; istep<N_step; istep++){
            for (int p=0; p<N_part; p++){
                int part_active = h_part_coords(p,2);
                num_push += (part_active+1 != 0);
            }
            Kokkos::parallel_for("particle-push-test-1", 1, KOKKOS_LAMBDA (int j) {
                for (int ipart=0; ipart<N_part; ipart++){
                    int part_active = part_coords(ipart,2);
                    if (part_active+1){
                        int isub, jsub, icell, jcell, bdry_hit, kcell_x1, kcell_x2;
                        bool in_domain;
                        int global_cell, topleft_node, bottomleft_node;
                        double q1 = part_coords(ipart,0);
                        double q2 = part_coords(ipart,1);
                        isub = part_coords(ipart,2);
                        jsub = part_coords(ipart,3);
                        icell = part_coords(ipart,4);
                        jcell = part_coords(ipart,5);

                        double Wgh2_x1, Wgh2_x2, Wgh1_x1, Wgh1_x2;
                        pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                        pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                        Wgh1_x1 = 1.0-Wgh2_x1;
                        Wgh1_x2 = 1.0-Wgh2_x2;
                        pumi::calc_global_cellID_and_nodeID(pumi_obj, isub, jsub, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        // double dq1 = -10.0 + 20.0*part_x1_disp(istep,ipart);
                        // double dq2 = -5.0 + 10.0*part_x2_disp(istep,ipart);
                        double dq1 = L_x1/dist_factor;
                        double dq2 = -L_x2/dist_factor;



                        pumi::push_particle(pumi_obj, q1, q2, dq1, dq2, &isub, &jsub, &icell, &jcell, &in_domain, &bdry_hit);
                        // part_activity(ipart) = in_domain;
                        if (!in_domain){
                            // printf("hit-edge=%2d -- is_bdry=%d\n",bdry_hit,pumi_obj.mesh(0).is_bdry(bdry_hit));
                            // part_activity(ipart) = false;
                            isub = -1;
                            jsub = -1;
                            icell = -1;
                            jcell = -1;
                        }
                        else{
                            pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                            pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                            Wgh1_x1 = 1.0-Wgh2_x1;
                            Wgh1_x2 = 1.0-Wgh2_x2;
                            pumi::calc_global_cellID_and_nodeID(pumi_obj, isub, jsub, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        }

                        part_coords(ipart,2) = isub;
                        part_coords(ipart,3) = jsub;
                        part_coords(ipart,4) = icell;
                        part_coords(ipart,5) = jcell;
                        part_coords(ipart,0) = q1+dq1;
                        part_coords(ipart,1) = q2+dq2;

                    }
                }
            });
            Kokkos::deep_copy(h_part_coords, part_coords);

            // write2file(h_part_coords, N_part, istep+2+N_step);
        }
        Kokkos::Profiling::popRegion();
        printf("Total number of particle pushes executed in Test-1 = %d\n",num_push );
    }
    if (test2){
        printf("push-test-2\n");
        for (int ipart=0; ipart<N_part; ipart++){
            bool part_set = false;
            while (!part_set){
                double rand_x1 = (double) rand()/RAND_MAX;
                double rand_x2 = (double) rand()/RAND_MAX;

                double q1 = x1_min + L_x1*rand_x1;
                double q2 = x2_min + L_x2*rand_x2;

                int isub, jsub;

                for (int i=1; i<=h_pumi_mesh(0).nsubmesh_x1; i++){
                    if (pumi_obj.host_submesh_x1[i].xmin < q1 && pumi_obj.host_submesh_x1[i].xmax > q1){
                        isub = i;
                        break;
                    }
                }
                for (int j=1; j<=h_pumi_mesh(0).nsubmesh_x2; j++){
                    if (pumi_obj.host_submesh_x2[j].xmin < q2 && pumi_obj.host_submesh_x2[j].xmax > q2){
                        jsub = j;
                        break;
                    }
                }

                if (h_pumi_mesh(0).host_isactive[isub][jsub]){
                    h_part_coords(ipart,0) = q1;
                    h_part_coords(ipart,1) = q2;
                    part_set = true;
                }
            }

        }

        Kokkos::deep_copy(part_coords, h_part_coords);

        Kokkos::parallel_for("particle-locate-2", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, jsub, icell, jcell;
            double q1 = part_coords(ipart,0);
            double q2 = part_coords(ipart,1);
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
            part_coords(ipart,2) = isub;
            part_coords(ipart,3) = jsub;
            part_coords(ipart,4) = icell;
            part_coords(ipart,5) = jcell;
        });
        Kokkos::deep_copy(h_part_coords, part_coords);
        // write2file(h_part_coords, N_part, 2*N_step+2);

        num_push = 0;
        Kokkos::Profiling::pushRegion("push_test_2");
        for (int istep=0; istep<N_step; istep++){
            for (int p=0; p<N_part; p++){
                int part_active = h_part_coords(p,2);
                num_push += (part_active+1 != 0);
            }
            Kokkos::parallel_for("particle-push-test-2", 1, KOKKOS_LAMBDA (int j) {
                for (int ipart=0; ipart<N_part; ipart++){
                    int part_active = part_coords(ipart,2);
                    if (part_active+1){
                        int isub, jsub, icell, jcell, bdry_hit, kcell_x1, kcell_x2;
                        bool in_domain;
                        int global_cell, topleft_node, bottomleft_node;
                        double q1 = part_coords(ipart,0);
                        double q2 = part_coords(ipart,1);
                        isub = part_coords(ipart,2);
                        jsub = part_coords(ipart,3);
                        icell = part_coords(ipart,4);
                        jcell = part_coords(ipart,5);

                        double Wgh2_x1, Wgh2_x2, Wgh1_x1, Wgh1_x2;
                        pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                        pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                        Wgh1_x1 = 1.0-Wgh2_x1;
                        Wgh1_x2 = 1.0-Wgh2_x2;
                        pumi::calc_global_cellID_and_nodeID(pumi_obj, isub, jsub, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        // double dq1 = -10.0 + 20.0*part_x1_disp(istep,ipart);
                        // double dq2 = -5.0 + 10.0*part_x2_disp(istep,ipart);
                        double dq1 = L_x1/dist_factor;
                        double dq2 = -L_x2/dist_factor;



                        pumi::push_particle_v2(pumi_obj, q1, q2, dq1, dq2, &isub, &jsub, &icell, &jcell, &in_domain, &bdry_hit);
                        // part_activity(ipart) = in_domain;
                        if (!in_domain){
                            // printf("hit-edge=%2d -- is_bdry=%d\n",bdry_hit,pumi_obj.mesh(0).is_bdry(bdry_hit));
                            // part_activity(ipart) = false;
                            isub = -1;
                            jsub = -1;
                            icell = -1;
                            jcell = -1;
                        }
                        else{
                            pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                            pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                            Wgh1_x1 = 1.0-Wgh2_x1;
                            Wgh1_x2 = 1.0-Wgh2_x2;
                            pumi::calc_global_cellID_and_nodeID(pumi_obj, isub, jsub, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                        }

                        part_coords(ipart,2) = isub;
                        part_coords(ipart,3) = jsub;
                        part_coords(ipart,4) = icell;
                        part_coords(ipart,5) = jcell;
                        part_coords(ipart,0) = q1+dq1;
                        part_coords(ipart,1) = q2+dq2;

                    }
                }
            });
            Kokkos::deep_copy(h_part_coords, part_coords);

            // write2file(h_part_coords, N_part, 2*N_step+3+istep);
        }
        Kokkos::Profiling::popRegion();
        printf("Total number of particle pushes executed in Test-2 = %d\n",num_push );
    }
    Kokkos::parallel_for("bdry-test-1", 1, KOKKOS_LAMBDA (int j) {
        int Nx = pumi_obj.mesh(0).nsubmesh_x1;
        int Ny = pumi_obj.mesh(0).nsubmesh_x2;
        for (int i=0; i<2*Nx*Ny+Nx+Ny; i++){
            printf("edge-%2d -- isbdry-%d\n",i,pumi_obj.mesh(0).is_bdry(i) );
        }
    });


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

    // reading submesh meshtypes
    char all_submesh_flag_x2[MAX_SUBMESHES*10];
    char each_submesh_flag_x2[MAX_SUBMESHES][10];
    strcpy(all_submesh_flag_x2, argv[7]);

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
    strcpy(all_p1_submesh_x2, argv[8]);

    tok = strtok(all_p1_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p1_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
        printf("ERROR: Number of p1_i_x2 arguments not equal to number of submeshes...\n");
        exit(0);
    }

    //reading submesh max elemsize
    char all_p2max_submesh_x2[MAX_SUBMESHES*10];
    char each_p2max_submesh_x2[MAX_SUBMESHES][10];
    strcpy(all_p2max_submesh_x2, argv[9]);

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
    strcpy(all_p2min_submesh_x2, argv[10]);

    tok = strtok(all_p2min_submesh_x2, ",");
    isubmesh=0;
    while (tok != NULL){
      strcpy (each_p2min_submesh_x2[isubmesh], tok);
      tok = strtok(NULL, ",");
      isubmesh++;
    }
    //print error if number of inputs do not match nsubmeshes
    if (isubmesh != pumi_inputs->nsubmesh_x2){
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
    if (isubmesh != pumi_inputs->nsubmesh_x2*pumi_inputs->nsubmesh_x1){
        printf("ERROR: Number of block_isactive arguments not equal to number of submeshes...\n");
        exit(0);
    }

    for (isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
        (pumi_inputs->meshtype).push_back(each_submesh_flag_x1[isubmesh]);
        *(pumi_inputs->p1_i_x1 + isubmesh) = atof(each_p1_submesh_x1[isubmesh]);
        *(pumi_inputs->p2max_i_x1 + isubmesh) = atof( each_p2max_submesh_x1[isubmesh]);
        *(pumi_inputs->p2min_i_x1 + isubmesh) = atof( each_p2min_submesh_x1[isubmesh]);
    }

    for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
        (pumi_inputs->meshtype).push_back(each_submesh_flag_x2[jsubmesh]);
        *(pumi_inputs->p1_i_x2 + jsubmesh) = atof(each_p1_submesh_x2[jsubmesh]);
        *(pumi_inputs->p2max_i_x2 + jsubmesh) = atof( each_p2max_submesh_x2[jsubmesh]);
        *(pumi_inputs->p2min_i_x2 + jsubmesh) = atof( each_p2min_submesh_x2[jsubmesh]);
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
}

void print_usage()
{
    printf("Execute the code with the following command line arguments -- \n\n" );
    printf("\t ./install/bin/pumiMBBL2D_Demo N_x1 \"typeflag_i_x1\" \"p1_i_x1\" \"Nel_i_x1\" \"p2min_i_x1\" N_x2 \"typeflag_i_x2\" \"p1_i_x2\" \"Nel_i_x2\" \"p2min_i_x2\" \"block_isactive\"\n\n\n");
    printf("\t N_x1     \t\t Total Number of submeshes along the x1-direction \n");
    printf("\t \"typeflag_i_x1\" \t Active mesh type segment in i-th submesh along the x1-direction \n" );
    printf("\t \"p1_i_x1\"  \t\t Number of Debye Lengths in i-th submesh along the x1-direction \n");
    printf("\t \"p2max_i_x1\" \t\t Maximum cell size in Debye lengths for i-th submesh along the x1-direction \n");
    printf("\t \"p2min_i_x1\"  \t\t For leftBL/rightBL, Minimum cell size in Debye lengths for i-th submesh for i-th submesh along the x1-direction \n");
    printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
    printf("\t N_x2     \t\t Total Number of submeshes along the x2-direction \n");
    printf("\t \"typeflag_i_x2\" \t Active mesh type segment in i-th submesh along the x2-direction \n" );
    printf("\t \"p1_i_x2\"  \t\t Number of Debye Lengths in i-th submesh along the x2-direction \n");
    printf("\t \"p2max_i_x2\" \t\t Maximum cell size in Debye lengths for i-th submesh along the x2-direction \n");
    printf("\t \"p2min_i_x2\"  \t\t For bottomBL/topBL, Minimum cell size in Debye lengths for i-th submesh for i-th submesh along the x2-direction \n");
    printf("\t \t  \t\t For uniform, the inputs will be ignored \n\n");
    printf("\t block_isactive \t Activity info of each submesh-block (N_x1*N_x2 inputs required)\n" );
    printf("\t \t  \t\t 0 is inactive \n\n");
    printf("\t ENSURE INPUTS FOR EACH SUBMESH ARE SEPARATED BY A COMMA AND WITHOUT ANY SPACES\n\n");
    // printf("  E.g.#1 [On-HOST]\n\n");
    // printf("    ./pumi-test.host 3 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" 3 \"maxBL,uniform,minBL\" \"50.0,20.0,50.0\" \"4.0,1.0,4.0\" \"1.0,2.0,1.0\" \"1,1,1,1,1,1,1,1,1\" \n\n");
    printf("  E.g.#1 [On-DEVICE]\n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 3 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" 3 \"maxBL,uniform,minBL\" \"50.0,20.0,50.0\" \"4.0,1.0,4.0\" \"1.0,2.0,1.0\" \"1,1,1,1,1,1,1,1,1\" \n\n");
    Kokkos::finalize();
    exit(0);
}

void write2file(Kokkos::View<double**>::HostMirror hp_coords, int N_part, int nstep){
    FILE *part_file;
    char part_filename[30];
    sprintf(part_filename,"part_coords_t%d.dat",nstep);
    part_file = fopen(part_filename,"w");

    for (int i=0; i<N_part; i++){
        int k;
        int part_active = hp_coords(i,2);
        if (part_active+1){
            k=1;
        }
        else{
            k=0;
        }
        fprintf(part_file, "%d %.5e %.5e\n", k, hp_coords(i,0), hp_coords(i,1));
    }

    fclose(part_file);
}
