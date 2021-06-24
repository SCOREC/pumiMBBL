#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_usage();
void write2file(Kokkos::View<double*[2]>::HostMirror hp_coords, Kokkos::View<bool*>::HostMirror hp_active, int N_part, int nstep);

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
    // std::cout << "Printing x1 grading ratio\n";
    // for (int i=1; i<h_pumi_mesh(0).Nel_tot_x1; i++){
    //     std::cout << "x1-r[" << i <<"] = " << pumi::return_gradingratio(pumi_obj, pumi::x1_dir, i) << "\n";
    // }
    // std::cout << "\n\n";
    // std::cout << "Printing x2 grading ratio\n";
    // for (int i=1; i<h_pumi_mesh(0).Nel_tot_x2; i++){
    //     std::cout << "x2-r[" << i <<"] = " << pumi::return_gradingratio(pumi_obj, pumi::x2_dir, i) << "\n";
    // }
    // std::cout << "\n\n";
    // std::cout << "Printing x1 cell lengths\n";
    // for (int i=0; i<h_pumi_mesh(0).Nel_tot_x1; i++){
    //     double cs = pumi::return_elemsize(pumi_obj, pumi::x1_dir, i, pumi::elem_input_offset);
    //     std::cout << "dx1[" << i << "] = " << cs << "\n";
    // }
    // std::cout << "\n\n";
    // std::cout << "Printing x2 cell lengths\n";
    // for (int i=0; i<h_pumi_mesh(0).Nel_tot_x2; i++){
    //     double cs = pumi::return_elemsize(pumi_obj, pumi::x2_dir, i, pumi::elem_input_offset);
    //     std::cout << "dx2[" << i << "] = " << cs << "\n";
    // }
    // std::cout << "\n\n";
    // int k=0;
    // for (int i=0; i<=h_pumi_mesh(0).Nel_tot_x2; i++){
    //     for (int j=0; j<=h_pumi_mesh(0).Nel_tot_x1; j++){
    //         double cv = pumi::return_covolume(pumi_obj, j, i);
    //         std::cout << "cv[" << j << "," << i << "] = " << cv << "\n";
    //         k++;
    //     }
    // }

    // std::srand((unsigned)(std::time(nullptr)));
    // double rand_val_x1 = (double) rand()/RAND_MAX;
    // double rand_val_x2 = (double) rand()/RAND_MAX;
    //
    // double x1_min = pumi::get_global_left_coord(pumi_obj);
    // double x1_max = pumi::get_global_right_coord(pumi_obj);
    // double L_x1 = x1_max-x1_min;
    // double x2_min = pumi::get_global_bottom_coord(pumi_obj);
    // double x2_max = pumi::get_global_top_coord(pumi_obj);
    // double L_x2 = x2_max-x2_min;
    //
    // double h_q1 = (0.9*x1_min + 0.1*x1_max) + 0.25*rand_val_x1*L_x1;
    // double h_q2 = (0.1*x2_min + 0.9*x2_max) - 0.25*rand_val_x2*L_x2;
    // double h_q1_new = h_q1 + rand_val_x1*L_x1;
    // double h_q2_new = h_q2 - rand_val_x2*L_x2;
    //
    // Kokkos::parallel_for("full-mesh-particle-ops-test", 1, KOKKOS_LAMBDA (const int) {
    //     int isub, jsub, icell, jcell, kcell_x1, kcell_x2, lb_node, lt_node, global_cell;
    //     double Wgh2;
    //     double q1 = h_q1;
    //     printf("q1=%2.4f\n",q1);
    //     pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
    //     printf("x1-submesh=%d\n",isub );
    //     printf("x1-cell=%d\n",icell );
    //     pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2);
    //     printf("x1-glob-cell=%d\n",kcell_x1 );
    //     printf("x1-W2=%2.4f\n\n",Wgh2 );
    //     double q2 = h_q2;
    //     printf("q2=%2.4f\n",q2);
    //     pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
    //     printf("x2-submesh=%d\n",jsub );
    //     printf("x2-cell=%d\n",jcell );
    //     pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2);
    //     printf("x2-glob-cell=%d\n",kcell_x2 );
    //     printf("x2-W2=%2.4f\n\n",Wgh2 );
    //
    //
    //     pumi::calc_global_cellID_and_nodeID_fullmesh(pumi_obj, kcell_x1, kcell_x2, &global_cell, &lb_node, &lt_node);
    //     printf("2D-global-cell=%d\n",global_cell);
    //     printf("leftbottom-node=%d    lefttop-node=%d\n",lb_node, lt_node);
    //
    //
    //     printf("\n\nAfter push...\n\n");
    //
    //     q1 = h_q1_new;
    //     printf("new q1=%2.4f\n",q1);
    //     pumi::update_submesh_and_cell_x1(pumi_obj, q1, isub, icell, &isub, &icell);
    //     printf("new-x1-submesh=%d\n", isub);
    //     printf("new-x1-cell=%d\n", icell);
    //     pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2);
    //     printf("new-x1-glob-cell=%d\n",kcell_x1 );
    //     printf("new-x1-W2=%2.4f\n\n",Wgh2 );
    //
    //     q2 = h_q2_new;
    //     printf("new q2=%2.4f\n",q2);
    //     pumi::update_submesh_and_cell_x2(pumi_obj, q2, jsub, jcell, &jsub, &jcell);
    //     printf("new-x2-submesh=%d\n", jsub);
    //     printf("new-x2-cell=%d\n", jcell);
    //     pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2);
    //     printf("new-x2-glob-cell=%d\n",kcell_x2 );
    //     printf("new-x2-W2=%2.4f\n\n",Wgh2 );
    //
    //
    //     pumi::calc_global_cellID_and_nodeID_fullmesh(pumi_obj, kcell_x1, kcell_x2, &global_cell, &lb_node, &lt_node);
    //     printf("new-2D-global-cell=%d\n",global_cell);
    //     printf("new-leftbottom-node=%d    new-lefttop-node=%d\n",lb_node, lt_node);
    //
    // });


    pumi::print_mesh_skeleton(pumi_obj);

    // Kokkos::parallel_for("bdry-test", 1, KOKKOS_LAMBDA (int j) {
    //     int Nx = pumi_obj.mesh(0).nsubmesh_x1;
    //     int Ny = pumi_obj.mesh(0).nsubmesh_x2;
    //     for (int i=0; i<2*Nx*Ny+Nx+Ny; i++){
    //         if (pumi_obj.mesh(0).is_bdry(i)){
    //             printf("edge-%2d is bdry\n",i );
    //         }
    //     }
    // });

    int N_part = 100;
    int N_step = 10;
    Kokkos::View<double*[2]> part_coords("particle-coordinates",N_part);
    // Kokkos::View<double**> part_x1_disp("part-x1-disp",N_step,N_part);
    // Kokkos::View<double**> part_x2_disp("part-x2-disp",N_step,N_part);
    Kokkos::View<bool*> part_activity("particle-activity",N_part);
    std::srand((unsigned)(std::time(nullptr)));

    double x1_min = pumi::get_global_left_coord(pumi_obj);
    double x1_max = pumi::get_global_right_coord(pumi_obj);
    double L_x1 = x1_max-x1_min;
    double x2_min = pumi::get_global_bottom_coord(pumi_obj);
    double x2_max = pumi::get_global_top_coord(pumi_obj);
    double L_x2 = x2_max-x2_min;

    Kokkos::View<double*[2]>::HostMirror h_part_coords = Kokkos::create_mirror_view(part_coords);
    Kokkos::View<bool*>::HostMirror h_part_activity = Kokkos::create_mirror_view(part_activity);
    // Kokkos::View<double**>::HostMirror h_part_x1_disp = Kokkos::create_mirror_view(part_x1_disp);
    // Kokkos::View<double**>::HostMirror h_part_x2_disp = Kokkos::create_mirror_view(part_x2_disp);

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
                h_part_activity(ipart) = true;
                part_set = true;
            }
        }

        // for (int i=0; i<N_step; i++){
        //     double x1_rand_val = (double) rand()/RAND_MAX;
        //     double x2_rand_val = (double) rand()/RAND_MAX;
        //     h_part_x1_disp(i,ipart) = x1_rand_val;
        //     h_part_x2_disp(i,ipart) = x2_rand_val;
        // }
    }

    Kokkos::deep_copy(part_coords, h_part_coords);
    Kokkos::deep_copy(part_activity, h_part_activity);
    // Kokkos::deep_copy(part_x1_disp, h_part_x1_disp);
    // Kokkos::deep_copy(part_x2_disp, h_part_x2_disp);

    write2file(h_part_coords, h_part_activity, N_part, 0);
    Kokkos::View<int*[2]> part_sub("particle-submesh-ID",N_part);
    Kokkos::View<int*[2]> part_cell("particle-cell-ID",N_part);
    Kokkos::parallel_for("particle-locate", N_part, KOKKOS_LAMBDA (int ipart) {
        int isub, jsub, icell, jcell;
        double q1 = part_coords(ipart,0);
        double q2 = part_coords(ipart,1);
        pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
        pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
        part_sub(ipart,0) = isub;
        part_sub(ipart,1) = jsub;
        part_cell(ipart,0) = icell;
        part_cell(ipart,1) = jcell;
    });

    int num_push = 0;
    Kokkos::Profiling::pushRegion("push_test_0");
    for (int istep=0; istep<N_step; istep++){
        Kokkos::deep_copy(h_part_activity, part_activity);
        for (int p=0; p<N_part; p++){
            num_push += h_part_activity(p);
        }

        Kokkos::parallel_for("particle-push-test-0", N_part, KOKKOS_LAMBDA (int ipart) {
            if (part_activity(ipart)){
                int isub, jsub, icell, jcell;
                double q1 = part_coords(ipart,0);
                double q2 = part_coords(ipart,1);
                // double dq1 = -10.0 + 20.0*part_x1_disp(istep,ipart);
                // double dq2 = -5.0 + 10.0*part_x2_disp(istep,ipart);
                double dq1 = 5.0;
                double dq2 = -3.0;

                isub = part_sub(ipart,0);
                jsub = part_sub(ipart,1);
                icell = part_cell(ipart,0);
                jcell = part_cell(ipart,1);

                q1 += dq1;
                q2 += dq2;

                if (q1 < x1_min || q1 > x1_max || q2 < x2_min || q2 > x2_max){
                    part_activity(ipart) = false;
                }
                else{
                    pumi::update_submesh_and_cell_x1(pumi_obj, q1, isub, icell, &isub, &icell);
                    pumi::update_submesh_and_cell_x2(pumi_obj, q2, jsub, jcell, &jsub, &jcell);
                }

                part_sub(ipart,0) = isub;
                part_sub(ipart,1) = jsub;
                part_cell(ipart,0) = icell;
                part_cell(ipart,1) = jcell;

                part_coords(ipart,0) = q1;
                part_coords(ipart,1) = q2;
            }
        });
        // Kokkos::deep_copy(h_part_coords, part_coords);
        // Kokkos::deep_copy(h_part_activity, part_activity);
        //
        // write2file(h_part_coords, h_part_activity, N_part, istep+1);
    }
    Kokkos::Profiling::popRegion();
    printf("Total number of particle pushes executed = %d\n",num_push );

    // int num_push = 0;
    // Kokkos::Profiling::pushRegion("push_test_1");
    // for (int istep=0; istep<N_step; istep++){
    //     Kokkos::deep_copy(h_part_activity, part_activity);
    //     for (int p=0; p<N_part; p++){
    //         num_push += h_part_activity(p);
    //     }
    //     Kokkos::parallel_for("particle-push-test-1", N_part, KOKKOS_LAMBDA (int ipart) {
    //         if (part_activity(ipart)){
    //             int isub, jsub, icell, jcell, bdry_hit;
    //             bool in_domain;
    //             double q1 = part_coords(ipart,0);
    //             double q2 = part_coords(ipart,1);
    //             // double dq1 = -10.0 + 20.0*part_x1_disp(istep,ipart);
    //             // double dq2 = -5.0 + 10.0*part_x2_disp(istep,ipart);
    //             double dq1 = 5.0;
    //             double dq2 = -3.0;
    //
    //             isub = part_sub(ipart,0);
    //             jsub = part_sub(ipart,1);
    //             icell = part_cell(ipart,0);
    //             jcell = part_cell(ipart,1);
    //
    //             pumi::push_particle(pumi_obj, q1, q2, dq1, dq2, &isub, &jsub, &icell, &jcell, &in_domain, &bdry_hit);
    //             part_activity(ipart) = in_domain;
    //             // if (!in_domain){
    //                 // printf("hit-edge=%2d -- is_bdry=%d\n",bdry_hit,pumi_obj.mesh(0).is_bdry(bdry_hit));
    //                 // part_activity(ipart) = false;
    //             // }
    //
    //             part_sub(ipart,0) = isub;
    //             part_sub(ipart,1) = jsub;
    //             part_cell(ipart,0) = icell;
    //             part_cell(ipart,1) = jcell;
    //
    //             part_coords(ipart,0) = q1+dq1;
    //             part_coords(ipart,1) = q2+dq2;
    //
    //         }
    //     });
    //     // Kokkos::deep_copy(h_part_coords, part_coords);
    //     // Kokkos::deep_copy(h_part_activity, part_activity);
    //
    //     // write2file(h_part_coords, h_part_activity, N_part, istep+1);
    // }
    // Kokkos::Profiling::popRegion();
    // printf("Total number of particle pushes executed = %d\n",num_push );

    // int num_push = 0;
    // Kokkos::Profiling::pushRegion("push_test_2");
    // for (int istep=0; istep<N_step; istep++){
    //     Kokkos::deep_copy(h_part_activity, part_activity);
    //     for (int p=0; p<N_part; p++){
    //         num_push += h_part_activity(p);
    //     }
    //     Kokkos::parallel_for("particle-push-test-2", N_part, KOKKOS_LAMBDA (int ipart) {
    //         if (part_activity(ipart)){
    //             int isub, jsub, icell, jcell, bdry_hit;
    //             bool in_domain;
    //             double q1 = part_coords(ipart,0);
    //             double q2 = part_coords(ipart,1);
    //             // double dq1 = -10.0 + 20.0*part_x1_disp(istep,ipart);
    //             // double dq2 = -5.0 + 10.0*part_x2_disp(istep,ipart);
    //             double dq1 = 5.0;
    //             double dq2 = -3.0;
    //
    //             isub = part_sub(ipart,0);
    //             jsub = part_sub(ipart,1);
    //             icell = part_cell(ipart,0);
    //             jcell = part_cell(ipart,1);
    //
    //             pumi::push_particle_v2(pumi_obj, q1, q2, dq1, dq2, &isub, &jsub, &icell, &jcell, &in_domain, &bdry_hit);
    //             part_activity(ipart) = in_domain;
    //             // if (!in_domain){
    //                 // printf("hit-edge=%2d -- is_bdry=%d\n",bdry_hit,pumi_obj.mesh(0).is_bdry(bdry_hit));
    //                 // part_activity(ipart) = false;
    //             // }
    //
    //             part_sub(ipart,0) = isub;
    //             part_sub(ipart,1) = jsub;
    //             part_cell(ipart,0) = icell;
    //             part_cell(ipart,1) = jcell;
    //
    //             part_coords(ipart,0) = q1+dq1;
    //             part_coords(ipart,1) = q2+dq2;
    //
    //         }
    //     });
    //     // Kokkos::deep_copy(h_part_coords, part_coords);
    //     // Kokkos::deep_copy(h_part_activity, part_activity);
    //
    //     // write2file(h_part_coords, h_part_activity, N_part, istep+1);
    // }
    // Kokkos::Profiling::popRegion();
    // printf("Total number of particle pushes executed = %d\n",num_push );

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

void write2file(Kokkos::View<double*[2]>::HostMirror hp_coords, Kokkos::View<bool*>::HostMirror hp_active, int N_part, int nstep){
    FILE *part_file;
    char part_filename[30];
    sprintf(part_filename,"part_coords_t%d.dat",nstep);
    part_file = fopen(part_filename,"w");

    for (int i=0; i<N_part; i++){
        fprintf(part_file, "%d %.5e %.5e\n", hp_active(i), hp_coords(i,0), hp_coords(i,1));
    }

    fclose(part_file);
}
