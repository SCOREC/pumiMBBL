#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_usage();

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

    // pumi::MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(mesh);
    // Kokkos::deep_copy(h_pumi_mesh, mesh);
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
    // double h_q1_new = h_q1 + rand_val_x1*L_x1*0.25;
    // double h_q2_new = h_q2 - rand_val_x2*L_x2*0.25;
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

    // Kokkos::parallel_for("node-ID-test", 1, KOKKOS_LAMBDA(const int){
    //
    //     int nodeID, isubmesh, jsubmesh;
    //
    //     for (int Jnp=pumi_obj.mesh(0).Nel_tot_x2; Jnp>=0; Jnp--){
    //
    //         for (jsubmesh=0; jsubmesh<pumi_obj.mesh(0).nsubmesh_x2; jsubmesh++){
    //             int submesh_min_node = pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    //             int submesh_max_node = pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative +
    //                                     pumi_obj.submesh_x2(jsubmesh)()->Nel;
    //             if (Jnp >= submesh_min_node && Jnp <= submesh_max_node){
    //                 break;
    //             }
    //         }
    //
    //         for (int Inp=0; Inp<=pumi_obj.mesh(0).Nel_tot_x1; Inp++ ){
    //
    //             for (isubmesh=0; isubmesh<pumi_obj.mesh(0).nsubmesh_x1; isubmesh++){
    //                 int submesh_min_node = pumi_obj.submesh_x1(isubmesh)()->Nel_cumulative;
    //                 int submesh_max_node = pumi_obj.submesh_x1(isubmesh)()->Nel_cumulative +
    //                                         pumi_obj.submesh_x1(isubmesh)()->Nel;
    //                 if (Inp >= submesh_min_node && Inp <= submesh_max_node){
    //                     break;
    //                 }
    //             }
    //
    //             nodeID = pumi::calc_global_nodeID(pumi_obj, isubmesh, jsubmesh, Inp, Jnp);
    //
    //             if (pumi_obj.mesh(0).isactive(isubmesh,jsubmesh)){
    //                 printf("%4d ", nodeID);
    //             }
    //             else{
    //                 pr// pumi::MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(mesh);
    // Kokkos::deep_copy(h_pumi_mesh, mesh);
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
    // std::cout << "\n\n";intf("%4d ", 0);
    //             }
    //         }
    //         printf("\n");
    //     }
    //     // for (int Jnp=pumi_obj.mesh(0).Nel_tot_x2; Jnp>=0; Jnp--){
    //     //     for (int isubmesh=0; isubmesh<pumi_obj.mesh(0).nsubmesh_x1; isubmesh++){
    //     //         printf("%3d            ", pumi_obj.mesh(0).nodeoffset(isubmesh,Jnp));
    //     //     }
    //     //     printf("\n");
    //     // }
    // });
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
