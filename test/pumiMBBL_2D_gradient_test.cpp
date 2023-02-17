#include "pumiMBBLGPU.hpp"
void parse_inputs(int argc, char* argv[], pumi::Mesh_Inputs *pumi_inputs);
void print_parsed_inputs(pumi::Mesh_Inputs *pumi_inputs);
void print_usage();
void write2file(Kokkos::View<pumi::ParticleData*>::HostMirror hp, int N_part, int nstep);
void write2file_2(pumi::DoubleView::HostMirror hp, int N_part);
void write2file_3(pumi::Vector3View::HostMirror hp, int N_part);
KOKKOS_INLINE_FUNCTION pumi::Vector3 lamb_oseen_vortex(pumi::Vector3 pos);
KOKKOS_INLINE_FUNCTION pumi::Vector3 trimmed_vortex(pumi::Vector3 pos);

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );
  {
    /*
    1 submesh block by 1 submesh block, 100x100 cells, unit square, uniform mesh
    ./pumiMBBL2D_gradient 1 -0.5 "uniform" "1.0" "100" "0.01" "NA" 1 -0.5 "uniform" "1.0" "100" "0.01" "NA" "1" 
    ./pumiMBBL2D_gradient 3 -1 "minBL,uniform,maxBL" "0.5,1.0,0.5" "0.05,0.01,0.05" "0.001,0.01,0.001" "NA,NA,NA" 3 -1 "minBL,uniform,maxBL" "0.5,1.0,0.5" "0.02,0.02,0.02" "0.01,0.02,0.01" "NA,NA,NA" "1,1,1,1,1,1,1,1,1"
    */

    if (argc != 16)
    {
        print_usage();
    }
    pumi::check_is_pumi_working();
    int nsubmesh_x1 = atoi( argv[1] );
    double domain_x1_min = atof( argv[2] );
    int nsubmesh_x2 = atoi( argv[8] );
    double domain_x2_min = atof( argv[9] );
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
    pumi_options.print_mesh_connectivity_option = pumi::print_mesh_connectivity_ON;

    pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);

    printf("Mesh volume = %2.2e\n",pumi::get_mesh_volume(pumi_obj) );
    pumi::inputs_deallocate(pumi_inputs);

    pumi::print_mesh_skeleton(pumi_obj);

    // pumi::print_blockwise_nodeIDs(pumi_obj);
    // pumi::print_node_submeshID(pumi_obj);
    // pumi::print_fullmesh_nodeIDs(pumi_obj);

    int Nnp_total = pumi_obj.mesh.Nnp_total;
    pumi::DoubleView phi = pumi::DoubleView("phi",Nnp_total);
    pumi::DoubleView::HostMirror h_phi = Kokkos::create_mirror_view(phi);
    pumi::Vector3View coords = pumi::Vector3View("coords", Nnp_total);
    pumi::Vector3View::HostMirror h_coords =Kokkos::create_mirror_view(coords);
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Ny = pumi_obj.mesh.nsubmesh_x2;
    for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
        for (int jsubmesh=1; jsubmesh<=Ny; jsubmesh++){
            if (is_block_active_host(pumi_obj,isubmesh,jsubmesh)){
                int Nnp_y = pumi_obj.host_submesh_x2[jsubmesh]->Nel+1;
                int Nnp_x = pumi_obj.host_submesh_x1[isubmesh]->Nel+1;
                for (int inp=0; inp<Nnp_x; inp++){
                    for (int jnp=0; jnp<Nnp_y; jnp++){
                        int Inp = inp + pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
                        int Jnp = jnp + pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
                        double x1_coord = pumi_obj.host_submesh_x1[isubmesh]->node_coords_host(inp);
                        double x2_coord = pumi_obj.host_submesh_x2[jsubmesh]->node_coords_host(jnp);
                        int nodeID = pumi::get_global_nodeID_2D(pumi_obj,Inp,Jnp);
                        h_coords(nodeID)[0] = x1_coord;
                        h_coords(nodeID)[1] = x2_coord;
                        h_phi(nodeID) = x1_coord*x2_coord;
                        // h_phi(nodeID) = 1.0;
                    }
                }
            }
        }
    }
    Kokkos::deep_copy(phi,h_phi);
    Kokkos::deep_copy(coords,h_coords);
    // write2file_2(h_phi,Nnp_total);

    pumi::Vector3View phi_grad_v2;
    pumi::Vector3View phi_grad_conv;

    int nstep=10000;
    // Kokkos::Profiling::pushRegion("grad_v1");
    // for (int istep=0; istep<nstep; istep++){
    //     // printf("t%d\n",istep );
    //     phi_grad = pumi::compute_2D_field_gradient(pumi_obj,phi);
    // }
    // Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("grad_convex");
    for (int istep=0; istep<nstep; istep++){
        // printf("t%d\n",istep );
        phi_grad_conv = pumi::compute_2D_field_gradient_convex(pumi_obj,phi);
    }
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("grad_v2");
    for (int istep=0; istep<nstep; istep++){
        // printf("t%d\n",istep );
        phi_grad_v2 = pumi::compute_2D_field_gradient_v2(pumi_obj,phi);
    }
    Kokkos::Profiling::popRegion();

    // Kokkos::parallel_for("verify solutions", Nnp_total, KOKKOS_LAMBDA (int i){
    //     double diff_conv_x1 = abs(coords(i)[0] - phi_grad_conv(i)[1]);
    //     double diff_conv_x2 = abs(coords(i)[1] - phi_grad_conv(i)[0]);
    //     double diff_v2_x1 = abs(coords(i)[0] - phi_grad_v2(i)[1]);
    //     double diff_v2_x2 = abs(coords(i)[1] - phi_grad_v2(i)[0]);

    //     if(diff_conv_x1 > 1e-8){
            
    //     }
    // });

    pumi::Vector3View::HostMirror h_phi_grad_conv = Kokkos::create_mirror_view(phi_grad_conv);
    Kokkos::deep_copy(h_phi_grad_conv,phi_grad_conv);
    pumi::Vector3View::HostMirror h_phi_grad_v2 = Kokkos::create_mirror_view(phi_grad_v2);
    Kokkos::deep_copy(h_phi_grad_v2,phi_grad_v2);
    
    write2file_3(h_phi_grad_conv,Nnp_total);
    // Kokkos::parallel_for("print-grad",1,KOKKOS_LAMBDA (const int){
    //     for (int i=0; i<Nnp_total; i++){
    //         printf("%.16e %.16e\n",phi_grad_new(i)[0]-phi_grad_uni(i)[0],phi_grad_new(i)[1]-phi_grad_uni(i)[1] );
    //     }
    // });

    // pumi::DoubleView phi_density = pumi::compute_2D_field_density(pumi_obj,phi);
    //
    // pumi::DoubleView phi_density_v2 = pumi::DoubleView("new-density-arr",pumi_obj.mesh.Nnp_total);
    // pumi::DoubleView::HostMirror h_phi_density_v2 = Kokkos::create_mirror_view(phi_density_v2);
    //
    // int knode=0;
    // for (int jnp=0; jnp<=pumi_obj.mesh.Nel_tot_x2; jnp++){
    //     for (int inp=0; inp<=pumi_obj.mesh.Nel_tot_x1; inp++){
    //         bool on_bdry, in_domain;
    //         int bdry_tag, bdry_dim;
    //         double cov;
    //         pumi::where_is_node(pumi_obj, inp, jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
    //         if (in_domain){
    //             cov = pumi::return_covolume(pumi_obj, inp, jnp);
    //             h_phi_density_v2(knode) = h_phi(knode)/cov;
    //             knode++;
    //         }
    //     }
    // }
    // Kokkos::deep_copy(phi_density_v2,h_phi_density_v2);
    //
    // Kokkos::parallel_for("print-density",1,KOKKOS_LAMBDA (const int){
    //     for (int i=0; i<Nnp_total; i++){
    //         printf("ID=%d  diff=%2.2e\n",i, fabs(phi_density(i)-phi_density_v2(i)));
    //     }
    // });
    // int Nel_total = pumi_obj.mesh.Nel_total;
    // Kokkos::parallel_for("print-density",1,KOKKOS_LAMBDA (const int){
    //     for (int i=0; i<Nel_total; i++){
    //         int kcell_x2 = i/pumi_obj.mesh.Nel_tot_x1;
    //         int kcell_x1 = i - kcell_x2*pumi_obj.mesh.Nel_tot_x1;
    //         int isub, jsub, icell, jcell;
    //         pumi::get_x1_submeshID_and_localcellID_of_x1_elem(pumi_obj, kcell_x1, &isub, &icell);
    //         pumi::get_x2_submeshID_and_localcellID_of_x2_elem(pumi_obj, kcell_x2, &jsub, &jcell);
    //         double inode, jnode, dx1, dx2;
    //         inode = pumi::get_x1_node_coord_in_submesh(pumi_obj, isub, icell);
    //         jnode = pumi::get_x2_node_coord_in_submesh(pumi_obj, jsub, jcell);
    //         dx1 = pumi::get_x1_elem_size_in_submesh(pumi_obj, isub, icell);
    //         dx2 = pumi::get_x2_elem_size_in_submesh(pumi_obj, jsub, jcell);
    //         printf("%1.4e %1.4e %1.4e %1.4e\n",inode,inode+dx1,jnode,jnode+dx2);
    //     }
    // });

    
    pumi::free_mbbl(pumi_obj);


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
    printf("    ./install/bin/pumiMBBL2D_Demo 3 0.0 \"minBL,uniform,maxBL\" \"20.0,10.0,20.0\" \"3.0,1.0,3.0\" \"1.0,1.0,1.0\" \"NA,NA,NA\" 3 1.0 \"maxBL,uniform,minBL\" \"50.0,20.0,50.0\" \"4.0,1.0,4.0\" \"1.0,2.0,1.0\" \"NA,NA,NA\" \"1,1,1,1,1,1,1,1,1\" \n\n");
    printf("  E.g.#2 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 4 0.0 \"minBL,uniform,uniform,maxBL\" \"10.0,5.0,5.0,10.0\" \"3.0,1.0,1.0,3.0\" \"1.0,1.0,1.0,1.0\" \"NA,NA,NA,NA\" 4 1.0 \"maxBL,uniform,uniform,minBL\" \"20.0,20.0,20.0,20.0\" \"4.0,1.0,1.0,4.0\" \"1.0,2.0,2.0,1.0\" \"NA,NA,NA,NA\" \"1,0,1,1,1,0,0,1,1,1,1,1,0,1,0,1\" \n\n");
    printf("  E.g.#3 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 1 1.0 \"uniform\" \"20.0\" \"1.0\" \"1.0\" \"NA\" 1 1.0 \"uniform\" \"20.0\" \"1.0\" \"1.0\" \"NA\" \"1\"\n\n");
    printf("  E.g.#4 \n\n");
    printf("    ./install/bin/pumiMBBL2D_Demo 4 0.0 \"minBL,uniform,arbitrary,maxBL\" \"10.0,5.0,5.0,10.0\" \"3.0,1.0,1.0,3.0\" \"1.0,1.0,1.0,1.0\" \"NA,NA,x1-size.dat,NA\" 4 1.0 \"maxBL,arbitrary,uniform,minBL\" \"20.0,20.0,20.0,20.0\" \"4.0,1.0,1.0,4.0\" \"1.0,2.0,2.0,1.0\" \"NA,x2-size.dat,NA,NA\" \"1,0,1,1,1,0,0,1,1,1,1,1,0,1,0,1\" \n\n");
    Kokkos::finalize();
    exit(0);
}

void write2file(Kokkos::View<pumi::ParticleData*>::HostMirror hp, int N_part, int nstep){
    FILE *part_file;
    char part_filename[30];
    sprintf(part_filename,"part_coords_t%d.dat",nstep);
    part_file = fopen(part_filename,"w");
    int skip = 1;
    for (int i=0; i<N_part; i=i+skip){
        int part_active = hp(i).part_active;
        int subID = hp(i).submeshID;
        int exit_faceID = hp(i).exit_faceID;
        fprintf(part_file, "%d %.5e %.5e %d %d\n", part_active, hp(i).x1, hp(i).x2, subID, exit_faceID);
    }

    fclose(part_file);
}

void write2file_2(pumi::DoubleView::HostMirror hp, int N_part){
    FILE *part_file;
    char part_filename[30];
    sprintf(part_filename,"phi_data.dat");
    part_file = fopen(part_filename,"w");
    for (int i=0; i<N_part; i=i+1){
        fprintf(part_file, "%.16e\n", hp(i));
    }

    fclose(part_file);
}

void write2file_3(pumi::Vector3View::HostMirror hp, int N_part){
    FILE *part_file;
    char part_filename[30];
    sprintf(part_filename,"phi_grad_data.dat");
    part_file = fopen(part_filename,"w");
    for (int i=0; i<N_part; i=i+1){
        fprintf(part_file, "%.16e %.16e\n", hp(i)[0], hp(i)[1]);
    }

    fclose(part_file);
}
