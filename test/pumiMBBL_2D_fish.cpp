#include "pumiMBBLGPU.hpp"
#include <stdio.h>
#include <Kokkos_Random.hpp>
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
    
    pumi::check_is_pumi_working();

    // Mesh specificatiaon - part is currently hard coded in parse_inputs()
    int nsubmesh_x1 = 3;
    double domain_x1_min = 0;
    int nsubmesh_x2 = 3;
    double domain_x2_min = 0;
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
    //end mesh specification/input

    /*
    'Bin' structure for particle inlet distribution
    Uses a nomralized cummulative distribution function for various distributions in order
    to introduce particles to the domain with position based on a distribution function
    Currently two options:
        0 - uniform
        1 - parabolic
    */
    int N_bins_inlet = 10;
    Kokkos::View<double*[2]> bins("bins",N_bins_inlet+1);
    Kokkos::View<double*[2]>::HostMirror h_bins = Kokkos::create_mirror_view(bins);
    for(int j = 0; j<N_bins_inlet+1;j++){
        double x = (j)/(double)N_bins_inlet; //normalized coordinate [0,1]
        h_bins(j,0) = x; //uniform distribution cdf
        h_bins(j,1) = 6*(x*x/2-x*x*x/3); //normalized parabolic distribution cdf
        
    }
    Kokkos::deep_copy(bins,h_bins);
    //end bins

    /*
    Properties of particle inlets for each inlet
    Indexed: inlets( <inlet number> , <property number> )
    */
    int num_inlets = 2;
    Kokkos::View<double*[6]> inlets("inlets",num_inlets);
    Kokkos::View<double*[6]>::HostMirror h_inlets = Kokkos::create_mirror_view(inlets);
    h_inlets(0,0) = 1;  h_inlets(1,0) = 0;  //orientation (0 = x1 algined, 1 = x2 aligned)
    h_inlets(0,1) = 1;  h_inlets(1,1) = 1;  //x_min for aligned direction
    h_inlets(0,2) = 3;  h_inlets(1,2) = 2;  //x_max for aligned direction
    h_inlets(0,3) = 0;  h_inlets(1,3) = 2;  //inlet coordinate for normal direction
    h_inlets(0,4) = 1;  h_inlets(1,4) = 0;  //distribution number
    h_inlets(0,5) = 10; h_inlets(1,5) = 10; //num particles to introduce per time step
    Kokkos::deep_copy(inlets,h_inlets);
    
    //Random number generator for particle location
    Kokkos::Random_XorShift64_Pool<> random_pool(12345);
    //end inlet setup

    /*
    Simulation parameters
    */
    int N_step = 300; //number of timesteps
    double dt = 0.01; //timestep size
    std::srand((unsigned)(std::time(nullptr)));

    double x1_min = pumi::get_global_x1_min_coord(pumi_obj);
    double x1_max = pumi::get_global_x1_max_coord(pumi_obj);
    double x2_min = pumi::get_global_x2_min_coord(pumi_obj);
    double x2_max = pumi::get_global_x2_max_coord(pumi_obj);
    
    // int num_push = 0;

    /*
    Setup particle data structure
    */
    int N_part_initial = 500; //initial number of ptcls in domain
    int N_part_current = N_part_initial; //total num ptcls in simulation so far (increments on time step)
    int n_part_in = 0; // number of ptcls to introduce each time step
    for(int i = 0; i < num_inlets; i++){
        n_part_in += h_inlets(i,5);
    }
    int m = 1; //buffering factor 
    int N_part_total = N_part_initial + N_step*n_part_in*m; //ptcl datastructure allocation
    Kokkos::View<pumi::ParticleData*> Partdata("particle-data",N_part_total);
    Kokkos::View<pumi::ParticleData*>::HostMirror h_Partdata = Kokkos::create_mirror_view(Partdata);
  
    /*
    Advection field setup
    velocity field import from CFD solve (python script orders data into global node index order)
    */
    Kokkos::View<double*[2]> velocity_field("velocity field", get_total_mesh_nodes(pumi_obj));
    Kokkos::View<double*[2]>::HostMirror h_velocity_field("velocity field host", get_total_mesh_nodes(pumi_obj));
    // txt file with velocity data for each mbbl node (interpolated from CFD mesh to mbbl mesh)
    // each row x1 and x2 advection values, ordered by mbbl global ordering
    std::ifstream vel_data("interp_velocity.txt");
    double val;
    for(int i = 0; i < get_total_mesh_nodes(pumi_obj); ++i){
        vel_data >> val;
        h_velocity_field(i,0) = val;
        vel_data >> val;
        h_velocity_field(i,1) = val;
    }
    vel_data.close();
    Kokkos::deep_copy(velocity_field,h_velocity_field);
    // end velocity field import

    /*
    Particle initialization
    */
    Kokkos::Profiling::pushRegion("particle initialization");
    printf("push-test\n");
    printf("prescribed velocity profile\n");

    /*
    Randomly introduce intial particles to domain
    */
    for (int ipart=0; ipart<N_part_initial; ipart++){
        pumi::Vector3 q = pumi::get_rand_point_in_mesh_host(pumi_obj);
        h_Partdata(ipart) = pumi::ParticleData(q[0],q[1]);
    }
    Kokkos::deep_copy(Partdata, h_Partdata);

    /*
    Locate particles 
    */
    Kokkos::parallel_for("particle-locate-0", N_part_initial, KOKKOS_LAMBDA (int ipart) {
        int isub, jsub, icell, jcell, submeshID, cellID;
        double q1 = Partdata(ipart).x1;
        double q2 = Partdata(ipart).x2;
        pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
        pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
        pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
        Partdata(ipart) = pumi::ParticleData(q1,q2,submeshID,cellID,true,-1);
        
    });
    Kokkos::deep_copy(h_Partdata,Partdata);
    Kokkos::Profiling::popRegion();

    /*
    Time loop with push
    */
    Kokkos::Profiling::pushRegion("push_test");
    for (int istep=0; istep<N_step; istep++){
        // for (int p=0; p<N_part_buffer; p++){
        //     bool part_active = h_Partdata(p).part_active;
        //     num_push += part_active;
        // }
        N_part_current = N_part_initial + n_part_in*istep;
        Kokkos::parallel_for("particle-push-test", N_part_current, KOKKOS_LAMBDA (int ipart) {
            // operate only on active particles
            bool part_active = Partdata(ipart).part_active;
            if (part_active){
                int isub, jsub, icell, jcell, kcell_x1, kcell_x2, submeshID, cellID,bdry_hit,bdry_faceID;
                bool in_domain;
                int global_cell, topleft_node, bottomleft_node;
                double fraction_done;
                double q1 = Partdata(ipart).x1;
                double q2 = Partdata(ipart).x2;
                submeshID = Partdata(ipart).submeshID;
                cellID = Partdata(ipart).cellID;

                //Calculate linear weighting for particles w.r.t its cell
                pumi::get_directional_submeshID_and_cellID(pumi_obj,submeshID,cellID,&isub,&icell,&jsub,&jcell);
                double Wgh2_x1, Wgh2_x2, Wgh1_x1, Wgh1_x2;
                pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                Wgh1_x1 = 1.0-Wgh2_x1;
                Wgh1_x2 = 1.0-Wgh2_x2;

                pumi::calc_global_cellID_and_nodeID(pumi_obj,isub,jsub,kcell_x1,kcell_x2,&global_cell,&bottomleft_node,&topleft_node);

                //velocity from CFD result and linear interpolation from cell nodes
                pumi::Vector3 velocity;
                velocity[0] = Wgh1_x1*Wgh1_x2*velocity_field(bottomleft_node,0) + Wgh2_x1*Wgh1_x2*velocity_field(bottomleft_node+1,0) +
                              Wgh1_x1*Wgh2_x2*velocity_field(topleft_node,0)    + Wgh2_x1*Wgh2_x2*velocity_field(topleft_node+1,0);
                velocity[1] = Wgh1_x1*Wgh1_x2*velocity_field(bottomleft_node,1) + Wgh2_x1*Wgh1_x2*velocity_field(bottomleft_node+1,1) +
                              Wgh1_x1*Wgh2_x2*velocity_field(topleft_node,1)    + Wgh2_x1*Wgh2_x2*velocity_field(topleft_node+1,1);

                double dq1 = dt*velocity[0];
                double dq2 = dt*velocity[1];

                pumi::Vector3 qnew = pumi::push_particle(pumi_obj, pumi::Vector3(q1,q2,0.0), pumi::Vector3(dq1,dq2,0.0), &isub, &jsub, &icell, &jcell,
                                    &in_domain, &bdry_hit, &fraction_done, &bdry_faceID);
                if (!in_domain){
                    //deactivate particles that have left the domain
                    q1 = qnew[0];
                    q2 = qnew[1];
                    pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
                    Partdata(ipart) = pumi::ParticleData(q1,q2,submeshID,cellID,false,bdry_faceID);
                }
                else{
                    //locate and update particles that remained in the domain
                    q1 = qnew[0];
                    q2 = qnew[1];
                    pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                    pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                    Wgh1_x1 = 1.0-Wgh2_x1;
                    Wgh1_x2 = 1.0-Wgh2_x2;
                    pumi::calc_global_cellID_and_nodeID(pumi_obj, isub, jsub, kcell_x1, kcell_x2, &global_cell, &bottomleft_node, &topleft_node);
                    pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
                    Partdata(ipart) = pumi::ParticleData(q1,q2,submeshID,cellID,true,-1,Wgh1_x1,Wgh2_x1,Wgh1_x2,Wgh2_x2);
                }
            }
        }); 

        Kokkos::parallel_for("particle initialization", n_part_in, KOKKOS_LAMBDA (int ipart) {
            // Introduce new particles here (introduce into their own memory range)
            int isub, jsub, icell, jcell, submeshID, cellID;
            double qx1 = 0;
            double qx2 = 0;
            auto generator = random_pool.get_state();
            //random number to select bin
            double bin_x = generator.drand(0.,1.);
            //random number to select position within bin
            double rand_x = generator.drand(0.,1.);
            random_pool.free_state(generator);

            int cum_num_in = 0;
            for(int k = 0; k < num_inlets; ++k){
                //find which inlet this particle should be added to 
                cum_num_in += inlets(k,5);
                if(ipart < cum_num_in){
                    int d = inlets(k,4);

                    // find which bin
                    int bin = 1;
                    for(int j = 1; j < N_bins_inlet; ++j){
                        bin = j;
                        if(bin_x < bins(j,d)) break;
                    }
                    
                    //global coordinate in algined direction
                    double p;
                    p = rand_x*(bins(bin,d)-bins(bin-1,d)) + bins(bin-1,d);
                    p = p*(inlets(k,2)-inlets(k,1)) + inlets(k,1);
                    //place on proper axis
                    int x2 = inlets(k,0);
                    int x1 = !x2;
                    // inlet(k,0) = 0 for x1 and = 1 for x2 axis aligned
                    qx1 = inlets(k,3) * x2 + p*x1;
                    qx2 = inlets(k,3) * x1 + p*x2;
                    break;
                }
            }

            //locate newly initialized particle
            pumi::locate_submesh_and_cell_x1(pumi_obj, qx1, &isub, &icell);
            pumi::locate_submesh_and_cell_x2(pumi_obj, qx2, &jsub, &jcell);
            pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
            Partdata(ipart + N_part_current) = pumi::ParticleData(qx1,qx2,submeshID,cellID,true,-1);
        });
        Kokkos::fence(); //synchronization step
        Kokkos::deep_copy(h_Partdata,Partdata);
        write2file(h_Partdata, N_part_current, istep); //write particle coords to file (not including newly introduce particles)
    }
    Kokkos::Profiling::popRegion();

    
    // printf("Total number of particle pushes executed in Test-0 = %d\n",num_push );
    //end push test

    pumi::free_mbbl(pumi_obj);
  }
  Kokkos::finalize();

  return 0;
}

KOKKOS_INLINE_FUNCTION
pumi::Vector3 trimmed_vortex(pumi::Vector3 pos){
    double x_max = 0.5;
    pumi::Vector3 lamb = lamb_oseen_vortex(pos);

    double v_x = lamb[0]*(1-(pos[0]/x_max)*(pos[0]/x_max));
    double v_y = lamb[1]*(1-(pos[1]/x_max)*(pos[1]/x_max));

    return pumi::Vector3(v_x,v_y,0);
}

KOKKOS_INLINE_FUNCTION
pumi::Vector3 lamb_oseen_vortex(pumi::Vector3 pos){
    /*
    v_r = 0
    v_{theta} = \frac{\Gamma}{2\pi r}(1-e^{\frac{-r^2}{4\nu t}})
    v_z = 0
    */

    //clockwise lamb oseen vortex
    double Gamma = 5;
    double nu = 3.5;
    double t = 1;
    double scaling = 100; 
    double pi = 3.14159265;
    
    double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    // scale radius to keep all interesting behavior of vortex in the unit square domain
    double r_scaled = scaling*r;
    double v_theta = Gamma/2/pi/r_scaled*(1-exp(-1*r_scaled*r_scaled/4/nu/t));

    double x_norm = pos[0]/r; //normalized coordinate of particle
    double y_norm = pos[1]/r;
    double v_x = y_norm*v_theta;
    double v_y = -1*x_norm*v_theta;

    return pumi::Vector3(v_x,v_y,0);
}

void parse_inputs(int , char* argv[], pumi::Mesh_Inputs *pumi_inputs)
{
    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int nsubmesh_x2 = pumi_inputs->nsubmesh_x2;
    // reading submesh meshtypes
    char all_submesh_flag_x1[MAX_SUBMESHES*10];
    char each_submesh_flag_x1[MAX_SUBMESHES][10];
    // strcpy(all_submesh_flag_x1, argv[3]);
    strcpy(all_submesh_flag_x1, "uniform,uniform,uniform");

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
    // strcpy(all_p1_submesh_x1, argv[4]);
    strcpy(all_p1_submesh_x1, "1.0,1.0,1.0");

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
    // strcpy(all_p2max_submesh_x1, argv[5]);
    strcpy(all_p2max_submesh_x1, "0.05,0.05,0.05");

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
    // strcpy(all_p2min_submesh_x1, argv[6]);
    strcpy(all_p2min_submesh_x1, "0.05,0.05,0.05");

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
    // strcpy(all_arb_submesh_x1, argv[7]);
    strcpy(all_arb_submesh_x1, "NA,NA,NA");

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
    // strcpy(all_submesh_flag_x2, argv[10]);
    strcpy(all_submesh_flag_x2, "uniform,uniform,uniform");

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
    // strcpy(all_p1_submesh_x2, argv[11]);
    strcpy(all_p1_submesh_x2, "1.0,1.0,1.0");

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
    // strcpy(all_p2max_submesh_x2, argv[12]);
    strcpy(all_p2max_submesh_x2, "0.05,0.05,0.05");

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
    // strcpy(all_p2min_submesh_x2, argv[13]);
    strcpy(all_p2min_submesh_x2, "0.05,0.05,0.05");

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
    // strcpy(all_arb_submesh_x2, argv[14]);
    strcpy(all_arb_submesh_x2, "NA,NA,NA");

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
    // strcpy(all_submesh_isactive, argv[15]);
    strcpy(all_submesh_isactive, "0,1,0,1,1,1,1,1,0");

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
